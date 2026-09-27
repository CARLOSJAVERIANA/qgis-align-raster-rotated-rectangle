# -*- coding: utf-8 -*-
# QGIS Align Raster to Rotated Rectangle
#
# Copyright (c) 2026 Carlos Acosta
# Licensed under the BSD 3-Clause License.
# See LICENSE in the project repository for full license information.
"""
Alinear raster a rectángulo rotado y exportar PNG
=============================================================

Script de Processing para QGIS 3.40+ (GDAL >= 3.8).

Crea una rejilla raster cuyos ejes coinciden con los lados de un rectángulo
rotado, remuestrea el GeoTIFF original sobre esa rejilla (GDAL Warp) y
exporta:

  * un GeoTIFF alineado que conserva todas las bandas, el tipo de dato,
    nodata, alfa, metadatos y tabla de colores; si el origen no tiene
    nodata ni alfa, se añade una máscara interna de validez;
  * un PNG rectangular (Gray/RGB + alfa) sin cuñas transparentes, con
    world file (.wld) y .aux.xml con el CRS.
"""

import os
from math import sqrt, acos, degrees, atan2, sin, radians, ceil
from xml.sax.saxutils import escape as xml_escape

import numpy as np
from osgeo import gdal, ogr

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterFileDestination,
    QgsProcessingOutputNumber,
    QgsCoordinateTransform,
    QgsProject,
    QgsWkbTypes,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingUtils,
)


# --------------------------------------------------------------------------
# Utilidades geométricas puras (trabajan con tuplas (x, y))
# --------------------------------------------------------------------------

MAX_SAMPLE_PIXELS = 4_000_000  # muestra máxima solo para percentiles 2–98 %
MASK_BLOCK_PIXELS = 1_048_576  # memoria acotada al construir máscaras


def _vsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _norm(v):
    return sqrt(v[0] * v[0] + v[1] * v[1])


def _dot(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1]


def _cross(v1, v2):
    return v1[0] * v2[1] - v1[1] * v2[0]


def _angle_deg(v1, v2):
    n1 = _norm(v1)
    n2 = _norm(v2)
    if n1 <= 0.0 or n2 <= 0.0:
        return 0.0
    c = max(-1.0, min(1.0, _dot(v1, v2) / (n1 * n2)))
    return degrees(acos(c))


def _clean_ring(ring, angle_tol_deg):
    """Elimina vértice de cierre, duplicados y vértices colineales."""
    pts = [(float(x), float(y)) for x, y in ring]
    if len(pts) < 3:
        raise QgsProcessingException("El anillo exterior tiene menos de 3 vértices.")

    scale = max(1.0, max(max(abs(x), abs(y)) for x, y in pts))
    eps = scale * 1e-9  # ~1 mm a 1000 km; robusto frente a redondeo float64

    out = []
    for p in pts:
        if not out or _norm(_vsub(p, out[-1])) > eps:
            out.append(p)
    while len(out) > 1 and _norm(_vsub(out[0], out[-1])) <= eps:
        out.pop()

    sin_tol = sin(radians(angle_tol_deg))
    changed = True
    while changed and len(out) > 3:
        changed = False
        n = len(out)
        for i in range(n):
            a = _vsub(out[i], out[i - 1])
            b = _vsub(out[(i + 1) % n], out[i])
            na, nb = _norm(a), _norm(b)
            if na <= eps or nb <= eps:
                del out[i]
                changed = True
                break
            if abs(_cross(a, b)) / (na * nb) <= sin_tol and _dot(a, b) > 0.0:
                del out[i]
                changed = True
                break
    return out


def _validate_rectangle(pts, angle_tol, length_tol_pct):
    """Devuelve (vecs, lengths, angles) o lanza QgsProcessingException."""
    if len(pts) != 4:
        raise QgsProcessingException(
            "El polígono debe tener exactamente 4 esquinas. Se encontraron {} "
            "tras eliminar vértices duplicados y colineales.".format(len(pts))
        )

    vecs, lengths = [], []
    for i in range(4):
        v = _vsub(pts[(i + 1) % 4], pts[i])
        L = _norm(v)
        if L <= 0.0:
            raise QgsProcessingException("El rectángulo contiene un lado de longitud cero.")
        vecs.append(v)
        lengths.append(L)

    angles = [_angle_deg(vecs[i], vecs[(i + 1) % 4]) for i in range(4)]

    bad = [(i, a) for i, a in enumerate(angles) if abs(a - 90.0) > angle_tol]
    if bad:
        detail = ", ".join("V{}={:.4f}°".format(i, a) for i, a in bad)
        raise QgsProcessingException(
            "El polígono no cumple la tolerancia rectangular de ±{:.3f}°: {}".format(angle_tol, detail)
        )

    sin_tol = sin(radians(angle_tol))
    for a_idx, b_idx in ((0, 2), (1, 3)):
        rel_cross = abs(_cross(vecs[a_idx], vecs[b_idx])) / (lengths[a_idx] * lengths[b_idx])
        if rel_cross > sin_tol:
            raise QgsProcessingException(
                "Los lados opuestos {} y {} no son suficientemente paralelos.".format(a_idx, b_idx)
            )
        denom = max(lengths[a_idx], lengths[b_idx])
        rel = abs(lengths[a_idx] - lengths[b_idx]) / denom * 100.0
        if rel > length_tol_pct:
            raise QgsProcessingException(
                "Las longitudes de los lados opuestos {} y {} difieren {:.4f} %, "
                "superando la tolerancia de {:.4f} %.".format(a_idx, b_idx, rel, length_tol_pct)
            )

    return vecs, lengths, angles


def _build_frame(pts, vecs, lengths, edge_mode):
    """
    Sistema de referencia local de la imagen:
      u = eje horizontal (columnas), v = eje vertical (filas hacia abajo).
    """
    pair0 = 0.5 * (lengths[0] + lengths[2])
    pair1 = 0.5 * (lengths[1] + lengths[3])

    if edge_mode == 0:
        edge_idx = 0 if pair0 >= pair1 else 1
    elif edge_mode == 1:
        edge_idx = 0 if pair0 <= pair1 else 1
    else:
        edge_idx = edge_mode - 2

    ux = vecs[edge_idx][0] / lengths[edge_idx]
    uy = vecs[edge_idx][1] / lengths[edge_idx]

    # El eje horizontal apunta siempre "hacia la derecha" (ux >= 0) para que
    # la imagen no salga invertida; el sentido del borde elegido no importa.
    if ux < 0.0 or (abs(ux) < 1e-15 and uy < 0.0):
        ux, uy = -ux, -uy

    vx, vy = uy, -ux  # u rotado -90°: filas crecientes hacia abajo (CRS con Y arriba)

    s_vals = sorted(x * ux + y * uy for x, y in pts)
    t_vals = sorted(x * vx + y * vy for x, y in pts)

    s_min = 0.5 * (s_vals[0] + s_vals[1])
    s_max = 0.5 * (s_vals[2] + s_vals[3])
    t_min = 0.5 * (t_vals[0] + t_vals[1])
    t_max = 0.5 * (t_vals[2] + t_vals[3])

    width = s_max - s_min
    height = t_max - t_min
    if width <= 0.0 or height <= 0.0:
        raise QgsProcessingException("Las dimensiones calculadas del rectángulo no son válidas.")

    theta = degrees(atan2(uy, ux))
    theta = ((theta + 90.0) % 180.0) - 90.0

    return {
        "edge_idx": edge_idx,
        "u": (ux, uy),
        "v": (vx, vy),
        "origin": (ux * s_min + vx * t_min, uy * s_min + vy * t_min),
        "width": width,
        "height": height,
        "theta": theta,
    }


def _source_pixel_step_in_direction(gt, ux, uy):
    """Tamaño de píxel del raster origen medido a lo largo de la dirección (ux, uy)."""
    a, b, c, d = gt[1], gt[2], gt[4], gt[5]
    det = a * d - b * c
    if abs(det) < 1e-30:
        raise QgsProcessingException("El GeoTransform del raster es singular o inválido.")
    p_col = (d * ux - b * uy) / det
    p_row = (-c * ux + a * uy) / det
    pix_per_unit = sqrt(p_col * p_col + p_row * p_row)
    if pix_per_unit <= 0.0:
        raise QgsProcessingException("No fue posible calcular la resolución equivalente.")
    return 1.0 / pix_per_unit


def _geotransform_polygon(gt, w, h):
    """Polígono (QgsGeometry) con la huella del raster, válido también con GT rotado."""
    corners = []
    for col, row in ((0, 0), (w, 0), (w, h), (0, h)):
        x = gt[0] + col * gt[1] + row * gt[2]
        y = gt[3] + col * gt[4] + row * gt[5]
        corners.append(QgsPointXY(x, y))
    corners.append(corners[0])
    return QgsGeometry.fromPolygonXY([corners])


# --------------------------------------------------------------------------
# Algoritmo
# --------------------------------------------------------------------------

class AlignRasterToRotatedRectangle(QgsProcessingAlgorithm):

    INPUT = "INPUT"
    RECTANGLE = "RECTANGLE"
    EDGE_MODE = "EDGE_MODE"
    FALLBACK_OMBB = "FALLBACK_OMBB"
    ANGLE_TOL = "ANGLE_TOL"
    LENGTH_TOL = "LENGTH_TOL"
    PIXEL_SIZE = "PIXEL_SIZE"
    RESAMPLING = "RESAMPLING"
    STRICT_SOURCE_CUTLINE = "STRICT_SOURCE_CUTLINE"
    WARP_MEMORY_MB = "WARP_MEMORY_MB"
    ALPHA_MAX = "ALPHA_MAX"
    PNG_MODE = "PNG_MODE"
    OUTPUT_TIF = "OUTPUT_TIF"
    OUTPUT_PNG = "OUTPUT_PNG"

    ANGLE_DEG = "ANGLE_DEG"
    WIDTH_PX = "WIDTH_PX"
    HEIGHT_PX = "HEIGHT_PX"

    RESAMPLING_NAMES = ["near", "bilinear", "cubic", "cubicspline", "lanczos"]

    PNG_EXACT = 0
    PNG_VISUAL_MINMAX = 1
    PNG_VISUAL_PCT = 2
    PNG_VISUAL_SHARED = 3
    PNG_VISUAL_UNIT = 4
    PNG_VISUAL_AUTO = 5

    # ---------------- Metadatos ----------------

    def name(self):
        return "align_raster_to_rotated_rectangle"

    def displayName(self):
        return "Alinear raster a rectángulo rotado y exportar PNG"

    def group(self):
        return "Raster personalizado"

    def groupId(self):
        return "raster_personalizado"

    def tags(self):
        return ["raster", "rotar", "rotate", "rectángulo", "png", "warp", "alinear"]

    def shortHelpString(self):
        return (
            "Crea una rejilla raster alineada con un rectángulo rotado, remuestrea el "
            "GeoTIFF original sobre esa rejilla (GDAL Warp) y exporta un GeoTIFF alineado "
            "con todas las bandas más un PNG rectangular (Gray/RGB + alfa) sin cuñas "
            "transparentes, con world file.\n\n"
            "<b>Rectángulo:</b> capa con exactamente 1 polígono (use «solo entidades "
            "seleccionadas» si la capa tiene varios). Se admiten vértice de cierre, "
            "vértices duplicados o colineales, multipartes de una parte y geometrías "
            "curvas. Si el polígono no pasa las tolerancias y está activada la opción "
            "OMBB, se usa su rectángulo envolvente mínimo orientado.\n\n"
            "<b>Borde:</b> los índices P0…P3 se refieren a las esquinas tras la limpieza "
            "(en el orden del anillo). El borde elegido pasa a ser horizontal; el sentido "
            "se normaliza para que la imagen no quede invertida.\n\n"
            "<b>Validez:</b> se combina por píxel la validez de todas las bandas y el "
            "alfa existente. El PNG usa un alfa construido con esa misma máscara; así, "
            "NoData distinto entre bandas nunca se muestra como un píxel opaco.\n\n"
            "<b>Recorte estricto de fuente:</b> por defecto el rectángulo también actúa "
            "como máscara antes del remuestreo. Desactívelo solamente si quiere que los "
            "bordes interpolados usen vecinos inmediatamente exteriores.\n\n"
            "<b>PNG exacto</b> requiere Byte o UInt16; para otros tipos se aplica "
            "automáticamente un modo visual que preserva mejor el balance RGB. "
            "El modo 0–1 sirve para TIFF Float32 cuyos colores ya son valores "
            "RGB normalizados. Los rasters con tabla de colores "
            "se expanden a RGB en el PNG."
        )

    def createInstance(self):
        return AlignRasterToRotatedRectangle()

    # ---------------- Parámetros ----------------

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT, "GeoTIFF de entrada"))

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.RECTANGLE,
                "Rectángulo rotado (exactamente 1 polígono)",
                [QgsProcessing.TypeVectorPolygon],
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.EDGE_MODE,
                "Borde / dirección que se convertirá en horizontal",
                options=[
                    "Automático: lado más largo",
                    "Automático: lado más corto",
                    "Borde 0: P0 → P1",
                    "Borde 1: P1 → P2",
                    "Borde 2: P2 → P3",
                    "Borde 3: P3 → P0",
                ],
                defaultValue=0,
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.STRICT_SOURCE_CUTLINE,
                "Limitar también la fuente al interior del rectángulo antes de interpolar",
                defaultValue=True,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.WARP_MEMORY_MB,
                "Memoria máxima para GDAL Warp (MB)",
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=256,
                minValue=64,
                maxValue=4096,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.ALPHA_MAX,
                "Máximo del alfa de entrada (0 = automático: NBITS/tipo de dato)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.FALLBACK_OMBB,
                "Si el polígono no es rectangular, usar su rectángulo envolvente mínimo orientado",
                defaultValue=False,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.ANGLE_TOL,
                "Tolerancia angular respecto a 90° (grados)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=2.0,
                minValue=0.01,
                maxValue=20.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.LENGTH_TOL,
                "Tolerancia de longitud entre lados opuestos (%)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=2.0,
                minValue=0.01,
                maxValue=50.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.PIXEL_SIZE,
                "Tamaño de píxel de salida (0 = automático, conserva la densidad del raster)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.RESAMPLING,
                "Remuestreo",
                options=["Nearest neighbour", "Bilinear", "Cubic", "Cubic spline", "Lanczos"],
                defaultValue=1,
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.PNG_MODE,
                "Conversión de valores para PNG",
                options=[
                    "Exacta: conservar Byte/UInt16",
                    "Visual: autoescalar cada banda a Byte (min–max de píxeles válidos)",
                    "Visual: autoescalar cada banda a Byte (corte acumulativo 2–98 %)",
                    "Visual RGB: escala común 2–98 % (conserva balance entre canales)",
                    "Visual RGB: valores 0–1 directamente a 0–255 (sin autoescalado)",
                    "Visual RGB: automático (0–1 si procede; si no, escala común)",
                ],
                defaultValue=5,
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT_TIF,
                "GeoTIFF alineado (conserva todas las bandas)",
            )
        )

        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT_PNG,
                "PNG rectangular",
                "PNG (*.png)",
            )
        )

        self.addOutput(QgsProcessingOutputNumber(self.ANGLE_DEG, "θ del eje horizontal (grados)"))
        self.addOutput(QgsProcessingOutputNumber(self.WIDTH_PX, "Ancho de salida (px)"))
        self.addOutput(QgsProcessingOutputNumber(self.HEIGHT_PX, "Alto de salida (px)"))

    # ---------------- Utilidades de E/S ----------------

    @staticmethod
    def _remove_existing(path):
        if not path:
            return
        base = os.path.splitext(path)[0]
        for p in (
            path, path + ".aux.xml", path + ".msk", path + ".ovr", path + ".msk.ovr",
            base + ".wld", base + ".pgw", base + ".tfw", path + "w",
        ):
            try:
                if os.path.exists(p):
                    os.remove(p)
            except Exception:
                pass

    @staticmethod
    def _unlink_vsimem(path):
        try:
            gdal.Unlink(path)
        except Exception:
            pass

    @staticmethod
    def _vrt_source_path(path):
        """Ruta para <SourceFilename>: absoluta en disco; intacta si es virtual (/vsimem/, /vsizip/...)."""
        if path.startswith("/vsi"):
            return path
        return os.path.abspath(path)

    @staticmethod
    def _copy_band_metadata(src_band, dst_band):
        md = src_band.GetMetadata()
        if md:
            dst_band.SetMetadata(md)
        desc = src_band.GetDescription()
        if desc:
            dst_band.SetDescription(desc)
        unit = src_band.GetUnitType()
        if unit:
            dst_band.SetUnitType(unit)
        scale = src_band.GetScale()
        if scale is not None:
            dst_band.SetScale(scale)
        offset = src_band.GetOffset()
        if offset is not None:
            dst_band.SetOffset(offset)
        ci = src_band.GetColorInterpretation()
        if ci is not None:
            dst_band.SetColorInterpretation(ci)
        ct = src_band.GetColorTable()
        if ct is not None:
            dst_band.SetColorTable(ct)
        cats = src_band.GetCategoryNames()
        if cats:
            dst_band.SetCategoryNames(cats)
        nd = src_band.GetNoDataValue()
        if nd is not None:
            dst_band.SetNoDataValue(nd)

    @staticmethod
    def _exterior_ring(geom):
        """Anillo exterior como lista de (x, y). Segmenta curvas y acepta multiparte de 1 parte."""
        if QgsWkbTypes.isCurvedType(geom.wkbType()):
            geom = QgsGeometry(geom.constGet().segmentize())

        if geom.isMultipart():
            parts = geom.asMultiPolygon()
            if len(parts) != 1:
                raise QgsProcessingException(
                    "El rectángulo debe ser un polígono simple; la geometría tiene {} partes.".format(len(parts))
                )
            poly = parts[0]
        else:
            poly = geom.asPolygon()

        if not poly or not poly[0]:
            raise QgsProcessingException("No fue posible obtener el anillo exterior del polígono.")
        if len(poly) != 1:
            raise QgsProcessingException("El rectángulo no debe contener huecos interiores.")
        return [(p.x(), p.y()) for p in poly[0]]

    @staticmethod
    def _get_visible_bands(ds, feedback):
        """Devuelve (índices de bandas visibles, índice alfa o None, 'RGB'|'GRAY')."""
        red = green = blue = alpha = None
        for i in range(1, ds.RasterCount + 1):
            ci = ds.GetRasterBand(i).GetColorInterpretation()
            if ci == gdal.GCI_RedBand and red is None:
                red = i
            elif ci == gdal.GCI_GreenBand and green is None:
                green = i
            elif ci == gdal.GCI_BlueBand and blue is None:
                blue = i
            elif ci == gdal.GCI_AlphaBand and alpha is None:
                alpha = i

        if red and green and blue:
            return [red, green, blue], alpha, "RGB"

        data_bands = [i for i in range(1, ds.RasterCount + 1) if i != alpha]
        if len(data_bands) >= 3:
            feedback.pushInfo(
                "El raster no declara RGB; el PNG usará las tres primeras bandas de datos. "
                "El GeoTIFF conserva todas."
            )
            return data_bands[:3], alpha, "RGB"
        if len(data_bands) == 2:
            feedback.pushInfo("Dos bandas de datos no RGB: el PNG usará la primera como gris.")
        if len(data_bands) >= 1:
            return [data_bands[0]], alpha, "GRAY"
        raise QgsProcessingException("No se encontraron bandas de datos utilizables.")

    @staticmethod
    def _block_shape(width, height):
        block_w = min(width, 2048)
        block_h = min(height, max(1, MASK_BLOCK_PIXELS // block_w))
        return block_w, block_h

    @staticmethod
    def _band_range(ds, valid_ds, band_idx, method):
        """Extremos exactos (min–max) o percentiles muestreados sobre la máscara común."""
        band = ds.GetRasterBand(band_idx)
        valid_band = valid_ds.GetRasterBand(1)
        w, h = ds.RasterXSize, ds.RasterYSize

        if method == AlignRasterToRotatedRectangle.PNG_VISUAL_MINMAX:
            lo = hi = None
            bw, bh = AlignRasterToRotatedRectangle._block_shape(w, h)
            for y in range(0, h, bh):
                rows = min(bh, h - y)
                for x in range(0, w, bw):
                    cols = min(bw, w - x)
                    arr = np.asarray(band.ReadAsArray(x, y, cols, rows), dtype=np.float64)
                    mask = np.asarray(valid_band.ReadAsArray(x, y, cols, rows)) > 0
                    vals = arr[np.isfinite(arr) & mask]
                    if vals.size:
                        mn, mx = float(vals.min()), float(vals.max())
                        lo = mn if lo is None else min(lo, mn)
                        hi = mx if hi is None else max(hi, mx)
            return None if lo is None else (lo, hi)

        factor = max(1, int(ceil(sqrt((w * h) / float(MAX_SAMPLE_PIXELS)))))
        bx = max(1, w // factor)
        by = max(1, h // factor)
        arr = np.asarray(band.ReadAsArray(0, 0, w, h, buf_xsize=bx, buf_ysize=by), dtype=np.float64)
        mask = np.asarray(valid_band.ReadAsArray(0, 0, w, h, buf_xsize=bx, buf_ysize=by)) > 0
        vals = arr[np.isfinite(arr) & mask]
        if vals.size == 0:
            return None
        lo, hi = np.percentile(vals, (2.0, 98.0))
        return float(lo), float(hi)

    @staticmethod
    def _alpha_max(band, manual_value):
        """Máximo nominal de alfa, priorizando valor explícito y NBITS."""
        if manual_value > 0.0:
            return float(manual_value), "valor indicado por el usuario"
        nbits_text = band.GetMetadataItem("NBITS", "IMAGE_STRUCTURE")
        try:
            nbits = int(nbits_text) if nbits_text else None
        except (TypeError, ValueError):
            nbits = None
        if nbits is not None and 1 <= nbits <= 63:
            return float((1 << nbits) - 1), "NBITS={}".format(nbits)
        dtype = band.DataType
        signed_max = {
            gdal.GDT_Byte: 255.0,
            gdal.GDT_UInt16: 65535.0,
            gdal.GDT_Int16: 32767.0,
            gdal.GDT_UInt32: 4294967295.0,
            gdal.GDT_Int32: 2147483647.0,
        }
        if hasattr(gdal, "GDT_UInt64"):
            signed_max[gdal.GDT_UInt64] = float((1 << 64) - 1)
        if hasattr(gdal, "GDT_Int64"):
            signed_max[gdal.GDT_Int64] = float((1 << 63) - 1)
        if dtype in signed_max:
            return signed_max[dtype], "tipo {}".format(gdal.GetDataTypeName(dtype))
        if dtype in (gdal.GDT_Float32, gdal.GDT_Float64):
            return 1.0, "alfa flotante normalizado (0–1)"
        return 255.0, "suposición conservadora para {}".format(gdal.GetDataTypeName(dtype))

    @staticmethod
    def _write_png_alpha(dst_ds, valid_ds, alpha_idx, source_alpha_max, png_alpha_ds, png_alpha_max):
        """Combina alfa (si existe) con la máscara común y lo normaliza para PNG."""
        w, h = dst_ds.RasterXSize, dst_ds.RasterYSize
        valid_band = valid_ds.GetRasterBand(1)
        png_band = png_alpha_ds.GetRasterBand(1)
        alpha_band = dst_ds.GetRasterBand(alpha_idx) if alpha_idx is not None else None
        out_dtype = np.uint16 if png_alpha_max > 255.0 else np.uint8
        bw, bh = AlignRasterToRotatedRectangle._block_shape(w, h)
        for y in range(0, h, bh):
            rows = min(bh, h - y)
            for x in range(0, w, bw):
                cols = min(bw, w - x)
                valid = np.asarray(valid_band.ReadAsArray(x, y, cols, rows)) > 0
                if alpha_band is None:
                    out = np.where(valid, png_alpha_max, 0.0)
                else:
                    raw = np.asarray(alpha_band.ReadAsArray(x, y, cols, rows), dtype=np.float64)
                    out = np.clip(raw * (png_alpha_max / source_alpha_max), 0.0, png_alpha_max)
                    out[~np.isfinite(out) | ~valid] = 0.0
                png_band.WriteArray(np.rint(out).astype(out_dtype), x, y)

    @staticmethod
    def _build_masked_source_vrt_xml(src_ds, src_path, valid_mask_path):
        """VRT de todas las bandas de origen con MaskBand común para GDAL Warp."""
        w, h = src_ds.RasterXSize, src_ds.RasterYSize
        gt_text = ", ".join("{:.17g}".format(v) for v in src_ds.GetGeoTransform())
        src_path_xml = xml_escape(AlignRasterToRotatedRectangle._vrt_source_path(src_path))
        mask_path_xml = xml_escape(AlignRasterToRotatedRectangle._vrt_source_path(valid_mask_path))
        src_rect = '      <SrcRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(w, h)
        dst_rect = '      <DstRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(w, h)
        parts = [
            '<VRTDataset rasterXSize="{}" rasterYSize="{}">'.format(w, h),
            "  <SRS>{}</SRS>".format(xml_escape(src_ds.GetProjection() or "")),
            "  <GeoTransform>{}</GeoTransform>".format(gt_text),
        ]
        for i in range(1, src_ds.RasterCount + 1):
            band = src_ds.GetRasterBand(i)
            nodata = band.GetNoDataValue()
            ci_name = gdal.GetColorInterpretationName(band.GetColorInterpretation())
            parts += [
                '  <VRTRasterBand dataType="{}" band="{}">'.format(gdal.GetDataTypeName(band.DataType), i),
                "    <ColorInterp>{}</ColorInterp>".format(xml_escape(ci_name)),
                "    <SimpleSource>",
                '      <SourceFilename relativeToVRT="0">{}</SourceFilename>'.format(src_path_xml),
                "      <SourceBand>{}</SourceBand>".format(i), src_rect, dst_rect,
                "    </SimpleSource>",
            ]
            if nodata is not None:
                parts.append("    <NoDataValue>{:.17g}</NoDataValue>".format(float(nodata)))
            parts.append("  </VRTRasterBand>")
        parts += [
            "  <MaskBand>",
            '    <VRTRasterBand dataType="Byte">',
            "      <SimpleSource>",
            '        <SourceFilename relativeToVRT="0">{}</SourceFilename>'.format(mask_path_xml),
            "        <SourceBand>1</SourceBand>",
            '        <SrcRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(w, h),
            '        <DstRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(w, h),
            "      </SimpleSource>",
            "    </VRTRasterBand>",
            "  </MaskBand>",
            "</VRTDataset>",
        ]
        return "\n".join(parts)

    @staticmethod
    def _build_png_vrt_xml(width, height, projection, geotransform, dtype_name, color_specs, alpha_spec):
        """
        color_specs: lista de dict(path, band, interp, nodata, scale=(offset, ratio)|None)
        alpha_spec : dict(path, band (int o 'mask,1'), ratio)
        """
        gt_text = ", ".join("{:.17g}".format(v) for v in geotransform)
        src_rect = '      <SrcRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(width, height)
        dst_rect = '      <DstRect xOff="0" yOff="0" xSize="{}" ySize="{}"/>'.format(width, height)

        parts = [
            '<VRTDataset rasterXSize="{}" rasterYSize="{}">'.format(width, height),
            "  <SRS>{}</SRS>".format(xml_escape(projection or "")),
            "  <GeoTransform>{}</GeoTransform>".format(gt_text),
        ]

        band_no = 1
        for spec in color_specs:
            path_xml = xml_escape(AlignRasterToRotatedRectangle._vrt_source_path(spec["path"]))
            parts.append(
                '  <VRTRasterBand dataType="{}" band="{}" colorInterp="{}">'.format(
                    dtype_name, band_no, spec["interp"]
                )
            )
            if spec.get("scale") is None:
                parts += [
                    "    <SimpleSource>",
                    '      <SourceFilename relativeToVRT="0">{}</SourceFilename>'.format(path_xml),
                    "      <SourceBand>{}</SourceBand>".format(spec["band"]),
                    src_rect,
                    dst_rect,
                    "    </SimpleSource>",
                ]
            else:
                offset, ratio = spec["scale"]
                parts += [
                    "    <ComplexSource>",
                    '      <SourceFilename relativeToVRT="0">{}</SourceFilename>'.format(path_xml),
                    "      <SourceBand>{}</SourceBand>".format(spec["band"]),
                    src_rect,
                    dst_rect,
                    "      <ScaleOffset>{:.17g}</ScaleOffset>".format(offset),
                    "      <ScaleRatio>{:.17g}</ScaleRatio>".format(ratio),
                ]
                if spec.get("nodata") is not None:
                    parts.append("      <NODATA>{:.17g}</NODATA>".format(float(spec["nodata"])))
                parts.append("    </ComplexSource>")
            parts.append("  </VRTRasterBand>")
            band_no += 1

        alpha_path_xml = xml_escape(AlignRasterToRotatedRectangle._vrt_source_path(alpha_spec["path"]))
        parts += [
            '  <VRTRasterBand dataType="{}" band="{}" colorInterp="Alpha">'.format(dtype_name, band_no),
            "    <ComplexSource>",
            '      <SourceFilename relativeToVRT="0">{}</SourceFilename>'.format(alpha_path_xml),
            "      <SourceBand>{}</SourceBand>".format(alpha_spec["band"]),
            src_rect,
            dst_rect,
            "      <ScaleOffset>0</ScaleOffset>",
            "      <ScaleRatio>{:.17g}</ScaleRatio>".format(alpha_spec["ratio"]),
            "    </ComplexSource>",
            "  </VRTRasterBand>",
            "</VRTDataset>",
        ]
        return "\n".join(parts)

    # ---------------- Ejecución ----------------

    def processAlgorithm(self, parameters, context, feedback):
        prev_exceptions = gdal.GetUseExceptions()
        gdal.UseExceptions()

        ds_refs = {}      # datasets abiertos, para cerrarlos siempre
        temp_files = []   # archivos temporales en disco
        vsimem = []       # rutas /vsimem/

        try:
            return self._run(parameters, context, feedback, ds_refs, temp_files, vsimem)
        finally:
            for k in list(ds_refs.keys()):
                ds_refs[k] = None
            ds_refs.clear()
            for p in vsimem:
                self._unlink_vsimem(p)
            for p in temp_files:
                self._remove_existing(p)
            if not prev_exceptions:
                gdal.DontUseExceptions()

    def _run(self, parameters, context, feedback, ds_refs, temp_files, vsimem):
        # ---------- Parámetros ----------
        raster_layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        rect_source = self.parameterAsSource(parameters, self.RECTANGLE, context)
        edge_mode = self.parameterAsEnum(parameters, self.EDGE_MODE, context)
        fallback_ombb = self.parameterAsBool(parameters, self.FALLBACK_OMBB, context)
        angle_tol = self.parameterAsDouble(parameters, self.ANGLE_TOL, context)
        length_tol_pct = self.parameterAsDouble(parameters, self.LENGTH_TOL, context)
        requested_pixel_size = self.parameterAsDouble(parameters, self.PIXEL_SIZE, context)
        resampling_idx = self.parameterAsEnum(parameters, self.RESAMPLING, context)
        strict_source_cutline = self.parameterAsBool(parameters, self.STRICT_SOURCE_CUTLINE, context)
        warp_memory_mb = self.parameterAsInt(parameters, self.WARP_MEMORY_MB, context)
        manual_alpha_max = self.parameterAsDouble(parameters, self.ALPHA_MAX, context)
        png_mode = self.parameterAsEnum(parameters, self.PNG_MODE, context)
        out_tif = self.parameterAsOutputLayer(parameters, self.OUTPUT_TIF, context)
        out_png = self.parameterAsFileOutput(parameters, self.OUTPUT_PNG, context)

        if raster_layer is None:
            raise QgsProcessingException("Raster de entrada inválido.")
        if rect_source is None:
            raise QgsProcessingException("Capa de rectángulo inválida.")

        if not out_tif.lower().endswith((".tif", ".tiff")):
            raise QgsProcessingException(
                "La salida GeoTIFF debe terminar en .tif o .tiff (ruta recibida: {}).".format(out_tif)
            )
        if not out_png.lower().endswith(".png"):
            raise QgsProcessingException("La salida PNG debe terminar en .png.")

        provider = raster_layer.dataProvider()
        if provider is None or provider.name() != "gdal":
            raise QgsProcessingException(
                "El raster debe ser una capa GDAL basada en archivo (proveedor actual: {}).".format(
                    provider.name() if provider else "desconocido"
                )
            )
        raster_path = provider.dataSourceUri().split("|")[0]

        # ---------- Rectángulo: exactamente 1 entidad ----------
        feats = []
        for f in rect_source.getFeatures():
            feats.append(f)
            if len(feats) > 1:
                break
        if len(feats) != 1:
            raise QgsProcessingException(
                "La entrada RECTANGLE debe contener exactamente 1 polígono ({} encontrados). "
                "Si la capa contiene varios, active «solo entidades seleccionadas».".format(
                    "0" if not feats else "≥2"
                )
            )
        geom = feats[0].geometry()
        if geom is None or geom.isEmpty():
            raise QgsProcessingException("El polígono está vacío.")
        if QgsWkbTypes.geometryType(geom.wkbType()) != QgsWkbTypes.PolygonGeometry:
            raise QgsProcessingException("La geometría debe ser poligonal.")
        if not geom.isGeosValid():
            raise QgsProcessingException("El polígono del rectángulo no es geométricamente válido.")
        if geom.area() <= 0.0:
            raise QgsProcessingException("El polígono del rectángulo tiene área nula.")

        # ---------- Raster origen ----------
        try:
            src_ds = gdal.Open(raster_path, gdal.GA_ReadOnly)
        except RuntimeError as exc:
            raise QgsProcessingException("GDAL no pudo abrir el raster {}: {}".format(raster_path, exc))
        if src_ds is None:
            raise QgsProcessingException("GDAL no pudo abrir el raster: {}".format(raster_path))
        ds_refs["src"] = src_ds

        n_bands = src_ds.RasterCount
        if n_bands < 1:
            raise QgsProcessingException("El raster no contiene bandas.")

        src_gt = src_ds.GetGeoTransform(can_return_null=True)
        if src_gt is None:
            raise QgsProcessingException("El raster no tiene GeoTransform válido.")

        raster_crs = raster_layer.crs()
        if not raster_crs.isValid():
            raise QgsProcessingException("El raster no tiene un CRS válido.")
        src_wkt = src_ds.GetProjection()
        if not src_wkt:
            # El archivo no lleva CRS embebido (p. ej. solo .prj/.aux.xml): usar el de la capa.
            src_wkt = raster_crs.toWkt()
            feedback.pushInfo("El raster no declara CRS internamente; se usa el CRS de la capa QGIS.")

        if raster_crs.isGeographic():
            feedback.pushWarning(
                "El CRS del raster es geográfico: ángulos y tamaños se calculan en el plano "
                "de sus coordenadas. Para fidelidad métrica conviene un CRS proyectado."
            )

        src_types = [src_ds.GetRasterBand(i).DataType for i in range(1, n_bands + 1)]
        if len(set(src_types)) != 1:
            raise QgsProcessingException(
                "Las bandas del TIFF tienen tipos de dato diferentes; se requiere un tipo homogéneo."
            )
        src_dtype = src_types[0]

        if png_mode == self.PNG_VISUAL_AUTO and src_dtype in (gdal.GDT_Byte, gdal.GDT_UInt16):
            png_mode = self.PNG_EXACT
            feedback.pushInfo("PNG automático: Byte/UInt16, se conservan los valores de los canales.")

        if png_mode == self.PNG_EXACT and src_dtype not in (gdal.GDT_Byte, gdal.GDT_UInt16):
            feedback.pushWarning(
                "PNG exacto solo admite Byte o UInt16; el raster es {}. Se usará el modo visual "
                "automático para el PNG. El GeoTIFF conserva el tipo original.".format(
                    gdal.GetDataTypeName(src_dtype)
                )
            )
            png_mode = self.PNG_VISUAL_AUTO

        # ---------- Geometría en el CRS del raster ----------
        rect_crs = rect_source.sourceCrs()
        if not rect_crs.isValid():
            raise QgsProcessingException("La capa del rectángulo no tiene CRS válido.")
        if rect_crs != raster_crs:
            geom = QgsGeometry(geom)
            ct = QgsCoordinateTransform(rect_crs, raster_crs, QgsProject.instance())
            geom.transform(ct)
            feedback.pushInfo("Rectángulo reproyectado de {} a {}.".format(rect_crs.authid(), raster_crs.authid()))

        ring = self._exterior_ring(geom)
        pts = _clean_ring(ring, angle_tol)

        used_ombb = False
        try:
            vecs, lengths, angles = _validate_rectangle(pts, angle_tol, length_tol_pct)
        except QgsProcessingException as exc:
            if not fallback_ombb:
                raise
            feedback.pushWarning(
                "El polígono no supera la validación ({}). Se usará su rectángulo "
                "envolvente mínimo orientado.".format(exc)
            )
            ombb = geom.orientedMinimumBoundingBox()
            ombb_geom = ombb[0] if isinstance(ombb, (tuple, list)) else ombb
            if ombb_geom is None or ombb_geom.isEmpty():
                raise QgsProcessingException("No fue posible calcular el rectángulo envolvente mínimo orientado.")
            pts = _clean_ring(self._exterior_ring(ombb_geom), angle_tol)
            vecs, lengths, angles = _validate_rectangle(pts, angle_tol, length_tol_pct)
            used_ombb = True

        frame = _build_frame(pts, vecs, lengths, edge_mode)
        ux, uy = frame["u"]
        vx, vy = frame["v"]
        origin_x, origin_y = frame["origin"]
        rect_width, rect_height = frame["width"], frame["height"]
        theta = frame["theta"]

        # ---------- Rejilla de salida ----------
        if requested_pixel_size > 0.0:
            target_du = target_dv = requested_pixel_size
        else:
            target_du = _source_pixel_step_in_direction(src_gt, ux, uy)
            target_dv = _source_pixel_step_in_direction(src_gt, vx, vy)

        width_px = max(1, int(ceil(rect_width / target_du)))
        height_px = max(1, int(ceil(rect_height / target_dv)))
        du = rect_width / width_px
        dv = rect_height / height_px

        dst_gt = (origin_x, ux * du, vx * dv, origin_y, uy * du, vy * dv)

        # Aviso si el rectángulo sobresale del raster
        rect_geom = QgsGeometry.fromPolygonXY([[QgsPointXY(x, y) for x, y in pts] + [QgsPointXY(*pts[0])]])
        footprint = _geotransform_polygon(src_gt, src_ds.RasterXSize, src_ds.RasterYSize)
        if not footprint.buffer(0.5 * max(du, dv), 1).contains(rect_geom):
            feedback.pushWarning(
                "El rectángulo sobresale de la huella del raster: la parte exterior quedará "
                "transparente/nodata en las salidas."
            )

        feedback.pushInfo("Rectángulo validado{}.".format(" (OMBB)" if used_ombb else ""))
        feedback.pushInfo("Esquinas: " + "; ".join("P{}=({:.4f}, {:.4f})".format(i, x, y) for i, (x, y) in enumerate(pts)))
        feedback.pushInfo("Ángulos: " + ", ".join("{:.6f}°".format(a) for a in angles))
        feedback.pushInfo("Borde elegido: {}".format(frame["edge_idx"]))
        feedback.pushInfo("θ = {:.8f}°".format(theta))
        feedback.pushInfo("Dimensiones: {:.6f} × {:.6f} unidades CRS".format(rect_width, rect_height))
        feedback.pushInfo(
            "Raster destino: {} × {} px; resolución efectiva {:.9g} × {:.9g}".format(width_px, height_px, du, dv)
        )

        # ---------- Estrategia de validez ----------
        source_alpha_idx = next(
            (i for i in range(1, n_bands + 1)
             if src_ds.GetRasterBand(i).GetColorInterpretation() == gdal.GCI_AlphaBand),
            None,
        )
        alpha_exists = source_alpha_idx is not None
        data_band_indices = [i for i in range(1, n_bands + 1) if i != source_alpha_idx]
        if not data_band_indices:
            raise QgsProcessingException("El raster no contiene bandas de datos además del alfa.")
        feedback.pushInfo(
            "Validez común: todas las {} bandas de datos{}{}.".format(
                len(data_band_indices),
                " y alfa" if alpha_exists else "",
                "; el rectángulo limita también la fuente" if strict_source_cutline else "",
            )
        )

        # ---------- GeoTIFF de salida ----------
        self._remove_existing(out_tif)
        self._remove_existing(out_png)
        os.makedirs(os.path.dirname(os.path.abspath(out_tif)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)

        gtiff = gdal.GetDriverByName("GTiff")
        if gtiff is None:
            raise QgsProcessingException("El driver GDAL GeoTIFF no está disponible.")

        predictor = "3" if src_dtype in (gdal.GDT_Float32, gdal.GDT_Float64) else "2"
        try:
            dst_ds = gtiff.Create(
                out_tif, width_px, height_px, n_bands, src_dtype,
                options=[
                    "TILED=YES",
                    "COMPRESS=DEFLATE",
                    "PREDICTOR={}".format(predictor),
                    "BIGTIFF=IF_SAFER",
                    "NUM_THREADS=ALL_CPUS",
                ],
            )
        except RuntimeError as exc:
            raise QgsProcessingException("No fue posible crear el GeoTIFF de salida: {}".format(exc))
        if dst_ds is None:
            raise QgsProcessingException("No fue posible crear el GeoTIFF de salida.")
        ds_refs["dst"] = dst_ds

        dst_ds.SetGeoTransform(dst_gt)
        dst_ds.SetProjection(src_wkt)
        ds_md = src_ds.GetMetadata()
        if ds_md:
            dst_ds.SetMetadata(ds_md)
        for i in range(1, n_bands + 1):
            self._copy_band_metadata(src_ds.GetRasterBand(i), dst_ds.GetRasterBand(i))

        # La máscara interna debe crearse ANTES de escribir píxeles. Si existe
        # una banda alfa se conserva esa semántica y el PNG emplea su máscara común.
        internal_mask_ok = False
        if not alpha_exists:
            prev_cfg = gdal.GetConfigOption("GDAL_TIFF_INTERNAL_MASK")
            gdal.SetConfigOption("GDAL_TIFF_INTERNAL_MASK", "YES")
            try:
                err = dst_ds.CreateMaskBand(gdal.GMF_PER_DATASET)
                internal_mask_ok = (err == gdal.CE_None)
            except RuntimeError as exc:
                feedback.pushWarning("No se pudo crear la máscara interna del GeoTIFF: {}".format(exc))
            finally:
                gdal.SetConfigOption("GDAL_TIFF_INTERNAL_MASK", prev_cfg)
            if not internal_mask_ok:
                feedback.pushWarning(
                    "El GeoTIFF no tendrá máscara de validez; el PNG sí tendrá alfa a partir de la máscara temporal."
                )

        # ---------- Máscara común de la fuente ----------
        # Se crea por bloques, evitando matrices del tamaño total del raster.
        src_valid_path = QgsProcessingUtils.generateTempFilename("rot_rect_source_valid.tif")
        self._remove_existing(src_valid_path)
        temp_files.append(src_valid_path)
        src_valid_ds = gtiff.Create(
            src_valid_path, src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GDT_Byte,
            options=["TILED=YES", "COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS"],
        )
        if src_valid_ds is None:
            raise QgsProcessingException("No fue posible crear la máscara temporal de validez.")
        ds_refs["src_valid"] = src_valid_ds
        src_valid_ds.SetGeoTransform(src_gt)
        src_valid_ds.SetProjection(src_wkt)
        src_valid_band = src_valid_ds.GetRasterBand(1)
        src_valid_band.Fill(0 if strict_source_cutline else 255)

        if strict_source_cutline:
            ogr_rect = ogr.CreateGeometryFromWkt(rect_geom.asWkt())
            if ogr_rect is None:
                raise QgsProcessingException("No fue posible convertir el rectángulo para la máscara de recorte.")
            vector_driver = ogr.GetDriverByName("Memory")
            if vector_driver is None:
                raise QgsProcessingException("El driver vectorial OGR Memory no está disponible.")
            vector_ds = vector_driver.CreateDataSource("")
            if vector_ds is None:
                raise QgsProcessingException("No fue posible crear la capa temporal de recorte.")
            vector_layer = vector_ds.CreateLayer("rectangle", geom_type=ogr.wkbPolygon)
            if vector_layer is None:
                raise QgsProcessingException("No fue posible crear la capa temporal del rectángulo.")
            vector_feature = ogr.Feature(vector_layer.GetLayerDefn())
            vector_feature.SetGeometry(ogr_rect)
            if vector_layer.CreateFeature(vector_feature) != ogr.OGRERR_NONE:
                raise QgsProcessingException("No fue posible insertar el rectángulo en la capa temporal.")
            vector_feature = None
            err = gdal.RasterizeLayer(src_valid_ds, [1], vector_layer, burn_values=[255])
            vector_layer = None
            vector_ds = None
            if err != gdal.CE_None:
                raise QgsProcessingException("No fue posible rasterizar el rectángulo sobre la máscara de fuente.")

        bw, bh = self._block_shape(src_ds.RasterXSize, src_ds.RasterYSize)
        feedback.pushInfo("Construyendo máscara común de bandas por bloques...")
        for y in range(0, src_ds.RasterYSize, bh):
            rows = min(bh, src_ds.RasterYSize - y)
            for x in range(0, src_ds.RasterXSize, bw):
                cols = min(bw, src_ds.RasterXSize - x)
                valid = np.asarray(src_valid_band.ReadAsArray(x, y, cols, rows)) > 0
                for band_idx in data_band_indices:
                    band_mask = np.asarray(
                        src_ds.GetRasterBand(band_idx).GetMaskBand().ReadAsArray(x, y, cols, rows)
                    ) > 0
                    valid &= band_mask
                if source_alpha_idx is not None:
                    alpha_values = np.asarray(
                        src_ds.GetRasterBand(source_alpha_idx).ReadAsArray(x, y, cols, rows)
                    )
                    valid &= np.isfinite(alpha_values) & (alpha_values > 0)
                src_valid_band.WriteArray(np.where(valid, 255, 0).astype(np.uint8), x, y)
            feedback.setProgress(2.0 + 3.0 * min(1.0, float(y + rows) / src_ds.RasterYSize))
            if feedback.isCanceled():
                raise QgsProcessingException("Proceso cancelado.")
        src_valid_ds.FlushCache()

        source_vrt_path = "/vsimem/rot_rect_masked_source.vrt"
        self._unlink_vsimem(source_vrt_path)
        vsimem.append(source_vrt_path)
        gdal.FileFromMemBuffer(
            source_vrt_path,
            self._build_masked_source_vrt_xml(src_ds, raster_path, src_valid_path).encode("utf-8"),
        )
        source_vrt_ds = gdal.Open(source_vrt_path, gdal.GA_ReadOnly)
        if source_vrt_ds is None:
            raise QgsProcessingException("No fue posible crear el VRT fuente con máscara común.")
        ds_refs["source_vrt"] = source_vrt_ds

        # ---------- Máscara común en la rejilla de salida ----------
        mask_path = QgsProcessingUtils.generateTempFilename("rot_rect_valid.tif")
        self._remove_existing(mask_path)
        temp_files.append(mask_path)
        mask_ds = gtiff.Create(
            mask_path, width_px, height_px, 1, gdal.GDT_Byte,
            options=["TILED=YES", "COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS"],
        )
        if mask_ds is None:
            raise QgsProcessingException("No fue posible crear la máscara temporal alineada.")
        ds_refs["mask"] = mask_ds
        mask_ds.SetGeoTransform(dst_gt)
        mask_ds.SetProjection(src_wkt)
        mask_ds.GetRasterBand(1).Fill(0)

        def mask_progress(complete, message, data):
            feedback.setProgress(5.0 + complete * 12.0)
            return 0 if feedback.isCanceled() else 1

        # La cobertura del color bilinear depende de vecinos; una máscara
        # nearest puede declarar opacos píxeles con canales incompletos.
        mask_resampling = "near" if self.RESAMPLING_NAMES[resampling_idx] == "near" else "bilinear"
        mask_opts = gdal.WarpOptions(
            srcSRS=src_wkt, dstSRS=src_wkt, resampleAlg=mask_resampling, multithread=True,
            errorThreshold=0.0, warpMemoryLimit=warp_memory_mb * 1024 * 1024,
            warpOptions=["INIT_DEST=0", "NUM_THREADS=ALL_CPUS"], callback=mask_progress,
        )
        try:
            mres = gdal.Warp(mask_ds, src_valid_ds, options=mask_opts)
        except RuntimeError as exc:
            raise QgsProcessingException("No fue posible alinear la máscara común: {}".format(exc))
        if mres is None or feedback.isCanceled():
            raise QgsProcessingException("Proceso cancelado." if feedback.isCanceled() else "Falló el Warp de la máscara común.")
        mres = None
        mask_ds.FlushCache()

        # Solo se muestra el borde cuando el soporte de la interpolación
        # tiene cobertura completa. Se valida además la finitud/NoData de
        # todas las bandas ya remuestreadas antes de crear el alfa del PNG.
        out_bw, out_bh = self._block_shape(width_px, height_px)
        out_mask_band = mask_ds.GetRasterBand(1)

        # ---------- Warp de los datos ----------
        def warp_progress(complete, message, data):
            feedback.setProgress(17.0 + complete * 55.0)
            return 0 if feedback.isCanceled() else 1

        feedback.pushInfo("Remuestreando sobre la rejilla rotada ({})...".format(self.RESAMPLING_NAMES[resampling_idx]))
        warp_options = gdal.WarpOptions(
            srcSRS=src_wkt,
            dstSRS=src_wkt,
            resampleAlg=self.RESAMPLING_NAMES[resampling_idx],
            multithread=True,
            errorThreshold=0.0,
            warpMemoryLimit=warp_memory_mb * 1024 * 1024,
            warpOptions=["INIT_DEST=NO_DATA", "NUM_THREADS=ALL_CPUS"],
            callback=warp_progress,
        )
        try:
            warped = gdal.Warp(dst_ds, source_vrt_ds, options=warp_options)
        except RuntimeError as exc:
            if feedback.isCanceled():
                raise QgsProcessingException("Proceso cancelado.")
            raise QgsProcessingException("GDAL Warp falló al generar el raster alineado: {}".format(exc))
        if warped is None or feedback.isCanceled():
            raise QgsProcessingException("Proceso cancelado." if feedback.isCanceled() else "GDAL Warp falló.")
        warped = None
        dst_ds.FlushCache()

        src_nodata_by_band = {
            idx: src_ds.GetRasterBand(idx).GetNoDataValue() for idx in data_band_indices
        }
        for y in range(0, height_px, out_bh):
            rows = min(out_bh, height_px - y)
            for x in range(0, width_px, out_bw):
                cols = min(out_bw, width_px - x)
                mask_array = np.asarray(out_mask_band.ReadAsArray(x, y, cols, rows))
                valid = mask_array >= (254 if mask_resampling == "bilinear" else 1)
                for idx in data_band_indices:
                    values = np.asarray(dst_ds.GetRasterBand(idx).ReadAsArray(x, y, cols, rows))
                    valid &= np.isfinite(values)
                    nd = src_nodata_by_band[idx]
                    if nd is not None and np.isfinite(nd):
                        valid &= values != nd
                out_mask_band.WriteArray(np.where(valid, 255, 0).astype(np.uint8), x, y)
            if feedback.isCanceled():
                raise QgsProcessingException("Proceso cancelado.")
        mask_ds.FlushCache()
        feedback.pushInfo("Borde: máscara coherente con el remuestreo {} y todas las bandas.".format(
            self.RESAMPLING_NAMES[resampling_idx]
        ))

        if internal_mask_ok:
            dst_mask_band = dst_ds.GetRasterBand(1).GetMaskBand()
            src_mask_band = mask_ds.GetRasterBand(1)
            for y in range(0, height_px, out_bh):
                rows = min(out_bh, height_px - y)
                for x in range(0, width_px, out_bw):
                    cols = min(out_bw, width_px - x)
                    buf = src_mask_band.ReadRaster(x, y, cols, rows, buf_type=gdal.GDT_Byte)
                    dst_mask_band.WriteRaster(x, y, cols, rows, buf, buf_type=gdal.GDT_Byte)
            dst_ds.FlushCache()

        feedback.setProgress(74.0)

        # ---------- Preparación del PNG ----------
        visible_bands, alpha_idx, png_layout = self._get_visible_bands(dst_ds, feedback)
        interps = ["Red", "Green", "Blue"] if png_layout == "RGB" else ["Gray"]
        feedback.pushInfo("Canales del PNG ({}): bandas {} del TIFF alineado.".format(
            png_layout, ", ".join(str(i) for i in visible_bands)
        ))

        is_palette = (
            png_layout == "GRAY"
            and dst_ds.GetRasterBand(visible_bands[0]).GetColorTable() is not None
        )

        ranges = {}
        if png_mode != self.PNG_EXACT and not is_palette:
            sampled_ranges = {}
            if png_mode in (self.PNG_VISUAL_AUTO, self.PNG_VISUAL_SHARED):
                for idx in visible_bands:
                    sampled_ranges[idx] = self._band_range(dst_ds, mask_ds, idx, self.PNG_VISUAL_PCT)

            if png_mode == self.PNG_VISUAL_AUTO:
                nonempty = [r for r in sampled_ranges.values() if r is not None]
                looks_normalized = (
                    len(nonempty) == len(visible_bands)
                    and all(lo >= -0.02 and hi <= 1.05 for lo, hi in nonempty)
                )
                if looks_normalized:
                    png_mode = self.PNG_VISUAL_UNIT
                    feedback.pushInfo("PNG automático: canales Float normalizados; se aplica 0–1 → 0–255.")
                else:
                    png_mode = self.PNG_VISUAL_SHARED if png_layout == "RGB" else self.PNG_VISUAL_PCT
                    feedback.pushInfo("PNG automático: se usa escala {}.".format(
                        "común para RGB" if png_layout == "RGB" else "2–98 % para gris"
                    ))

            if png_mode == self.PNG_VISUAL_UNIT:
                ranges = {idx: (0.0, 1.0) for idx in visible_bands}
                feedback.pushInfo("PNG: mapeo idéntico en todos los canales, 0–1 → 0–255.")
            elif png_mode == self.PNG_VISUAL_SHARED and png_layout == "RGB":
                if not sampled_ranges:
                    for idx in visible_bands:
                        sampled_ranges[idx] = self._band_range(dst_ds, mask_ds, idx, self.PNG_VISUAL_PCT)
                nonempty = [r for r in sampled_ranges.values() if r is not None]
                common = (min(r[0] for r in nonempty), max(r[1] for r in nonempty)) if nonempty else (0.0, 1.0)
                ranges = {idx: common for idx in visible_bands}
                feedback.pushInfo("PNG: escala RGB común {:.6g} – {:.6g}; se conserva el balance entre canales.".format(*common))
            else:
                if png_mode == self.PNG_VISUAL_SHARED:
                    png_mode = self.PNG_VISUAL_PCT
                    feedback.pushInfo("Raster gris: se aplica el corte 2–98 %.")
                range_label = "min–max exacto" if png_mode == self.PNG_VISUAL_MINMAX else "percentiles 2–98 %"
                feedback.pushInfo("Calculando rangos visuales ({})...".format(range_label))
                for idx in visible_bands:
                    rng = self._band_range(dst_ds, mask_ds, idx, png_mode)
                    if rng is None:
                        feedback.pushWarning("La banda {} no tiene píxeles válidos; se escalará 0–1.".format(idx))
                        rng = (0.0, 1.0)
                    ranges[idx] = rng
                    feedback.pushInfo("  Banda {}: {:.6g} – {:.6g}".format(idx, rng[0], rng[1]))

        nodata_by_band = {i: dst_ds.GetRasterBand(i).GetNoDataValue() for i in visible_bands}

        # ---------- Alfa PNG: alfa original × máscara común ----------
        if alpha_exists and alpha_idx is None:
            raise QgsProcessingException("La banda alfa del origen no se conservó en el GeoTIFF alineado.")
        if is_palette or png_mode != self.PNG_EXACT:
            png_dtype = gdal.GDT_Byte
            dst_alpha_max = 255.0
        else:
            png_dtype = gdal.GDT_UInt16 if src_dtype == gdal.GDT_UInt16 else gdal.GDT_Byte
            dst_alpha_max = 65535.0 if png_dtype == gdal.GDT_UInt16 else 255.0

        src_alpha_max = 1.0
        if alpha_exists:
            src_alpha_max, alpha_basis = self._alpha_max(src_ds.GetRasterBand(source_alpha_idx), manual_alpha_max)
            feedback.pushInfo("Normalización del alfa: máximo {:.9g} ({})".format(src_alpha_max, alpha_basis))

        png_alpha_path = QgsProcessingUtils.generateTempFilename("rot_rect_png_alpha.tif")
        self._remove_existing(png_alpha_path)
        temp_files.append(png_alpha_path)
        png_alpha_ds = gtiff.Create(
            png_alpha_path, width_px, height_px, 1, png_dtype,
            options=["TILED=YES", "COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS"],
        )
        if png_alpha_ds is None:
            raise QgsProcessingException("No fue posible crear el alfa temporal para el PNG.")
        ds_refs["png_alpha"] = png_alpha_ds
        png_alpha_ds.SetGeoTransform(dst_gt)
        png_alpha_ds.SetProjection(src_wkt)
        self._write_png_alpha(dst_ds, mask_ds, alpha_idx if alpha_exists else None,
                              src_alpha_max, png_alpha_ds, dst_alpha_max)
        png_alpha_ds.FlushCache()

        # Cerrar las fuentes antes de referenciarlas desde el VRT PNG.
        dst_ds.FlushCache()
        ds_refs["dst"] = None
        dst_ds = None
        ds_refs["mask"] = None
        mask_ds = None
        ds_refs["png_alpha"] = None
        png_alpha_ds = None

        feedback.setProgress(80.0)

        # Fuente de color: el GeoTIFF, o un VRT expandido a RGB si es paleta
        color_source_path = out_tif
        if is_palette:
            expand_vrt = "/vsimem/rot_rect_expand.vrt"
            self._unlink_vsimem(expand_vrt)
            vsimem.append(expand_vrt)
            try:
                expand_ds = gdal.Translate(expand_vrt, out_tif, options="-of VRT -b {} -expand rgb".format(visible_bands[0]))
            except RuntimeError as exc:
                raise QgsProcessingException("No fue posible expandir la tabla de colores a RGB: {}".format(exc))
            ds_refs["expand"] = expand_ds
            expand_ds.FlushCache()
            ds_refs["expand"] = None
            expand_ds = None
            color_source_path = expand_vrt
            visible_bands = [1, 2, 3]
            interps = ["Red", "Green", "Blue"]
            feedback.pushInfo("Raster con tabla de colores: el PNG se expande a RGB.")

        if is_palette or png_mode != self.PNG_EXACT:
            dtype_name = "Byte"
        else:
            dtype_name = gdal.GetDataTypeName(src_dtype)

        color_specs = []
        for src_idx, interp in zip(visible_bands, interps):
            spec = {"path": color_source_path, "band": src_idx, "interp": interp, "nodata": None, "scale": None}
            if not is_palette and png_mode != self.PNG_EXACT:
                vmin, vmax = ranges[src_idx]
                if vmax > vmin:
                    ratio = 255.0 / (vmax - vmin)
                    offset = -vmin * ratio
                else:
                    ratio, offset = 0.0, 0.0
                spec["scale"] = (offset, ratio)
                spec["nodata"] = nodata_by_band.get(src_idx)
            color_specs.append(spec)

        alpha_spec = {"path": png_alpha_path, "band": 1, "ratio": 1.0}

        vrt_xml = self._build_png_vrt_xml(
            width_px, height_px, src_wkt, dst_gt, dtype_name, color_specs, alpha_spec
        )

        png_vrt_path = "/vsimem/rot_rect_png.vrt"
        self._unlink_vsimem(png_vrt_path)
        vsimem.append(png_vrt_path)
        gdal.FileFromMemBuffer(png_vrt_path, vrt_xml.encode("utf-8"))

        try:
            png_vrt_ds = gdal.Open(png_vrt_path, gdal.GA_ReadOnly)
        except RuntimeError as exc:
            raise QgsProcessingException("No se pudo construir el VRT para exportar el PNG: {}".format(exc))
        ds_refs["png_vrt"] = png_vrt_ds

        png_driver = gdal.GetDriverByName("PNG")
        if png_driver is None:
            raise QgsProcessingException("El driver PNG de GDAL no está disponible.")

        def png_progress(complete, message, data):
            feedback.setProgress(82.0 + complete * 18.0)
            return 0 if feedback.isCanceled() else 1

        feedback.pushInfo("Exportando PNG rectangular ({} + alfa, {})...".format(png_layout, dtype_name))
        try:
            png_ds = png_driver.CreateCopy(
                out_png, png_vrt_ds, strict=1,
                options=["WORLDFILE=YES", "ZLEVEL=6"],
                callback=png_progress,
            )
        except RuntimeError as exc:
            if feedback.isCanceled():
                raise QgsProcessingException("Proceso cancelado.")
            raise QgsProcessingException("GDAL no pudo crear el PNG: {}".format(exc))
        if png_ds is None:
            raise QgsProcessingException("GDAL no pudo crear el PNG.")
        ds_refs["png"] = png_ds
        png_ds.FlushCache()
        ds_refs["png"] = None
        png_ds = None
        ds_refs["png_vrt"] = None
        png_vrt_ds = None

        feedback.setProgress(100.0)
        feedback.pushInfo("Proceso terminado correctamente.")
        feedback.pushInfo("GeoTIFF: {}".format(out_tif))
        feedback.pushInfo("PNG: {}".format(out_png))

        return {
            self.OUTPUT_TIF: out_tif,
            self.OUTPUT_PNG: out_png,
            self.ANGLE_DEG: theta,
            self.WIDTH_PX: width_px,
            self.HEIGHT_PX: height_px,
        }
