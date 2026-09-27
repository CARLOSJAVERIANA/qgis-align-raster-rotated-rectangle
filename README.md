# QGIS Align Raster to Rotated Rectangle

A Python Processing script for **QGIS** that aligns a GeoTIFF raster to the orientation of a rotated rectangle, resamples the original raster onto the resulting grid using **GDAL Warp**, and exports both an aligned **GeoTIFF** and a georeferenced **PNG**.

The script is designed for cases where a raster needs to be exported according to a specific non-axis-aligned rectangular area while preserving its spatial reference and raster information.

## Features

- Aligns a raster grid with a rotated rectangular polygon.
- Automatically determines the raster orientation from the rectangle.
- Allows the user to choose which rectangle edge becomes the horizontal axis.
- Preserves the original raster resolution automatically or allows a custom output pixel size.
- Uses GDAL Warp for raster resampling.
- Supports:
  - Nearest neighbour
  - Bilinear
  - Cubic
  - Cubic spline
  - Lanczos
- Preserves all raster bands in the GeoTIFF output.
- Preserves:
  - data type
  - NoData values
  - alpha bands
  - metadata
  - color tables
- Builds a common validity mask across raster bands.
- Supports rasters with different NoData configurations between bands.
- Generates an internal validity mask when appropriate.
- Supports grayscale and RGB PNG export.
- Preserves or reconstructs alpha transparency.
- Generates a georeferenced PNG with:
  - World file (`.wld`)
  - Auxiliary CRS file (`.aux.xml`)
- Supports Byte and UInt16 exact PNG output.
- Provides several visualization modes for other raster data types, including Float32.
- Supports RGB data normalized between 0 and 1.
- Uses a common RGB percentile range when requested to preserve color balance.
- Supports multithreaded GDAL processing.
- Allows the maximum memory assigned to GDAL Warp to be controlled.
- Processes large masks and raster statistics in blocks to reduce memory consumption.

## Requirements

The script was developed for:

- **QGIS 3.40 or newer**
- **GDAL 3.8 or newer**
- Python environment included with QGIS
- NumPy
- GDAL/OGR Python bindings

The script imports functionality from:

- `qgis.core`
- `osgeo.gdal`
- `osgeo.ogr`
- `numpy`

It is intended to run inside the **QGIS Processing framework**, rather than as a standalone Python program.

## Inputs

### Input GeoTIFF

A file-based raster accessible through the QGIS GDAL provider.

The raster must:

- contain at least one data band;
- have a valid GeoTransform;
- have a valid CRS;
- use a homogeneous data type across its bands.

A projected CRS is recommended when accurate metric dimensions and rotation angles are important.

### Rotated rectangle

A polygon layer containing **exactly one polygon feature**.

If a layer contains multiple polygons, select the desired feature and use the QGIS option to process selected features only.

The script accepts and cleans:

- closing vertices;
- duplicated vertices;
- collinear vertices;
- single-part multipart geometries;
- curved polygon geometries.

The resulting geometry is validated as a rectangle according to configurable angular and opposite-side length tolerances.

If the geometry does not satisfy the rectangular tolerances, the optional **Oriented Minimum Bounding Box (OMBB)** fallback can be used.

## Orientation

The user can choose which side of the rectangle becomes the horizontal axis of the output raster.

Available modes are:

- Automatic: longest side
- Automatic: shortest side
- Edge 0: P0 → P1
- Edge 1: P1 → P2
- Edge 2: P2 → P3
- Edge 3: P3 → P0

The horizontal direction is normalized automatically to prevent the exported image from being unintentionally inverted.

## Output resolution

The output pixel size can be defined manually.

When the value is set to `0`, the script automatically estimates the source raster pixel spacing along the two axes of the rotated output grid and preserves approximately the original raster density.

## Resampling methods

The following GDAL resampling algorithms are available:

- Nearest neighbour
- Bilinear
- Cubic
- Cubic spline
- Lanczos

Bilinear is used as the default.

## Raster validity and transparency

The script constructs a common validity mask before resampling.

A source pixel is considered valid only when all data bands are valid and, when an alpha band exists, its alpha value is greater than zero.

This prevents pixels containing invalid values in individual RGB channels from appearing as opaque valid pixels.

The rectangle can also be used as a strict source mask before interpolation. This prevents interpolation kernels such as Bilinear, Cubic, and Lanczos from incorporating samples outside the requested rectangular area.

## PNG conversion modes

The script provides several methods for converting raster values to PNG:

1. **Exact**
   - Preserves Byte or UInt16 values.

2. **Visual Min–Max**
   - Scales each band independently using the exact minimum and maximum of valid pixels.

3. **Visual 2–98%**
   - Scales each band independently using the 2nd and 98th percentiles.

4. **Visual RGB shared 2–98%**
   - Uses a common range across RGB channels to better preserve the original color balance.

5. **Visual RGB 0–1**
   - Maps normalized Float RGB values directly from `0–1` to `0–255`.

6. **Automatic RGB**
   - Automatically uses the `0–1` method when appropriate; otherwise it uses a shared RGB scaling strategy.

For unsupported exact PNG data types, the script automatically falls back to a visual conversion while preserving the original data type in the GeoTIFF.

## Outputs

The algorithm generates two principal files.

### Aligned GeoTIFF

The GeoTIFF preserves the raster data as faithfully as possible, including:

- all bands;
- original data type;
- NoData information;
- alpha information;
- raster metadata;
- color interpretation;
- color tables when available;
- CRS;
- rotated GeoTransform.

The GeoTIFF is created using tiled DEFLATE compression and multithreaded processing.

### Georeferenced PNG

The PNG is exported as grayscale or RGB with an alpha channel.

The export also creates georeferencing information through:

- a World File (`.wld`);
- a GDAL auxiliary file (`.aux.xml`) containing CRS information.

This allows the PNG to retain its spatial relationship when used in GIS-compatible workflows.

## Additional outputs

The QGIS Processing algorithm also reports:

- rotation angle of the horizontal axis;
- output raster width in pixels;
- output raster height in pixels.

## Installation in QGIS

1. Download `align_raster_to_rotated_rectangle.py`.

2. Open **QGIS**.

3. Open:

   `Processing → Toolbox`

4. In the Processing Toolbox, locate:

   `Scripts`

5. Choose:

   `Create New Script from File...`

   or add the script to your QGIS Processing scripts directory.

6. Once installed, the algorithm should appear in the Processing Toolbox as:

   **Alinear raster a rectángulo rotado y exportar PNG**

   under the group:

   **Raster personalizado**

## Basic usage

1. Load the source GeoTIFF into QGIS.

2. Create or load a polygon representing the desired rotated rectangular area.

3. Make sure the rectangle and raster have valid CRS information.

4. Run:

   **Alinear raster a rectángulo rotado y exportar PNG**

5. Select:
   - input GeoTIFF;
   - rotated rectangle;
   - orientation mode;
   - resampling method;
   - output pixel size, if required;
   - PNG conversion mode;
   - GeoTIFF output path;
   - PNG output path.

6. Run the algorithm.

The resulting GeoTIFF and PNG will be aligned with the selected rectangle.

## Notes

For geometrically meaningful distances and rotation calculations, a **projected CRS** is recommended.

When using a geographic CRS, the script operates directly in the planar coordinate values of that CRS, which may not provide metrically accurate distances.

If the selected rectangle extends beyond the source raster footprint, the exterior area will be exported as transparent or NoData.

## Version

Current script version:

**3.2**

Main improvements in version 3.2 include:

- improved common validity-mask handling;
- conservative border masking for interpolated rasters;
- RGB `0–1` mapping for Float32 imagery;
- common percentile scaling across RGB channels;
- improved preservation of RGB color balance;
- GDAL 3.10-compatible rectangle rasterization.

## Author

**Carlos Acosta**

If you use, modify, or redistribute this software, please retain the attribution to the original author as specified in the license.

### Suggested citation

Acosta, C. (2026). *QGIS Align Raster to Rotated Rectangle* (Version 3.2) [Python software].

Once a permanent GitHub repository URL and/or DOI is available, it can be added to this citation.

## License

This software is distributed under the **BSD 3-Clause License**.

Copyright © 2026 Carlos Acosta.

Redistribution and modification are permitted under the conditions described in the `LICENSE` file. Redistributions must retain the copyright notice, license conditions, and disclaimer.

See `LICENSE` for the complete terms.

## Acknowledgements

This script relies on the QGIS Python API, GDAL/OGR, and NumPy.

QGIS and GDAL are independent open-source projects and are subject to their respective licenses.# qgis-align-raster-rotated-rectangle
A Python Processing script for QGIS that aligns GeoTIFF rasters to a rotated rectangle and exports georeferenced GeoTIFF and PNG outputs.
