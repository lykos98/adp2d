# ADP2D API Documentation

## Overview

ADP2D is a Python library for detecting and analyzing clusters in 2D data using the Advanced Density Peak (ADP) method. The library interfaces with a C shared library (`libadp2d.so`) for performance-critical operations.

## Main Class: `Data`

The primary interface to the library is through the `Data` class.

### Initialization

```python
data = Data(array)
```

**Parameters:**
- `array` (numpy.ndarray): A 2D array of data values to analyze

**Raises:**
- `TypeError`: If the input is not a 2D numpy array

### Methods

#### `computeDensityFromImg`

Computes the density of data points from an image.

```python
data.computeDensityFromImg(img, mask=None, r=15, algorithm="MEAN", use_log=True, use_adaptive_radius=False, param=None)
```

**Parameters:**
- `img` (numpy.ndarray): 2D image array
- `mask` (numpy.ndarray, optional): Integer mask array (defaults to ones_like(img))
- `r` (int, optional): Neighborhood radius (default: 15)
- `algorithm` (str, optional): Density algorithm ("MEAN", "MEDIAN", "GAUSSIAN", or "SPLINE") (default: "MEAN")
- `use_log` (bool, optional): Whether to use logarithmic density (default: True)
- `use_adaptive_radius` (bool, optional): Whether to use adaptive radius (default: False)
- `param` (int, optional): Smoothing parameter. If None, defaults to r // 2 for "SPLINE", or r // 3 for other algorithms (default: None)

**Notes:**
- This method must be called before clustering operations
- Internally calls the C function `computeDensityFromImg`
- Supported algorithms: MEAN (adaptive kernel), MEDIAN (fixed kernel with median), GAUSSIAN (gaussian kernel), SPLINE (cubic spline kernel)

#### `computeClusteringADP`

Performs clustering using the Advanced Density Peak method.

```python
data.computeClusteringADP(Z, halo=False, min_area=10, useSparse="auto", splitPerThread=False)
```

**Parameters:**
- `Z` (float): Z parameter for the ADP method
- `halo` (bool, optional): Whether to compute halo properties (default: False)
- `min_area` (int, optional): Minimum area for cluster consideration (default: 10)
- `useSparse` (str, optional): Border storage method ("auto", "yes", or "no") (default: "auto")
- `splitPerThread` (bool, optional): Whether to split processing per thread (default: False)

**Raises:**
- `ValueError`: If density has not been computed yet

**Notes:**
- Automatically determines sparse storage usage based on data size when `useSparse="auto"`
- Internally calls the C function `adpWrapper`

#### `getClusterAssignment`

Retrieves the cluster label for each data point.

```python
labels = data.getClusterAssignment()
```

**Returns:**
- `numpy.ndarray`: 2D array of cluster labels with same shape as input data

**Raises:**
- `ValueError`: If clustering has not been computed yet

**Notes:**
- Internally calls the C function `export_cluster_assignment`

#### `getDensity`

Retrieves the density values for each data point.

```python
density = data.getDensity()
```

**Returns:**
- `numpy.ndarray`: Array of density values (logarithmic if `use_log=True` was set)

**Raises:**
- `ValueError`: If density has not been computed yet

#### `getKstar`

Retrieves the k-star values for each data point.

```python
kstar = data.getKstar()
```

**Returns:**
- `numpy.ndarray`: Array of k-star values (nearest higher density neighbor)

**Notes:**
- No explicit check for density computation (relies on internal state)

#### `getDensityError`

Retrieves the density error values for each data point.

```python
density_error = data.getDensityError()
```

**Returns:**
- `numpy.ndarray`: Array of density error values

**Raises:**
- `ValueError`: If density has not been computed yet

#### `setRhoErrK`

Sets density error and k values (typically used internally).

```python
data.setRhoErrK(rho, err, k)
```

**Parameters:**
- `rho` (numpy.ndarray): Density values
- `err` (numpy.ndarray): Density error values
- `k` (numpy.ndarray): K-star values

**Notes:**
- Internally calls the C function `setRhoErrK`
- Marks density computation as complete

#### `computeSourcesProperties`

Computes properties of detected sources (clusters).

```python
properties = data.computeSourcesProperties(min_area=10, header=None)
```

**Parameters:**
- `min_area` (int, optional): Minimum area for source consideration (default: 10)
- `header` (astropy.io.fits.Header, optional): FITS header for world coordinate conversion

**Returns:**
- `dict`: Dictionary containing source properties:
  - `centers_of_mass`: Array of [x, y] centroid coordinates
  - `world_coord`: Array of [RA, Dec] world coordinates (if header provided)
  - `areas`: Array of pixel counts per source
  - `ellipticities`: Array of ellipticity values
  - `major_axes`: Array of semi-major axis lengths
  - `minor_axes`: Array of semi-minor axis lengths
  - `position_angle`: Array of position angles (degrees)
  - `r_max`: Array of maximum radius values
  - `semi_major_angles`: Array of semi-major axis angles (degrees)
  - `x_limits`: Array of [xmin, xmax] bounds per source
  - `y_limits`: Array of [ymin, ymax] bounds per source
  - `flux`: Array of total flux values per source
  - `parent_id`: Array of parent IDs per source

**Notes:**
- Internally calls C functions `compute_covs` and `compute_eigensystems`
- Requires cluster assignment to be computed first

#### `exportSourcesCatalogue`

Exports source properties to a FITS file.

```python
data.exportSourcesCatalogue(fname="catalogue.fits")
```

**Parameters:**
- `fname` (str, optional): Output filename (default: "catalogue.fits")

**Raises:**
- `ValueError`: If source properties have not been computed

**Notes:**
- Creates a binary table FITS file with columns for all source properties
- Uses astropy.io.fits for FITS file handling

#### `writePNG`

Writes a PNG visualization of the clustering results.

```python
data.writePNG(fname, scale=0.5)
```

**Parameters:**
- `fname` (str): Output filename
- `scale` (float, optional): Scale factor for output image (default: 0.5)

**Notes:**
- Automatically adjusts scale if dimensions exceed 4000 pixels
- Internally calls the C function `tiny_colorize`

## Internal Data Structures

The library uses several ctypes structures to interface with the C library:

- `HeapNode`: Value and array index for heap operations
- `Heap`: Heap data structure for nearest neighbor search
- `luDynamicArray`: Dynamic array of indices
- `DatapointInfo`: Information about each data point including density, neighbors, and cluster assignment
- `Border_t`: Border information between clusters
- `SparseBorder_t`: Sparse representation of cluster borders
- `AdjList`: Adjacency list for sparse border storage
- `Clusters`: Cluster collection containing centers, borders, and metadata

## Usage Example

```python
import numpy as np
from adp2d import Data

# Create or load 2D data
data = np.random.rand(100, 100)  # Example data

# Initialize Data object
d = Data(data)

# Compute density
d.computeDensityFromImg(data, r=20)

# Perform clustering
d.computeClusteringADP(Z=2.0, min_area=15)

# Get results
labels = d.getClusterAssignment()
density = d.getDensity()

# Compute and export source properties
props = d.computeSourcesProperties(min_area=15)
d.exportSourcesCatalogue("sources.fits")

# Create visualization
d.writePNG("clustering.png", scale=0.8)
```

## Requirements

- Python 3.x
- NumPy
- SciPy
- Astropy
- Numba
- Shared library: `libadp2d.so` (must be in the same directory as adp2d.py)

## Notes

- The library automatically detects whether to use float32 or float64 based on the shared library compilation
- The library automatically detects whether to use int32 or uint64 indices based on the shared library compilation
- Memory management is handled automatically through the `__del__` method
- For large datasets (>20,000 points), sparse border storage is automatically enabled