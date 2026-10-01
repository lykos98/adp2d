import ctypes as ct
import numpy as np  
import os
import time
import numba
from astropy.io import fits
from astropy.wcs import WCS

# Default ctypes scalar types for the float and index ABI widths
# These are overridden at runtime in __init__ based on FloatAndUintSize()
# (32-bit builds switch them to ct.c_float / ct.c_uint32)
ctFloatType = ct.c_double
ctIdxType   = ct.c_uint64

class HeapNode(ct.Structure):
    """A single node in the nearest-neighbour heap"""
    _fields_ = [
        ("value", ctFloatType),     # density of the neighbour
        ("array_idx", ctIdxType)    # flat index of the neighbour pixel
    ]

class Heap(ct.Structure):
    """Fixed-capacity binary heap used to keep the k nearest neighbours per pixel"""
    _fields_ = [
        ("N", ctIdxType),          # capacity
        ("count", ctIdxType),      # current number of elements
        ("data", ct.POINTER(HeapNode))
    ]

class luDynamicArray(ct.Structure):
    """Growable array of indices (e.g. cluster center / member lists)"""
    _fields_ = [
        ("data", ct.POINTER(ctIdxType)),  # elements
        ("size", ctIdxType),              # allocated capacity
        ("count", ctIdxType)              # number of elements in use
    ]

class DatapointInfo(ct.Structure):
    """Per-pixel info produced by the density estimation, consumed by clustering"""
    _fields_ = [
        ("g", ctFloatType),           # density minus error (used for border/max-g selection)
        ("ngbh", Heap),               # nearest-neighbour heap
        ("array_idx", ctIdxType),     # flat pixel index
        ("log_rho", ctFloatType),     # log-density estimate
        ("log_rho_c", ctFloatType),   # corrected/shifted log-density (>=0)
        ("log_rho_err", ctFloatType), # density error
        ("kstar", ctIdxType),         # adaptive kernel radius used for this pixel
        ("is_center", ct.c_int),      # 1 if this pixel is a cluster center (peak)
        ("cluster_idx", ct.c_int)     # assigned cluster label (-1 if none)
    ]

    def __repr__(self):
        return f"{{ g: {self.g}, ngbh: {self.ngbh}, array_idx: {self.array_idx}, log_rho: {self.log_rho}, log_rho_c: {self.log_rho_c}, log_rho_err: {self.log_rho_err}, kstar: {self.kstar}, isCenter: {self.is_center}, clusterIdx: {self.cluster_idx} }}"

class Border_t(ct.Structure):
    """Dense per-pair border record between two clusters"""
    _fields_ = [
        ("idx", ctIdxType),        # flat index of the chosen border pixel
        ("density", ctFloatType),  # saddle density (log_rho_c) for this pair
        ("error", ctFloatType)     # uncertainty on the saddle density
    ]

class SparseBorder_t(ct.Structure):
    """Sparse border entry (used when dense nclus x nclus storage is too big)"""
    _fields_ = [
        ("i", ctIdxType),         # first cluster
        ("j", ctIdxType),         # second cluster
        ("idx", ctIdxType),       # flat index of the border pixel
        ("density", ctFloatType), # saddle density
        ("error", ctFloatType)    # saddle error
    ]

class AdjList(ct.Structure):
    """Adjacency list of borders incident to a cluster (sparse storage)"""
    _fields_ = [
        ("count", ctIdxType),              # number of borders
        ("size", ctIdxType),               # allocated capacity
        ("data", ct.POINTER(SparseBorder_t))
    ]

class Clusters(ct.Structure):
    """Result of the clustering step: cluster centers and their borders"""
    _fields_ = [
        ("UseSparseBorders", ct.c_int),                # 1 if using sparse border storage
        ("SparseBorders", ct.POINTER(AdjList)),        # per-cluster border adjacency lists
        ("centers", luDynamicArray),                   # indices of cluster centers
        ("borders", ct.POINTER(ct.POINTER(Border_t))), # dense nclus x nclus border matrix
        ("__borders_data", ct.POINTER(Border_t)),      # backing storage for 'borders'
        ("n", ctIdxType)                               # number of clusters
    ]

class Data():
    def __init__(self, data : np.array):
        """Data object for dadaC library

        Args:
            data (np.array): 2d array to use in searching for clusters 

        Raises:
            TypeError: Raises TypeError if a type different from a matrix is passed 

        """
        # Load the compiled shared library. This must be rebuilt (make) after any C change
        path     = os.path.join(os.path.dirname(__file__), "bin/libadp2d.so")
        self.lib = ct.CDLL(path)

        # Query the C library for the runtime width of float/idx types so the ctypes
        # field types (ctFloatType/ctIdxType) match the compiled ABI
        global ctFloatType, ctIdxType
        s = self.lib.FloatAndUintSize()

        self.__useFloat32 = s < 1
        self.__useInt32   = (s == 0) or (s == 2)

        # Per-step state flags used to guard the order of the ADP pipeline
        # (ngbh -> density -> clustering)
        self.state = {
                    "ngbh"        : False,
                    "id"          : False,
                    "density"     : False,
                    "clustering"  : False,
                    "useFloat32"  : self.__useFloat32,
                    "useInt32"    : self.__useInt32,
                    "useSparse"   : None,
                    "computeHalo" : None 
                    }

        # Convert the input to the same float width as the compiled library
        if self.__useFloat32:
            self.data  = np.ascontigousarray(data.astype(np.float32))
            ctFloatType = ct.c_float
        else:
            self.data = np.ascontiguousarray(data.astype(np.float64))
        
        # Match the index type to the compiled ABI
        if self.__useInt32:
            ctIdxType = ct.c_uint32

        # ---------------------------------------------------------------------------
        # Bind the C functions to ctypes signatures (function pointers from the .so).
        # ---------------------------------------------------------------------------
        self.__computeDensityFromImg = self.lib.computeDensityFromImg
        self.__computeDensityFromImg.argtypes = [   np.ctypeslib.ndpointer(ctFloatType),
                                                    np.ctypeslib.ndpointer(np.int32),
                                                    ct.c_int32,
                                                    ct.c_int32,
                                                    ct.c_int32,
                                                    ct.c_int32,
                                                    ct.c_bool,
                                                    ct.c_bool,
                                                    ct.c_int32]
        self.__computeDensityFromImg.restype  = ct.POINTER(DatapointInfo)

        # overwrite rho/err/kstar from external arrays
        self.__setRhoErrK = self.lib.setRhoErrK  
        self.__setRhoErrK.argtypes = [  ct.POINTER(DatapointInfo),    
                                        np.ctypeslib.ndpointer(ctFloatType), 
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctIdxType),
                                        ct.c_uint64]
        
        # shift log_rho -> log_rho_c (>=0)
        self.__computeCorrection          = self.lib.computeCorrection
        self.__computeCorrection.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(ct.c_int32), ctIdxType, ct.c_double]

        # density-peak assignment (cluster centers)
        self.__H1          = self.lib.Heuristic1
        self.__H1.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(ct.c_int32), ct.c_uint64, ct.c_uint64, ct.c_int32, ct.c_bool]
        self.__H1.restype  = Clusters
        
        # Runtime border-mode setter (mode 0 = max-g, 3 = percentile).
        self.lib.set_adp_border.argtypes = [ct.c_int, ct.c_float]
        self.lib.set_adp_border.restype  = None

        # main entry: H1+H2+H3 + area filter
        self.__adpWrapper          = self.lib.adpWrapper
        self.__adpWrapper.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(ct.c_int32), 
                                      ct.c_uint64, ct.c_uint64, ct.c_int32, ct.c_float, ct.c_bool, ct.c_bool]
        self.__adpWrapper.restype  = Clusters
        
        # Allocate the Clusters struct; arg2 (s) selects 1=sparse / 0=dense border storage
        self.__ClustersAllocate          = self.lib.Clusters_allocate
        self.__ClustersAllocate.argtypes = [ct.POINTER(Clusters), ct.c_int]

        # find borders between clusters
        self.__H2          = self.lib.Heuristic2
        self.__H2.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), 
                              np.ctypeslib.ndpointer(ct.c_int32), ct.c_uint64, ct.c_uint64, ct.c_int32, ct.c_bool]

        # merge clusters / compute halo
        self.__H3          = self.lib.Heuristic3
        self.__H3.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), ct.c_double, ct.c_int, ct.c_int, ct.c_bool]
        
        # Free various Datapoints and Clusters
        self.__freeDatapoints          = self.lib.freeDatapointArray
        self.__freeDatapoints.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]
        self.__freeClusters            = self.lib.Clusters_free
        self.__freeClusters.argtypes   = [ct.POINTER(Clusters)]
        
        # segmap output
        self._export_cluster_assignment          = self.lib.export_cluster_assignment
        self._export_cluster_assignment.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(np.int32), ct.c_uint64]

        # shape (cov) eigendecomposition
        self.__compute_eigensystems          = self.lib.compute_eigensystems
        self.__compute_eigensystems.argtypes = [np.ctypeslib.ndpointer(ct.c_double),
                                                np.ctypeslib.ndpointer(ct.c_double),
                                                np.ctypeslib.ndpointer(ct.c_double),
                                                ct.c_int32]
        
        # per-object moments/covariances
        self.__compute_covs          = self.lib.compute_covs
        self.__compute_covs.argtypes = [np.ctypeslib.ndpointer(ct.c_double),
                                        np.ctypeslib.ndpointer(ct.c_int32),
                                        np.ctypeslib.ndpointer(ct.c_int32),
                                        ct.c_int32, ct.c_int32, ct.c_int32,
                                        np.ctypeslib.ndpointer(ct.c_double),
                                        np.ctypeslib.ndpointer(ct.c_double),
                                        np.ctypeslib.ndpointer(ct.c_double),
                                        np.ctypeslib.ndpointer(ct.c_int32),
                                        np.ctypeslib.ndpointer(ct.c_double),
                                        np.ctypeslib.ndpointer(ct.c_int32),
                                        np.ctypeslib.ndpointer(ct.c_int32),
                                        np.ctypeslib.ndpointer(ct.c_int32)]

        # PNG export of the segmentation
        self.__write_png          = self.lib.tiny_colorize
        self.__write_png.argtypes = [ct.c_char_p,
                                     ct.POINTER(DatapointInfo),
                                     np.ctypeslib.ndpointer(ctFloatType),
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,]
        
        # Validate input shape and initialise per-object cached outputs
        if len(self.data.shape) != 2:
            raise TypeError("Please provide a 2d numpy array")

        # Cached results, filled in as the pipeline runs (density -> clustering)
        self.__datapoints   = None # C-side Datapoint_info array (density output)
        self.__clusters     = None # C-side Clusters (H1/H2/H3 output)

        # Dimensions of the input data
        self.n              = self.data.shape[0] # number of rows (or total points)
        self.dims           = self.data.shape[1] # number of columns (or feature dims)

        # Optional per-step results (cached on first access)
        self.k              = None # adaptive kernel radius map
        self.id             = None

        # Cached output arrays, populated by the respective get* methods / export steps
        self.clusterAssignment = None # 2D segmentation map (labels per pixel), from getClusterAssignment()
        self.neighbors         = None # neighbour info (cached)
        self.borders           = None # per-cluster border data (cached)
        self.density           = None # per-pixel density array, from getDensity()
        self.densityError      = None # per-pixel density error, from getDensityError()

    def computeDensityFromImg(self, img, mask=None, r=15, algorithm="MEAN", use_log=True, use_adaptive_radius=False, param=None):
        """Estimate the density field from the image.

        Builds a per-pixel log-density (and its error) via the selected algorithm
        and stores it for the subsequent clustering steps.

        Args:
            img:      2D image array.
            mask:     detection mask (0 = background/excluded). If None, all pixels used.
            r:        base density radius (pixels).
            algorithm: "MEAN" (adaptive), "MEDIAN", "GAUSSIAN", or "SPLINE".
            use_log:  apply log to the density.
            use_adaptive_radius: grow the kernel radius adaptively until stable.
            param:    extra smoothing parameter; if None, derived from r.

        """
        if mask is None:
            mask = np.ones_like(img, dtype=np.int32)
        self.n = np.prod(img.shape)
        mask = mask.astype(np.int32)
        self.nrows, self.ncols = img.shape
        self.img = img
        self.mask = mask

        # Map the density algorithm keyword to its integer code (mirrors density_alg_t in C)
        if algorithm == "MEAN":
            alg_val = 0
        elif algorithm == "MEDIAN":
            alg_val = 1
        elif algorithm == "GAUSSIAN":
            alg_val = 2
        elif algorithm == "SPLINE":
            alg_val = 3
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}. Use 'MEAN', 'MEDIAN', 'GAUSSIAN', or 'SPLINE'.")

        if param is None:  # default smoothing param per algorithm
            if algorithm == "SPLINE":
                param = r // 2
            else:
                param = r // 3

        # Run the C density estimation; returns the Datapoint_info array.
        self.__datapoints = self.__computeDensityFromImg(img, mask, self.nrows, self.ncols, r, alg_val, use_log, use_adaptive_radius, param)
        self.state["density"] = True

    def computeClusteringADP(self,Z : float, halo = False, min_area = 10, useSparse = "auto", splitPerThread = False,  border="maxg", border_perc=0.8):
        """Compute clustering via the Advanced Density Peak method

        Args:
            Z:            significance level for the merge test.
            halo:         compute halo assignment.
            min_area:     minimum pixel area to keep a cluster.
            useSparse:    "auto"/"yes"/"no" — sparse border storage.
            splitPerThread: process each detection patch independently.
            border:       "maxg" (merge-aggressive) or "percentile" (split-capable).
            border_perc:  percentile used by the "percentile" border mode. 

        Raises:
            ValueError: Raises value error if density is not computed, use `Data.computeDensity()` method

        """
        if not self.state["density"]:
            raise ValueError("Please compute density before calling this function")
        if useSparse == "auto":
            if self.n > 2e4:
                self.state["useSparse"] = True
            else: 
                self.state["useSparse"] = False 
        elif useSparse == "y":
            self.state["useSparse"] = True
        else:
            self.state["useSparse"] = False
        
        # Border mode is chosen at RUNTIME so we can switch between
        # "maxg" and "percentile" runs without recompiling the C library
        border_map = {"maxg":0, "percentile":3}
        bstat      = border_map[border] 
        self.lib.set_adp_border(ct.c_int(bstat), ct.c_float(border_perc))

        # change import of libs, import only the wrapper 
        self.state["computeHalo"] = halo 
        self.Z                    = Z
        self.n                    = np.prod(self.img.shape)
        self.min_area             = min_area
        self.halo                 = halo

        # Call the main C clustering routine (H1/H2/H3 wrapped in adpWrapper)
        self.__clusters = self.__adpWrapper(self.__datapoints, self.mask, self.nrows, self.ncols, self.min_area, self.Z, self.halo, splitPerThread)

        self.state["clustering"] = True
        self.clusterAssignment   = None

    def getClusterAssignment(self) -> list:
        """Retrieve cluster assignment

        Raises:
            ValueError: Raises error if clustering is not computed, use `Data.computeClusteringADP(Z)` 

        Returns:
            List of cluster labels
            
        """
        self.clusterAssignment = np.ascontiguousarray(np.zeros((self.nrows, self.ncols), np.int32))
        self._export_cluster_assignment(self.__datapoints, self.clusterAssignment, self.n)
        return self.clusterAssignment

    def getBorders(self):
        raise NotImplemented("It's difficult I have to think about it")

    def getDensity(self) -> list:
        """Retrieve list of density values

        Raises:
            ValueError: Raise error if density is not computed, use `Data.computeDensity()` 

        Returns:
            List of density values
            
        """
        if self.density is None:
            if self.state["density"]:
                self.density = np.array([float(self.__datapoints[j].log_rho) for j in range(self.n)])
                return self.density
            else:
                raise ValueError("Density is not computed yet")
        else:
            return self.density

    def getKstar(self) -> list:
        """Retrieve list of density values

        Raises:
            ValueError: Raise error if density is not computed, use `Data.computeDensity()` 

        Returns:
            List of density values

        """
        self.kstar = np.array([float(self.__datapoints[j].kstar) for j in range(self.n)]) 
        return self.kstar

    def getDensityError(self):
        """Retrieve list of density error values

        Raises:
            ValueError: Raise error if density is not computed, use `Data.computeDensity()` 

        Returns: List of density error values
        """
        if self.densityError is None:
            if self.state["density"]:
                self.densityError = np.array([float(self.__datapoints[j].log_rho_err) for j in range(self.n)])
                return self.densityError
            else:
                raise ValueError("Density Error is not computed yet")
        else:
            return self.densityError

    def setRhoErrK(self, rho, err, k):
        """Overwrite the per-pixel density, error and kernel radius from external arrays."""
        self.__setRhoErrK(self.__datapoints, rho, err, k, self.n)
        self.state["density"] = True
        return

    def exportSourcesCatalogue(self, fname = "catalogue.fits" ):
        """Write the per-source measurements to a FITS binary table
        
        One row per deblended source. PARENT_ID is the detection-patch label the
        source was deblended from (-1 for non-deblended sources); it is the key used
        by the two-pass over-split detection.
        """
        columns = [
            fits.Column(name = 'SOURCE_ID', format = 'K'),
            fits.Column(name = 'PARENT_ID', format = 'J'),
            fits.Column(name = 'X_CENTER', format = 'D'),
            fits.Column(name = 'Y_CENTER', format = 'D'),
            fits.Column(name = 'XWIN_WORLD', format = 'D'),
            fits.Column(name = 'YWIN_WORLD', format = 'D'),
            fits.Column(name = 'X_MIN', format = 'J'),
            fits.Column(name = 'X_MAX', format = 'J'),
            fits.Column(name = 'Y_MIN', format = 'J'),
            fits.Column(name = 'Y_MAX', format = 'J'),
            fits.Column(name = 'ALPHA', format = 'D'),
            fits.Column(name = 'BETA', format = 'D'),
            fits.Column(name = 'ELLIPTICITY', format = 'D'),
            fits.Column(name = 'R_MAX', format = 'D'),
            fits.Column(name = 'POSITION_ANGLE', format = 'D'),
            fits.Column(name = 'SEMI_MAJOR_ANGLE', format = 'D'),
            fits.Column(name = 'FLUX_TOT', format = 'D'),
            fits.Column(name = 'ISOAREA', format = 'J'),
            fits.Column(name = 'SKIPPED', format = 'J'),]
        if self.sources_properties is None:
            raise ValueError("Sources properties are not computed")

        nsources =  self.sources_properties["areas"].shape[0]

        hdu       = fits.BinTableHDU.from_columns(columns)
        rec_array = np.zeros(nsources, dtype = hdu.data.dtype)
         
        if nsources > 1:
            for i in range(nsources):
                rec_array[i]['SOURCE_ID']        = i + 1
                rec_array[i]['PARENT_ID']        = self.sources_properties["parent_id"][i] 
                rec_array[i]['X_CENTER']         = self.sources_properties["centers_of_mass"][i][0]
                rec_array[i]['Y_CENTER']         = self.sources_properties["centers_of_mass"][i][1]
                rec_array[i]['XWIN_WORLD']       = self.sources_properties["world_coord"][i][0]
                rec_array[i]['YWIN_WORLD']       = self.sources_properties["world_coord"][i][1]
                rec_array[i]['X_MIN']            = self.sources_properties["x_limits"][i][0]
                rec_array[i]['X_MAX']            = self.sources_properties["x_limits"][i][1]
                rec_array[i]['Y_MIN']            = self.sources_properties["y_limits"][i][0]
                rec_array[i]['Y_MAX']            = self.sources_properties["y_limits"][i][1]
                rec_array[i]['ALPHA']            = self.sources_properties["minor_axes"][i]
                rec_array[i]['BETA']             = self.sources_properties["major_axes"][i]
                rec_array[i]['ELLIPTICITY']      = self.sources_properties["ellipticities"][i]
                rec_array[i]['R_MAX']            = self.sources_properties["r_max"][i]
                rec_array[i]['POSITION_ANGLE']   = self.sources_properties["position_angle"][i]
                rec_array[i]['SEMI_MAJOR_ANGLE'] = self.sources_properties["semi_major_angles"][i]
                rec_array[i]['FLUX_TOT']         = self.sources_properties["flux"][i]
                rec_array[i]['ISOAREA']          = self.sources_properties["areas"][i]
                rec_array[i]['SKIPPED']          = 0
        else:
            # Single-source case: scalar assignment (no per-row loop)
            rec_array['SOURCE_ID']        = 1
            rec_array['PARENT_ID']        = self.sources_properties["parent_id"] 
            rec_array['X_CENTER']         = self.sources_properties["centers_of_mass"][0][0]
            rec_array['Y_CENTER']         = self.sources_properties["centers_of_mass"][0][1]
            rec_array['XWIN_WORLD']       = self.sources_properties["world_coord"][0][0]
            rec_array['YWIN_WORLD']       = self.sources_properties["world_coord"][0][1]
            rec_array['X_MIN']            = self.sources_properties["x_limits"][0][0]
            rec_array['X_MAX']            = self.sources_properties["x_limits"][0][1]
            rec_array['Y_MIN']            = self.sources_properties["y_limits"][0][0]
            rec_array['Y_MAX']            = self.sources_properties["y_limits"][0][1]
            rec_array['ALPHA']            = self.sources_properties["minor_axes"]
            rec_array['BETA']             = self.sources_properties["major_axes"]
            rec_array['ELLIPTICITY']      = self.sources_properties["ellipticities"]
            rec_array['R_MAX']            = self.sources_properties["r_max"]
            rec_array['POSITION_ANGLE']   = self.sources_properties["position_angle"]
            rec_array['SEMI_MAJOR_ANGLE'] = self.sources_properties["semi_major_angles"]
            rec_array['FLUX_TOT']         = self.sources_properties["flux"]
            rec_array['ISOAREA']          = self.sources_properties["areas"]
            rec_array['SKIPPED']          = 0

        hdu.data = rec_array
        hdu.writeto(fname, overwrite=True)

    def computeSourcesProperties(self, min_area = 10, header=None):
        """Compute per-source measurements from the cluster assignment (segmap)

        Fills self.sources_properties (com, flux, shape, parent_id, limits, ...),
        then applies the min_area filter. PARENT_ID comes from the detection mask
        """
        wcs = WCS(header) if header is not None else None

        start = time.time()
        print("Computing sources properties")
        self.getClusterAssignment()

        segmentation_map = self.clusterAssignment.reshape((self.nrows, self.ncols)).astype(np.int32)
        nrows, ncols     = self.data.shape
        nlabs            = self.__clusters.centers.count
        
        # Per-cluster output buffers filled by the C compute_covs call
        coms                = np.zeros((nlabs,2), dtype = np.float64)
        ra_dec              = np.zeros((nlabs, 2), dtype=np.float64)
        covariance_matrices = np.zeros((nlabs,2,2), dtype = np.float64)
        areas               = np.zeros((nlabs), dtype = np.int32)
        a_images            = np.zeros((nlabs), dtype = np.float64)
        b_images            = np.zeros((nlabs), dtype = np.float64)
        ellipticities       = np.zeros((nlabs), dtype = np.float64)
        position_angle      = np.zeros((nlabs), dtype = np.float64)
        r_max               = np.zeros((nlabs), dtype = np.float64)
        semi_major_angles   = np.zeros((nlabs), dtype = np.float64)
        flux                = np.zeros((nlabs), dtype = np.float64)
        parent_id           = np.zeros((nlabs), dtype = np.int32)
        x_limits            = np.zeros((nlabs, 2), dtype = np.int32)
        y_limits            = np.zeros((nlabs, 2), dtype = np.int32)

        print("Computing cov matrices")
        start_sub = time.monotonic()
        
        # C: per-object pixel moments, flux, area, parent_id and bounding boxes
        self.__compute_covs(self.data, segmentation_map, self.mask, nrows, ncols, nlabs, coms, covariance_matrices,
                            flux, areas, r_max, parent_id, x_limits, y_limits)

        stop_sub = time.monotonic()
        print(f"Time: {stop_sub - start_sub: .2f}")

        print("Computing cov features")
        start_sub = time.monotonic()
        
        # Derive shape parameters (ellipticity, PA, axis ratio) from the covariance eigendecomposition
        eigenvals = np.zeros((nlabs,2), dtype = np.float64)
        eigenvecs = np.zeros((nlabs,4), dtype = np.float64)
        self.__compute_eigensystems(covariance_matrices, eigenvals, eigenvecs, nlabs)

        sig_x             = np.sqrt(eigenvals[:, 0])
        sig_y             = np.sqrt(eigenvals[:, 1])
        a_images          = sig_x
        b_images          = sig_y
        ellipticities     = (a_images - b_images) / a_images
        semi_major_angles = np.rad2deg(np.arctan2(eigenvecs[:, 1], eigenvecs[:, 0]))
        position_angle    = ((semi_major_angles + 90.) % 180.) - 90. 
        ra_dec            = wcs.all_pix2world(coms, 0)  # pixel -> world (WCS) coordinates

        stop_sub = time.monotonic()
        print(f"    Time: {stop_sub - start_sub: .2f}")

        # Apply the minimum-area cut and store the filtered per-source measurements
        f = np.where(areas > min_area)

        self.sources_properties = {}
        self.sources_properties["centers_of_mass"]   = coms[f]
        self.sources_properties["world_coord"]       = ra_dec[f]
        self.sources_properties["areas"]             = areas[f]
        self.sources_properties["ellipticities"]     = ellipticities[f]
        self.sources_properties["major_axes"]        = b_images[f]
        self.sources_properties["minor_axes"]        = a_images[f]
        self.sources_properties["position_angle"]    = position_angle[f]
        self.sources_properties["r_max"]             = r_max[f]
        self.sources_properties["semi_major_angles"] = semi_major_angles[f]
        self.sources_properties["x_limits"]          = x_limits[f]
        self.sources_properties["y_limits"]          = y_limits[f]
        self.sources_properties["flux"]              = flux[f]
        self.sources_properties["parent_id"]         = parent_id[f] 

        stop = time.time()
        print(f"\tElapsed time: {stop - start:.2f}s")

        return self.sources_properties

    def writePNG(self,fname, scale = 0.5):
        """Render the segmentation to a PNG, downscaling if above max_allowed_dim."""
        max_allowed_dim = 4000
        if self.nrows * scale > max_allowed_dim or self.ncols * scale > max_allowed_dim:
            dim = max(self.nrows, self.ncols)
            scale = max_allowed_dim / dim

        self.__write_png(fname.encode("utf-8"), self.__datapoints, self.data, self.__clusters.centers.count, self.ncols, self.nrows, np.uint64(self.ncols * scale), np.uint64(self.nrows * scale) )
    def __del__(self):
        # Release the C-side allocations to avoid leaks when the Data object is freed.
        if not self.__datapoints is None:
            self.__freeDatapoints(self.__datapoints, self.n)
        if not self.__clusters is None:
            self.__freeClusters(ct.pointer(self.__clusters))
