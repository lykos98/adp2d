import ctypes as ct
import numpy as np  
import os
import time
from scipy.ndimage import center_of_mass
from concurrent.futures import ThreadPoolExecutor
import numba


ctFloatType = ct.c_double
ctIdxType = ct.c_uint64

class HeapNode(ct.Structure):
    _fields_ = [
        ("value", ctFloatType),
        ("array_idx", ctIdxType)
    ]


class Heap(ct.Structure):
    _fields_ = [
        ("N", ctIdxType),
        ("count", ctIdxType),
        ("data", ct.POINTER(HeapNode))
    ]

class luDynamicArray(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ctIdxType)),
        ("size", ctIdxType),
        ("count", ctIdxType)
    ]

class DatapointInfo(ct.Structure):
    _fields_ = [
        ("g", ctFloatType),
        ("ngbh", Heap),
        ("array_idx", ctIdxType),
        ("log_rho", ctFloatType),
        ("log_rho_c", ctFloatType),
        ("log_rho_err", ctFloatType),
        ("kstar", ctIdxType),
        ("is_center", ct.c_int),
        ("cluster_idx", ct.c_int)
    ]

    def __repr__(self):
        return f"{{ g: {self.g}, ngbh: {self.ngbh}, array_idx: {self.array_idx}, log_rho: {self.log_rho}, log_rho_c: {self.log_rho_c}, log_rho_err: {self.log_rho_err}, kstar: {self.kstar}, isCenter: {self.is_center}, clusterIdx: {self.cluster_idx} }}"

class Border_t(ct.Structure):
    _fields_ = [
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class SparseBorder_t(ct.Structure):
    _fields_ = [
        ("i", ctIdxType),
        ("j", ctIdxType),
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class AdjList(ct.Structure):
    _fields_ = [
        ("count", ctIdxType),
        ("size", ctIdxType),
        ("data", ct.POINTER(SparseBorder_t))
    ]

class Clusters(ct.Structure):
    _fields_ = [
        ("UseSparseBorders", ct.c_int),
        ("SparseBorders", ct.POINTER(AdjList)),
        ("centers", luDynamicArray),
        ("borders", ct.POINTER(ct.POINTER(Border_t))),
        ("__borders_data", ct.POINTER(Border_t)),
        ("n", ctIdxType)
    ]



@numba.njit
def _compute_coms_and_covs(segmentation_map: np.array, label_idx: np.array, coms: np.array, areas: np.array, covariance_matrices: np.array):
    """
        segmentation_map: np.array( nrows, ncols )
        label_idx: map from label to the correct index
        coms: np.array( nlabels, 2 )
        areas: np.array( nlabels )
        covariance_matrices: np.array( nlabels, 2, 2)
    """
    
    xrange  = segmentation_map.shape[1]
    yrange  = segmentation_map.shape[0]
    nlabels = coms.shape[0] 

    for xx in range(xrange):
        for yy in range(yrange):
            lab = segmentation_map[yy,xx] 
            if lab != -1:
                lab_idx = label_idx[lab]
                coms[lab_idx,0] += float(xx)
                coms[lab_idx,1] += float(yy)
                areas[lab_idx]  += 1.
                #numba.cuda.atomic.add(coms,(lab,0),xx)
                #numba.cuda.atomic.add(coms,(lab,1),yy)
                #numba.cuda.atomic.add(areas,lab,xx)

    for lab_idx in range(nlabels):
        if areas[lab_idx] > 1:
            coms[lab_idx] = coms[lab_idx]/areas[lab_idx]

    for xx in range(xrange):
        for yy in range(yrange):
            lab = segmentation_map[yy,xx] 
            if lab != -1:
                lab_idx = label_idx[lab]
                
                x_n = float(xx) - coms[lab_idx,0]
                y_n = float(yy) - coms[lab_idx,1]
                covariance_matrices[lab_idx, 0, 0] += x_n ** 2
                covariance_matrices[lab_idx, 0, 1] += x_n * y_n
                covariance_matrices[lab_idx, 1, 0] += x_n * y_n
                covariance_matrices[lab_idx, 1, 1] += y_n ** 2
                
    for lab_idx in range(nlabels):
        if areas[lab_idx] > 1:
            covariance_matrices[lab_idx,:,:] = covariance_matrices[lab_idx,:,:]/(areas[lab_idx] - 1)

@numba.njit(parallel = True)
def _compute_cov_properties(covariance_matrices, areas, ellipticities, b_images, a_images, semi_major_angles):
    nlabels = covariance_matrices.shape[0]
    for i in numba.prange(nlabels):
        if areas[i] > 1:
            cov = covariance_matrices[i]
            eigenvectors, eigenvalues, V = np.linalg.svd(cov, full_matrices=True)
            ind                          = np.argsort(eigenvalues)[::-1]
            eigenvalues                  = eigenvalues[ind]
            eigenvectors                 = eigenvectors[:, ind]
        
            # Asterism like ellipticity
            sig_x       = np.sqrt(eigenvalues[0])
            sig_y       = np.sqrt(eigenvalues[1])
            a_image     = np.sqrt(sig_x * sig_x + sig_y * sig_y)
            b_image     = sig_y / sig_x * a_image
            ellipticity = (a_image - b_image) / a_image
    
            # Asterism like PA
            semi_major_angle = np.rad2deg(np.arctan2(eigenvectors[0][1], eigenvectors[0][0]))
    
            a_images[i] = a_image
            b_images[i] = b_image
            ellipticities[i] = ellipticity
            semi_major_angles[i] = semi_major_angle

class Data():
    def __init__(self, data : np.array):
        """Data object for dadaC library

        Args:
            data (np.array): 2d array to use in searching for clusters 

        Raises:
            TypeError: Raises TypeError if a type different from a matrix is passed 
        """
        path = os.path.join(os.path.dirname(__file__), "bin/libadp2d.so")
        self.lib = ct.CDLL(path)

        global ctFloatType, ctIdxType
        s = self.lib.FloatAndUintSize()

        self.__useFloat32 = s < 1
        self.__useInt32   = (s == 0) or (s == 2)

        #initialize class
        self.state = {
                    "ngbh" : False,
                    "id" : False,
                    "density" : False,
                    "clustering" : False,
                    "useFloat32": self.__useFloat32,
                    "useInt32": self.__useInt32,
                    "useSparse": None,
                    "computeHalo": None 
                    }
        if self.__useFloat32:
            self.data = np.ascontigousarray(data.astype(np.float32))
            ctFloatType = ct.c_float
        else:
            self.data = np.ascontiguousarray(data.astype(np.float64))

        if self.__useInt32:
            ctIdxType = ct.c_uint32

        #retrieve function pointers form .so file

        #Datapoint_info* computeDensityFromImg(FLOAT_TYPE* vals, int* mask, size_t nrows, size_t ncols)
        self.__computeDensityFromImg = self.lib.computeDensityFromImg
        self.__computeDensityFromImg.argtypes = [   np.ctypeslib.ndpointer(ctFloatType), 
                                                    np.ctypeslib.ndpointer(np.int32),
                                                    ct.c_int32,
                                                    ct.c_int32,
                                                    ct.c_int32]
        self.__computeDensityFromImg.restype  = ct.POINTER(DatapointInfo)

        #void setRhoErrK(Datapoint_info* points, FLOAT_TYPE* rho, FLOAT_TYPE* rhoErr, idx_t* k, size_t n)
        self.__setRhoErrK = self.lib.setRhoErrK
        self.__setRhoErrK.argtypes = [  ct.POINTER(DatapointInfo),    
                                        np.ctypeslib.ndpointer(ctFloatType), 
                                        np.ctypeslib.ndpointer(ctFloatType),
                                        np.ctypeslib.ndpointer(ctIdxType),
                                        ct.c_uint64 ]


        self.__computeCorrection = self.lib.computeCorrection
        self.__computeCorrection.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(np.int32), ctIdxType, ct.c_double]

        self.__H1 = self.lib.Heuristic1
        #Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols);
        self.__H1.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(np.int32), ct.c_int32, ct.c_int32]
        self.__H1.restype = Clusters

        self.__ClustersAllocate = self.lib.Clusters_allocate
        self.__ClustersAllocate.argtypes = [ct.POINTER(Clusters), ct.c_int]

        self.__H2 = self.lib.Heuristic2
        #void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols);
        self.__H2.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), 
                              np.ctypeslib.ndpointer(np.int32), ct.c_uint64, ct.c_uint64]

        self.__H3 = self.lib.Heuristic3
        self.__H3.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), ct.c_double, ct.c_int]
        
        self.__freeDatapoints = self.lib.freeDatapointArray
        self.__freeDatapoints.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]

        self.__freeClusters = self.lib.Clusters_free
        self.__freeClusters.argtypes = [ct.POINTER(Clusters)]

        self._export_cluster_assignment = self.lib.export_cluster_assignment
        self._export_cluster_assignment.argtypes = [
            ct.POINTER(DatapointInfo),
            np.ctypeslib.ndpointer(np.int32),
            ct.c_uint64,
        ]


        #void tiny_colorize(
        #        const char* fname, 
        #        Datapoint_info* dp, 
        #        FLOAT_TYPE* data, 
        #        uint32_t n_clusters, 
        #        uint32_t og_width, 
        #        uint32_t og_height, 
        #        uint32_t target_width,
        #        uint32_t target_height)
        #{
        self.__write_png = self.lib.tiny_colorize
        self.__write_png.argtypes = [ct.c_char_p,
                                     ct.POINTER(DatapointInfo),
                                     np.ctypeslib.ndpointer(ctFloatType),
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,
                                     ct.c_uint32,]


        if len(self.data.shape) != 2:
            raise TypeError("Please provide a 2d numpy array")


        self.__datapoints     = None
        self.__clusters       = None
        self.n              = self.data.shape[0]
        self.dims           = self.data.shape[1]
        self.k              = None
        self.id             = None

        self.clusterAssignment = None
        self.neighbors         = None
        self.borders           = None
        self.density           = None
        self.densityError      = None

    def computeDensityFromImg(self,img, mask = None, r = 15):
        if mask is None:
            mask = np.ones_like(img, dtype = np.int32)
        self.n = np.prod(img.shape)
        mask = mask.astype(np.int32)
        self.nrows, self.ncols = img.shape
        self.img = img
        self.mask = mask
        self.__datapoints = self.__computeDensityFromImg(img,mask, img.shape[0], img.shape[1], r)
        self.state["density"] = True


    def computeClusteringADP(self,Z : float, halo = True, useSparse = "auto"):

        """Compute clustering via the Advanced Density Peak method

        Args:
            Z (float): Z value for the method 
            useSparse (str): optional [``auto``,`yes`,`no`], use sparse implementation of border storage between clusters. Memory usage for big datsets is significant. 

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

        self.state["computeHalo"] = halo 
        self.Z = Z
        self.n = np.prod(self.img.shape)
        self.__computeCorrection(self.__datapoints, self.mask, self.n, self.Z)
        self.__clusters = self.__H1(self.__datapoints, self.mask, self.nrows, self.ncols)
        self.__ClustersAllocate(ct.pointer(self.__clusters), 1)
        self.__H2(ct.pointer(self.__clusters), self.__datapoints, self.mask, self.nrows, self.ncols)
        self.__H3(ct.pointer(self.__clusters), self.__datapoints, self.Z, 1 if halo else 0 )
        self.state["clustering"] = True
        self.clusterAssignment = None

    def getClusterAssignment(self) -> list:


        """Retrieve cluster assignment

        Raises:
            ValueError: Raises error if clustering is not computed, use `Data.computeClusteringADP(Z)` 

        Returns:
            List of cluster labels
            
        """
        #if self.clusterAssignment is None:
        #    if self.state["clustering"]:
        #        self.clusterAssignment = np.array([int(self.__datapoints[j].cluster_idx) for j in range(self.n)])
        #        return self.clusterAssignment
        #    else:
        #        raise ValueError("Clustering is not computed yet")
        #else:
        #    return self.clusterAssignment
        self.clusterAssignment = np.ascontiguousarray(np.zeros((self.nrows, self.ncols), np.int32))
        self._export_cluster_assignment(self.__datapoints, self.clusterAssignment, self.n)
        #self.clusterAssignment = np.array([int(self.__datapoints[j].cluster_idx) for j in range(self.n)], dtype = np.int32)
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
        self.__setRhoErrK(self.__datapoints, rho, err, k, self.n)
        self.state["density"] = True
        return



    def computeSourcesProperties(self, min_area = 10):
        start = time.time()
        print("Computing sources properties")
        self.getClusterAssignment()
        # Get unique labels in dadaC segmentation map
        unique_labels = np.unique(self.clusterAssignment)
        segmentation_map = self.clusterAssignment.reshape((self.nrows, self.ncols)).astype(np.int32)
        


        label_idx = np.array([0 for _ in range(max(unique_labels + 1))], dtype = np.int32)

        for i,l in enumerate(unique_labels[1:]):
            label_idx[l] = i


        nlabs = len(unique_labels) - 1

        coms = np.zeros((nlabs,2), dtype = np.float32)
        covariance_matrices = np.zeros((nlabs,2,2), dtype = np.float32)
        areas         = np.zeros((nlabs), dtype = np.float32)
        a_images      = np.zeros((nlabs), dtype = np.float32)
        b_images      = np.zeros((nlabs), dtype = np.float32)
        ellipticities = np.zeros((nlabs), dtype = np.float32)
        semi_major_angles = np.zeros((nlabs), dtype = np.float32)

        # Use ThreadPoolExecutor to parallelize the computation

        _compute_coms_and_covs(segmentation_map, label_idx, coms, areas, covariance_matrices)
        _compute_cov_properties(covariance_matrices, areas, ellipticities, b_images, a_images, semi_major_angles)
        #print(coms)
        f = np.where(areas > min_area)
        self.sources_properties = {}
        self.sources_properties["centers_of_mass"]   = coms[f].T
        self.sources_properties["areas"]             = areas[f]
        self.sources_properties["ellipticities"]     = ellipticities[f]
        self.sources_properties["major_axes"]        = b_images[f]
        self.sources_properties["minor_axes"]        = a_images[f]
        self.sources_properties["position_angles"]   = semi_major_angles[f]

        stop = time.time()
        print(f"\tElapsed time: {stop - start:.2f}s")

        return self.sources_properties



    def writePNG(self,fname, scale = 0.5):
        max_allowed_dim = 4000
        if self.nrows * scale > max_allowed_dim or self.ncols * scale > max_allowed_dim:
            dim = max(self.nrows, self.ncols)
            scale = max_allowed_dim / dim

        self.__write_png(fname.encode("utf-8"), self.__datapoints, self.data, self.__clusters.centers.count, self.ncols, self.nrows, np.uint64(self.ncols * scale), np.uint64(self.nrows * scale) )

    def __del__(self):
        if not self.__datapoints is None:
            self.__freeDatapoints(self.__datapoints, self.n)
        if not self.__clusters is None:
            self.__freeClusters(ct.pointer(self.__clusters))
