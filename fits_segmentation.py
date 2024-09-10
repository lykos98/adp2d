import numpy as np 
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import time


#------------------------------------
#-- SET DADAC PATH ------------------
ADP_PATH = "." 

sys.path.append(ADP_PATH)
import adp2d

#-- Read parameters from stdin

Z = 15
nbkg = 3.5
radius = 10
halo = False
use_asterism_bg_reject = True


fitsPath = "../../robavaria/euclid/Parallel_sci_cut_f200w.fits"
#fitsPath = "/home/francesco/Downloads/mosaic_elgordo_cluster_nircam_f444w_30mas_20230718_drz.fits"

mask_path = "../../robavaria/euclid/binaryMask.fits"

outPath = "Mask_2.fits"
outImg  = "img.png"




#-- Get fits file 

print(f"Opening fits file -> {fitsPath}\n")
ref_file = fits.open(fitsPath)
ref_frame = ref_file[0] 

if use_asterism_bg_reject:

	#-- Set a threshold to filter out bg from the image 
	#-- Using strategy used in paper DBSCAN+DENCLUE

	rf_adj = ref_frame.data + np.abs(ref_frame.data.min())
	submats = []
	nrows = rf_adj.data.shape[0]
	for i in range(10):
	    low = nrows//10 * i 
	    high  = min(nrows//10 * (i+1), nrows)
	    sm = rf_adj[low:high,:]
	    submats.append(sm)

	integratedFlux = [np.sum(m) for m in submats]
	bgMatIdx = np.argmin(integratedFlux)
	bgMat = submats[bgMatIdx]

	bb = bgMat.ravel()
	bb = bb[np.where(bb < 1)]
	#count, vals = np.histogram(bgMat,bins=5000)
	count, vals = np.histogram(bb,bins=5000)
	mode = vals[np.argmax(count)]
	sigma = np.std(bgMat)

	ns = [0.1*i for i in range(5,21)]
	from scipy.stats import skew
	skews = []
	for n in ns:
	    ff = np.where((bgMat > mode - n*sigma) & (bgMat < mode + n*sigma)) 
	    #print(ff)
	    skews.append(skew(bgMat[ff]))

	#print(skews)
	#print(mode)
	nstar = ns[np.argmin(np.abs(skews))]

	ffstar = np.where((bgMat > mode - nstar*sigma) & (bgMat < mode + nstar*sigma)) 
	bgMatN = bgMat[ffstar]
	count, vals = np.histogram(bgMatN,bins=5000)
	mode = vals[np.argmax(count)]
	sigma = np.std(bgMatN)

	bkg = nbkg*sigma + mode 

	print(f"Using as threshold -> {bkg:.2f}\n")

	#-- Filter out image
	rr = ref_frame.data
	rr = rr + np.abs(rr.min())

	mask = rr > bkg 
	mask = mask.astype(np.int32)
else:
	rr = ref_frame.data

	rr = rr + np.abs(rr.min())
	mask = fits.open(mask_path)[0].data

#-- Let DadaC cook 
print(f"Warming up the clustering engine\n")
data = adp2d.Data(rr)

data.computeDensityFromImg(rr.astype(np.float64),mask,radius)
data.computeClusteringADP(Z, halo = halo)

clusterLabels = data.getClusterAssignment()

print("Exporting to png")
t1 = time.monotonic()
data.writePNG(outImg)

t2 = time.monotonic()
print(f" -> took {(t2 - t1): .3f}s")


# %%
print("Exporting to fits")
t1 = time.monotonic()

hdu = fits.PrimaryHDU(clusterLabels.reshape(ref_frame.data.shape))
hdul = fits.HDUList([hdu])
hdul.writeto(outPath, overwrite=True)

t2 = time.monotonic()
print(f" -> took {(t2 - t1): .3f}s")

# %%



