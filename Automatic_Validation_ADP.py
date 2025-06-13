#!/home/francesco/py_envs/nope/bin/python3
import os
import sys
import gc
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D
from sklearn.metrics import confusion_matrix, accuracy_score, jaccard_score
from matplotlib.colors import ListedColormap


# ADP Path
ADP_PATH = "/euclid_data/mlepinzan/ADP_vs_Asterism_Q_1/ADP_RUN/adp2d" 

sys.path.append(ADP_PATH)
import adp2d

# Convert S/N to Flux
def flux_from_snr(snr, area, background):
    # Solve: SNR = F / sqrt(F + N * B)
    # Positive solution:
    return (snr**2 / 2) * (1 + np.sqrt(1 + (4 * area * background) / snr**2))

# Segmentation map creation based on criteria according to which a pixels goes to the max flux source in the overlapping region
def build_segmentation_map_from_snr(source_images, background_level, snr_threshold, min_pixels=10):
    """
    Build a segmentation map from source images using per-pixel S/N and max assignment.
    
    Parameters
    ----------
    source_images : list of 2D np.ndarray
        Flux images for each source.
    background_level : float
        Background flux level per pixel (used for S/N estimation).
    snr_threshold : float
        Minimum S/N to consider a pixel as part of a source.
    min_pixels : int
        Minimum number of pixels to retain a labeled source.
    
    Returns
    -------
    segmap : 2D np.ndarray (int)
        Segmentation map (0 = background, 1, 2, ... = source labels).
    snr_mask : 2D np.ndarray (bool)
        Binary mask of pixels above threshold in any source.
    """
    if len(source_images) == 0:
        raise ValueError("source_images must contain at least one image.")

    # Ensure all source images have the same shape
    shape = source_images[0].shape
    if not all(img.shape == shape for img in source_images):
        raise ValueError("All source images must have the same shape.")
    
    # Compute per-pixel S/N maps
    snr_maps = [img / np.sqrt(img + background_level) for img in source_images]

    # Build combined mask (pixels where any source has S/N > threshold)
    snr_mask = np.any([s > snr_threshold for s in snr_maps], axis=0)

    # Stack all S/N maps 
    snr_stack = np.stack(snr_maps, axis=0)  # shape: (N_sources, H, W)

    # Find index of the maximum S/N per pixel (among sources)
    dominant_idx = np.argmax(snr_stack, axis=0)  # values in [0, N_sources-1]

    # Initial segmentation map
    segmap           = np.zeros_like(snr_mask, dtype=int)
    segmap[snr_mask] = dominant_idx[snr_mask] + 1  # 1-based indexing

    # Remove small areas
    final_segmap = np.zeros_like(segmap)
    label_id     = 1
    for src_label in range(1, len(snr_maps) + 1):
        mask = segmap == src_label
        if np.sum(mask) >= min_pixels:
            final_segmap[mask] = label_id
            label_id += 1

    return final_segmap, snr_mask

def build_segmentation_map_with_importance(source_images, background_level, n=0.1, min_pixels=10):
    """
    Builds a segmentation map from multiple source images using importance scoring.

    Parameters
    ----------
    source_images : list of 2D np.ndarray
        Each element is the flux image of one source (same shape).
    background_level : float
        Mean background level per pixel (Poisson mean).
    n : float
        Threshold coefficient, e.g., 0.1 → t = bg + n * sqrt(bg)
    min_pixels : int
        Minimum number of pixels required to keep a labeled source

    Returns
    -------
    segmap : 2D np.ndarray (int)
        Final segmentation map (0 = background, 1, 2, ..., N = source labels)
    """
    num_sources = len(source_images)
    h, w        = source_images[0].shape

    # Stack sources: shape (N_sources, H, W)
    sources = np.stack(source_images, axis=0)

    # Compute total flux for each source (Fs)
    Fs = sources.sum(axis=(1, 2)) + 1e-12  # prevent divide-by-zero

    # Compute detectability threshold t
    t = background_level + n * np.sqrt(background_level)

    # Compute importance score: (fs,p)^2 / Fs
    score = (sources ** 2) / Fs[:, None, None]

    # Mask out undetectable pixels (below flux threshold)
    detectable = sources >= t
    score[~detectable]    = 0
    detection_mask        = np.any(detectable, axis=0)

    # Assign each pixel to the source with the highest score
    best_idx   = np.argmax(score, axis=0)
    best_score = np.max(score, axis=0)
    assigned   = best_score > 0

    segmap           = np.zeros((h, w), dtype=int)
    segmap[assigned] = best_idx[assigned] + 1  # use 1-based labels

    # Filter out small areas
    final_segmap = np.zeros_like(segmap)
    new_label = 1
    for label_val in range(1, num_sources + 1):
        mask = segmap == label_val
        if np.sum(mask) >= min_pixels:
            final_segmap[mask] = new_label
            new_label += 1

    return final_segmap, detection_mask


def automatic_run(snrs, positions_list, sigmas, background_levels, gt_type="importance", output_root="runs", n=0.1, snr_threshold=3, min_pixels=10):
    """
    Run a grid of simulations over different parameter combinations.

    Parameters
    ----------
    snrs : list of float
        Signal-to-noise ratios for sources.
    positions_list : list of list of tuples
        Each sublist defines the (x, y) positions of sources for a run.
    sigmas : list of float
        Standard deviations (sizes) of sources.
    background_levels : list of float
        Background levels (mean of Poisson noise).
    gt_type : str
        Ground truth type: 'importance' or 'snr'.
    output_root : str
        Directory to save results.
    snr_threshold : float
        Threshold for S/N based ground truth segmentation.
    min_pixels : int
        Minimum number of pixels to keep a labeled region.

    Returns
    -------
    run_results : list of dict
        List of results (metrics + run tag) for each combination.
    """
    
    # Set image size and create output root directory if doesn't exits. Inside this each directory will contain a sub directory indicating the parameter combinations of the run
    y, x = np.mgrid[0:500, 0:500]
    os.makedirs(output_root, exist_ok=True)

    run_results = []

    # Loop over the parameter input
    for snr in snrs:
        for bg in background_levels:
            for sigma in sigmas:
                for positions in positions_list:
                    
                    # Directory and file naming tag
                    pos_str = "_".join([f"x{x}_y{y}" for (x, y) in positions])
                    tag = f"snr{snr}_bg{bg}_sigma{sigma}_{pos_str}"
                    run_dir = os.path.join(output_root, tag)
                    os.makedirs(run_dir, exist_ok=True)

                    # Compute source properties
                    source_area = 2 * np.pi * sigma**2
                    fluxes      = [flux_from_snr(snr, source_area, bg)] * len(positions)
                    amplitudes  = [f / source_area for f in fluxes]
                    
                    # Generate individual source images
                    images = []
                    for amp, (x0, y0) in zip(amplitudes, positions):
                        g = Gaussian2D(amplitude=amp, x_mean=x0, y_mean=y0, x_stddev=sigma, y_stddev=sigma)
                        images.append(g(x, y))

                    # Combine sources + background
                    poisson_background = np.random.poisson(lam=bg, size=(500, 500))
                    image              = poisson_background + sum(images)

                    # Create segmentation map (ground truth)
                    if gt_type == "importance":
                        segmap, mask = build_segmentation_map_with_importance(images, background_level=bg, n=0.1, min_pixels=10)
                    else:
                        segmap, mask = build_segmentation_map_from_snr(images, background_level=bg, snr_threshold=3, min_pixels=10)
                    
                    # Save image, ground truth segmentation map and mask
                    fits.writeto(f"{run_dir}/Simulated_image.fits", image, overwrite=True)
                    fits.writeto(f"{run_dir}/Ground_truth_segmap.fits", segmap.astype(np.int32), overwrite=True)
                    fits.writeto(f"{run_dir}/ADP_mask.fits", mask.astype(np.int32), overwrite=True)

                    # RUN ADP
                    # ADP parameters
                    Z      = 3.5
                    nbkg   = 3.5
                    radius = 5
                    halo   = False

                    # File paths
                    image_path       = f"{run_dir}/Simulated_image.fits"
                    mask_path        = f"{run_dir}/ADP_mask.fits"
                    adp_segmap_path  = f"{run_dir}/ADP_segmentation_map.fits"

                    # Load image and mask
                    image_data = np.ascontiguousarray(fits.open(image_path)[0].data)
                    mask_data  = np.ascontiguousarray(fits.open(mask_path)[0].data)
                                        
                    # Open simulated image
                    print(f"Opening fits file -> {image_path}\n")

                    # Let ADP cook 
                    print(f"Warming up the clustering engine\n")
                    
                    data = adp2d.Data(image_data)

                    data.computeDensityFromImg(image_data.astype(np.float64), mask_data,radius)

                    data.computeClusteringADP(Z, halo = halo)

                    clusterLabels = data.getClusterAssignment()

                    print("Exporting ADP segmentation map to fits")
                    t1 = time.monotonic()

                    hdu  = fits.PrimaryHDU(clusterLabels.reshape(image_data.data.shape))
                    hdul = fits.HDUList([hdu])
                    hdul.writeto(adp_segmap_path, overwrite=True)

                    t2 = time.monotonic()
                    print(f" -> took {(t2 - t1): .3f}s")

                    # Cleanup
                    del image_data, mask_data, data, clusterLabels
                    gc.collect()

                    # Compare prediction vs ground truth
                    # Load ADP result
                    adp_segmap_data = fits.getdata(adp_segmap_path) + 1  # Shift -1 → 0, 0+ → 1+
                    
                    # Evaluation mask
                    mask_eval = (segmap > 0) | (adp_segmap_data > 0)
                    gt_flat   = segmap[mask_eval].flatten()
                    pred_flat = adp_segmap_data[mask_eval].flatten()

                    # Compute metrics
                    conf_mat = confusion_matrix(gt_flat, pred_flat)
                    accuracy = accuracy_score(gt_flat, pred_flat)
                    labels   = np.unique(np.concatenate((gt_flat, pred_flat)))
                    ious     = jaccard_score(gt_flat, pred_flat, labels=labels, average=None)

                    metrics = {
                        "labels": labels.tolist(),
                        "confusion_matrix": conf_mat.tolist(),
                        "pixel_accuracy": float(accuracy),
                        "iou_scores_per_class": ious.tolist(),
                    }

                    with open(f"{run_dir}/metrics.json", "w") as f:
                        json.dump(metrics, f, indent=2)
                    with open(f"{run_dir}/params.json", "w") as f:
                        json.dump({"snr": snr, "positions": positions, "sigma": sigma, "background": bg, "gt_type": gt_type}, f, indent=2)

                    run_results.append({"run": tag, **metrics})

    return run_results


if __name__ == "__main__":
    
    # Example configuration
    snrs              = [150,150]
    sigmas            = [15]
    background_levels = [10]
    positions_list    = [[(220, 240), (270, 240)]]
    
    # positions_list    = [[(220, 240), (270, 240)], [(240, 240), (245, 280)]]

    # Launch the run
    results = automatic_run(
              snrs=snrs,
              positions_list=positions_list,
              sigmas=sigmas,
              background_levels=background_levels,
              gt_type="importance",  # or "snr"
              output_root="runs",
              n=0.1, 
              snr_threshold=3, 
              min_pixels=10
        )

    # Optional: Print summary
    for r in results:
        print(f"[{r['run']}]")
        print(f"  Accuracy: {r['pixel_accuracy']:.3f}")
        print(f"  IoUs: {r['iou_scores_per_class']}")
        print(f"  Confusion Matrix:\n{np.array(r['confusion_matrix'])}")
