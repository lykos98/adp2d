import time
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.ndimage import gaussian_filter
import json

ADP_PATH  = "/euclid_data/mlepinzan/ADP_last_version"
sys.path.append(ADP_PATH)

import adp2d

# ============================================================
# CONFIGURATION
# ============================================================

config_path = sys.argv[1]

with open(config_path, "r") as f:
    config = json.load(f)

OUTPUT_DIR = Path(config["output"]["directory"])
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# First ADP run: full detection segmentation map
# ------------------------------------------------------------

RUN1 = {
    "R": 2,
    "Z": 1.2,
    "density_algorithm": "MEAN",
    "param": None,
    "border": "percentile",
    "border_perc": 0.8,
}


# ------------------------------------------------------------
# Second ADP run: moderate over-deblending
# ------------------------------------------------------------

RUN2 = {
    "R": 2,
    "Z": 1.2,
    "density_algorithm": "GAUSSIAN",
    "param": 2,
    "border": "maxg",
}

# ------------------------------------------------------------
# Third ADP run: extreme over-deblending
# ------------------------------------------------------------

RUN3 = {
    "R": 3,
    "Z": 2.0,
    "density_algorithm": "GAUSSIAN",
    "param": 2,
    "border": "maxg",
}


# Percentiles used to classify problematic parents
MODERATE_PERCENTILE = 95
EXTREME_PERCENTILE  = 99


# ============================================================
# ADP RUN
# ============================================================

def run_adp(
    image_data,
    detection_segmap,
    header,
    output_dir,
    run_name,
    R,
    Z,
    density_algorithm,
    param,
    border,
    border_perc=None,
):
    """
    Perform one complete ADP run.

    Returns
    -------
    segmap : ndarray
        ADP segmentation map.

    cat : astropy.table.Table
        ADP source catalogue.
    """

    print("\n" + "=" * 60)
    print(f"Starting {run_name}")
    print("=" * 60)

    print(f"Density estimator : {density_algorithm}")
    print(f"R                 : {R}")
    print(f"param             : {param}")
    print(f"Z                 : {Z}")
    print(f"Border            : {border}")

    if border == "percentile":
        print(f"Border percentile : {border_perc}")

    t_run = time.perf_counter()

    # --------------------------------------------------------
    # Create ADP object
    # --------------------------------------------------------

    data = adp2d.Data(image_data)

    # --------------------------------------------------------
    # Density estimation
    # --------------------------------------------------------

    t0 = time.perf_counter()

    data.computeDensityFromImg(
        image_data,
        detection_segmap,
        R,
        algorithm=density_algorithm,
        param=param,
        use_log=True,
        use_adaptive_radius=True,
    )

    density_time = time.perf_counter() - t0

    # --------------------------------------------------------
    # Clustering
    # --------------------------------------------------------

    clustering_kwargs = {
        "Z": Z,
        "halo": False,
        "splitPerThread": True,
        "border": border,
    }

    if border == "percentile":

        if border_perc is None:
            raise ValueError(
                "border_perc must be provided when"
                "border='percentile'."
            )

        clustering_kwargs["border_perc"] = border_perc

    t0 = time.perf_counter()

    data.computeClusteringADP(
        **clustering_kwargs
    )

    clustering_time = time.perf_counter() - t0

    # --------------------------------------------------------
    # Create ADP segmentation map
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cluster_labels = data.getClusterAssignment() + 1

    segmap = cluster_labels.reshape(
        image_data.shape
    )

    segmap_path = (
        output_dir /
        f"{run_name}_segmentation.fits"
    )

    fits.writeto(
        segmap_path,
        segmap,
        overwrite=True,
        header=header,
    )

    segmentation_time = time.perf_counter() - t0


    # --------------------------------------------------------
    # Source properties
    # --------------------------------------------------------

    t0 = time.perf_counter()

    data.computeSourcesProperties(
        header=header
    )

    properties_time = time.perf_counter() - t0


    # --------------------------------------------------------
    # Catalogue
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cat_path = (
        output_dir /
        f"{run_name}_catalogue.fits"
    )

    data.exportSourcesCatalogue(
        str(cat_path)
    )

    cat = Table.read(cat_path)

    catalogue_time = time.perf_counter() - t0

    total_time = time.perf_counter() - t_run


    print(f"\n{run_name} completed")
    print(f"Sources          : {len(cat)}")
    print(f"Density          : {density_time:.2f} s")
    print(f"Clustering       : {clustering_time:.2f} s")
    print(f"Segmentation     : {segmentation_time:.2f} s")
    print(f"Properties       : {properties_time:.2f} s")
    print(f"Catalogue        : {catalogue_time:.2f} s")
    print(f"Total            : {total_time:.2f} s")

    return segmap, cat, total_time


# ============================================================
# CLASSIFY PARENTS FROM FIRST RUN
# ============================================================

def classify_parents(
    catalogue,
    moderate_percentile=95,
    extreme_percentile=99,
):
    """
    Classify over-deblended parents using the child-count
    distribution of the first ADP run.
    """

    parent_ids = np.asarray(
        catalogue["PARENT_ID"],
        dtype=np.int64,
    )

    # PARENT_ID = -1 corresponds to sources that were not split
    parent_ids = parent_ids[parent_ids != -1]

    ids, counts = np.unique(
        parent_ids,
        return_counts=True,
    )

    p_moderate = np.percentile(
        counts,
        moderate_percentile,
    )

    p_extreme = np.percentile(
        counts,
        extreme_percentile,
    )

    # Same logic that gave us 13 and 69:
    #
    # P95 = 12 -> threshold = 13
    # P99 = 68 -> threshold = 69

    moderate_threshold = (
        int(np.floor(p_moderate)) # Ricordati di aggiungere il + 1 dopo il debugging
    )

    extreme_threshold = (
        int(np.floor(p_extreme)) # Ricordati di aggiungere il + 1 dopo il debugging
    )

    moderate_ids = ids[
        (counts >= moderate_threshold)
        & (counts < extreme_threshold)
    ]

    extreme_ids = ids[
        counts >= extreme_threshold
    ]

    print("\n" + "=" * 60)
    print("FIRST-RUN CLASSIFICATION")
    print("=" * 60)

    print(
        f"P{moderate_percentile} = "
        f"{p_moderate:.2f}"
    )

    print(
        f"P{extreme_percentile} = "
        f"{p_extreme:.2f}"
    )

    print(
        f"Moderate threshold: "
        f">= {moderate_threshold} and "
        f"< {extreme_threshold}"
    )

    print(
        f"Extreme threshold: "
        f">= {extreme_threshold}"
    )

    print(
        f"Moderate parents: {len(moderate_ids)}"
    )

    print(
        f"Extreme parents: {len(extreme_ids)}"
    )

    return {
        "ids": ids,
        "counts": counts,
        "moderate_ids": moderate_ids,
        "extreme_ids": extreme_ids,
        "moderate_threshold": moderate_threshold,
        "extreme_threshold": extreme_threshold,
        "p_moderate": p_moderate,
        "p_extreme": p_extreme,
    }


# ============================================================
# CREATE DETECTION SUBSET
# ============================================================

# This is working but np.isin is quite slow on a large image. Change directly in the run_iterative_adp()
def create_detection_subset(
    detection_segmap,
    selected_ids,
):
    """
    Keep only selected original detection patches.
    """

    return np.where(
        np.isin(
            detection_segmap,
            selected_ids,
        ),
        detection_segmap,
        0,
    ).astype(
        detection_segmap.dtype
    )


# ============================================================
# MERGE THE THREE RUNS
# ============================================================

def merge_runs(
    segmentation_map_detection,
    segmap_first,
    cat_first,
    segmap_second,
    cat_second,
    moderate_mask,
    segmap_third,
    cat_third,
    extreme_mask,
):
    """
    Replace first-run sources in moderate/extreme detection
    patches with the corresponding second/third-run sources.
    """

    print("\n" + "=" * 60)
    print("MERGING RUNS")
    print("=" * 60)

    t_merge = time.perf_counter()

    replacement_runs = [
        (
            "Moderate",
            cat_second,
            segmap_second,
            moderate_mask,
        ),
        (
            "Extreme",
            cat_third,
            segmap_third,
            extreme_mask,
        ),
    ]

    # Start from the complete first run
    final_cat    = cat_first.copy()
    final_segmap = segmap_first.copy()

    first_ids = np.asarray(
        cat_first["SOURCE_ID"],
        dtype=np.int64,
    )

    if len(first_ids) != len(
        np.unique(first_ids)
    ):
        raise ValueError(
            "First-run catalogue contains "
            "duplicate SOURCE_ID values."
        )

    # New source IDs begin after the maximum ID
    # used by the first run
    next_id = int(first_ids.max()) + 1

    for (
        name,
        cat_new,
        segmap_new,
        replace_mask,
    ) in replacement_runs:

        print(
            f"\nProcessing {name} group ..."
        )

        # ----------------------------------------------------
        # Find first-run sources that must disappear
        # ----------------------------------------------------

        old_source_ids = np.unique(
            final_segmap[replace_mask]
        )

        old_source_ids = old_source_ids[
            old_source_ids > 0
        ]


        # ----------------------------------------------------
        # Remove them from catalogue
        # ----------------------------------------------------

        current_ids = np.asarray(
            final_cat["SOURCE_ID"],
            dtype=np.int64,
        )

        keep = ~np.isin(
            current_ids,
            old_source_ids,
        )

        final_cat = final_cat[keep]


        # ----------------------------------------------------
        # IDs produced by this new ADP run
        # ----------------------------------------------------

        new_ids = np.asarray(
            cat_new["SOURCE_ID"],
            dtype=np.int64,
        )

        if len(new_ids) != len(
            np.unique(new_ids)
        ):
            raise ValueError(
                f"{name} catalogue contains "
                "duplicate SOURCE_ID values."
            )


        # ----------------------------------------------------
        # Assign globally unique IDs
        # ----------------------------------------------------

        global_ids = np.arange(
            next_id,
            next_id + len(new_ids),
            dtype=np.int64,
        )


        # ----------------------------------------------------
        # Extract replacement segmentation pixels
        # ----------------------------------------------------

        selected_labels = segmap_new[
            replace_mask
        ]

        seg_ids = np.unique(
            selected_labels
        )

        seg_ids = seg_ids[
            seg_ids > 0
        ]


        # Catalogue and segmentation must contain
        # exactly the same source IDs
        if not np.array_equal(
            np.sort(new_ids),
            seg_ids,
        ):
            raise ValueError(
                f"{name}: catalogue and "
                "segmentation IDs differ."
            )


        # ----------------------------------------------------
        # Update catalogue IDs
        # ----------------------------------------------------

        cat_relabelled = cat_new.copy()

        cat_relabelled[
            "SOURCE_ID"
        ] = global_ids

        final_cat = vstack(
            [
                final_cat,
                cat_relabelled,
            ],
            metadata_conflicts="silent",
        )


        # ----------------------------------------------------
        # Vectorized segmentation relabeling
        # ----------------------------------------------------

        order = np.argsort(new_ids)

        sorted_old_ids = new_ids[order]

        sorted_global_ids = (
            global_ids[order]
        )

        mapped = np.zeros(
            selected_labels.shape,
            dtype=np.int64,
        )

        valid = selected_labels > 0

        positions = np.searchsorted(
            sorted_old_ids,
            selected_labels[valid],
        )

        if np.any(
            positions >= len(sorted_old_ids)
        ):
            raise ValueError(
                f"{name}: segmentation contains "
                "unknown source IDs."
            )

        if not np.array_equal(
            sorted_old_ids[positions],
            selected_labels[valid],
        ):
            raise ValueError(
                f"{name}: segmentation/catalogue "
                "ID mismatch."
            )

        mapped[valid] = (
            sorted_global_ids[positions]
        )


        # ----------------------------------------------------
        # Replace pixels
        # ----------------------------------------------------

        final_segmap[
            replace_mask
        ] = mapped


        next_id += len(new_ids)


        # ----------------------------------------------------
        # Immediate duplicate check
        # ----------------------------------------------------

        current_ids = np.asarray(
            final_cat["SOURCE_ID"],
            dtype=np.int64,
        )

        if len(current_ids) != len(
            np.unique(current_ids)
        ):
            raise ValueError(
                f"Duplicate SOURCE_ID after "
                f"{name} merge."
            )


        print(
            f"Replaced original sources: "
            f"{len(old_source_ids)}"
        )

        print(
            f"Inserted new sources: "
            f"{len(new_ids)}"
        )

        print(
            f"Next available SOURCE_ID: "
            f"{next_id}"
        )


    merge_time = (
        time.perf_counter() - t_merge
    )

    print(
        f"\nMerging completed in "
        f"{merge_time:.2f} s"
    )

    return final_segmap, final_cat, merge_time


# ============================================================
# FINAL VALIDATION
# ============================================================

def validate_final_products(
    final_segmap,
    final_cat,
):
    """
    Check catalogue/segmentation consistency.
    """

    print("\n" + "=" * 60)
    print("FINAL VALIDATION")
    print("=" * 60)

    catalogue_ids = np.asarray(
        final_cat["SOURCE_ID"],
        dtype=np.int64,
    )

    segmentation_ids = np.unique(
        final_segmap
    )

    segmentation_ids = (
        segmentation_ids[
            segmentation_ids > 0
        ]
    )

    duplicate_count = (
        len(catalogue_ids)
        - len(np.unique(catalogue_ids))
    )

    cat_missing = np.setdiff1d(
        catalogue_ids,
        segmentation_ids,
    )

    seg_missing = np.setdiff1d(
        segmentation_ids,
        catalogue_ids,
    )

    print(
        f"Final catalogue rows: "
        f"{len(final_cat)}"
    )

    print(
        f"Final segmentation labels: "
        f"{len(segmentation_ids)}"
    )

    print(
        f"Duplicate catalogue IDs: "
        f"{duplicate_count}"
    )

    print(
        f"Catalogue IDs absent from "
        f"segmentation: {len(cat_missing)}"
    )

    print(
        f"Segmentation IDs absent from "
        f"catalogue: {len(seg_missing)}"
    )

    if duplicate_count != 0:
        raise ValueError(
            "Duplicate SOURCE_ID values "
            "in final catalogue."
        )

    if len(cat_missing) != 0:
        raise ValueError(
            "Some catalogue sources are absent "
            "from the segmentation."
        )

    if len(seg_missing) != 0:
        raise ValueError(
            "Some segmentation labels are absent "
            "from the catalogue."
        )

    print("\nAll final checks passed.")


# ============================================================
# MAIN ITERATIVE PIPELINE
# ============================================================

def run_iterative_adp(
    image_data,
    segmentation_map_detection,
    header,
):
    """
    Complete three-stage ADP pipeline.
    """

    pipeline_start = time.perf_counter()

    # ========================================================
    # RUN 1
    # ========================================================

    segmap_first, cat_first, time_run1 = run_adp(
        image_data       = image_data,
        detection_segmap = segmentation_map_detection,
        header           = header,
        output_dir       = OUTPUT_DIR,
        run_name         = "RUN1",
        **RUN1,
    )

    # ========================================================
    # ANALYZE FIRST-RUN MULTIPLICITY
    # ========================================================

    classification = classify_parents(
        cat_first,
        moderate_percentile = MODERATE_PERCENTILE,
        extreme_percentile  = EXTREME_PERCENTILE,
    )

    moderate_ids = classification[
        "moderate_ids"
    ]

    extreme_ids = classification[
        "extreme_ids"
    ]

    # ========================================================
    # CREATE INPUT MAPS FOR RUNS 2 AND 3
    # ========================================================

    print(
        "\nCreating moderate detection map ..."
    )

    t0 = time.perf_counter()

    # --------------------------------------------------------
    # Build lookup table
    # --------------------------------------------------------

    max_parent_id = int(segmentation_map_detection.max())

    # 0 = normal
    # 1 = moderate
    # 2 = extreme
    group_lut = np.zeros(
        max_parent_id + 1,
        dtype=np.uint8
    )

    group_lut[moderate_ids] = 1
    group_lut[extreme_ids]  = 2

    # --------------------------------------------------------
    # Classify all pixels
    # --------------------------------------------------------

    group_map = group_lut[
        segmentation_map_detection
    ]

    # --------------------------------------------------------
    # Build masks
    # --------------------------------------------------------

    moderate_mask = (group_map == 1)
    extreme_mask  = (group_map == 2)

    # --------------------------------------------------------
    # Build detection maps for Run 2 / Run 3
    # --------------------------------------------------------

    segmap_moderate_detection = np.where(
        moderate_mask,
        segmentation_map_detection,
        0
    ).astype(
        segmentation_map_detection.dtype
    )

    segmap_extreme_detection = np.where(
        extreme_mask,
        segmentation_map_detection,
        0
    ).astype(
        segmentation_map_detection.dtype
    )

    subset_time = (time.perf_counter() - t0)

    # group_map and LUT are no longer needed
    del group_map
    del group_lut

    # Optional: save these for diagnostics/reproducibility.

    t0 = time.perf_counter()

    fits.writeto(
        OUTPUT_DIR /
        "moderate_detection_segmap.fits",
        segmap_moderate_detection,
        header    = header,
        overwrite = True,
    )

    fits.writeto(
        OUTPUT_DIR /
        "extreme_detection_segmap.fits",
        segmap_extreme_detection,
        header    = header,
        overwrite = True,
    )

    subset_write_time = (time.perf_counter() - t0)

    # ========================================================
    # RUN 2 — MODERATE
    # ========================================================

    segmap_second, cat_second, time_run2 = run_adp(
        image_data       = image_data,
        detection_segmap = segmap_moderate_detection,
        header           = header,
        output_dir       = OUTPUT_DIR,
        run_name         = "RUN2_MODERATE",
        **RUN2,
    )

    # We no longer need the moderate input segmentation map.
    del segmap_moderate_detection

    # ========================================================
    # RUN 3 — EXTREME
    # ========================================================

    segmap_third, cat_third, time_run3 = run_adp(
        image_data       = image_data,
        detection_segmap = segmap_extreme_detection,
        header           = header,
        output_dir       = OUTPUT_DIR,
        run_name         = "RUN3_EXTREME",
        **RUN3,
    )

    del segmap_extreme_detection

    # ========================================================
    # MERGE
    # ========================================================

    final_segmap, final_cat, merge_time = merge_runs(
        segmentation_map_detection = segmentation_map_detection,

        segmap_first = segmap_first,
        cat_first    = cat_first,

        segmap_second = segmap_second,
        cat_second    = cat_second,
        moderate_mask = moderate_mask,

        segmap_third = segmap_third,
        cat_third    = cat_third,
        extreme_mask = extreme_mask,
    )
    
    del moderate_mask
    del extreme_mask

    # ========================================================
    # VALIDATION
    # ========================================================

    validate_final_products(
        final_segmap,
        final_cat,
    )

    # ========================================================
    # SAVE FINAL PRODUCTS
    # ========================================================

    final_segmap_path = (
        OUTPUT_DIR /
        "ADP_final_segmentation.fits"
    )

    final_cat_path = (
        OUTPUT_DIR /
        "ADP_final_catalogue.fits"
    )

    fits.writeto(
        final_segmap_path,
        final_segmap,
        header=header,
        overwrite=True,
    )

    final_cat.write(
        final_cat_path,
        format="fits",
        overwrite=True,
    )


    # ========================================================
    # TIMING SUMMARY
    # ========================================================

    total_time     = (time.perf_counter() - pipeline_start)
    accounted_time = (time_run1 + subset_time + subset_write_time + time_run2 + time_run3 + merge_time)
    overhead_time  = total_time - accounted_time

    print("\n" + "=" * 60)
    print("TIMING SUMMARY")
    print("=" * 60)

    print(
        f"Run 1              : {time_run1:8.2f} s"
    )

    print(
        f"Subset maps        : {subset_time:8.2f} s"
    )

    print(
        f"Subset FITS writes : {subset_write_time:8.2f} s"
    )

    print(
        f"Run 2 Moderate     : {time_run2:8.2f} s"
    )

    print(
        f"Run 3 Extreme      : {time_run3:8.2f} s"
    )

    print(
        f"Merge              : {merge_time:8.2f} s"
    )

    print(
        f"Other / overhead   : {overhead_time:8.2f} s"
    )

    print("-" * 60)

    print(
        f"TOTAL PIPELINE      : {total_time:8.2f} s"
    )


    return {
        "segmap": final_segmap,
        "catalogue": final_cat,
        "classification": classification,
        "timings": {
            "run1": time_run1,
            "run2": time_run2,
            "run3": time_run3,
            "merge": merge_time,
            "total": total_time,
        },
    }

def load_adp_inputs(
    image_path,
    detection_path,
    sigma_pix=0.5,
):
    # ---------------------------------------------------------
    # Load science image
    # ---------------------------------------------------------
    with fits.open(image_path, memmap=False) as hdul:
        image_data = np.asarray(
            hdul[0].data,
            dtype=np.float64,
            order="C",
        )
        header = hdul[0].header.copy()

    # ---------------------------------------------------------
    # Load detection segmentation map
    # ---------------------------------------------------------
    with fits.open(detection_path, memmap=False) as hdul:
        segmentation_map_detection = np.asarray(
            hdul[0].data,
            dtype=np.int32,
            order="C",
        )

    # ---------------------------------------------------------
    # Basic consistency check
    # ---------------------------------------------------------
    if image_data.shape != segmentation_map_detection.shape:
        raise ValueError(
            "Image and detection segmentation map have "
            f"different shapes: {image_data.shape} vs "
            f"{segmentation_map_detection.shape}"
        )

    # ---------------------------------------------------------
    # Asterism-like smoothing
    # ---------------------------------------------------------
    if sigma_pix is not None and sigma_pix > 0:
        image_data = gaussian_filter(
            image_data,
            sigma=sigma_pix,
        )

    return image_data, segmentation_map_detection, header


# ---------------------------------------------------------
# Load inputs
# ---------------------------------------------------------
image_data, segmentation_map_detection, header = load_adp_inputs(
    image_path=config["image_path"],
    detection_path=config["detection_path"],
    sigma_pix=config["preprocessing"]["sigma_pix"],
)

# Iterative Run
results = run_iterative_adp(
    image_data,
    segmentation_map_detection,
    header,
)
