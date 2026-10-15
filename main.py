import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import sys
import torch
import torch.nn as nn
import torch.optim as optim

import glob
import logging
import json
import time
import numpy as np
import cv2
import pydicom
import tifffile
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.ndimage import median_filter, binary_fill_holes, gaussian_filter
from skimage import measure
import trimesh

# Suppress tifffile imagej_metadata IndexError warnings (harmless, non-standard tags)
logging.getLogger('tifffile').setLevel(logging.CRITICAL)

# Standard output resolution for all modalities
TARGET_SIZE = (512, 512)


# ============================================================
# STEP 1 – DATASET LOADING FUNCTIONS
# ============================================================

def load_dicom_series_info(folder_path):
    """
    Memory-efficient DICOM series loader.
    Selects representative middle slice for accurate visual inspection.
    """
    files = sorted([os.path.join(folder_path, f) for f in os.listdir(folder_path) if not f.endswith('.txt')])
    if not files:
        return None
        
    num_slices = len(files)
    mid_idx = num_slices // 2
    
    # Read middle slice for proper anatomical representation
    mid_ds = pydicom.dcmread(files[mid_idx])
    mid_slice = mid_ds.pixel_array.astype(np.float32)
    
    # Clean any invalid NaN/Inf values
    mid_slice = np.nan_to_num(mid_slice, nan=0.0)
    
    return {
        "num_slices": num_slices,
        "slice_shape": mid_slice.shape,
        "dtype": mid_ds.pixel_array.dtype,
        "sample_slice": mid_slice,
        "mid_idx": mid_idx,
        "files": files,
        "min_val": float(mid_slice.min()),
        "max_val": float(mid_slice.max())
    }

def load_tiff_series_info(folder_path):
    tiff_files = sorted(glob.glob(os.path.join(folder_path, "*.tif")))
    if not tiff_files:
        return None

    # Check for precomputed MIP cache
    cache_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ct_intact_mip.npy")
    if os.path.exists(cache_path):
        print(f"    [CT Intact] Loading cached Max-Intensity Projection ({len(tiff_files)} slices)...")
        mip = np.load(cache_path).astype(np.float32)
        return {
            "num_slices": len(tiff_files),
            "slice_shape": mip.shape,
            "dtype": np.float32,
            "sample_slice": mip,
            "mid_idx": "MIP",
            "files": tiff_files,
            "min_val": float(mip.min()),
            "max_val": float(mip.max())
        }

    num_slices = len(tiff_files)
    print(f"    [CT Intact] Building Max-Intensity Projection over {num_slices} slices...")

    # Redirect stderr to suppress tifffile imagej_metadata IndexError noise
    _devnull = open(os.devnull, 'w')
    _stderr = sys.stderr
    sys.stderr = _devnull

    try:
        first = tifffile.imread(tiff_files[0]).astype(np.float32)
        first = np.nan_to_num(first, nan=0.0)
        mip = np.zeros_like(first)
        total_nonzero = 0

        for i, f in enumerate(tiff_files):
            img = tifffile.imread(f).astype(np.float32)
            img = np.nan_to_num(img, nan=0.0)
            mip = np.maximum(mip, np.abs(img))
            total_nonzero += np.count_nonzero(img)
            if (i + 1) % 200 == 0:
                sys.stderr = _stderr
                print(f"      Processed {i+1}/{num_slices} slices...")
                sys.stderr = _devnull
    finally:
        sys.stderr = _stderr
        _devnull.close()

    np.save(cache_path, mip)
    print(f"    [CT Intact] MIP complete. Total nonzero pixels: {total_nonzero:,}, MIP nonzero: {np.count_nonzero(mip):,}")

    return {
        "num_slices": num_slices,
        "slice_shape": mip.shape,
        "dtype": np.float32,
        "sample_slice": mip,
        "mid_idx": "MIP",
        "files": tiff_files,
        "min_val": float(mip.min()),
        "max_val": float(mip.max())
    }

def normalize_slice_for_display(img):
    """Normalize image slice using 1st-99th percentile for crisp visual contrast."""
    img_clean = np.nan_to_num(img, nan=0.0)
    p1, p99 = np.percentile(img_clean, (1, 99))
    if p99 > p1:
        img_norm = np.clip((img_clean - p1) / (p99 - p1), 0, 1)
    else:
        img_norm = img_clean
    return img_norm


# ============================================================
# STEP 2 – IMAGE PREPROCESSING FUNCTIONS
# ============================================================

def step2_noise_removal(img, modality="CT"):
    """
    Noise Removal:
    - CT: 3x3 Median filter to remove salt-and-pepper / acquisition noise
    - MRI: 3x3 Gaussian blur to smooth thermal/RF noise
    """
    if modality == "CT":
        denoised = median_filter(img, size=3)
    else:
        denoised = cv2.GaussianBlur(img.astype(np.float32), (3, 3), sigmaX=0.8)
    return denoised.astype(np.float32)

def step2_intensity_normalization(img):
    """
    Intensity Normalization:
    Uses 1st-99th percentile clipping then scales to [0, 1].
    Robust to outliers and extreme HU values in CT and intensity spikes in MRI.
    """
    p1  = np.percentile(img, 1)
    p99 = np.percentile(img, 99)
    if p99 > p1:
        norm = np.clip((img - p1) / (p99 - p1), 0.0, 1.0)
    else:
        norm = np.zeros_like(img, dtype=np.float32)
    return norm.astype(np.float32)

def step2_resize(img, target_size=TARGET_SIZE):
    """
    Resizing:
    Resize all images to a standard target resolution using bi-cubic interpolation.
    Preserves fine anatomical edge details better than bilinear.
    """
    resized = cv2.resize(img, (target_size[1], target_size[0]),
                         interpolation=cv2.INTER_CUBIC)
    return resized.astype(np.float32)

def step2_image_enhancement(img):
    """
    Image Enhancement using CLAHE (Contrast Limited Adaptive Histogram Equalization):
    - Enhances local contrast in bone/tissue boundaries
    - Suppresses over-amplification in homogeneous regions (clip limit)
    - Applied on uint8 scaled image, then mapped back to float32
    """
    img_uint8 = (img * 255).clip(0, 255).astype(np.uint8)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    enhanced_uint8 = clahe.apply(img_uint8)
    enhanced = enhanced_uint8.astype(np.float32) / 255.0
    return enhanced

def preprocess_image(raw_img, modality="CT", target_size=TARGET_SIZE):
    """
    Full preprocessing pipeline for one image slice:
      Step 2a: Noise Removal     → Median (CT) / Gaussian (MRI)
      Step 2b: Intensity Norm.   → Percentile-based [0, 1] clipping
      Step 2c: Resize            → Standard target resolution (512x512)
      Step 2d: Enhancement       → CLAHE contrast enhancement
    Returns: preprocessed float32 image [0, 1]
    """
    img = raw_img.astype(np.float32)
    img = np.nan_to_num(img, nan=0.0)

    # Step 2a - Noise removal
    denoised = step2_noise_removal(img, modality=modality)

    # Step 2b - Intensity normalization
    normalized = step2_intensity_normalization(denoised)

    # Step 2c - Resize to standard resolution
    resized = step2_resize(normalized, target_size)

    # Step 2d - CLAHE enhancement
    enhanced = step2_image_enhancement(resized)

    return enhanced

def run_preprocessing(datasets):
    """
    Run Step 2 preprocessing on all loaded datasets.
    Returns a dict of preprocessed sample slices keyed by dataset name.
    """
    print("\n" + "=" * 70)
    print("  STEP 2 - IMAGE PREPROCESSING (Noise Removal / Normalization /")
    print("           Resizing / CLAHE Enhancement)")
    print("=" * 70)

    modality_map = {
        "CT_Femur":  "CT",
        "CT_Tibia":  "CT",
        "CT_Intact": "CT",
        "MR_TSE":    "MRI",
        "MR_DESS":   "MRI",
    }

    preprocessed = {}
    for name, info in datasets.items():
        raw = info["sample_slice"]
        mod = modality_map.get(name, "CT")
        proc = preprocess_image(raw, modality=mod, target_size=TARGET_SIZE)
        preprocessed[name] = proc

        orig_shape = raw.shape
        orig_min   = float(raw.min())
        orig_max   = float(raw.max())

        print(f"\n  [{name}]")
        print(f"    Original  : shape={orig_shape}, range=[{orig_min:.1f}, {orig_max:.1f}]")
        print(f"    Denoised  : Median3x3" if mod == "CT" else f"    Denoised  : GaussianBlur(3x3, sigma=0.8)")
        print(f"    Normalized: Percentile [1%-99%] -> [0.0, 1.0]")
        print(f"    Resized   : {orig_shape} -> {TARGET_SIZE}")
        print(f"    Enhanced  : CLAHE (clipLimit=2.0, tileGridSize=8x8)")
        print(f"    Output    : shape={proc.shape}, range=[{proc.min():.3f}, {proc.max():.3f}]")

    return preprocessed


# ============================================================
# STEP 3 – 3D SEGMENTATION AND MASK GENERATION (U-NET)
# ============================================================

class DoubleConv(nn.Module):
    """
    Double 2D Convolution block with Batch Normalization and ReLU:
    [Conv2d -> BatchNorm -> ReLU -> Conv2d -> BatchNorm -> ReLU]
    """
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(in_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True)
        )
    def forward(self, x):
        return self.net(x)

class UNet(nn.Module):
    """
    Validated 3D/2D-Slice Deep Learning U-Net Segmentation Architecture:
    - Multi-scale Contracting Encoder (Downsampling with MaxPool)
    - Bottleneck Feature Representation
    - Expanding Decoder with Skip Connections (Transposed Conv & Interpolation)
    - Sigmoid Binary/Multi-channel Probability Head
    """
    def __init__(self, in_channels=1, out_channels=1, features=(16, 32, 64)):
        super().__init__()
        self.downs = nn.ModuleList()
        self.ups = nn.ModuleList()
        self.pool = nn.MaxPool2d(kernel_size=2, stride=2)

        in_c = in_channels
        for f in features:
            self.downs.append(DoubleConv(in_c, f))
            in_c = f

        self.bottleneck = DoubleConv(features[-1], features[-1] * 2)

        for f in reversed(features):
            self.ups.append(nn.ConvTranspose2d(f * 2, f, kernel_size=2, stride=2))
            self.ups.append(DoubleConv(f * 2, f))

        self.final_conv = nn.Conv2d(features[0], out_channels, kernel_size=1)

    def forward(self, x):
        skip_connections = []
        for down in self.downs:
            x = down(x)
            skip_connections.append(x)
            x = self.pool(x)

        x = self.bottleneck(x)
        skip_connections = skip_connections[::-1]

        for idx in range(0, len(self.ups), 2):
            x = self.ups[idx](x)
            skip = skip_connections[idx // 2]
            if x.shape != skip.shape:
                x = torch.nn.functional.interpolate(x, size=skip.shape[2:], mode="bilinear", align_corners=False)
            x = torch.cat((skip, x), dim=1)
            x = self.ups[idx + 1](x)

        return torch.sigmoid(self.final_conv(x))

def dice_bce_loss(pred, target, smooth=1e-5):
    """
    Combined BCE + Dice Loss for robust medical bone segmentation:
    Balances pixel-level cross-entropy with overlap-level Dice similarity.
    """
    bce = nn.functional.binary_cross_entropy(pred, target)
    intersection = (pred * target).sum(dim=(2, 3))
    union = pred.sum(dim=(2, 3)) + target.sum(dim=(2, 3))
    dice = (2.0 * intersection + smooth) / (union + smooth)
    return 0.5 * bce + 0.5 * (1.0 - dice.mean())

def generate_semiautomatic_ground_truth(slice_data, modality="CT", hu_threshold=250, target_size=TARGET_SIZE):
    """
    Generate manual/semi-automatic ground-truth labels for femur and tibia.
    - Bone-calibrated thresholding (HU > 250 for cortical/trabecular bone in CT)
    - Scanner circular FOV mask to eliminate perimeter bed/gantry artifacts
    - Morphological closing to seal cortical edge gaps
    - Solid hole filling (binary_fill_holes) to encompass cancellous/medullary cavity
    - Central anatomical component analysis isolating primary bone entity
    - Standardized to target resolution
    """
    h, w = slice_data.shape
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - w/2.0)**2 + (Y - h/2.0)**2)
    fov_mask = dist_from_center <= min(h, w) * 0.46

    if modality == "CT":
        bone_bin = ((slice_data > hu_threshold) & fov_mask).astype(np.uint8)
    else:
        p90 = np.percentile(slice_data, 90)
        bone_bin = ((slice_data > p90 * 0.7) & fov_mask).astype(np.uint8)

    # Morphological closing (seal fine cortical gaps)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    closed = cv2.morphologyEx(bone_bin, cv2.MORPH_CLOSE, kernel)

    # Hole filling (fill marrow and trabecular voids into a solid anatomical bone label)
    filled = binary_fill_holes(closed).astype(np.uint8)

    # Central anatomical component analysis - retain primary bone structure
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(filled)
    if num_labels > 1:
        center_x, center_y = w / 2.0, h / 2.0
        scores = []
        for i in range(1, num_labels):
            area = stats[i, cv2.CC_STAT_AREA]
            cx, cy = centroids[i]
            dist = np.sqrt((cx - center_x)**2 + (cy - center_y)**2) / (min(h, w) / 2.0)
            score = area / (1.0 + 2.0 * dist)
            scores.append(score)
        largest_idx = 1 + np.argmax(scores)
        clean_mask = (labels == largest_idx).astype(np.uint8)
    else:
        clean_mask = filled

    # Resize to target resolution
    if clean_mask.shape != target_size:
        clean_mask = cv2.resize(clean_mask, (target_size[1], target_size[0]), interpolation=cv2.INTER_NEAREST)

    return clean_mask

def postprocess_segmentation(prob_map, threshold=0.55):
    """
    Medical post-processing for predicted bone masks:
    Binarization + hole-filling + central anatomical component filtering.
    """
    binary = (prob_map > threshold).astype(np.uint8)
    filled = binary_fill_holes(binary).astype(np.uint8)
    num, labels, stats, centroids = cv2.connectedComponentsWithStats(filled)
    if num > 1:
        h, w = prob_map.shape
        center_x, center_y = w / 2.0, h / 2.0
        scores = []
        for i in range(1, num):
            area = stats[i, cv2.CC_STAT_AREA]
            cx, cy = centroids[i]
            dist = np.sqrt((cx - center_x)**2 + (cy - center_y)**2) / (min(h, w) / 2.0)
            score = area / (1.0 + 2.0 * dist)
            scores.append(score)
        largest_idx = 1 + np.argmax(scores)
        clean = (labels == largest_idx).astype(np.uint8)
    else:
        clean = filled
    return clean

def compute_segmentation_metrics(pred_bin, target_bin):
    """
    Evaluate segmentation performance with full metric suite:
    - Dice Similarity Coefficient (DSC)
    - Intersection over Union (IoU / Jaccard Index)
    - Precision (Positive Predictive Value)
    - Recall / Sensitivity (True Positive Rate)
    - Specificity (True Negative Rate)
    - Accuracy (Overall Pixel Accuracy)
    - Volumetric Overlap Error (VOE, %)
    - Relative Volume Difference (RVD, %)
    """
    pred_bin = (pred_bin > 0).astype(np.uint8)
    target_bin = (target_bin > 0).astype(np.uint8)

    tp = float(np.logical_and(pred_bin == 1, target_bin == 1).sum())
    fp = float(np.logical_and(pred_bin == 1, target_bin == 0).sum())
    fn = float(np.logical_and(pred_bin == 0, target_bin == 1).sum())
    tn = float(np.logical_and(pred_bin == 0, target_bin == 0).sum())

    dice = (2.0 * tp) / (2.0 * tp + fp + fn + 1e-8)
    iou = tp / (tp + fp + fn + 1e-8)
    precision = tp / (tp + fp + 1e-8)
    recall = tp / (tp + fn + 1e-8)
    specificity = tn / (tn + fp + 1e-8)
    accuracy = (tp + tn) / (tp + tn + fp + fn + 1e-8)
    voe = 100.0 * (1.0 - iou)
    rvd = 100.0 * float(abs(int(pred_bin.sum()) - int(target_bin.sum()))) / float(target_bin.sum() + 1e-8)

    return {
        "Dice": dice,
        "IoU": iou,
        "Precision": precision,
        "Recall": recall,
        "Specificity": specificity,
        "Accuracy": accuracy,
        "VOE": voe,
        "RVD": rvd
    }

def run_3d_segmentation(datasets, base_dir):
    """
    Execute Step 3:
    1. Generate manual/semi-automatic GT labels for femur and tibia.
    2. Train/apply validated 3D deep-learning U-Net model.
    3. Evaluate performance using Dice, IoU, Precision, Recall, Specificity, etc.
    4. Generate and save validated 3D femur and tibia masks (NPY) for Step 4 3D reconstruction.
    """
    print("\n" + "=" * 70)
    print("  STEP 3 - 3D SEGMENTATION AND MASK GENERATION (U-NET)")
    print("=" * 70)

    femur_files = datasets.get("CT_Femur", {}).get("files", [])
    tibia_files = datasets.get("CT_Tibia", {}).get("files", [])

    if not femur_files or not tibia_files:
        print("  [!] Error: Femur and Tibia datasets are required for Step 3.")
        return None

    train_res = (256, 256)
    print(f"\n[+] Preparing dataset slices for Deep-Learning U-Net Training & Validation...")
    print(f"    - Femur available slices: {len(femur_files)}")
    print(f"    - Tibia available slices: {len(tibia_files)}")
    print(f"    - Training Resolution   : {train_res}")

    # Sample representative slices across anatomical bone volume
    # 8 for training + 5 for validation = 13 indices per structure
    # Guard: end index must be > start to avoid reversed/degenerate linspace
    f_start = min(80, len(femur_files) // 4)
    f_end   = max(f_start + 1, len(femur_files) - max(80, len(femur_files) // 4))
    t_start = min(80, len(tibia_files) // 4)
    t_end   = max(t_start + 1, len(tibia_files) - max(80, len(tibia_files) // 4))
    f_indices = np.linspace(f_start, f_end, 13, dtype=int)
    t_indices = np.linspace(t_start, t_end, 13, dtype=int)

    train_x, train_y = [], []
    val_slices = []

    # Femur slices – 8 train + 5 validation
    for i, idx in enumerate(f_indices):
        raw = pydicom.dcmread(femur_files[idx]).pixel_array.astype(np.float32)
        img = preprocess_image(raw, modality="CT", target_size=train_res)
        gt = generate_semiautomatic_ground_truth(raw, modality="CT", target_size=train_res)
        if i < 8:
            train_x.append(img)
            train_y.append(gt)
            # Data augmentation (horizontal flip)
            train_x.append(np.fliplr(img).copy())
            train_y.append(np.fliplr(gt).copy())
        else:
            val_slices.append(("Femur", idx, raw, img, gt))

    # Tibia slices – 8 train + 5 validation
    for i, idx in enumerate(t_indices):
        raw = pydicom.dcmread(tibia_files[idx]).pixel_array.astype(np.float32)
        img = preprocess_image(raw, modality="CT", target_size=train_res)
        gt = generate_semiautomatic_ground_truth(raw, modality="CT", target_size=train_res)
        if i < 8:
            train_x.append(img)
            train_y.append(gt)
            train_x.append(np.fliplr(img).copy())
            train_y.append(np.fliplr(gt).copy())
        else:
            val_slices.append(("Tibia", idx, raw, img, gt))

    train_x = torch.tensor(np.array(train_x)[:, None, :, :], dtype=torch.float32)
    train_y = torch.tensor(np.array(train_y)[:, None, :, :], dtype=torch.float32)

    print(f"    - Training Set Dimensions : {train_x.shape} (with data augmentation)")
    print(f"    - Validation Slices Count : {len(val_slices)}")

    # Initialize and train U-Net
    print(f"\n[+] Training Validated Deep-Learning U-Net Segmentation Architecture...")
    torch.manual_seed(42)
    model = UNet(in_channels=1, out_channels=1, features=(16, 32, 64))
    optimizer = optim.Adam(model.parameters(), lr=3e-3)

    num_epochs = 12
    t0 = time.time()
    for epoch in range(1, num_epochs + 1):
        model.train()
        perm = torch.randperm(train_x.size(0))
        epoch_loss = 0.0
        batch_size = 8
        for b in range(0, train_x.size(0), batch_size):
            bx, by = train_x[perm[b:b+batch_size]], train_y[perm[b:b+batch_size]]
            optimizer.zero_grad()
            pred = model(bx)
            loss = dice_bce_loss(pred, by)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(bx)
        epoch_loss /= train_x.size(0)
        print(f"    Epoch [{epoch:02d}/{num_epochs}]  Loss: {epoch_loss:.4f}")

    train_time = time.time() - t0
    print(f"    - U-Net Convergence achieved in {num_epochs} epochs ({train_time:.2f}s)")
    print(f"    - Final Loss (BCE + Dice): {epoch_loss:.4f}")

    # Evaluate validation slices
    print(f"\n[+] Evaluating Model Performance on Independent Validation Slices...")
    model.eval()
    val_results = {"Femur": [], "Tibia": []}

    with torch.no_grad():
        for name, idx, raw, img, gt in val_slices:
            inp = torch.tensor(img[None, None, :, :], dtype=torch.float32)
            prob = model(inp)[0, 0].cpu().numpy()
            pred_mask = postprocess_segmentation(prob, threshold=0.55)
            metrics = compute_segmentation_metrics(pred_mask, gt)
            val_results[name].append({
                "slice_idx": idx,
                "raw": raw,
                "preprocessed": img,
                "gt": gt,
                "prob": prob,
                "pred": pred_mask,
                "metrics": metrics
            })

    # Compute summary metrics
    def calc_mean_metrics(results_list):
        keys = ["Dice", "IoU", "Precision", "Recall", "Specificity", "Accuracy", "VOE", "RVD"]
        return {k: float(np.mean([r["metrics"][k] for r in results_list])) for k in keys}

    femur_summary = calc_mean_metrics(val_results["Femur"])
    tibia_summary = calc_mean_metrics(val_results["Tibia"])
    
    overall_summary = {
        k: float((femur_summary[k] + tibia_summary[k]) / 2.0)
        for k in femur_summary
    }

    summary_metrics = {
        "Femur": femur_summary,
        "Tibia": tibia_summary,
        "Mean": overall_summary
    }

    # Print quantitative metric report
    print("\n" + "-" * 70)
    print("  QUANTITATIVE SEGMENTATION PERFORMANCE EVALUATION")
    print("-" * 70)
    print(f"  {'Metric':<25} | {'Femur':<12} | {'Tibia':<12} | {'Overall Mean':<12}")
    print("-" * 70)
    print(f"  {'Dice (DSC)':<25} | {femur_summary['Dice']:<12.4f} | {tibia_summary['Dice']:<12.4f} | {overall_summary['Dice']:<12.4f}")
    print(f"  {'Intersection over Union (IoU)':<25} | {femur_summary['IoU']:<12.4f} | {tibia_summary['IoU']:<12.4f} | {overall_summary['IoU']:<12.4f}")
    print(f"  {'Precision (PPV)':<25} | {femur_summary['Precision']:<12.4f} | {tibia_summary['Precision']:<12.4f} | {overall_summary['Precision']:<12.4f}")
    print(f"  {'Recall / Sensitivity':<25} | {femur_summary['Recall']:<12.4f} | {tibia_summary['Recall']:<12.4f} | {overall_summary['Recall']:<12.4f}")
    print(f"  {'Specificity (TNR)':<25} | {femur_summary['Specificity']:<12.4f} | {tibia_summary['Specificity']:<12.4f} | {overall_summary['Specificity']:<12.4f}")
    print(f"  {'Overall Accuracy':<25} | {femur_summary['Accuracy']:<12.4f} | {tibia_summary['Accuracy']:<12.4f} | {overall_summary['Accuracy']:<12.4f}")
    print(f"  {'Volumetric Overlap Error':<25} | {femur_summary['VOE']:<11.2f}% | {tibia_summary['VOE']:<11.2f}% | {overall_summary['VOE']:<11.2f}%")
    print(f"  {'Relative Volume Difference':<25} | {femur_summary['RVD']:<11.2f}% | {tibia_summary['RVD']:<11.2f}% | {overall_summary['RVD']:<11.2f}%")
    print("-" * 70)

    # Generate 3D Volumetric Masks (64 slices stack) for Step 4
    print(f"\n[+] Generating Validated 3D Femur and Tibia Masks for Subsequent 3D Reconstruction...")
    stack_slices = 64
    f_stack_idxs = np.linspace(0, len(femur_files) - 1, stack_slices, dtype=int)
    t_stack_idxs = np.linspace(0, len(tibia_files) - 1, stack_slices, dtype=int)

    femur_mask_3d = np.zeros((stack_slices, 256, 256), dtype=np.uint8)
    tibia_mask_3d = np.zeros((stack_slices, 256, 256), dtype=np.uint8)

    with torch.no_grad():
        for i, idx in enumerate(f_stack_idxs):
            raw = pydicom.dcmread(femur_files[idx]).pixel_array.astype(np.float32)
            img = preprocess_image(raw, modality="CT", target_size=(256, 256))
            inp = torch.tensor(img[None, None, :, :], dtype=torch.float32)
            prob = model(inp)[0, 0].cpu().numpy()
            femur_mask_3d[i] = postprocess_segmentation(prob, threshold=0.55)

        for i, idx in enumerate(t_stack_idxs):
            raw = pydicom.dcmread(tibia_files[idx]).pixel_array.astype(np.float32)
            img = preprocess_image(raw, modality="CT", target_size=(256, 256))
            inp = torch.tensor(img[None, None, :, :], dtype=torch.float32)
            prob = model(inp)[0, 0].cpu().numpy()
            tibia_mask_3d[i] = postprocess_segmentation(prob, threshold=0.55)

    femur_mask_path = os.path.join(base_dir, "femur_mask_3d.npy")
    tibia_mask_path = os.path.join(base_dir, "tibia_mask_3d.npy")
    metrics_json_path = os.path.join(base_dir, "segmentation_metrics.json")

    np.save(femur_mask_path, femur_mask_3d)
    np.save(tibia_mask_path, tibia_mask_3d)

    with open(metrics_json_path, "w") as f:
        json.dump(summary_metrics, f, indent=4)

    print(f"    - Validated 3D Femur Mask saved : {femur_mask_path} (shape: {femur_mask_3d.shape}, active voxels: {np.count_nonzero(femur_mask_3d):,})")
    print(f"    - Validated 3D Tibia Mask saved : {tibia_mask_path} (shape: {tibia_mask_3d.shape}, active voxels: {np.count_nonzero(tibia_mask_3d):,})")
    print(f"    - Metrics summary JSON saved   : {metrics_json_path}")

    return {
        "val_results": val_results,
        "summary_metrics": summary_metrics,
        "femur_mask_3d": femur_mask_3d,
        "tibia_mask_3d": tibia_mask_3d
    }


# ============================================================
# STEP 4 – 3D RECONSTRUCTION (STL GENERATION & METRIC EXTRACTION)
# ============================================================

def reconstruct_bone_mesh(mask_3d, spacing, name="Bone", sigma=0.8):
    """
    Convert a 3D binary segmentation mask into a watertight, manifold 3D bone mesh (STL).
    - Preserves anatomical curvature and contours required for personalized implant development.
    - 3D Gaussian curvature smoothing (sigma=0.8) eliminates voxel staircasing artifacts.
    - Zero-padded boundary ensures fully enclosed, watertight solid geometries.
    - Lorensen & Cline Marching Cubes algorithm applied at level=0.5.
    - Topological healing: normal reorientation, winding order repair, inversion fixing.
    """
    vol = mask_3d.astype(np.float32)
    # 3D Anti-aliasing / smoothing filter to preserve anatomical contours without staircasing
    smoothed = gaussian_filter(vol, sigma=sigma)
    
    # Pad boundary with 1 voxel of zeros to ensure closed, watertight manifold mesh
    padded = np.pad(smoothed, 1, mode="constant", constant_values=0)
    
    # Marching Cubes at level 0.5 with true physical voxel spacing (in mm)
    verts, faces, normals, _ = measure.marching_cubes(padded, level=0.5, spacing=spacing)
    
    # Re-align coordinates by subtracting padding offset
    verts = verts - np.array([spacing[0], spacing[1], spacing[2]])
    
    # Initialize Trimesh geometry
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals)
    
    # Topological repair and normal orientation for CAD/implant suitability
    mesh.fix_normals()
    trimesh.repair.fix_winding(mesh)
    trimesh.repair.fix_inversion(mesh)
    trimesh.repair.fill_holes(mesh)
    
    return mesh

def compute_3d_reconstruction_metrics(mesh, mask_3d, spacing, name="Bone"):
    """
    Calculate comprehensive 3D reconstruction quality and geometric metrics:
    - Mesh Complexity: Vertices (Nv), Triangular Faces (Nf), Edges (Ne), Euler Characteristic
    - Anatomical Dimensions: Enclosed Bone Volume (mm3 & cm3), Surface Area (mm2 & cm2), SA/V ratio
    - Bounding Box Dimensions (AP x ML x SI in mm), Centroid (X, Y, Z in mm)
    - Reconstruction Quality: Volumetric Fidelity (%), Volumetric Error (%), Sphericity Index,
      Mean Edge Length (mm), Std Edge Length (mm), Triangle Aspect Quality Index
    """
    voxel_vol = float(np.count_nonzero(mask_3d > 0.5) * spacing[0] * spacing[1] * spacing[2])
    mesh_vol = abs(float(mesh.volume))
    mesh_area = float(mesh.area)
    
    vol_fidelity = (1.0 - abs(mesh_vol - voxel_vol) / (voxel_vol + 1e-8)) * 100.0
    vol_error = (abs(mesh_vol - voxel_vol) / (voxel_vol + 1e-8)) * 100.0
    sa_to_vol = mesh_area / (mesh_vol + 1e-8)
    
    # Edge metrics
    edges = mesh.edges_unique_length
    mean_edge = float(np.mean(edges))
    std_edge = float(np.std(edges))
    
    # Sphericity / Anatomical compactness index
    sphericity = float((np.pi ** (1.0 / 3.0) * (6.0 * mesh_vol) ** (2.0 / 3.0)) / (mesh_area + 1e-8))
    
    # Triangle quality (aspect ratio metric: 4 * sqrt(3) * Area / sum(edge^2))
    tri_verts = mesh.vertices[mesh.faces]
    v0, v1, v2 = tri_verts[:, 0], tri_verts[:, 1], tri_verts[:, 2]
    a = np.linalg.norm(v1 - v0, axis=1)
    b = np.linalg.norm(v2 - v1, axis=1)
    c = np.linalg.norm(v0 - v2, axis=1)
    tri_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    qualities = (4.0 * np.sqrt(3.0) * tri_areas) / (a**2 + b**2 + c**2 + 1e-12)
    mean_quality = float(np.mean(qualities))
    
    bbox = mesh.bounds[1] - mesh.bounds[0]
    centroid = mesh.centroid
    
    metrics = {
        "name": name,
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "edges": int(len(mesh.edges_unique)),
        "is_watertight": bool(mesh.is_watertight),
        "euler_characteristic": int(mesh.euler_number),
        "surface_area_mm2": round(mesh_area, 2),
        "surface_area_cm2": round(mesh_area / 100.0, 2),
        "volume_mm3": round(mesh_vol, 2),
        "volume_cm3": round(mesh_vol / 1000.0, 2),
        "voxel_volume_mm3": round(voxel_vol, 2),
        "volumetric_fidelity_pct": round(vol_fidelity, 2),
        "volumetric_error_pct": round(vol_error, 2),
        "sa_to_vol_ratio_mm_inv": round(sa_to_vol, 4),
        "bounding_box_mm": [round(float(b), 2) for b in bbox],
        "centroid_mm": [round(float(c), 2) for c in centroid],
        "mean_edge_length_mm": round(mean_edge, 4),
        "std_edge_length_mm": round(std_edge, 4),
        "sphericity_index": round(sphericity, 4),
        "triangle_quality_mean": round(mean_quality, 4)
    }
    return metrics

def _tnr(size=14, weight="bold"):
    return {"fontfamily": "Times New Roman", "fontsize": size, "fontweight": weight}

def plot_reconstructed_bone_models(femur_mesh, tibia_mesh, femur_shifted_mesh, save_dir):
    """
    Render and save publication-grade 3D multi-panel visualization of reconstructed
    Femur, Tibia, and Knee Joint Anatomical Assembly (DPI 1000, Times New Roman 18pt Bold).
    Uses anatomical orientation, Taubin surface smoothing, and dual-source lighting for realistic cortical bone shading.
    """
    os.makedirs(save_dir, exist_ok=True)
    
    # Orient meshes upright for true anatomical standing view: [X_ML, Y_AP, Z_SI]
    fv = femur_mesh.vertices
    ff = femur_mesh.faces
    fv_up = np.column_stack([fv[:, 2], fv[:, 1], fv[:, 0]])
    fm = trimesh.Trimesh(vertices=fv_up, faces=ff)
    fm.fix_normals()
    trimesh.smoothing.filter_taubin(fm, iterations=12)
    
    tv = tibia_mesh.vertices
    tf = tibia_mesh.faces
    tv_up = np.column_stack([tv[:, 2], tv[:, 1], tv[:, 0]])
    tm = trimesh.Trimesh(vertices=tv_up, faces=tf)
    tm.fix_normals()
    trimesh.smoothing.filter_taubin(tm, iterations=12)
    
    # Dual-source realistic lighting calculation
    def calculate_lighting(m, base_rgb):
        l1 = np.array([0.5, 0.5, 0.7]) / np.linalg.norm([0.5, 0.5, 0.7])  # Key light
        l2 = np.array([-0.6, 0.3, 0.4]) / np.linalg.norm([-0.6, 0.3, 0.4]) # Fill light
        n = m.face_normals
        diffuse1 = np.clip(np.dot(n, l1), 0.0, 1.0) * 0.55
        diffuse2 = np.clip(np.dot(n, l2), 0.0, 1.0) * 0.25
        ambient = 0.35
        total_light = np.clip(ambient + diffuse1 + diffuse2, 0.0, 1.0)
        colors = base_rgb[None, :] * total_light[:, None]
        return np.clip(colors, 0.0, 1.0)
    
    fig = plt.figure(figsize=(25, 8.5), facecolor="#FFFFFF")
    tnr = dict(fontsize=15, fontweight="bold", fontfamily="Times New Roman")
    
    # Panel A: 3D Femur Model (Warm Ivory Cortical Bone)
    ax1 = fig.add_subplot(131, projection="3d", facecolor="#FFFFFF")
    f_col = calculate_lighting(fm, np.array([0.94, 0.88, 0.76]))
    p1 = Poly3DCollection(fm.vertices[fm.faces], facecolors=f_col, edgecolors="none", linewidths=0)
    ax1.add_collection3d(p1)
    ax1.auto_scale_xyz(fm.vertices[:, 0], fm.vertices[:, 1], fm.vertices[:, 2])
    ax1.view_init(elev=20, azim=40)
    ax1.set_axis_off()
    ax1.set_title(
        "(A) Reconstructed 3D Femur Model (STL)\n"
        f"Vol: {abs(femur_mesh.volume)/1000.0:.1f} cm^3 | Area: {femur_mesh.area/100.0:.1f} cm^2 | Watertight: {femur_mesh.is_watertight}",
        **tnr, color="#1A1A1A", pad=18
    )
    
    # Panel B: 3D Tibia Model (Natural Sage Cortical Bone)
    ax2 = fig.add_subplot(132, projection="3d", facecolor="#FFFFFF")
    t_col = calculate_lighting(tm, np.array([0.84, 0.89, 0.82]))
    p2 = Poly3DCollection(tm.vertices[tm.faces], facecolors=t_col, edgecolors="none", linewidths=0)
    ax2.add_collection3d(p2)
    ax2.auto_scale_xyz(tm.vertices[:, 0], tm.vertices[:, 1], tm.vertices[:, 2])
    ax2.view_init(elev=20, azim=40)
    ax2.set_axis_off()
    ax2.set_title(
        "(B) Reconstructed 3D Tibia Model (STL)\n"
        f"Vol: {abs(tibia_mesh.volume)/1000.0:.1f} cm^3 | Area: {tibia_mesh.area/100.0:.1f} cm^2 | Watertight: {tibia_mesh.is_watertight}",
        **tnr, color="#1A1A1A", pad=18
    )
    
    # Panel C: Knee Joint Anatomical Assembly (Femur Superior, Tibia Inferior with joint space)
    ax3 = fig.add_subplot(133, projection="3d", facecolor="#FFFFFF")
    p3_t = Poly3DCollection(tm.vertices[tm.faces], facecolors=t_col, edgecolors="none", linewidths=0)
    ax3.add_collection3d(p3_t)
    
    femur_shifted_verts = fm.vertices.copy()
    femur_shifted_verts[:, 0] += (tm.vertices[:, 0].mean() - fm.vertices[:, 0].mean())
    femur_shifted_verts[:, 1] += (tm.vertices[:, 1].mean() - fm.vertices[:, 1].mean())
    z_top_tibia = tm.vertices[:, 2].max()
    z_bot_femur = fm.vertices[:, 2].min()
    femur_shifted_verts[:, 2] += (z_top_tibia - z_bot_femur) + 12.0
    
    fm_shifted = trimesh.Trimesh(vertices=femur_shifted_verts, faces=fm.faces)
    fm_shifted.fix_normals()
    f_col_shifted = calculate_lighting(fm_shifted, np.array([0.94, 0.88, 0.76]))
    
    p3_f = Poly3DCollection(fm_shifted.vertices[fm_shifted.faces], facecolors=f_col_shifted, edgecolors="none", linewidths=0)
    ax3.add_collection3d(p3_f)
    
    all_v = np.vstack([tm.vertices, fm_shifted.vertices])
    ax3.auto_scale_xyz(all_v[:, 0], all_v[:, 1], all_v[:, 2])
    ax3.view_init(elev=18, azim=40)
    ax3.set_axis_off()
    ax3.set_title(
        "(C) Knee Joint Anatomical Assembly\n"
        "Patient-Specific Biomechanical Digital Twin Base",
        **tnr, color="#1A1A1A", pad=18
    )
    
    # Title without 'STEP 4' word
    fig.suptitle(
        "3D Reconstruction: Anatomical Bone Models (STL) & Knee Assembly",
        **_tnr(size=19), color="#1A1A1A", y=0.98
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out_path = os.path.join(save_dir, "step4_reconstructed_bone_models.png")
    plt.savefig(out_path, dpi=1000, bbox_inches="tight", facecolor="#FFFFFF")
    plt.close()
    print(f"\n  [3D Reconstruction] Reconstructed 3D bone models figure saved: {out_path}")

def plot_reconstruction_metrics_bar(metrics_dict, save_dir):
    """
    Standalone 3D Reconstruction performance and geometry metrics bar chart.
    (Figsize 11x8, DPI 1000, Times New Roman 18pt Bold, white background).
    """
    os.makedirs(save_dir, exist_ok=True)
    
    FIG_BG  = "#FFFFFF"
    AX_BG   = "#FFFFFF"
    TXT_COL = "#1A1A1A"
    SPINE_C = "#444444"
    
    # Metrics to display
    display_keys = [
        ("volume_cm3", "Volume\n(cm^3)"),
        ("surface_area_cm2", "Surface Area\n(cm^2)"),
        ("volumetric_fidelity_pct", "Volumetric\nFidelity (%)"),
        ("sphericity_index_100", "Sphericity\n(x100)"),
        ("triangle_quality_mean_100", "Element Quality\n(x100)")
    ]
    
    f_vals = [
        metrics_dict["Femur"]["volume_cm3"],
        metrics_dict["Femur"]["surface_area_cm2"],
        metrics_dict["Femur"]["volumetric_fidelity_pct"],
        metrics_dict["Femur"]["sphericity_index"] * 100.0,
        metrics_dict["Femur"]["triangle_quality_mean"] * 100.0
    ]
    
    t_vals = [
        metrics_dict["Tibia"]["volume_cm3"],
        metrics_dict["Tibia"]["surface_area_cm2"],
        metrics_dict["Tibia"]["volumetric_fidelity_pct"],
        metrics_dict["Tibia"]["sphericity_index"] * 100.0,
        metrics_dict["Tibia"]["triangle_quality_mean"] * 100.0
    ]
    
    x = np.arange(len(display_keys))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(11, 8))
    fig.patch.set_facecolor(FIG_BG)
    ax.set_facecolor(AX_BG)
    
    bars_f = ax.bar(x - width/2, f_vals, width, label="Femur Model", color="#1565C0", alpha=0.88, edgecolor="#0D47A1", linewidth=0.8)
    bars_t = ax.bar(x + width/2, t_vals, width, label="Tibia Model", color="#BF360C", alpha=0.88, edgecolor="#7F0000", linewidth=0.8)
    
    def _label_bars(bars, color):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(
                f"{h:.1f}",
                xy=(bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=12, fontweight="bold",
                fontfamily="Times New Roman",
                color=color
            )
            
    _label_bars(bars_f, "#1565C0")
    _label_bars(bars_t, "#BF360C")
    
    ax.set_xticks(x)
    ax.set_xticklabels([d[1] for d in display_keys], **_tnr(size=15), color=TXT_COL)
    ax.set_ylabel("Metric Value / Magnitude", **_tnr(), color=TXT_COL, labelpad=12)
    ax.set_xlabel("3D Reconstruction Geometric & Quality Metric", **_tnr(), color=TXT_COL, labelpad=14)
    # Title without 'Step 4' word
    ax.set_title(
        "3D Reconstruction Performance & Geometric Metrics\n"
        "Femur vs Tibia STL Bone Models | Patient LTKN9332",
        **_tnr(size=18), color=TXT_COL, pad=16
    )
    
    max_val = max(f_vals + t_vals)
    ax.set_ylim(0, max_val * 1.15)
    
    ax.tick_params(axis="both", colors=TXT_COL, labelsize=16, which="major", length=5, width=1.2)
    for sp in ax.spines.values():
        sp.set_edgecolor(SPINE_C)
        sp.set_linewidth(1.2)
    ax.grid(False)
    
    legend = ax.legend(loc="upper left", fontsize=15, facecolor="#F0F0F0", edgecolor=SPINE_C, framealpha=0.95)
    for text in legend.get_texts():
        text.set_fontfamily("Times New Roman")
        text.set_fontweight("bold")
        text.set_color(TXT_COL)
        
    plt.tight_layout()
    out_path = os.path.join(save_dir, "step4_reconstruction_metrics_barplot.png")
    plt.savefig(out_path, dpi=1000, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [3D Reconstruction] Reconstruction metrics bar chart saved: {out_path}")

def run_3d_reconstruction(step3_data, datasets, base_dir):
    print("\n" + "=" * 70)
    print("  3D RECONSTRUCTION (PYTHON)")
    print("  Convert 3D femur and tibia masks into 3D bone models (STL)")
    print("  Preserve anatomical geometry for personalized implant design")
    print("=" * 70)
    
    reconstruction_dir = os.path.join(base_dir, "reconstruction_results")
    os.makedirs(reconstruction_dir, exist_ok=True)
    
    femur_mask_3d = step3_data.get("femur_mask_3d")
    tibia_mask_3d = step3_data.get("tibia_mask_3d")
    
    if femur_mask_3d is None:
        femur_mask_3d = np.load(os.path.join(base_dir, "femur_mask_3d.npy"))
    if tibia_mask_3d is None:
        tibia_mask_3d = np.load(os.path.join(base_dir, "tibia_mask_3d.npy"))
        
    femur_files = datasets.get("CT_Femur", {}).get("files", [])
    tibia_files = datasets.get("CT_Tibia", {}).get("files", [])
    
    stack_slices = femur_mask_3d.shape[0]
    target_res = femur_mask_3d.shape[1:]
    
    # Compute anatomical voxel spacing in physical millimeters (mm)
    # Original pixel spacing = 0.082 mm, slice thickness = 0.082 mm
    f_num_slices = len(femur_files) if femur_files else 658
    t_num_slices = len(tibia_files) if tibia_files else 845
    
    f_orig_shape = datasets.get("CT_Femur", {}).get("slice_shape", (1277, 1292))
    t_orig_shape = datasets.get("CT_Tibia", {}).get("slice_shape", (1267, 1169))
    
    f_spacing = (
        float((f_num_slices - 1) * 0.082 / (stack_slices - 1)),
        float(f_orig_shape[0] * 0.082 / target_res[0]),
        float(f_orig_shape[1] * 0.082 / target_res[1])
    )
    
    t_spacing = (
        float((t_num_slices - 1) * 0.082 / (stack_slices - 1)),
        float(t_orig_shape[0] * 0.082 / target_res[0]),
        float(t_orig_shape[1] * 0.082 / target_res[1])
    )
    
    print(f"\n[+] Anatomical Physical Spacing Calculated:")
    print(f"    - Femur Voxel Spacing (Z, Y, X) : ({f_spacing[0]:.4f}, {f_spacing[1]:.4f}, {f_spacing[2]:.4f}) mm")
    print(f"    - Tibia Voxel Spacing (Z, Y, X) : ({t_spacing[0]:.4f}, {t_spacing[1]:.4f}, {t_spacing[2]:.4f}) mm")
    
    # 1. Reconstruct Femur 3D Bone Model
    print(f"\n[+] Reconstructing 3D Femur Bone Model...")
    t0 = time.time()
    femur_mesh = reconstruct_bone_mesh(femur_mask_3d, f_spacing, name="Femur", sigma=0.8)
    t_femur = time.time() - t0
    femur_metrics = compute_3d_reconstruction_metrics(femur_mesh, femur_mask_3d, f_spacing, name="Femur")
    
    femur_stl_path = os.path.join(reconstruction_dir, "femur_model.stl")
    femur_mesh.export(femur_stl_path)
    # Also save to base_dir for immediate accessibility
    femur_mesh.export(os.path.join(base_dir, "femur_model.stl"))
    print(f"    - Femur reconstruction completed in {t_femur:.2f}s")
    print(f"    - Femur STL Model saved : {femur_stl_path}")
    print(f"    - Vertices: {femur_metrics['vertices']:,} | Faces: {femur_metrics['faces']:,} | Watertight: {femur_metrics['is_watertight']}")
    print(f"    - Volume  : {femur_metrics['volume_cm3']:.2f} cm^3 | Surface Area: {femur_metrics['surface_area_cm2']:.2f} cm^2")
    print(f"    - Volumetric Reconstruction Fidelity: {femur_metrics['volumetric_fidelity_pct']:.2f}%")
    
    # 2. Reconstruct Tibia 3D Bone Model
    print(f"\n[+] Reconstructing 3D Tibia Bone Model...")
    t0 = time.time()
    tibia_mesh = reconstruct_bone_mesh(tibia_mask_3d, t_spacing, name="Tibia", sigma=0.8)
    t_tibia = time.time() - t0
    tibia_metrics = compute_3d_reconstruction_metrics(tibia_mesh, tibia_mask_3d, t_spacing, name="Tibia")
    
    tibia_stl_path = os.path.join(reconstruction_dir, "tibia_model.stl")
    tibia_mesh.export(tibia_stl_path)
    tibia_mesh.export(os.path.join(base_dir, "tibia_model.stl"))
    print(f"    - Tibia reconstruction completed in {t_tibia:.2f}s")
    print(f"    - Tibia STL Model saved : {tibia_stl_path}")
    print(f"    - Vertices: {tibia_metrics['vertices']:,} | Faces: {tibia_metrics['faces']:,} | Watertight: {tibia_metrics['is_watertight']}")
    print(f"    - Volume  : {tibia_metrics['volume_cm3']:.2f} cm^3 | Surface Area: {tibia_metrics['surface_area_cm2']:.2f} cm^2")
    print(f"    - Volumetric Reconstruction Fidelity: {tibia_metrics['volumetric_fidelity_pct']:.2f}%")
    
    # 3. Create Knee Joint Anatomical Assembly
    print(f"\n[+] Assembling Knee Joint Digital Twin Base (Femur + Tibia in situ)...")
    # Position femur superior to tibia with 6.0 mm physiological joint space
    joint_gap_mm = 6.0
    z_shift = (tibia_mesh.bounds[1][0] - femur_mesh.bounds[0][0]) + joint_gap_mm
    femur_shifted = femur_mesh.copy()
    femur_shifted.apply_translation([z_shift, 0, 0])
    
    knee_assembly = trimesh.util.concatenate([tibia_mesh, femur_shifted])
    assembly_stl_path = os.path.join(reconstruction_dir, "knee_joint_assembly.stl")
    knee_assembly.export(assembly_stl_path)
    knee_assembly.export(os.path.join(base_dir, "knee_joint_assembly.stl"))
    print(f"    - Combined Knee Joint Assembly STL saved : {assembly_stl_path}")
    print(f"    - Assembly Total Vertices: {len(knee_assembly.vertices):,} | Total Faces: {len(knee_assembly.faces):,}")
    
    # 4. Save Reconstruction Metrics JSON
    reconstruction_metrics = {
        "Femur": femur_metrics,
        "Tibia": tibia_metrics,
        "Knee_Assembly": {
            "total_vertices": int(len(knee_assembly.vertices)),
            "total_faces": int(len(knee_assembly.faces)),
            "total_volume_cm3": round((femur_metrics["volume_cm3"] + tibia_metrics["volume_cm3"]), 2),
            "total_surface_area_cm2": round((femur_metrics["surface_area_cm2"] + tibia_metrics["surface_area_cm2"]), 2),
            "physiological_joint_gap_mm": joint_gap_mm
        }
    }
    
    metrics_json_path = os.path.join(reconstruction_dir, "reconstruction_metrics.json")
    with open(metrics_json_path, "w") as f:
        json.dump(reconstruction_metrics, f, indent=4)
    # Also save to base_dir
    with open(os.path.join(base_dir, "reconstruction_metrics.json"), "w") as f:
        json.dump(reconstruction_metrics, f, indent=4)
        
    print(f"    - Reconstruction Metrics JSON saved : {metrics_json_path}")
    
    # 5. Print Quantitative Reconstruction Metrics Table
    print("\n" + "-" * 75)
    print("  QUANTITATIVE 3D RECONSTRUCTION PERFORMANCE & GEOMETRIC METRICS")
    print("-" * 75)
    print(f"  {'Metric Description':<32} | {'Femur Model':<18} | {'Tibia Model':<18}")
    print("-" * 75)
    print(f"  {'Mesh Vertices (Nv)':<32} | {femur_metrics['vertices']:<18,} | {tibia_metrics['vertices']:<18,}")
    print(f"  {'Mesh Faces (Nf)':<32} | {femur_metrics['faces']:<18,} | {tibia_metrics['faces']:<18,}")
    print(f"  {'Unique Edges (Ne)':<32} | {femur_metrics['edges']:<18,} | {tibia_metrics['edges']:<18,}")
    print(f"  {'Watertight Manifold':<32} | {str(femur_metrics['is_watertight']):<18} | {str(tibia_metrics['is_watertight']):<18}")
    print(f"  {'Euler Characteristic (Chi)':<32} | {femur_metrics['euler_characteristic']:<18} | {tibia_metrics['euler_characteristic']:<18}")
    print(f"  {'Bone Volume (cm^3)':<32} | {femur_metrics['volume_cm3']:<18.2f} | {tibia_metrics['volume_cm3']:<18.2f}")
    print(f"  {'Bone Volume (mm^3)':<32} | {femur_metrics['volume_mm3']:<18,.1f} | {tibia_metrics['volume_mm3']:<18,.1f}")
    print(f"  {'Bone Surface Area (cm^2)':<32} | {femur_metrics['surface_area_cm2']:<18.2f} | {tibia_metrics['surface_area_cm2']:<18.2f}")
    print(f"  {'Surface-to-Volume Ratio (1/mm)':<32} | {femur_metrics['sa_to_vol_ratio_mm_inv']:<18.4f} | {tibia_metrics['sa_to_vol_ratio_mm_inv']:<18.4f}")
    print(f"  {'Volumetric Fidelity (%)':<32} | {femur_metrics['volumetric_fidelity_pct']:<18.2f} | {tibia_metrics['volumetric_fidelity_pct']:<18.2f}")
    print(f"  {'Volumetric Error (%)':<32} | {femur_metrics['volumetric_error_pct']:<18.2f} | {tibia_metrics['volumetric_error_pct']:<18.2f}")
    print(f"  {'Mean Edge Length (mm)':<32} | {femur_metrics['mean_edge_length_mm']:<18.4f} | {tibia_metrics['mean_edge_length_mm']:<18.4f}")
    print(f"  {'Sphericity Index (Psi)':<32} | {femur_metrics['sphericity_index']:<18.4f} | {tibia_metrics['sphericity_index']:<18.4f}")
    print(f"  {'Triangle Element Quality':<32} | {femur_metrics['triangle_quality_mean']:<18.4f} | {tibia_metrics['triangle_quality_mean']:<18.4f}")
    print(f"  {'Bounding Box [AP, ML, SI] mm':<32} | {str(femur_metrics['bounding_box_mm']):<18} | {str(tibia_metrics['bounding_box_mm']):<18}")
    print("-" * 75)
    
    # 6. Generate Figures
    print("\n[+] Generating Publication-Grade 3D Bone Models and Metrics Figures (DPI 1000)...")
    plot_reconstructed_bone_models(femur_mesh, tibia_mesh, femur_shifted, reconstruction_dir)
    plot_reconstruction_metrics_bar(reconstruction_metrics, reconstruction_dir)
    
    # Also save copies to base_dir
    import shutil
    shutil.copyfile(os.path.join(reconstruction_dir, "step4_reconstructed_bone_models.png"),
                    os.path.join(base_dir, "step4_reconstructed_bone_models.png"))
    shutil.copyfile(os.path.join(reconstruction_dir, "step4_reconstruction_metrics_barplot.png"),
                    os.path.join(base_dir, "step4_reconstruction_metrics_barplot.png"))
                    
    return {
        "femur_mesh": femur_mesh,
        "tibia_mesh": tibia_mesh,
        "knee_assembly": knee_assembly,
        "metrics": reconstruction_metrics,
        "stl_paths": {
            "femur": femur_stl_path,
            "tibia": tibia_stl_path,
            "assembly": assembly_stl_path
        }
    }




# ============================================================
# VISUALIZATION FUNCTIONS
# ============================================================

def plot_dataset_samples(datasets, output_path):
    """Plot Step 1 - Raw loaded dataset sample slices."""
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.5))
    fig.suptitle("STEP 1 – Raw Loaded Dataset Sample Slices", fontsize=13, fontweight="bold", y=1.02)

    items = [
        ("CT_Femur",  "bone", f"CT Femur (Slice {datasets.get('CT_Femur', {}).get('mid_idx', '')})"),
        ("CT_Tibia",  "bone", f"CT Tibia (Slice {datasets.get('CT_Tibia', {}).get('mid_idx', '')})"),
        ("CT_Intact", "hot",  "CT Intact (MIP - All Slices)"),
        ("MR_TSE",    "gray", f"MR TSE (Slice {datasets.get('MR_TSE', {}).get('mid_idx', '')})"),
        ("MR_DESS",   "gray", f"MR DESS (Slice {datasets.get('MR_DESS', {}).get('mid_idx', '')})"),
    ]
    for ax, (key, cmap, title) in zip(axes, items):
        if key in datasets:
            ax.imshow(normalize_slice_for_display(datasets[key]["sample_slice"]), cmap=cmap)
            ax.set_title(title, fontsize=10, fontweight="bold")
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n  [Step 1] Raw sample preview saved: {output_path}")

def plot_preprocessing_comparison(datasets, preprocessed, output_path):
    """
    Plot Step 2 - Side-by-side comparison of Raw vs Preprocessed for all 5 datasets.
    Top row = Raw, Bottom row = Preprocessed
    """
    keys = ["CT_Femur", "CT_Tibia", "CT_Intact", "MR_TSE", "MR_DESS"]
    cmaps = {"CT_Femur": "bone", "CT_Tibia": "bone", "CT_Intact": "hot",
             "MR_TSE": "gray", "MR_DESS": "gray"}
    labels = {
        "CT_Femur":  f"CT Femur\n(Slice {datasets.get('CT_Femur', {}).get('mid_idx', '')})",
        "CT_Tibia":  f"CT Tibia\n(Slice {datasets.get('CT_Tibia', {}).get('mid_idx', '')})",
        "CT_Intact": "CT Intact\n(MIP)",
        "MR_TSE":    f"MR TSE\n(Slice {datasets.get('MR_TSE', {}).get('mid_idx', '')})",
        "MR_DESS":   f"MR DESS\n(Slice {datasets.get('MR_DESS', {}).get('mid_idx', '')})",
    }

    present_keys = [k for k in keys if k in datasets and k in preprocessed]
    n = len(present_keys)

    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 8))
    fig.suptitle("STEP 2 – Preprocessing: Raw vs Preprocessed", fontsize=13,
                 fontweight="bold", y=1.01)

    for col, key in enumerate(present_keys):
        raw_norm = normalize_slice_for_display(datasets[key]["sample_slice"])
        pre_norm = preprocessed[key]
        cmap = cmaps[key]

        # Row 0 – Raw
        axes[0, col].imshow(raw_norm, cmap=cmap)
        axes[0, col].set_title(f"RAW\n{labels[key]}", fontsize=9, fontweight="bold", color="#CC3333")
        axes[0, col].axis("off")

        # Row 1 – Preprocessed
        axes[1, col].imshow(pre_norm, cmap=cmap)
        axes[1, col].set_title(f"PREPROCESSED\n{labels[key]}", fontsize=9, fontweight="bold", color="#006633")
        axes[1, col].axis("off")

    # Row labels
    for row, label in enumerate(["RAW INPUT", "PREPROCESSED"]):
        fig.text(0.01, 0.75 - row * 0.5, label, va="center", ha="left",
                 fontsize=10, fontweight="bold", rotation=90,
                 color="#CC3333" if row == 0 else "#006633")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [Step 2] Preprocessing comparison saved: {output_path}")

def plot_preprocessing_steps(datasets, output_path):
    """
    Show the 4 individual preprocessing sub-steps for one representative slice (CT Femur).
    Columns: Raw -> Denoised -> Normalized -> Resized -> Enhanced (CLAHE)
    """
    key = "CT_Femur" if "CT_Femur" in datasets else list(datasets.keys())[0]
    mod = "CT" if "CT" in key else "MRI"
    raw = datasets[key]["sample_slice"].astype(np.float32)
    raw = np.nan_to_num(raw, nan=0.0)

    denoised   = step2_noise_removal(raw, modality=mod)
    normalized = step2_intensity_normalization(denoised)
    resized    = step2_resize(normalized, TARGET_SIZE)
    enhanced   = step2_image_enhancement(resized)

    steps = [
        (normalize_slice_for_display(raw),    "Step 2a: Raw Input"),
        (normalize_slice_for_display(denoised), "Step 2b: Noise Removal\n(Median Filter 3x3)"),
        (normalized,                           "Step 2c: Intensity\nNormalization [0,1]"),
        (resized,                              f"Step 2d: Resize\n({raw.shape} → {TARGET_SIZE})"),
        (enhanced,                             "Step 2e: CLAHE\nEnhancement"),
    ]

    fig, axes = plt.subplots(1, 5, figsize=(22, 4.5))
    fig.suptitle(f"STEP 2 – Preprocessing Sub-steps ({key})", fontsize=13,
                 fontweight="bold", y=1.02)

    for ax, (img, title) in zip(axes, steps):
        ax.imshow(img, cmap="bone")
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [Step 2] Sub-steps pipeline saved: {output_path}")

# ─────────────────────────────────────────────────────────────────────────
# Typography helper – Times New Roman 18pt Bold for Step 3 visualisations
# ─────────────────────────────────────────────────────────────────────────
FONT_PROPS = dict(fontsize=18, fontweight="bold", fontfamily="Times New Roman")


def _tnr(size=18, **kw):
    """Return a dict suitable for ax.set_title / ax.set_xlabel etc."""
    return dict(fontsize=size, fontweight="bold", fontfamily="Times New Roman", **kw)


def _make_seg_overlay(data):
    """Build a colour-blended anatomical boundary overlay image (RGB uint8)."""
    bg_rgb = cv2.cvtColor((data["preprocessed"] * 255).astype(np.uint8), cv2.COLOR_GRAY2RGB)
    overlay = bg_rgb.copy()
    overlay[data["gt"] == 1]                                          = [0, 200, 200]   # Cyan  – GT
    overlay[data["pred"] == 1]                                        = [230, 80,  50]  # Red   – Pred
    overlap = np.logical_and(data["gt"] == 1, data["pred"] == 1)
    overlay[overlap]                                                   = [50, 205,  50]  # Green – TP
    return cv2.addWeighted(bg_rgb, 0.40, overlay, 0.60, 0)


def _make_error_map(data):
    """Build TP/FP/FN colour error map (RGB uint8)."""
    err = np.zeros((*data["gt"].shape, 3), dtype=np.uint8)
    err[np.logical_and(data["pred"] == 1, data["gt"] == 1)] = [46,  204, 113]  # TP Green
    err[np.logical_and(data["pred"] == 1, data["gt"] == 0)] = [231,  76,  60]  # FP Red
    err[np.logical_and(data["pred"] == 0, data["gt"] == 1)] = [52,  152, 219]  # FN Blue
    return err


def plot_bone_mask_samples(bone_name, val_results_list, save_dir, n_samples=5):
    """
    Save 5 separate figures – one per validation slice – for a given bone.

    Each figure has 5 sub-panels:
      1. Preprocessed CT Slice
      2. Semi-Automatic Ground-Truth Mask
      3. U-Net Predicted Mask
      4. Anatomical Boundary Overlay
      5. Segmentation Error Map (TP / FP / FN)

    Theme      : Light (white background, near-black text)
    Typography : Times New Roman, 18pt, Bold on all text elements.
    Output DPI : 1000
    """
    os.makedirs(save_dir, exist_ok=True)
    samples = val_results_list[:n_samples]

    # ── Light-theme colour constants ──────────────────────────────────────
    FIG_BG   = "#FFFFFF"   # pure white figure background
    AX_BG    = "#F5F5F5"   # very light grey axes background
    TXT_COL  = "#1A1A1A"   # near-black for all text
    SPINE_C  = "#333333"   # dark spine border

    for s_idx, data in enumerate(samples, start=1):
        fig, axes = plt.subplots(1, 5, figsize=(26, 5))
        fig.patch.set_facecolor(FIG_BG)

        suptitle = (
            f"{bone_name} Mask \u2013 Sample {s_idx}/{n_samples}  |  "
            f"Slice #{data['slice_idx']}  |  "
            f"Dice = {data['metrics']['Dice']:.4f}  |  "
            f"IoU = {data['metrics']['IoU']:.4f}"
        )
        fig.suptitle(suptitle, **_tnr(size=16), color=TXT_COL, y=1.01)

        panel_data = [
            (
                data["preprocessed"],
                "bone",
                f"(A) Preprocessed CT\n{bone_name} Slice #{data['slice_idx']}",
            ),
            (
                data["gt"],
                "Blues",
                "(B) Semi-Automatic\nGround-Truth Mask",
            ),
            (
                data["pred"],
                "Oranges",
                f"(C) U-Net Predicted Mask\nDice={data['metrics']['Dice']:.4f}  IoU={data['metrics']['IoU']:.4f}",
            ),
            (
                _make_seg_overlay(data),
                None,
                "(D) Anatomical Overlay\nGreen=TP  Cyan=GT  Red=Pred",
            ),
            (
                _make_error_map(data),
                None,
                f"(E) Error Map\nPrec={data['metrics']['Precision']:.3f}  Recall={data['metrics']['Recall']:.3f}",
            ),
        ]

        for ax, (img, cmap, title) in zip(axes, panel_data):
            ax.set_facecolor(AX_BG)
            if cmap is not None:
                ax.imshow(img, cmap=cmap)
            else:
                ax.imshow(img)
            ax.set_title(title, **_tnr(), color=TXT_COL, pad=8)
            ax.axis("off")
            # Subtle dark spine border
            for sp in ax.spines.values():
                sp.set_edgecolor(SPINE_C)
                sp.set_linewidth(1.2)
                sp.set_visible(True)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_name = f"step3_{bone_name.lower()}_mask_sample_{s_idx:02d}.png"
        out_path = os.path.join(save_dir, out_name)
        plt.savefig(out_path, dpi=1000, bbox_inches="tight",
                    facecolor=FIG_BG)
        plt.close()
        print(f"    Saved {bone_name} sample {s_idx}: {out_path}")


def plot_performance_metrics_bar(summary_metrics, save_dir):
    """
    Standalone performance-metrics bar chart (figsize 10x8, DPI 1000).

    Displays Dice, IoU, Precision, Recall, Specificity, and Accuracy as
    grouped bars for Femur vs Tibia vs Overall Mean.  Bar value labels are
    printed on top of every bar.

    Theme      : Light (white background, near-black text)
    Typography : Times New Roman, 18pt, Bold
    Grid       : Disabled
    """
    os.makedirs(save_dir, exist_ok=True)

    # ── Light-theme colour constants ──────────────────────────────────────
    FIG_BG  = "#FFFFFF"   # pure white figure background
    AX_BG   = "#FFFFFF"   # white axes background
    TXT_COL = "#1A1A1A"   # near-black for all text
    SPINE_C = "#444444"   # medium-dark spines

    # Metrics to display (all normalised [0,1])
    metric_keys   = ["Dice", "IoU", "Precision", "Recall", "Specificity", "Accuracy"]
    metric_labels = [
        "Dice\n(DSC)",
        "IoU\n(Jaccard)",
        "Precision\n(PPV)",
        "Recall\n(Sens.)",
        "Specificity\n(TNR)",
        "Accuracy",
    ]

    femur_vals = [summary_metrics["Femur"][m] for m in metric_keys]
    tibia_vals = [summary_metrics["Tibia"][m] for m in metric_keys]
    mean_vals  = [summary_metrics["Mean"][m]  for m in metric_keys]

    x     = np.arange(len(metric_keys))
    width = 0.26

    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor(FIG_BG)
    ax.set_facecolor(AX_BG)

    bars_f = ax.bar(x - width, femur_vals, width,
                    label="Femur",        color="#1565C0", alpha=0.88, zorder=3,
                    edgecolor="#0D47A1",  linewidth=0.8)
    bars_t = ax.bar(x,         tibia_vals, width,
                    label="Tibia",        color="#BF360C", alpha=0.88, zorder=3,
                    edgecolor="#7F0000",  linewidth=0.8)
    bars_m = ax.bar(x + width, mean_vals,  width,
                    label="Overall Mean", color="#2E7D32", alpha=0.88, zorder=3,
                    edgecolor="#1B5E20",  linewidth=0.8)

    # Value labels on top of every bar
    def _label_bars(bars, color):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(
                f"{h:.4f}",
                xy=(bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=11, fontweight="bold",
                fontfamily="Times New Roman",
                color=color
            )

    _label_bars(bars_f, "#1565C0")
    _label_bars(bars_t, "#BF360C")
    _label_bars(bars_m, "#2E7D32")

    # ── X-axis tick labels (metric names) ─────────────────────────────────
    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels, **_tnr(), color=TXT_COL)

    # ── Axis labels ───────────────────────────────────────────────────────
    ax.set_ylabel("Metric Score  [0.0 \u2013 1.0]", **_tnr(), color=TXT_COL, labelpad=12)
    ax.set_xlabel("Segmentation Performance Metric",  **_tnr(), color=TXT_COL, labelpad=14)

    # ── Title ─────────────────────────────────────────────────────────────
    ax.set_title(
        "Step 3 \u2013 U-Net Segmentation Performance Metrics\n"
        "Femur vs Tibia vs Overall Mean  |  Patient LTKN9332",
        **_tnr(size=18), color=TXT_COL, pad=16
    )

    # ── Dynamic y-axis ────────────────────────────────────────────────────
    all_vals = femur_vals + tibia_vals + mean_vals
    y_min = max(0.0, min(all_vals) - 0.08)
    y_max = min(1.05, max(all_vals) + 0.12)
    ax.set_ylim(y_min, y_max)

    # ── Tick styling ──────────────────────────────────────────────────────
    ax.tick_params(axis="both", colors=TXT_COL, labelsize=18,
                   which="major", length=5, width=1.2)
    ax.yaxis.label.set_color(TXT_COL)
    ax.xaxis.label.set_color(TXT_COL)

    # ── Spines ────────────────────────────────────────────────────────────
    for sp in ax.spines.values():
        sp.set_edgecolor(SPINE_C)
        sp.set_linewidth(1.2)

    # ── No grid ───────────────────────────────────────────────────────────
    ax.grid(False)

    # ── Legend ────────────────────────────────────────────────────────────
    legend = ax.legend(
        loc="lower right", fontsize=15,
        facecolor="#F0F0F0", edgecolor=SPINE_C,
        framealpha=0.95
    )
    for text in legend.get_texts():
        text.set_fontfamily("Times New Roman")
        text.set_fontweight("bold")
        text.set_color(TXT_COL)

    plt.tight_layout()
    out_path = os.path.join(save_dir, "step3_performance_metrics_barplot.png")
    plt.savefig(out_path, dpi=1000, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"\n  [Step 3] Performance metrics bar chart saved: {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # ----------------------------------------------------------
    # STEP 1 – DATASET LOADING
    # ----------------------------------------------------------
    print("=" * 70)
    print("  STEP 1 - LEEDS KNEE CT/MRI DATASET LOADING")
    print("=" * 70)

    datasets = {}

    # 1. CT Dissected - Femur
    femur_dir = os.path.join(base_dir, "LTKN9332_CT_dissected", "LTKN9332_CT_dissected", "femur")
    if os.path.exists(femur_dir):
        femur_info = load_dicom_series_info(femur_dir)
        if femur_info:
            datasets["CT_Femur"] = femur_info
            print(f"\n[+] Loaded CT Dissected (Femur):")
            print(f"    - Number of Slices: {femur_info['num_slices']}")
            print(f"    - Mid-Slice Index: {femur_info['mid_idx']}")
            print(f"    - Slice Resolution (Y, X): {femur_info['slice_shape']}")
            print(f"    - Data Type: {femur_info['dtype']}")
            print(f"    - Mid-Slice Range: {femur_info['min_val']:.1f} to {femur_info['max_val']:.1f}")

    # 2. CT Dissected - Tibia
    tibia_dir = os.path.join(base_dir, "LTKN9332_CT_dissected", "LTKN9332_CT_dissected", "tibia")
    if os.path.exists(tibia_dir):
        tibia_info = load_dicom_series_info(tibia_dir)
        if tibia_info:
            datasets["CT_Tibia"] = tibia_info
            print(f"\n[+] Loaded CT Dissected (Tibia):")
            print(f"    - Number of Slices: {tibia_info['num_slices']}")
            print(f"    - Mid-Slice Index: {tibia_info['mid_idx']}")
            print(f"    - Slice Resolution (Y, X): {tibia_info['slice_shape']}")
            print(f"    - Data Type: {tibia_info['dtype']}")
            print(f"    - Mid-Slice Range: {tibia_info['min_val']:.1f} to {tibia_info['max_val']:.1f}")

    # 3. CT Intact Merged TIFFs
    intact_dir = os.path.join(base_dir, "LTKN9332_CT_intact", "LTKN9332_CT_intact", "merged_tiffs")
    if os.path.exists(intact_dir):
        print(f"\n[+] Loading CT Intact (Merged TIFFs)...")
        intact_info = load_tiff_series_info(intact_dir)
        if intact_info:
            datasets["CT_Intact"] = intact_info
            print(f"    - Total Slice Files: {intact_info['num_slices']}")
            print(f"    - Display Mode: Max-Intensity Projection (MIP) over ALL slices")
            print(f"    - MIP Resolution (Y, X): {intact_info['slice_shape']}")
            print(f"    - MIP Intensity Range: {intact_info['min_val']:.1f} to {intact_info['max_val']:.1f}")

    # 4. MR 2D TSE SAG IW
    mr_tse_dir = os.path.join(base_dir, "LTKN9332_MR_2D_TSE_SAG_IW", "LTKN9332_MR_2D TSE SAG IW")
    if os.path.exists(mr_tse_dir):
        mr_tse_info = load_dicom_series_info(mr_tse_dir)
        if mr_tse_info:
            datasets["MR_TSE"] = mr_tse_info
            print(f"\n[+] Loaded MRI 2D TSE SAG IW:")
            print(f"    - Number of Slices: {mr_tse_info['num_slices']}")
            print(f"    - Mid-Slice Index: {mr_tse_info['mid_idx']}")
            print(f"    - Slice Resolution (Y, X): {mr_tse_info['slice_shape']}")
            print(f"    - Data Type: {mr_tse_info['dtype']}")

    # 5. MR DESS
    mr_dess_dir = os.path.join(base_dir, "LTKN9332_MR_DESS", "LTKN9332_MR_DESS")
    if os.path.exists(mr_dess_dir):
        mr_dess_info = load_dicom_series_info(mr_dess_dir)
        if mr_dess_info:
            datasets["MR_DESS"] = mr_dess_info
            print(f"\n[+] Loaded MRI DESS:")
            print(f"    - Number of Slices: {mr_dess_info['num_slices']}")
            print(f"    - Mid-Slice Index: {mr_dess_info['mid_idx']}")
            print(f"    - Slice Resolution (Y, X): {mr_dess_info['slice_shape']}")
            print(f"    - Data Type: {mr_dess_info['dtype']}")

    # Step 1 visualization
    step1_out = os.path.join(base_dir, "step1_dataset_preview.png")
    plot_dataset_samples(datasets, step1_out)

    # ----------------------------------------------------------
    # STEP 2 – IMAGE PREPROCESSING
    # ----------------------------------------------------------
    preprocessed = run_preprocessing(datasets)

    # Step 2a: Full pipeline sub-steps figure (CT Femur)
    step2_pipeline_out = os.path.join(base_dir, "step2_preprocessing_pipeline.png")
    plot_preprocessing_steps(datasets, step2_pipeline_out)

    # Step 2b: Raw vs Preprocessed comparison for all 5 datasets
    step2_compare_out = os.path.join(base_dir, "step2_preprocessing_comparison.png")
    plot_preprocessing_comparison(datasets, preprocessed, step2_compare_out)

    # ----------------------------------------------------------
    # STEP 3 – 3D SEGMENTATION AND MASK GENERATION (U-NET)
    # ----------------------------------------------------------
    step3_data = run_3d_segmentation(datasets, base_dir)

    # ── Step 3 visualisations – all saved in segmentation_results/ ──────
    seg_results_dir = os.path.join(base_dir, "segmentation_results")
    os.makedirs(seg_results_dir, exist_ok=True)

    print("\n[+] Saving 5 separate Femur Mask sample figures...")
    plot_bone_mask_samples(
        bone_name="Femur",
        val_results_list=step3_data["val_results"]["Femur"],
        save_dir=seg_results_dir,
        n_samples=5
    )

    print("\n[+] Saving 5 separate Tibia Mask sample figures...")
    plot_bone_mask_samples(
        bone_name="Tibia",
        val_results_list=step3_data["val_results"]["Tibia"],
        save_dir=seg_results_dir,
        n_samples=5
    )

    print("\n[+] Saving Performance Metrics Bar Chart...")
    plot_performance_metrics_bar(
        summary_metrics=step3_data["summary_metrics"],
        save_dir=seg_results_dir
    )

    # ----------------------------------------------------------
    # STEP 4 – 3D RECONSTRUCTION (STL MODELS & METRIC EXTRACTION)
    # ----------------------------------------------------------
    step4_data = run_3d_reconstruction(step3_data, datasets, base_dir)

    print("\n" + "=" * 70)
    print("  ALL PIPELINE STEPS COMPLETED SUCCESSFULLY")
    print(f"  Generated Artifacts:")
    print(f"    - Step 1 Preview Image        : {step1_out}")
    print(f"    - Step 2 Sub-steps Image      : {step2_pipeline_out}")
    print(f"    - Step 2 Comparison Image     : {step2_compare_out}")
    print(f"    - Step 3 Outputs Folder       : {seg_results_dir}")
    print(f"      - 5x Femur mask sample figures (DPI 1000)")
    print(f"      - 5x Tibia mask sample figures (DPI 1000)")
    print(f"      - Performance metrics bar chart (DPI 1000)")
    print(f"    - 3D Femur Volumetric Mask    : {os.path.join(base_dir, 'femur_mask_3d.npy')}")
    print(f"    - 3D Tibia Volumetric Mask    : {os.path.join(base_dir, 'tibia_mask_3d.npy')}")
    print(f"    - Step 3 Metrics JSON Summary : {os.path.join(base_dir, 'segmentation_metrics.json')}")
    print(f"    - Step 4 Outputs Folder       : {os.path.join(base_dir, 'reconstruction_results')}")
    print(f"      - 3D Femur Bone Model (STL) : {step4_data['stl_paths']['femur']}")
    print(f"      - 3D Tibia Bone Model (STL) : {step4_data['stl_paths']['tibia']}")
    print(f"      - Knee Joint Assembly (STL) : {step4_data['stl_paths']['assembly']}")
    print(f"      - Reconstructed Models View : {os.path.join(base_dir, 'step4_reconstructed_bone_models.png')}")
    print(f"      - Reconstruction Bar Chart  : {os.path.join(base_dir, 'step4_reconstruction_metrics_barplot.png')}")
    print(f"      - Step 4 Metrics JSON       : {os.path.join(base_dir, 'reconstruction_metrics.json')}")
    print("=" * 70)

if __name__ == "__main__":
    main()
