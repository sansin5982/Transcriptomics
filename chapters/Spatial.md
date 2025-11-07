# Spatial Transcriptomics

> **“In biology, where a gene is expressed is as important as what is
> expressed.”**

## 1. Introduction to Spatial Transcriptomics

### 1.1 Definition and Historical Context

-   **Spatial transcriptomics (ST)**: Technology that measures gene
    expression while preserving the native spatial coordinates of RNA
    molecules in tissue.
-   **First breakthrough: 2016 – Ståhl et al.** introduced the **first
    array-based ST** using barcoded spots \[1\].
-   **Current Status**: Over 50 commercial platforms, **&gt;10,000
    publications**, and **routine clinical** use in oncology.

### 1.2 Importance in Preserving Tissue Architecture

-   Bulk RNA-Seq → **averages**
-   scRNA-Seq → **dissociates**
-   **ST → reconstructs tissue context**
-   Enables study of **cell-cell interactions, gradients, niches**, and
    **microenvironments**.

### 1.3 Comparison with Bulk and Single-Cell RNA-Seq

<table>
<thead>
<tr>
<th style="text-align: left;">Feature</th>
<th style="text-align: left;">Bulk</th>
<th style="text-align: left;">scRNA-Seq</th>
<th style="text-align: left;">Spatial</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Resolution</td>
<td style="text-align: left;">Tissue</td>
<td style="text-align: left;">scRNA-Seq</td>
<td style="text-align: left;">Spatial</td>
</tr>
<tr>
<td style="text-align: left;">Spatial Info</td>
<td style="text-align: left;">Lost</td>
<td style="text-align: left;">Lost</td>
<td style="text-align: left;">Preserved</td>
</tr>
<tr>
<td style="text-align: left;">Throughput</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">Variable</td>
</tr>
<tr>
<td style="text-align: left;">Cost</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">Medium-High</td>
</tr>
</tbody>
</table>

## 2. Principles of Spatially Resolved Gene Expression

### 2.1 Spatial Barcoding and Coordinate Systems

-   Each **location** in tissue gets a **unique barcode** (DNA
    sequence).
-   RNA captured → tagged with barcode → sequenced → **mapped back to
    (x,y)**.

### 2.2 In Situ vs. Capture-Based Approaches

<table>
<thead>
<tr>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">Example</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">In situ</td>
<td style="text-align: left;">Direct imaging in tissue</td>
<td style="text-align: left;">MERFISH, seqFISH</td>
</tr>
<tr>
<td style="text-align: left;">Capture-based</td>
<td style="text-align: left;">RNA diffuses to barcoded array</td>
<td style="text-align: left;">Visium, Slide-seq</td>
</tr>
</tbody>
</table>

### 2.3 Resolution Scales

<table>
<thead>
<tr>
<th style="text-align: left;">Scale</th>
<th style="text-align: left;">Size</th>
<th style="text-align: left;">Technology</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Subcellular</td>
<td style="text-align: left;">&lt;1 µm</td>
<td style="text-align: left;">Stereo-seq, NanoString</td>
</tr>
<tr>
<td style="text-align: left;">Cellular</td>
<td style="text-align: left;">1–10 µm</td>
<td style="text-align: left;">MERFISH, Xenium</td>
</tr>
<tr>
<td style="text-align: left;">Multicellular</td>
<td style="text-align: left;">50–100 µm</td>
<td style="text-align: left;">Visium</td>
</tr>
</tbody>
</table>

## Imaging-Based Spatial Technologies

### 3.1 MERFISH (Multiplexed Error-Robust FISH)

-   **Principle**: Sequential hybridization + error-correcting barcodes
-   **Genes**: Up to 1,000
-   **Resolution**: Single-molecule
-   MERFISH v2 → **10,000 genes, 3D brain mapping \[3\]**

### 3.2 seqFISH and seqFISH+

-   **seqFISH+**: 10,000 genes in **single cells**
-   Used in **mouse embryo atlas** → revealed **Hox code gradients**
    \[4\]

### 3.3 STARmap and ExSeq

-   **STARmap**: 3D expansion + sequencing
-   **ExSeq**: Expansion sequencing for **proteins + RNA**

### 3.4 HCR Imaging

-   Signal amplification for **low-expressed genes**

## 4. Sequencing-Based Spatial Technologies

### 4.1 10x Genomics Visium Platform

-   **Spot size**: 55 µm (~10–20 cells)
-   **Genes**: Genome-wide
-   **2025: Visium HD → 2 µm spots**
-   **Clinical use**: Intraoperative tumor margin assessment

### 4.2 Slide-seq and Slide-seqV2

-   **Bead-based, 10 µm resolution**
-   **2025**: Used in **mouse brain connectomics**

### 4.3 Stereo-seq (BGI)

-   **DNA nanoball (DNB)** array
-   **Resolution: 500 nm (subcellular)**
-   **2025: Full mouse brain in 3D \[5\]**

### 4.4 HDST (High-Density Spatial Transcriptomics)

-   **2 µm beads, high throughput**

## 5. Emerging and Hybrid Technologies

### 5.1 Xenium (10x In Situ Sequencing)

-   **In situ sequencing** → no capture loss
-   **Panel**: 5,000 genes
-   **2025**: FFPE-compatible

### 5.2 GeoMX DSP (NanoString)

-   **ROI selection** + NGS or nCounter
-   **Proteins + RNA**

### 5.3 NanoString CosMx

-   **Single-cell, 6,000-plex**
-   **2025: Human cancer atlas**

### 5.4 Light-seq and DBiT

-   **Light-activated barcoding**
-   **DNA-barcoded imaging + sequencing**

## 6. Sample Preparation and Tissue Handling

### 6.1 Fresh-Frozen vs. FFPE

<table>
<thead>
<tr>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Pros</th>
<th style="text-align: left;">Cons</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Fresh-frozen</td>
<td style="text-align: left;">High RNA quality</td>
<td style="text-align: left;">Needs −80°C</td>
</tr>
<tr>
<td style="text-align: left;">FFPE</td>
<td style="text-align: left;">Archival</td>
<td style="text-align: left;">RNA degradation</td>
</tr>
</tbody>
</table>

### 6.2 Sectioning and Mounting

-   Cryosection: **10 µm**
-   Visium: **Capture area 6.5 × 6.5 mm**

### 6.3 Permeabilization

-   **Visium**: 12–18 min tissue optimization
-   **MERFISH**: Cell fixation + probe penetration

### 6.4 Quality Control

-   **H&E staining**
-   **RNA integrity (RIN)**
-   **Permeabilization test**

## 7. Data Generation and Preprocessing

### 7.1 Image Acquisition

-   **Brightfield + fluorescence**
-   **Fiducial markers** for alignment

### 7.2 Barcode Demultiplexing

    spaceranger count --id=sample --image=img.tif --transcriptome=refdata-gex-GRCh38-2020-A

### 7.3 Spatial Registration

-   Align **spots to histology**
-   **Louve Browser** (10x)

### 7.4 Noise Filtering

-   Remove **ambient RNA, low-quality spots**

## 8. Computational Analysis Pipeline

### 8.1 Spot-Based vs. Pixel-Based

<table>
<thead>
<tr>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Use</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Spot</td>
<td style="text-align: left;">Visium</td>
</tr>
<tr>
<td style="text-align: left;">Pixel</td>
<td style="text-align: left;">Stereo-seq, MERFISH</td>
</tr>
</tbody>
</table>

### 8.2 Cell Segmentation and Deconvolution

-   **BayesSpace, STdeconvolve**
-   **RCTD**: Map scRNA reference to spots

### 8.3 Spatial Clustering

-   **SpaGCN, SEDR**
-   Identify **domains** (tumor core vs. stroma)

### 8.4 Integration with scRNA

-   **Tangram, Cell2location**

## 9. Visualization of Spatial Data

### 9.1 2D and 3D Heatmaps

-   **Seurat, Giotto, Squidpy**

### 9.2 Interactive Tools

-   **Loupe Browser** (10x)
-   **Vitessce** (web-based)

### 9.3 Overlay with Histology

-   **H&E + gene expression**

### 9.4 VR/AR Interfaces

-   **2025: VR brain atlas** using Stereo-seq

## 10. Biological Applications

### 10.1 Developmental Biology

-   Mouse gastrulation → anterior-posterior gradients

### 10.2 Tumor Microenvironment

-   Breast cancer → immune exclusion zone at invasive front

### 10.3 Neuroanatomy

-   Human cortex → layer-specific markers

### 10.4 Host-Pathogen

-   TB lung → granuloma center vs. rim

## 11. Clinical and Translational Applications

### 11.1 Biomarker Discovery

-   PD-L1 spatial pattern predicts immunotherapy response

### 11.2 Digital Pathology

-   AI + Visium → automated tumor grading

### 11.3 Drug Response Mapping

-   Organoids → drug penetration gradients

### 11.4 Precision Medicine

-   2025: Intraoperative ST in glioma surgery → real-time margin

## 12. Challenges and Limitations

<table>
<thead>
<tr>
<th style="text-align: left;">Challenge</th>
<th style="text-align: left;">Current Solution</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Resolution vs. throughput</td>
<td style="text-align: left;">Hybrid platforms</td>
</tr>
<tr>
<td style="text-align: left;">FFPE RNA degradation</td>
<td style="text-align: left;">Xenium, DSP</td>
</tr>
<tr>
<td style="text-align: left;">Batch effects</td>
<td style="text-align: left;">ComBat-seq, Harmony</td>
</tr>
<tr>
<td style="text-align: left;">Data size (TB-scale)</td>
<td style="text-align: left;">Cloud computing</td>
</tr>
</tbody>
</table>

### Summary Table: Spatial Transcriptomics

<table>
<thead>
<tr>
<th style="text-align: left;">Technology</th>
<th style="text-align: left;">Resoultion</th>
<th style="text-align: left;">Genes</th>
<th style="text-align: left;">Tissue</th>
<th style="text-align: left;">Clinical Use</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Visium HD</td>
<td style="text-align: left;">2 µm</td>
<td style="text-align: left;">Genome-wide</td>
<td style="text-align: left;">Fresh/FFPE</td>
<td style="text-align: left;">Yes</td>
</tr>
<tr>
<td style="text-align: left;">Stereo-seq</td>
<td style="text-align: left;">500 nm</td>
<td style="text-align: left;">Genome-wide</td>
<td style="text-align: left;">Fresh</td>
<td style="text-align: left;">Research</td>
</tr>
<tr>
<td style="text-align: left;">MERFISH</td>
<td style="text-align: left;">Single-molecule</td>
<td style="text-align: left;">10,000</td>
<td style="text-align: left;">Fresh</td>
<td style="text-align: left;">Brain</td>
</tr>
<tr>
<td style="text-align: left;">Xenium</td>
<td style="text-align: left;">Single-cell</td>
<td style="text-align: left;">5,000</td>
<td style="text-align: left;">FFPE</td>
<td style="text-align: left;">Pathology</td>
</tr>
</tbody>
</table>
