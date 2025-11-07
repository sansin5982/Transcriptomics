# Spatial Transcriptomics

> **“In biology, where a gene is expressed is as important as what is
> expressed.”**

## 1. Introduction to Spatial Transcriptomics

Spatial transcriptomics (ST) is a technique that maps gene expression
profiles within the intact context of a tissue, preserving the spatial
organization of cells and their molecular states. Unlike traditional
transcriptomics methods, ST combines high-throughput RNA sequencing or
imaging with spatial coordinates, allowing researchers to visualize how
genes are expressed in specific locations and how cells interact within
their microenvironment. This is achieved by slicing tissue samples,
treating them to release RNA that binds to barcoded spots on a slide,
sequencing the barcoded RNA, and integrating this data with imaging to
associate transcripts with precise positions. The core principle builds
on single-cell genomics, using histological staining and microscopy to
link RNA data to cellular and tissue architecture.

#### Importance in Biology

ST is transformative for understanding multicellular organisms, where
spatial organization drives organ function, development, and disease. It
reveals cellular heterogeneity, cell-cell interactions, and
microenvironmental dynamics that are lost in dissociated samples. By
providing context to gene expression—such as how immune cells infiltrate
tumors or how injury propagates in tissues—ST bridges historical
histology with modern omics, accelerating discoveries in developmental
biology, immunology, and pathology. It is particularly vital for
studying complex diseases like cancer and neurodegeneration, informing
drug targeting, precision medicine, and tissue engineering.

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Key Concept</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Biological Relevance</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Spatial Coordinates</td>
<td style="text-align: left;">GPS-like mapping of transcripts to tissue
locations.</td>
<td style="text-align: left;">Enables neighborhood analysis (e.g.,
tumor-immune interfaces).</td>
</tr>
<tr>
<td style="text-align: left;">Transcript Capture</td>
<td style="text-align: left;">Barcoded probes or arrays bind RNA in
situ.</td>
<td style="text-align: left;">Captures whole transcriptome or targeted
genes without dissociation.</td>
</tr>
<tr>
<td style="text-align: left;">Multi-Modal Integration</td>
<td style="text-align: left;">Combines RNA with proteins,
epigenomics.</td>
<td style="text-align: left;">Provides holistic views of cellular
states.</td>
</tr>
</tbody>
</table>

## 2. Technologies and Methods

ST technologies are categorized into imaging-based (using fluorescence
in situ hybridization for targeted detection) and sequencing-based
(using barcoded arrays for unbiased profiling). Imaging-based methods
excel in high-resolution, targeted studies, while sequencing-based offer
broader coverage but lower resolution.

#### Imaging-Based Technologies

These rely on single-molecule FISH (smFISH) for subcellular precision,
suitable for hypothesis-driven research.

<table>
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Technology</th>
<th style="text-align: left;">Resolution</th>
<th style="text-align: left;">Throughput</th>
<th style="text-align: left;">Pros</th>
<th style="text-align: left;">Cons</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Xenium (10x Genomics)</td>
<td style="text-align: left;">~1 μm (subcellular)</td>
<td style="text-align: left;">Low (hours-days)</td>
<td style="text-align: left;">High sensitivity/specificity;
FFPE-compatible; up to 6,000 genes; custom panels.</td>
<td style="text-align: left;">Targeted only; long imaging; no protein
co-detection; high custom cost.</td>
</tr>
<tr>
<td style="text-align: left;">MERSCOPE (Akoya Biosciences</td>
<td style="text-align: left;">~1 μm (subcellular)</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">Customizable (up to 1,000 genes);
RNA/protein co-detection; FFPE/cells.</td>
<td style="text-align: left;">Gene limit; long times; lower FFPE
sensitivity.</td>
</tr>
<tr>
<td style="text-align: left;">CosMx SMI (NanoString)</td>
<td style="text-align: left;">~1 μm (subcellular)</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">Up to 6,000 genes; RNA/protein
co-detection; fast cycles; FFPE.</td>
<td style="text-align: left;">Lower specificity; sequential imaging;
costly.</td>
</tr>
</tbody>
</table>

#### Sequencing-Based Technologies

These use next-generation sequencing for whole-transcriptome analysis,
ideal for discovery.

<table>
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Technology</th>
<th style="text-align: left;">Resolution</th>
<th style="text-align: left;">Throughput</th>
<th style="text-align: left;">Pros</th>
<th style="text-align: left;">Cons</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Visium (10x Genomics)</td>
<td style="text-align: left;">~55 μm (multi-cell)</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">Whole transcriptome; FFPE/fresh; large
area; protein panels.</td>
<td style="text-align: left;">No single-cell; low sensitivity for rare
transcripts; needs deconvolution.</td>
</tr>
<tr>
<td style="text-align: left;">Visium HD (10x Genomics)</td>
<td style="text-align: left;">~2 μm (near single-cell)</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">High resolution; whole transcriptome;
human/mouse FFPE.</td>
<td style="text-align: left;">Reduced genes at high res;
species-limited; smaller area.</td>
</tr>
<tr>
<td style="text-align: left;">Stereo-seq</td>
<td style="text-align: left;">~0.2 μm (near single-cell)</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">Nanoscale res; whole transcriptome; all
species/FFPE; large area.</td>
<td style="text-align: left;">Diffusion artifacts; specialized
sequencer; lower sensitivity.</td>
</tr>
<tr>
<td style="text-align: left;">GeoMx DSP (NanoString)</td>
<td style="text-align: left;">~50 μm (ROI-based)</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">Targeted ROI; RNA/protein; FFPE
human/mouse; multi-sample.</td>
<td style="text-align: left;">Predefined ROIs; no single-cell;
species-limited.</td>
</tr>
</tbody>
</table>
