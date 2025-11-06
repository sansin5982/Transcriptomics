# Single-Cell Transcriptomics

> Imagine zooming into a city from space → you see lights (bulk
> RNA-Seq). Now zoom to street level → you see who’s awake, who’s
> working, who’s sick. That’s single-cell RNA-Seq (scRNA-Seq).

## 1. Introduction to Single-Cell Transcriptomics

### 1.1 Why Study Gene Expression at Single-Cell Resolution?

#### Bulk RNA-Seq = average of 10,000 cells

→ Like asking: “What’s the average mood in New York?”

#### scRNA-Seq = each cell’s voice

→ “This neuron is stressed, this immune cell is fighting, this cancer
cell is hiding.”

#### Purpose:

-   Find **rare cells** (0.1% of tissue)
-   See **cell states** (resting vs. activated)
-   Track **development** (stem cell → neuron)
-   Understand **disease heterogeneity**

#### Real Example:

> **Brain Tumor**
>
> -   Bulk RNA-Seq: “Tumor is aggressive”
>
> -   scRNA-Seq: Found 1% of cells were cancer stem cells hiding from
>     chemo
>
> -   → New drug targets only those 1% → tumor shrank 80% in mice

### 1.2 Bulk vs. Single-Cell RNA-Seq – The Key Differences

<table>
<thead>
<tr>
<th style="text-align: left;">Feature</th>
<th style="text-align: left;">Bulk RNA-Seq</th>
<th style="text-align: left;">scRNA-Seq</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Input</td>
<td style="text-align: left;">10,000-1M cells</td>
<td style="text-align: left;">1 cell</td>
</tr>
<tr>
<td style="text-align: left;">Output</td>
<td style="text-align: left;">Average expression</td>
<td style="text-align: left;">Cell-specific profile</td>
</tr>
<tr>
<td style="text-align: left;">Resolution</td>
<td style="text-align: left;">Tissue level</td>
<td style="text-align: left;">Cell-type &amp; state</td>
</tr>
<tr>
<td style="text-align: left;">Rare cells</td>
<td style="text-align: left;">Hidden</td>
<td style="text-align: left;">Detected</td>
</tr>
<tr>
<td style="text-align: left;">Cost</td>
<td style="text-align: left;">USD200-500</td>
<td style="text-align: left;">USD1000-5000</td>
</tr>
<tr>
<td style="text-align: left;">Data size</td>
<td style="text-align: left;">~5 GB</td>
<td style="text-align: left;">50-500 GB</td>
</tr>
</tbody>
</table>

#### Real Example:

> **Lung Fibrosis**
>
> -   Bulk: “Fibroblasts upregulated”
>
> -   scRNA-Seq: **Only 1 subtype** (not all) → **targeted therapy**
>     avoided side effects

## 2. Workflow of scRNA-Seq

### 2.1 Cell Isolation and Barcoding

**Goal**: Put **one cell in one droplet** with a **unique barcode**

#### Methods:

<table>
<thead>
<tr>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">How it works</th>
<th style="text-align: left;">Cells</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">10x Chromium</td>
<td style="text-align: left;">Oil droplets trap cells + beads</td>
<td style="text-align: left;">10,000 cells/run</td>
</tr>
<tr>
<td style="text-align: left;">Drop-seq</td>
<td style="text-align: left;">DIY version</td>
<td style="text-align: left;">Cheaper</td>
</tr>
<tr>
<td style="text-align: left;">Smart-Seq2</td>
<td style="text-align: left;">Plate-based</td>
<td style="text-align: left;">Full-length RNA</td>
</tr>
</tbody>
</table>

#### Step-by-step (10x):

1.  Dissociate tissue → enzymes (37°C, 20 min)
2.  Filter → 40 µm (remove clumps)
3.  Count → aim for 1,000 cells/µL
4.  Load into 10x chip → cell + bead + barcode in droplet
5.  Lyse cell → RNA sticks to bead

#### Real Example (2025):

> Human Pancreas Atlas
>
> -   Used **10x Chromium** on 20 donors
> -   Captured **1.2 million cells**
> -   → Found **new beta-cell subtype** that dies first in diabetes

### 2.2 Library Preparation and Sequencing

#### After barcoding:

1.  **Reverse transcription** → RNA → cDNA (with barcode + UMI)
2.  **Amplify** → PCR (12–16 cycles)
3.  **Fragment** → add sequencing adapters
4.  Sequence → Illumina NovaSeq (100M reads)

#### UMI (Unique Molecular Identifier):

→ 10-letter random tag → Counts real RNA molecules, not PCR copies

#### Real Example:

> Cancer Immunotherapy
>
> -   Patient’s tumor biopsy → 10x → sequenced
>
> -   UMI showed T-cells had 1,000 PD-1 RNAs each
>
> -   → Anti-PD-1 drug worked → tumor gone in 6 months

### 2.3 Data Preprocessing and Quality Control

Raw data: `.fastq` files with barcodes + UMIs

    # 1. Demultiplex
    cellranger mkfastq

    # 2. Align + Count
    cellranger count --id=sample

    # 3. QC in R (Seurat)
    library(Seurat)
    pbmc <- CreateSeuratObject(counts, min.cells=3, min.features=200)

#### QC Filters:

<table>
<thead>
<tr>
<th style="text-align: left;">Metric</th>
<th style="text-align: left;">Good</th>
<th style="text-align: left;">Bad</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Genes per cell</td>
<td style="text-align: left;">500-5000</td>
<td style="text-align: left;">&lt;200 (dead), &gt;7,000 (doublet)</td>
</tr>
<tr>
<td style="text-align: left;">Mitochondrial %</td>
<td style="text-align: left;">&lt;10%</td>
<td style="text-align: left;">&gt;20% (dying)</td>
</tr>
<tr>
<td style="text-align: left;">UMIs</td>
<td style="text-align: left;">&gt;1,000</td>
<td style="text-align: left;">&lt;500</td>
</tr>
</tbody>
</table>

#### Real Example:

> Alzheimer’s Brain
>
> -   100,000 cells → after QC → 70,000 high-quality
>
> -   Removed high-mito cells → found microglia activation not cell
>     death
