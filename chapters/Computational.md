# Computational Methods in Transcriptomics

## 1. Data Preprocessing

### 1.1 Quality Control of Sequencing Reads

#### Why?

Sequencers make mistakes — **bad reads = wrong answers**

Tools: `FastQC`, `MultiQC`

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">MultiQC</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;">Run <code>FastQC</code> on
<code>.fastq</code> files</td>
</tr>
<tr>
<td style="text-align: left;">2</td>
<td style="text-align: left;">Check <strong>Per Base Quality</strong>
(green = good, red = bad)</td>
</tr>
<tr>
<td style="text-align: left;">3</td>
<td style="text-align: left;">Trim <strong>low-quality ends</strong> and
<strong>adapters</strong> with <code>Trim Galore!</code> or any other
tool</td>
</tr>
</tbody>
</table>

#### Real Example:

> **Lung Cancer RNA-Seq**
>
> -   100 patient samples
>
> -   `FastQC` showed **adapter contamination** in 30 samples
>
> -   Trimmed → **recovered 15% more reads** → found EGFR mutation
>     missed before

### 1.2 Read Alignment and Mapping

**Goal**: Match each 150-letter read to the **human genome map**.

**Tools**:

<table>
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Speed</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><strong>STAR</strong></td>
<td style="text-align: left;">Spliced Alignment</td>
<td style="text-align: left;">Fast</td>
</tr>
<tr>
<td style="text-align: left;"><strong>HISTA2</strong></td>
<td style="text-align: left;">Graph-based</td>
<td style="text-align: left;">Memory-light</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Bowtie2</strong></td>
<td style="text-align: left;">Short reads</td>
<td style="text-align: left;">Ultra-fast</td>
</tr>
</tbody>
</table>

**Step-by-step (STAR)**:

    STAR --genomeDir hg38_index --readFilesIn sample_R1.fastq sample_R2.fastq --outSAMtype BAM

#### Real Example:

> **Alzheimer’s Brain Study**
>
> -   Used **STAR** to align 1 billion reads
>
> -   Found **90% mapped** → **10% unmapped = novel RNAs**
>
> -   → Discovered **new lncRNA** near APP gene

### 1.3 Handling Alternative Splicing & Novel Transcripts

Problem: One gene → multiple RNA versions (isoforms)

**Tools**:

<table>
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">What it does</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">StringTie</td>
<td style="text-align: left;">Assemble transcripts</td>
</tr>
<tr>
<td style="text-align: left;">Salmon (quasi-mapping)</td>
<td style="text-align: left;">Quantifies isoforms</td>
</tr>
</tbody>
</table>

#### Real Example:

> **Muscular Dystrophy**
>
> -   Patient had normal DMD gene DNA but weak muscles
>
> -   StringTie found exon 45 skipped → short protein
>
> -   → Exon-skipping therapy (Exondys 51) saved walking ability

## 2. Quantification of Gene Expression

### 2.1 RPKM, FPKM, TPM – What’s the Difference?

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">Formula</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><strong>RPKM/FPKM</strong></td>
<td style="text-align: left;"><span class="math inline">$\frac{reads *
10^9}{\text {gene length} * \text {total reads}}$</span></td>
</tr>
<tr>
<td style="text-align: left;"><strong>TPM</strong></td>
<td style="text-align: left;"><span
class="math inline">$\frac{\frac{reads}{\text{gene length}}} *
{\frac{10^6}{sum(\text{all TPM)}}$</span></td>
</tr>
</tbody>
</table>

#### Why TPM wins:

-   **Gene A**: 1000 reads, 1000 bp
-   **Gene B**: 1000 reads, 2000 bp → **RPKM says equal**, but **TPM
    says Gene A is 2x higher**

#### Real Example:

> **Liver vs. Brain Comparison**
>
> -   Used **TPM** → ALB (liver protein) **100 x higher in liver**
>
> -   RPKM gave **wrong ratio** due to library size

### 2.2 Tools for Quantification

<table>
<thead>
<tr>
<th style="text-align: left;">Tools</th>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">Speed</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Salmon</td>
<td style="text-align: left;">Quasi-mapping (no alignment)</td>
<td style="text-align: left;">&lt;10 min</td>
</tr>
<tr>
<td style="text-align: left;">Kallisto</td>
<td style="text-align: left;">Pseudoalignment</td>
<td style="text-align: left;">&lt;5 min</td>
</tr>
<tr>
<td style="text-align: left;">featureCounts</td>
<td style="text-align: left;">After STAR</td>
<td style="text-align: left;">Accurate</td>
</tr>
</tbody>
</table>

#### Real Example:

> **COVID-19 Drug Screen**
>
> -   10,000 compounds, 96-well plates
>
> -   Used **Kallisto** → quantified **20,000 genes in 3 minutes per
>     plate**
>
> -   → Found **remdesivir blocks viral RNA polymerase**

## 3. Differential Expression Analysis

### 3.1 Statistical Tools: DESeq2, edgeR

<table>
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">Input</th>
<th style="text-align: left;">Output</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">DESeq2</td>
<td style="text-align: left;">Raw counts</td>
<td style="text-align: left;">log2FoldChange, p-value, padj</td>
</tr>
<tr>
<td style="text-align: left;">edgeR</td>
<td style="text-align: left;">Raw counts</td>
<td style="text-align: left;">Same + better for small n</td>
</tr>
</tbody>
</table>

#### DESeq2 Formula (simplified):

    dds <- DESeqDataSetFromMatrix(countData, colData, design = ~condition)
    dds <- DESeq(dds)
    results <- results(dds, contrast=c("condition","disease","control"))

### 3.2 Fold Change, p-value, FDR

<table>
<thead>
<tr>
<th style="text-align: left;">Term</th>
<th style="text-align: left;">Meaning</th>
<th style="text-align: left;">Threshold</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">log2FC &gt; 1</td>
<td style="text-align: left;">2x upregulated</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">p &lt; 0.05</td>
<td style="text-align: left;">Statistically significant</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">FDR &lt; 0.05</td>
<td style="text-align: left;">Controls false positives</td>
<td style="text-align: left;">Use this!</td>
</tr>
</tbody>
</table>

#### Volcano plot

    x = log2FC
    y = -log10(padj)

→ Top-right = **up in disease, significant** \#### Real Example:

> **Parkinson’s vs. Healthy Brain**
>
> -   DESeq2 → 1,200 genes changed
>
> -   Top hit: SNCA (alpha-synuclein) ↑ 4-fold, FDR = 10<sup>−</sup>15
>
> -   → New drug targeting SNCA RNA in clinical trial

## 4. Functional Annotation & Pathway Analysis

### 4.1 Gene Ontology (GO) & KEGG

<table>
<thead>
<tr>
<th style="text-align: left;">Database</th>
<th style="text-align: left;">What it tells you</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">GO</td>
<td style="text-align: left;">Biological Process, Cellular Component,
Molecular Function</td>
</tr>
<tr>
<td style="text-align: left;">KEGG</td>
<td style="text-align: left;">Metabolic pathways, diseases</td>
</tr>
</tbody>
</table>

#### Real Example:

> **Obesity RNA-Seq**
>
> -   500 upregulated genes → input to Enrichr
>
> -   Top GO: “Lipid metabolic process” (*p* = 10<sup>−</sup>25)
>
> -   KEGG: “PPAR signaling pathway”
>
> -   → Drug pioglitazone activates PPAR → weight loss in trial

### 4.2 Tools: DAVID, Enrichr, GSEA

<table>
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">Best For</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Enrichr</td>
<td style="text-align: left;">Fast, web-based</td>
</tr>
<tr>
<td style="text-align: left;">GSEA</td>
<td style="text-align: left;">Pathway ranking (uses all genes)</td>
</tr>
<tr>
<td style="text-align: left;">clusterProfiler (R)</td>
<td style="text-align: left;">Custom plots</td>
</tr>
</tbody>
</table>

#### GSEA Example (2025):

> **Cancer vs. Normal**
>
> -   GSEA showed **“Hypoxia pathway”** enriched (NES = 2.8)
>
> → Tumors **resist chemotherapy** in low oxygen
>
> → New **hypoxia-activated** prodrug in Phase II

## 5. Visualization of Transcriptomic Data

### 5.1 Heatmaps, Volcano, MA Plots

<table>
<thead>
<tr>
<th style="text-align: left;">Plot</th>
<th style="text-align: left;">What it shows</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><strong>Heatmap</strong></td>
<td style="text-align: left;">Expression patterns across samples</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Volcano</strong></td>
<td style="text-align: left;">Significance vs. fold change</td>
</tr>
<tr>
<td style="text-align: left;"><strong>MA PLot</strong></td>
<td style="text-align: left;">log ratio vs. mean expression</td>
</tr>
</tbody>
</table>

Real Example:

> **5 Cancer Types**
>
> -   Heatmap → **clear clusters**: breast, lung, colon
>
> -   → AI model **99% accurate** in diagnosing cancer type from RNA

### 5.2 Dimensionality Reduction: PCA, t-SNE, UMAP

<table>
<thead>
<tr>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">Use</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">PCA</td>
<td style="text-align: left;">Find main sources of variation</td>
</tr>
<tr>
<td style="text-align: left;">UMAP</td>
<td style="text-align: left;">Visualize clusters (better than
t-SNE)</td>
</tr>
</tbody>
</table>

#### PCA Example:

> **scRNA-Seq of Immune Cells**
>
> -   PCA → **PC1 = T cells vs. B cells**
>
> -   PC2 = **activated vs. resting**
>
> → Found **new exhausted T-cell state** in chronic infection

### 5.3 Tools: R, Python, Tableau

<table>
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">Best For</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">R (ggplot2, pheatmap)</td>
<td style="text-align: left;">Publication-quality</td>
</tr>
<tr>
<td style="text-align: left;">Python (scanpy, seaborn)</td>
<td style="text-align: left;">scRNA-Seq</td>
</tr>
<tr>
<td style="text-align: left;">Tableau</td>
<td style="text-align: left;">Interactive dashboards</td>
</tr>
</tbody>
</table>

#### Real Example:

> **Hospital Dashboard**
>
> -   Tableau + RNA-Seq → **live cancer subtype predictor**
>
> -   Surgeon uploads biopsy → **result in 2 hours**

#### First Analysis Pipeline

    # 1. QC
    fastqc *.fastq
    multiqc .

    # 2. Trim
    trim_galore --paired R1.fastq R2.fastq

    # 3. Align & Quantify
    salmon quant -i hg38_index -l A -1 R1_trimmed.fq -2 R2_trimmed.fq -o quant

    # 4. Differential Expression (R)
    library(DESeq2)
    counts <- tximport("quant.sf", type="salmon")
    dds <- DESeqDataSetFromTximport(counts, colData, ~condition)
    dds <- DESeq(dds)
    res <- results(dds)

    # 5. Visualize
    volcanoPlot(res)

#### Summary

<table>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">Tool</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">QC</td>
<td style="text-align: left;">FastQC</td>
</tr>
<tr>
<td style="text-align: left;">Align</td>
<td style="text-align: left;">STAR</td>
</tr>
<tr>
<td style="text-align: left;">Quantify</td>
<td style="text-align: left;">Salmon</td>
</tr>
<tr>
<td style="text-align: left;">DE</td>
<td style="text-align: left;">DESeq2</td>
</tr>
<tr>
<td style="text-align: left;">Pathway</td>
<td style="text-align: left;">GSEA</td>
</tr>
<tr>
<td style="text-align: left;">Plot</td>
<td style="text-align: left;">UMAP</td>
</tr>
</tbody>
</table>
