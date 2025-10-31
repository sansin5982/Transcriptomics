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

> **Lung Cancer RNA-Seq** \* 100 patient samples `FastQC` showed
> **adapter contamination** in 30 samples Trimmed → **recovered 15% more
> reads** → found EGFR mutation missed before

### .2 Read Alignment and Mapping

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

> **Alzheimer’s Brain Study** \* Used **STAR** to align 1 billion reads
> \* Found **90% mapped** → **10% unmapped = novel RNAs** \* →
> Discovered **new lncRNA** near APP gene

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

> **Muscular Dystrophy** \* Patient had normal DMD gene DNA but weak
> muscles \* StringTie found exon 45 skipped → short protein \* →
> Exon-skipping therapy (Exondys 51) saved walking ability

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

#### Real Example (2025):

> **Liver vs. Brain Comparison** \* Used **TPM** → ALB (liver protein)
> **100 x higher in liver** \* RPKM gave **wrong ratio** due to library
> size
