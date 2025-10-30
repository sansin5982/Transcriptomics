# Experimental Design in Transcriptomics

This chapter is the blueprint for turning curiosity into discovery.
Whether you’re a novice or refining your skills, experimental design is
where science meets strategy. A well-designed study ensures your
transcriptomic data is reproducible, interpretable, and biologically
meaningful.

## Key Considerations

### 1. Defining Biological Questions and Hypotheses

**Why it matters**: Your question drives everything—sample choice,
method, analysis. A vague question yields vague results.

#### How to do it (step-by-step):

-   **1. Start with a biological problem**: “How does obesity affect
    liver gene expression?”
-   **2. Formulate a testable hypothesis**: “Obesity upregulates lipid
    synthesis genes in hepatocytes.”
-   **3. Make it specific and measurable:** “We hypothesize that FASN
    and SCD1 are &gt;2-fold upregulated in obese vs. lean human liver
    biopsies.”

**Novice tip**: Use the **PICO framework**:

-   **P**opulation: Obese vs. lean adults
-   **I**ntervention/Exposure: Obesity
-   **C**omparison: Lean controls
-   **O**utcome: Differential gene expression

**Example**: A 2025 study asked: “Does Alzheimer’s disease alter
microglial activation in the hippocampus?” → Hypothesis: TREM2 and
CX3CR1 are upregulated in AD patient microglia. → This directed
scRNA-Seq on hippocampal tissue, leading to targeted biomarker
discovery.

### 2. Selection of Appropriate Samples

**Why it matters**: The transcriptome is **tissue-, cell-, and
condition-specific**. Wrong samples = wrong answers.

**Step-by-step guide**: |Step|Action|Example| |:—|:—|:—| |1|Define the
biological context|Liver for metabolic disease, not blood| |2|Choose
sample type|Fresh-frozen tissue, FFPE, single cells| |3|Match
conditions|Age, sex, BMI, medication| |4|Avoid confounders|Smoking,
inflammation|

**Novice trap**: Using whole blood for brain-specific questions →
diluted signal.

**2025 Example**: A lung cancer study compared **tumor core**
vs. **tumor edge** vs. **healthy lung** tissue using spatial
transcriptomics. They found PD-L1 highly expressed only at the invasive
edge—critical for immunotherapy targeting.

### 3. Replication and Sample Size

**Why it matters**: Biology is variable. Replication distinguishes
signal from noise.

**Step-by-step**:

-   1.  **Biological replicates**: Different individuals (n ≥ 3 per
        group)
-   1.  **Technical replicates**: Same sample, different runs (rarely
        needed with NGS)
-   1.  **Power calculation**: Use tools like RNASeqPower or Scotty

#### **Rule of thumb**:

-   Bulk RNA-Seq: n = 3–6 per group
-   scRNA-Seq: n = 3 patients, &gt;5,000 cells per sample

**Example**: A 2025 diabetes study used n = 5 obese and n = 5 lean
patients. Power analysis showed 80% power to detect 2-fold changes at
FDR &lt; 0.05. Novice tip: Under-replication → false positives.
Over-replication → wasted resources.

## Sample Preparation

### 1. RNA Extraction and Quality Control

**Goal**: Isolate intact, pure RNA. **Step-by-step protocol** (for
tissue):

1.  **Homogenize tissue** in TRIzol or Qiagen RNeasy kit
2.  **Phase separation** (chloroform) → aqueous phase = RNA
3.  **Column purification** → DNase treatment (remove DNA)
4.  **Quantify**: NanoDrop (concentration), Qubit (accurate),
    Bioanalyzer (integrity)

#### Quality metrics:

<table>
<thead>
<tr>
<th style="text-align: left;">Metric</th>
<th style="text-align: left;">Good</th>
<th style="text-align: left;">Poor</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">RIN (RNA Integrity Number)</td>
<td style="text-align: left;">≥ 7.0</td>
<td style="text-align: left;">&lt; 6.0</td>
</tr>
<tr>
<td style="text-align: left;">260/280 ratio</td>
<td style="text-align: left;">1.8–2.1</td>
<td style="text-align: left;">&lt; 1.7 (protein)</td>
</tr>
<tr>
<td style="text-align: left;">260/230 ratio</td>
<td style="text-align: left;">&gt; 2.0</td>
<td style="text-align: left;">&lt; 1.8 (salts)</td>
</tr>
</tbody>
</table>

**Novice tip**: Degraded RNA → 3’ bias in sequencing. Always check
electropherogram.

**2025 Example**: A brain bank study rejected 30% of samples with RIN
&lt; 7, ensuring reliable Alzheimer’s transcriptomes.

### 2. Handling Low-Input or Degraded RNA Samples

**Challenge**: Biopsies, laser-capture microdissection, FFPE samples.
**Solutions**: |Sample Type|Method|Key Step| |:—|:—|:—| |**Low input
(&lt;10 ng)**|Smart-Seq2, CEL-Seq2|Full-length amplification|
|**FFPE**|TruSeq RNA Exome, Ovation FFPE|Deparaffinization + rRNA
depletion| |**Degraded**|RNAtag-Seq, BRB-Seq|3’ end capture|

**Example**: A 2025 pancreatic cancer study used **FFPE blocks from
2015** with a new **FFPE-optimized RNA-Seq kit**, recovering 80% of
transcripts despite degradation.

### 3. Considerations for Single-Cell Studies

#### Unique challenges:

-   Cell viability (&gt;80%)
-   Doublet removal
-   Ambient RNA contamination

#### Step-by-step:

**1. Tissue dissociation**: Enzymatic (collagenase), gentle **2. Cell
counting**: Trypan blue or automated (Countess) **3. Loading**: 10x
Chromium → target 5,000–10,000 cells **4. Library QC**: Check cDNA size
(Tapestation)

**Novice tip**: Dead cells release RNA → false ambient signal. Use
**CellBender** to remove. **2025 Example**: A kidney scRNA-Seq study
used **fresh biopsies within 30 minutes** to achieve 92% viability,
revealing a novel podocyte subtype in lupus nephritis.

## Controls and Normalization

### 1. Use of Housekeeping Genes and Spike-In Controls

**Housekeeping genes** (e.g., GAPDH, ACTB): Assumed stable → **not
always true**! **Better: ERCC spike-ins** (synthetic RNAs of known
concentration)

**Step-by-step**: 1. Add ERCC mix **before** library prep 2. Sequence →
count ERCC reads 3. Use for **absolute quantification** and **technical
normalization**

**Example**: A 2025 study found GAPDH upregulated in hypoxia →
misleading normalization. ERCC spike-ins corrected this, revealing true
hypoxia response genes.

### 2. Normalization Methods for Transcriptomic Data

#### Goal: Remove technical variation (sequencing depth, RNA input)

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">When to use</th>
<th style="text-align: left;">Formula</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><strong>TPM (Transcripts Per
Million)</strong></td>
<td style="text-align: left;">Compare genes within sample</td>
<td style="text-align: left;">(reads × 10⁶) / (gene length × total
reads)</td>
</tr>
<tr>
<td style="text-align: left;"><strong>RPKM/FPKM</strong></td>
<td style="text-align: left;">Older bulk RNA-Seq</td>
<td style="text-align: left;">Similar to TPM</td>
</tr>
<tr>
<td style="text-align: left;"><strong>TMM (Trimmed Mean of
M-values)</strong></td>
<td style="text-align: left;">DESeq2, edgeR</td>
<td style="text-align: left;">Scales by effective library size</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Spike-in
normalization</strong></td>
<td style="text-align: left;">Low-input, scRNA-Seq</td>
<td style="text-align: left;">ERCC-based scaling</td>
</tr>
</tbody>
</table>

**Tip: Never use raw counts for comparison across samples.** **2025
Example**: A multi-lab RNA-Seq study used **TMM + batch correction
(ComBat-seq)** to harmonize data from 12 centers, identifying robust
COVID-19 biomarkers.
