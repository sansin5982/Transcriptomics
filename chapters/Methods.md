# Methods for Transcriptome Profiling

Think of the **transcriptome** as a **library of messages** (RNA) that
cells send. These methods are **different ways to read those messages**.

## 1. Hybridization-Based Methods

**(“Sticky Notes” Methods – RNA sticks to DNA)**

### 1.1 Northern Blotting

**“The Old-School Detective”**

#### What it is:

A way to **detect ONE specific RNA** using a **sticky DNA probe**.

#### Step-by-Step (Like a Recipe):

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">What You Do</th>
<th style="text-align: left;">Picture</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;"><strong>Extract RNA</strong> from
cells</td>
<td style="text-align: left;"></td>
</tr>
<tr>
<td style="text-align: left;">2</td>
<td style="text-align: left;"><strong>Run RNA on a gel</strong> (like
sorting by size)</td>
<td style="text-align: left;">Small RNAs move far, big ones stay near
top</td>
</tr>
<tr>
<td style="text-align: left;">3</td>
<td style="text-align: left;"><strong>Transfer to membrane
(blot)</strong></td>
<td style="text-align: left;">Like pressing a stamp</td>
</tr>
<tr>
<td style="text-align: left;">4</td>
<td style="text-align: left;"><strong>Add radioactive/fluorescent DNA
probe</strong></td>
<td style="text-align: left;">Probe “sticks” only to matching RNA</td>
</tr>
<tr>
<td style="text-align: left;">5</td>
<td style="text-align: left;">Detect signal (X-ray or camera)</td>
<td style="text-align: left;">Bright band = lots of RNA</td>
</tr>
</tbody>
</table>

#### Real Example (2025):

> **Alzheimer’s Research** \* Wanted to know: “Is APP gene RNA increased
> in old brains?” \* Used **Northern blot** on 10 human brains \*
> Result: **APP RNA 3x higher in Alzheimer’s brains** → confirmed
> protein buildup

**Pros**: Simple, no fancy machine **Cons**: Only 1 gene at a time, slow

## 1.2 Microarrays

**“The Sticker Wall” – 20,000 genes at once!**

#### What it is:

A **glass slide with 20,000 tiny DNA spots** (probes). RNA from your
sample **sticks** to matching spots → glows.

#### Workflow:

<table>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">Action</th>
<th style="text-align: left;">Example</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;"><strong>Label RNA</strong> with green/red
dye</td>
<td style="text-align: left;">Cancer = green, Healthy = red</td>
</tr>
<tr>
<td style="text-align: left;">2</td>
<td style="text-align: left;"><strong>Pour on microarray
chip</strong></td>
<td style="text-align: left;">RNA sticks to matching DNA spots</td>
</tr>
<tr>
<td style="text-align: left;">3</td>
<td style="text-align: left;"><strong>Scan with laser</strong></td>
<td style="text-align: left;">Green spot = cancer gene ON</td>
</tr>
<tr>
<td style="text-align: left;">4</td>
<td style="text-align: left;"><strong>Analyze colors</strong></td>
<td style="text-align: left;">Yellow = equal, Red = healthy gene ON</td>
</tr>
</tbody>
</table>

\####Real Example (2025):

> **Breast Cancer Subtypes** \* 100 patients → tumor RNA (green), normal
> (red) \* Microarray showed **HER2 gene glowing bright green** in 20
> patients \* → Doctors gave **Herceptin drug** → 5-year survival ↑ 30%

**Pros**: Cheap, fast for known genes **Cons**: Can’t find new genes,
background noise

## 2. Sequence-Based Methods

**(“Reading the Letters” – Actual RNA sequence!)**

### 2.1 Sanger Sequencing

**“The First DNA Reader” – 1977 Nobel Prize**

#### How it works:

-   Copy RNA → cDNA
-   Add **color-coded DNA letters** (A=T=green, T=A=red, etc.)
-   Machine reads **one base at a time** → “ATCGGT…”

#### Real Example:

> **HIV Drug Resistance (2025)** \* Patient not responding to drug \*
> Sanger sequenced **HIV pol gene RNA** \* Found **M184V mutation** →
> switched to new drug → virus gone in 3 months

**Pros**: Gold standard for accuracy **Cons**: Only ~1,000 bases per run

### 2.2 RNA-Seq

**“The Google of Transcriptomics” – Reads EVERYTHING**

#### Library Preparation (Step-by-Step):

<table>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">WHat Happens</th>
<th style="text-align: left;">Why</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;"><strong>Extract RNA</strong></td>
<td style="text-align: left;">Get all messages</td>
</tr>
<tr>
<td style="text-align: left;">2</td>
<td style="text-align: left;"><strong>Break into tiny pieces</strong>
(~200 letters)</td>
<td style="text-align: left;">So sequencer can read</td>
</tr>
<tr>
<td style="text-align: left;">3</td>
<td style="text-align: left;"><strong>Add adapters</strong>
(barcodes)</td>
<td style="text-align: left;">Like address labels</td>
</tr>
<tr>
<td style="text-align: left;">4</td>
<td style="text-align: left;"><strong>Amplify</strong> (make
copies)</td>
<td style="text-align: left;">Need millions of reads</td>
</tr>
<tr>
<td style="text-align: left;">5</td>
<td style="text-align: left;"><strong>Sequence</strong> (Illumina)</td>
<td style="text-align: left;">Machine reads 150 letters × 50 million
times</td>
</tr>
</tbody>
</table>

#### Read Mapping:

-   Computer matches each 150-letter piece to **human genome map**
-   Count how many land on each gene → **expression level**

#### Real Example (2025):

> **Long COVID Brain Fog** \* 50 patients vs. 50 healthy \* RNA-Seq on
> blood → found **100 genes stuck “ON”** (inflammation) \* → New drug
> trial targeting IL-6 pathway

**Pros**: Finds new genes, isoforms, mutations **Cons**: Expensive,
needs bioinformatics

### 2.3 Tag-Based Methods (SAGE, CAGE)

**“Counting Tags” – Like counting book titles**

#### SAGE (Serial Analysis of Gene Expression):

-   Cut **10-letter tag** from each RNA
-   Glue tags together → sequence long chain
-   Count tags → gene expression

#### CAGE (Cap Analysis of Gene Expression):

-   Focus on **5’ end** (start of RNA)
-   Finds exact **transcription start sites (TSS)**

#### Real Example (2025):

> **FANTOM6 Project** \* Used **CAGE** on 1,800 human samples \* Found
> **200,000 new gene start points** \* → Discovered **lincRNAs**
> controlling heart development

### 3. Quantitative PCR (qPCR)

**“The Zoom Lens” – Validate RNA-Seq results**

#### Principles of Real-Time PCR

<table>
<thead>
<tr>
<th style="text-align: left;">Step</th>
<th style="text-align: left;">What Happens</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">1</td>
<td style="text-align: left;">Add <strong>primers</strong> (short DNA
hooks) + <strong>fluorescent dye</strong></td>
</tr>
<tr>
<td style="text-align: left;">2</td>
<td style="text-align: left;">Heat → DNA strands separate</td>
</tr>
<tr>
<td style="text-align: left;">3</td>
<td style="text-align: left;">Cool → primers stick</td>
</tr>
<tr>
<td style="text-align: left;">4</td>
<td style="text-align: left;">Polymerase copies → <strong>dye glows
more</strong></td>
</tr>
<tr>
<td style="text-align: left;">5</td>
<td style="text-align: left;">Measure glow <strong>every cycle → Ct
value</strong> (earlier = more RNA)</td>
</tr>
</tbody>
</table>

#### Formula

$$
\large 2^{-(\triangle Ct\_{Target} - \triangle Ct\_{housekeeping})}
$$

#### Real Example:

> **COVID-19 Testing** \* qPCR on nasal swab \* **Ct = 20** → high virus
> (very sick) \* **Ct = 35** → low virus (mild) \* Used in **billions of
> tests worldwide**

### High-Throughput qPCR

-   **96-well or 384-well plates**
-   **Fluidigm Biomark**: 9,000 reactions in one run!

#### Example:

> **Cancer Drug Screening** \* Test **96 drugs** on **96 patient
> tumors** \* qPCR for **5 key genes** → predict response in 4 hours

## 4. Emerging Technologies

**“The Future is Here”**

### 4.1 Long-Read Sequencing

**(Oxford Nanopore, PacBio)**

#### Why long reads?

-   Short reads (150 bp) → **can’t solve complex puzzles**
-   Long reads (10,000+ bp) → **full gene isoforms**

#### Nanopore: “DNA through a tiny hole”

-   RNA passes through protein pore
-   Current change = base identity
-   **Portable** (size of a USB!)

#### Real Example:

> **Rare Disease Diagnosis** \* Child with unknown muscle disease \*
> **Nanopore sequenced full DMD gene** (2.2 million bases) \* Found
> **huge deletion** missed by short-read → correct diagnosis in 48 hours

### 4.2 Spatial Transcriptomics

**(Visium, Slide-seq)**

**“Where in the tissue?”**

#### Visium (10x Genomics):

-   Tissue section on slide with **5,000 spots**
-   Each spot = **55 µm diameter** (10–20 cells)
-   RNA sticks to barcoded beads → sequence → **map genes to location**

#### Real Example

> **Brain Tumor Surgery** \* Surgeon needed to know **tumor border** \*
> Visium on biopsy → **red zone = cancer genes, blue = healthy** →
> Removed **only tumor**, saved speech center

#### 4.3 In Situ Sequencing

**“Sequence RNA inside the cell!”**

#### How?

-   Fix tissue → add probes → sequence in place → image

#### Real Example (2025):

> **Alzheimer’s Plaques** \* Used **in situ sequencing** on brain slice
> Found **APP RNA concentrated around amyloid plaques** → Proved plaques
> **trigger local gene changes**

#### Summary Table: Which Method When?

<table>
<thead>
<tr>
<th style="text-align: left;">Goal</th>
<th style="text-align: left;">Method</th>
<th style="text-align: left;">Example</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><strong>One gene, confirm</strong></td>
<td style="text-align: left;">Northern / qPCR</td>
<td style="text-align: left;">Validate BRCA1 in cancer</td>
</tr>
<tr>
<td style="text-align: left;"><strong>20,000 genes, fast</strong></td>
<td style="text-align: left;">Microarray</td>
<td style="text-align: left;">Classify leukemia subtypes</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Everything, discover</strong></td>
<td style="text-align: left;">RNA-seq</td>
<td style="text-align: left;">Find new COVID variants</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Full gene structure</strong></td>
<td style="text-align: left;">Nanopore/PacBio</td>
<td style="text-align: left;">Solve splicing in autism</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Where in tissue?</strong></td>
<td style="text-align: left;">Visium</td>
<td style="text-align: left;">Map tumor-immune border</td>
</tr>
<tr>
<td style="text-align: left;"><strong>Inside single cell?</strong></td>
<td style="text-align: left;">In situ seq</td>
<td style="text-align: left;">See RNA in neuron dendrites</td>
</tr>
</tbody>
</table>

#### Final Message

> **“Every method is a tool. Pick the right one for your question.”**

-   **Northern blot** → detective for 1 gene
-   **Microarray** → wall of stickers
-   **RNA-Seq** → Google search
-   **Nanopore** → full book
-   **Visium** → treasure map
