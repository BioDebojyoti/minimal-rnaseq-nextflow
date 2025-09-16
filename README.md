# minimal-rnaseq-nextflow


This is a modular **Nextflow DSL2 RNA-seq pipeline** for processing FASTQ data.  
It supports:

- QC: **FastQC**  
- Adapter trimming: **Trim Galore**  
- Quantification:
  - Alignment-free: **Salmon**
  - Alignment-based: **STAR + featureCounts**  
- Optional automatic indexing of references (STAR/Salmon) if indices are not provided  
- Aggregation of counts into a **gene count matrix**  
- Differential expression analysis using **DESeq2**  
- **MultiQC** report generation

The pipeline is designed to run in **GitHub Codespaces** or any Linux environment with Nextflow.

---

## Folder Structure

```bash
├── data
├── design.tsv
├── envs
│   └── rnaseq.yml
├── main.nf
├── modules
│   ├── deseq2.nf
│   ├── fastqc.nf
│   ├── featurecounts.nf
│   ├── merge_counts.nf
│   ├── multiqc.nf
│   ├── salmon.nf
│   ├── salmon_index.nf
│   ├── star.nf
│   ├── star_index.nf
│   └── trim_galore.nf
├── nextflow.config
├── results
└── samples.tsv
```


---

## Input Files

### 1. `samples.tsv` (required)

Tab-delimited table with sample IDs and FASTQ paths:


|sample_id&emsp;|&emsp;fastq_1&emsp;|&emsp;fastq_2|  
|---------------|-------------------|-------------|
|sampleA&emsp;|&emsp;data/sampleA_R1.fastq.gz&emsp;|&emsp;data/sampleA_R2.fastq.gz| 
|sampleB&emsp;|&emsp;data/sampleB_R1.fastq.gz&emsp;|&emsp;data/sampleB_R2.fastq.gz|
|sampleC&emsp;|&emsp;data/sampleC_R1.fastq.gz&emsp;|&emsp;data/sampleC_R2.fastq.gz| 
|sampleD&emsp;|&emsp;data/sampleD_R1.fastq.gz&emsp;|&emsp;data/sampleD_R2.fastq.gz|

### 2. `design.tsv` (required for DESeq2)

Tab-delimited table specifying sample conditions:

|sample_id&emsp;|&emsp;condition|  
|---------------|-------------------|
|sampleA&emsp;|&emsp;control| 
|sampleB&emsp;|&emsp;control|
|sampleC&emsp;|&emsp;treated| 
|sampleD&emsp;|&emsp;treated|


- `sample_id` must match `samples.tsv`.
- You can add additional covariates as needed.

### 3. Reference files

- Place your reference genome and annotation in `reference/`:

```bash
reference/genome.fa.gz
reference/annotation.gtf.gz
````


---

## Configuration (`nextflow.config`)

Set pipeline parameters:

```groovy
params {
  samples = 'samples.tsv'               // path to sample sheet
  design  = 'design.tsv'                // path to metadata for DESeq2
  outdir  = 'results'                   // output folder
  fasta   = 'reference/genome.fa.gz'    // reference genome
  gtf     = 'reference/annotation.gtf.gz'  // annotation
  salmon_index = null                    // optional prebuilt Salmon index
  star_index   = null                    // optional prebuilt STAR index
}
```

If salmon_index or star_index is null, the pipeline will build them automatically.
outdir is configurable; can point to 
```bash 
/mnt/data
``` 
or another mounted volume.

# Running the Pipeline

1. Using Salmon (alignment-free)

```bash
nextflow run . -profile conda \
  --samples samples.tsv \
  --fasta reference/genome.fa.gz \
  --gtf reference/annotation.gtf.gz \
  --aligner salmon
```

2. Using STAR + featureCounts (alignment-based)

```bash
nextflow run . -profile conda \
  --samples samples.tsv \
  --fasta reference/genome.fa.gz \
  --gtf reference/annotation.gtf.gz \
  --design design.tsv \
  --aligner star
```

#### Optional: specify prebuilt indices:

```bash
--salmon_index reference/salmon/index
--star_index reference/star/index
```

# Output
#### 1. QC reports

results/fastqc/ → per-sample FastQC reports

results/multiqc/ → combined MultiQC HTML summary

#### 2. Quantification

Salmon: results/salmon/<sample>_quant/

STAR: results/star/<sample>.bam + results/counts/<sample>_counts.txt

#### 3. Aggregated counts

results/counts/gene_counts_matrix.txt → merged gene × sample counts

#### 4. Differential expression (STAR branch only)

results/deseq2/deseq2_results.tsv → DESeq2 results table

results/deseq2/deseq2_plots.pdf → MA plot + PCA plot