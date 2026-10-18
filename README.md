# RNA-seq Initiation: From Raw Reads to a Gene Count Matrix

A step-by-step Google Colab notebook ([`RNAseq_initiation.ipynb`](RNAseq_initiation.ipynb)) that runs a complete bulk RNA-seq preprocessing pipeline. It downloads public sequencing data from NCBI SRA, checks read quality, aligns the reads to the human genome, and produces a gene-level count matrix that is ready for differential expression analysis.

Each step has a short explanation alongside the code, so the notebook also works as a learning resource for anyone new to RNA-seq.

## Biological Context

The data come from **Machino *et al.* (2023), *Experimental & Molecular Medicine*: "Integrative analysis reveals early epigenetic alterations in high-grade serous ovarian carcinomas"** ([doi:10.1038/s12276-023-01090-1](https://doi.org/10.1038/s12276-023-01090-1)).

High-grade serous ovarian carcinoma (HGSOC) is the most lethal gynaecological cancer. Many HGSOCs are now thought to start in the **secretory epithelial cells of the fallopian tube**, passing through a precursor lesion called serous tubal intraepithelial carcinoma (STIC). The mutations and copy-number changes in HGSOC are well described. Much less is known about the early **epigenetic** and transcription-factor changes that drive the disease.

To model the first steps of tumour formation, the authors started from immortalised human fallopian tube secretory epithelial cells (**HF1**). They then added oncogenic changes one at a time:

```
HFTSEC (primary, normal)
   └─ HF1 (immortalised)
        └─ HF1_TP53         (TP53 alteration)
             └─ HF1_TP53_KRAS
                  ├─ HF1_TP53_KRAS_AKT
                  └─ HF1_TP53_KRAS_MYC
```

They profiled each stage with **RNA-seq**, **H3K27ac ChIP-seq** (active enhancers) and **ATAC-seq** (open chromatin). They also used siRNA to knock down candidate regulators (**MAF**, **GATA6**, **DAB2**) in HF1 cells. The main findings were:

- During early tumour formation, the **AP-1 complex** (JUN, FOSL2) became more active and **GATA family** binding changed.
- **GATA6** and **DAB2** were lost in STIC lesions. This was linked to epithelial–mesenchymal transition (EMT) and lower proteasome expression.
- Oncogenic signalling epigenetically silenced a **cadherin gene cluster** near the FRA16D fragile site.
- Inhibiting the proteasome increased **chemokine** gene expression, which may help tumours evade the immune system.
- **MEK inhibition** reversed these oncogenic changes. This suggests MEK inhibitors as a possible treatment for a subgroup of chemotherapy-resistant HGSOC.

### What this notebook covers

This notebook reprocesses **part of the RNA-seq** from that study, starting from the raw reads. By default it covers two stages of the model:

| Run | Sample | Stage |
|-----|--------|-------|
| DRR416925 | HF1_TP53_KRAS_RNA_REP2 | TP53 + KRAS |
| DRR416926 | HF1_TP53_KRAS_RNA_REP3 | TP53 + KRAS |
| DRR416927 | HF1_TP53_KRAS_AKT_RNA_REP1 | TP53 + KRAS + AKT |
| DRR416928 | HF1_TP53_KRAS_AKT_RNA_REP2 | TP53 + KRAS + AKT |
| DRR416929 | HF1_TP53_KRAS_AKT_RNA_REP3 | TP53 + KRAS + AKT |

The resulting count matrix can be used to look at the effect of adding **AKT** on top of TP53 + KRAS. If you add replicate 1 of the TP53 + KRAS stage (DRR416924), you get a balanced 3 vs 3 design.

## Pipeline Overview

```
NCBI SRA (DRA015121)
        │  prefetch + fasterq-dump
        ▼
  Paired-end FASTQ ──► FastQC (quality control)
        │
        │  STAR (splice-aware alignment to GRCh38)
        ▼
  Sorted BAM files
        │
        │  featureCounts (Ensembl 109 GTF)
        ▼
  Gene count matrix ──► pandas / seaborn (exploration) ──► CSV export
```

| Step | Description | Tool |
|------|-------------|------|
| 1 | Install bioinformatics dependencies | `apt`, source build of STAR |
| 2 | Download the human reference genome (FASTA) and annotation (GTF) | Ensembl release 109, GRCh38 |
| 3 | Build the STAR genome index | STAR 2.7.10b |
| 4 | Query NCBI and download raw reads for project `DRA015121` | NCBI E-utilities, SRA Toolkit |
| 4.5 | Read quality control | FastQC |
| 5 | Align reads and count reads per gene | STAR, Subread `featureCounts` |
| 6 | Load the counts, plot them, and export to CSV | pandas, matplotlib, seaborn |

## Dataset

- **Submission:** [`DRA015121`](https://www.ncbi.nlm.nih.gov/sra/?term=DRA015121), BioProject [`PRJDB14752`](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB14752) (DDBJ/NCBI SRA)
- **Organism:** *Homo sapiens*, cultured fallopian tube secretory epithelial cells
- **Library:** RNA-seq, SMART-Seq v4 Ultra Low Input, Illumina NovaSeq 6000, paired-end
- **Runs:** The submission has 46 runs: 30 RNA-seq (`DRR416915` to `DRR416944`), 10 ChIP-seq and 6 ATAC-seq. The notebook processes **5** of the RNA-seq runs (listed above) to stay within Colab's disk and time limits.

> **Note:** The NCBI `esearch` call in Step 4 returns only 20 IDs by default. The notebook then takes the first 5 of those. The five samples you end up with therefore depend on the order NCBI returns them in, not on the experimental design. To pick samples on purpose, set the list of runs yourself:
>
> ```python
> run_accessions = ["DRR416924", "DRR416925", "DRR416926",   # HF1_TP53_KRAS  rep 1-3
>                   "DRR416927", "DRR416928", "DRR416929"]   # HF1_TP53_KRAS_AKT rep 1-3
> for run in run_accessions:
>     ...
> ```
>
> Alternatively, add `&retmax=100` to the `esearch` URL so that every run in the project is returned.

## Requirements

### Compute

The STAR index for the human genome needs **about 30–40 GB of RAM**. The notebook's metadata requests a **High-RAM** runtime (`machine_shape: hm`).

- **Google Colab:** Use Colab Pro/Pro+ with a High-RAM runtime (*Runtime → Change runtime type*).
- **Disk:** The reference genome, the index and the FASTQ files take up a lot of space. The pipeline deletes FASTQ and BAM files once it has finished with them, to save room.
- **Time:** Building the index and aligning the reads can take several hours.

### Software

The notebook installs everything automatically:

- [SRA Toolkit](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [STAR 2.7.10b](https://github.com/alexdobin/STAR) (compiled from source)
- [Subread / featureCounts](https://subread.sourceforge.net/)
- Python: `pandas`, `numpy`, `matplotlib`, `seaborn`

## How to Run

1. Open the notebook in Google Colab:
   - Upload `RNAseq_initiation.ipynb` at [colab.research.google.com](https://colab.research.google.com), **or**
   - Open it directly from GitHub: *File → Open notebook → GitHub* and paste this repository's URL.
2. Switch to a **High-RAM** runtime.
3. Run the cells in order, from top to bottom (*Runtime → Run all*).
4. When it finishes, the notebook downloads the count matrices to your computer.

> **Running outside Colab:** All paths are hard-coded to `/content/...`, and the export cells use `google.colab.files`. On a local Linux machine or an HPC cluster, replace `/content` with your working directory, and replace `files.download(...)` with a plain `to_csv(...)`.

## Outputs

| File | Description |
|------|-------------|
| `/content/fastqc_results/*.html` | FastQC quality report for each FASTQ file |
| `/content/gene_counts.txt` | Raw `featureCounts` output, plus a `.summary` file with assignment statistics |
| `counts_df.csv` | Full count table: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, plus one column per sample |
| `final_counts_df.csv` | Simplified matrix: `Geneid`, `Length`, plus one column per sample (e.g. `DRR416925`, …) |

Genes are identified by **Ensembl gene IDs** (e.g. `ENSG00000160072`).

The notebook also plots two figures for each sample:
- a histogram of log2(counts + 1) across all genes
- the 20 most highly expressed genes

## Next Steps

`final_counts_df.csv` contains **raw (unnormalised) counts**, which is the input that differential expression tools expect:

- [DESeq2](https://bioconductor.org/packages/DESeq2/) or [edgeR](https://bioconductor.org/packages/edgeR/) in R
- [PyDESeq2](https://github.com/owkin/PyDESeq2) in Python

You can use the `Length` column to calculate length-normalised values such as TPM or FPKM.

## Repository Structure

```
.
├── RNAseq_initiation.ipynb   # Main pipeline notebook
└── README.md                 # This file
```

## Author

**Walusimbi Uthman**

## Acknowledgements

- Data source: Machino H., Dozen A., Konaka M., *et al.* (2023). Integrative analysis reveals early epigenetic alterations in high-grade serous ovarian carcinomas. *Exp Mol Med* 55(10):2205–2219. [doi:10.1038/s12276-023-01090-1](https://doi.org/10.1038/s12276-023-01090-1)
- Reference genome and annotation: [Ensembl](https://www.ensembl.org) (GRCh38, release 109)
- Raw data: [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) / [DDBJ](https://www.ddbj.nig.ac.jp)
- Dobin A. *et al.* (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15–21.
- Liao Y., Smyth G.K., Shi W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30(7):923–930.
