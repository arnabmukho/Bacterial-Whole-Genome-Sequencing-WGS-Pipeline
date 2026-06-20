# 🦠 Bacterial Whole Genome Sequencing Pipeline

A robust, fully containerized, high-throughput bacterial whole-genome sequencing (WGS) pipeline for assembly, annotation, taxonomic profiling, antimicrobial resistance characterization, and genome visualization from paired-end sequencing data.

Designed for reproducible microbial genomics research, genomic surveillance programs, and large-scale bacterial sequencing.

---

# ✨ Features

* End-to-end bacterial genome analysis from raw FASTQ files
* Fully containerized workflow (Docker / Singularity compatible)
* Batch processing of multiple samples
* Thread-safe parallel execution
* Automated QC, assembly, annotation, typing, AMR profiling, and visualization
* HPC-friendly design
* Modular and customizable workflow
* Reproducible software environment
* Publication-quality outputs

---

# 🧬 Workflow Overview

```text
Raw FASTQ
    │
    ▼
FastQC
    │
    ▼
fastp
    │
    ▼
FastQC (trimmed)
    │
    ▼
SPAdes Assembly
    │
    ▼
Contig Filtering
    │
    ▼
CheckM
    │
    ▼
Prokka Annotation
    │
    ├────────► Kraken2
    │
    ├────────► MLST
    │
    ├────────► Abricate
    │           ├─ ResFinder
    │           ├─ CARD
    │           └─ PlasmidFinder
    │
    ├────────► MacSyFinder (TXSScan)
    │
    ├────────► antiSMASH
    │
    └────────► CGView
                    │
                    ▼
               MultiQC
```

---

# 🚀 Quick Start

## 1. Pull Docker Image

```bash
docker pull arnabmukho/bacterial_wgs:1.1.0-ubuntu
```

---

## 2. Clone Repository

```bash
https://github.com/arnabmukho/Bacterial-Whole-Genome-Sequencing-WGS-Pipeline.git

cd Bacterial-Whole-Genome-Sequencing-WGS-Pipeline
```

---

## 3. Prepare Input Files

Place paired-end Illumina reads in the working directory:

```text
Sample1_R1.fastq.gz
Sample1_R2.fastq.gz

Sample2_R1.fastq.gz
Sample2_R2.fastq.gz
```

---

## 4. Configure Databases

Update paths in:

```bash
bacterial_wgs_pipeline.sh
```

Example:

```bash
export KRAKEN_DB=/path/to/kraken2
export CHECKM_DB=/path/to/checkm
export ANTISMASH_DB=/path/to/antismash
```

---

# 📥 Database Setup (One-Time Installation)

The Docker image contains all required software for downloading and configuring databases. No local installation of Kraken2, CheckM, antiSMASH, Abricate, MLST, or associated dependencies is required.

Create a database directory:

```bash
mkdir -p db/{kraken2,checkm,antismash}
```

## Automated Database Installation (Recommended)

```text
setup_databases.sh
```
Run:

```bash
chmod +x setup_databases.sh

./setup_databases.sh
```

## 5. ▶️ Run Pipeline

### Docker

```bash
docker run --rm \
-v $(pwd):/data \
-v /path/to/databases:/databases \
-w /data \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
bash bacterial_wgs_pipeline.sh
```

### Singularity / Apptainer

```bash
singularity exec \
--bind $(pwd):$(pwd) \
--bind /path/to/databases:/databases \
bacterial_wgs.sif \
bash bacterial_wgs_pipeline.sh
```

---

# 🔬 Pipeline Components

## 1. Quality Control

### FastQC

Evaluates:

* Per-base sequence quality
* GC content
* Adapter contamination
* Sequence duplication

### fastp

Performs:

* Adapter trimming
* Quality filtering
* PolyG removal
* Read pruning

---

## 2. Genome Assembly

### SPAdes

Assembly performed using:

```bash
--isolate
```

Optimized for bacterial isolate genomes.

### Contig Filtering

Custom filtering removes:

```text
Contigs < 1000 bp
```

and normalizes FASTA headers for downstream compatibility.

---

## 3. Assembly Quality Assessment

### CheckM

Reports:

* Completeness
* Contamination
* Strain heterogeneity

using lineage-specific marker genes.

---

## 4. Genome Annotation

### Prokka

Rapid bacterial genome annotation including:

* CDS prediction
* rRNA annotation
* tRNA annotation
* Functional annotation

Outputs:

```text
.gff
.gbk
.faa
.ffn
.tsv
```

---

## 5. Taxonomic Classification

### Kraken2

Exact k-mer based taxonomic assignment.

Outputs:

* Classification report
* Taxonomic abundance profile

### Interactive Sankey Visualization

Custom Python + Plotly workflow generates:

* Full taxonomic Sankey diagram
* Filtered-depth Sankey diagram

Useful for:

* Contamination assessment
* Mixed culture detection
* Taxonomic exploration

---

## 6. Molecular Typing

### MLST

Assigns:

* Sequence Type (ST)
* Allelic profile

using PubMLST schemes.

---

## 7. Antimicrobial Resistance Analysis

### Abricate

Database screening against:

#### ResFinder

Detection of acquired resistance genes.

#### CARD

Comprehensive Antibiotic Resistance Database profiling.

#### PlasmidFinder

Plasmid replicon identification.

---

## 8. Secretion System Detection

### MacSyFinder (TXSScan)

Detection of:

* Type I Secretion System
* Type II Secretion System
* Type III Secretion System
* Type IV Secretion System
* Type V Secretion System
* Type VI Secretion System
* Type VII Secretion System
* Type IX Secretion System

and associated bacterial appendages.

---

## 9. Secondary Metabolite Discovery

### antiSMASH v7.1.0

Detection and annotation of:

* NRPS
* PKS
* RiPPs
* Terpenes
* Siderophores
* Hybrid clusters

Uses version-matched antiSMASH databases.

---

## 10. Circular Genome Visualization

### CGView

Publication-ready circular genome maps displaying:

* Coding sequences
* rRNA/tRNA genes
* GC content
* GC skew

Generated directly from Prokka annotations.

---

## 11. MultiQC Reporting

Aggregates results from:

* FastQC
* fastp
* CheckM
* Kraken2
* Prokka

into a single interactive HTML report.

---

# 📂 Output Structure

```text
SAMPLE_WGS/

├── QC/
├── TRIM/
├── ASSEMBLY/
├── CHECKM/
├── CHECKM_bins/
├── PROKKA/
├── KRAKEN/
├── MLST/
├── AMR/
├── PLASMID/
├── TXSSCAN/
├── ANTISMASH/
├── PLOTS/
├── TABLES/
├── LOGS/
└── multiqc_report.html
```

---

# ⚙️ Resource Configuration

Edit at top of:

```bash
bacterial_wgs_pipeline.sh
```

```bash
THREADS=6
MAX_PARALLEL=5
```

Default behavior:

```text
5 simultaneous samples
6 threads per sample
30 active threads total
```

Suitable for:

* Workstations
* HPC clusters
* Large surveillance studies

---

# 🐳 Docker Image

Image:

```text
arnabmukho/bacterial_wgs:1.1.0-ubuntu
```

Includes:

* FastQC
* fastp
* SPAdes
* CheckM
* Prokka
* Kraken2
* MLST
* Abricate
* MacSyFinder
* antiSMASH 7.1.0
* CGView
* MultiQC
* GNU Parallel
* BLAST+
* Java Runtime

Pre-configured under Ubuntu 22.04.

---

# 💻 HPC Compatibility

Tested with:

* Docker
* Singularity
* Apptainer
* SLURM environments
* Multi-user Linux clusters

Recommended for production deployment using Singularity/Apptainer.

---

# 📖 Citation

If you use this pipeline in your research, please cite the individual software packages used within the workflow as well as this repository.
