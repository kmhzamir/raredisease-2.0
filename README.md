# ![https://www.google.com/url?sa=i&url=https%3A%2F%2Fsynapselaboratory.com%2F&psig=AOvVaw2HswlH2hgLQjqqQcnxxWrZ&ust=1720234805818000&source=images&cd=vfe&opi=89978449&ved=0CBEQjRxqFwoTCNC-75j0jocDFQAAAAAdAAAAABAE](https://synapselaboratory.com/wp-content/uploads/2022/12/Synapse-logo.png)

## Introduction

**raredisease-2.0** is a Synapsys pipeline for rare disease analysis. It uses fastp for preprocessing, bwa-mem2 for alignment, MarkDuplicates for duplicate removal, qCBAM for BAM QC, and HaplotypeCaller for variant calling and annotation. It now includes the MT pipeline and supports SNVs and SMNCopyNumberCaller, excluding SVs and GermlineCNVCaller.

## Pipeline summary

<!-- prettier-ignore -->
<p align="center">
    <img title="nf-core/raredisease workflow" src="docs/images/Rare Disease Pipeline Modularization (1).jpg">
</p>



## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,lane,fastq_1,fastq_2,sex,phenotype,paternal_id,maternal_id,case_id
hugelymodelbat,1,reads_1.fastq.gz,reads_2.fastq.gz,1,2,,,justhusky
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Second, ensure that you have defined the path to reference files and parameters required for the type of analysis that you want to perform. More information about this can be found [here](https://github.com/nf-core/raredisease/blob/dev/docs/usage.md).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/raredisease \
   -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/raredisease/usage) and the [parameter documentation](https://nf-co.re/raredisease/parameters).
