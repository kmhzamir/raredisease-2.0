# ![https://www.google.com/url?sa=i&url=https%3A%2F%2Fsynapselaboratory.com%2F&psig=AOvVaw2HswlH2hgLQjqqQcnxxWrZ&ust=1720234805818000&source=images&cd=vfe&opi=89978449&ved=0CBEQjRxqFwoTCNC-75j0jocDFQAAAAAdAAAAABAE](https://synapselaboratory.com/wp-content/uploads/2022/12/Synapse-logo.png)

## Introduction

**raredisease-2.0** is a Synapsys pipeline for rare disease analysis. It uses fastp for preprocessing, bwa-mem2 for alignment, MarkDuplicates for duplicate removal, qCBAM for BAM QC, and HaplotypeCaller for variant calling and annotation. It now includes the MT pipeline and supports SNVs and SMNCopyNumberCaller, excluding SVs and GermlineCNVCaller..

## Pipeline summary

<!-- prettier-ignore -->
<p align="center">
    <img title="nf-core/raredisease workflow" src="docs/images/Rare Disease Pipeline Modularization - raredisease 2.0.jpg" width=40%>
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

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/raredisease/output).

## Credits

nf-core/raredisease was written in a collaboration between the Clinical Genomics nodes in Sweden, with major contributions from [Ramprasad Neethiraj](https://github.com/ramprasadn), [Anders Jemt](https://github.com/jemten), [Lucia Pena Perez](https://github.com/Lucpen), and [Mei Wu](https://github.com/projectoriented) at Clinical Genomics Stockholm.

Additional contributors were [Sima Rahimi](https://github.com/sima-r), [Gwenna Breton](https://github.com/Gwennid) and [Emma Västerviga](https://github.com/EmmaCAndersson) (Clinical Genomics Gothenburg); [Halfdan Rydbeck](https://github.com/hrydbeck) and [Lauri Mesilaakso](https://github.com/ljmesi) (Clinical Genomics Linköping); [Subazini Thankaswamy Kosalai](https://github.com/sysbiocoder) (Clinical Genomics Örebro); [Annick Renevey](https://github.com/rannick) and [Peter Pruisscher](https://github.com/peterpru) (Clinical Genomics Stockholm); [Ryan Kennedy](https://github.com/ryanjameskennedy) (Clinical Genomics Lund); [Anders Sune Pedersen](https://github.com/asp8200) (Danish National Genome Center) and [Lucas Taniguti](https://github.com/lmtani).
<iframe width="768" height="432" src="https://miro.com/app/live-embed/uXjVNBRfigU=/?moveToViewport=-5663,-6922,5529,3310&embedId=871456525688" frameborder="0" scrolling="no" allow="fullscreen; clipboard-read; clipboard-write" allowfullscreen></iframe>
We thank the nf-core community for their extensive assistance in the development of this pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#raredisease` channel](https://nfcore.slack.com/channels/raredisease) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/raredisease for your analysis, please cite it using the following doi: [10.5281/zenodo.7995798](https://doi.org/10.5281/zenodo.7995798)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

You can read more about MIP's use in healthcare in,

> Stranneheim H, Lagerstedt-Robinson K, Magnusson M, et al. Integration of whole genome sequencing into a healthcare setting: high diagnostic rates across multiple clinical entities in 3219 rare disease patients. Genome Med. 2021;13(1):40. doi:10.1186/s13073-021-00855-5
