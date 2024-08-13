process VCFANNO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.5--h9ee0642_0':
        'biocontainers/vcfanno:0.3.5--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(specific_resources)
    path toml
    path lua
    path resources

    output:
    tuple val(meta), path("*.vcf")     , emit: vcf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def lua_cmd = lua ? "--lua ${lua}" : ""

    # Download files from S3 if they are referenced in the paths
    if (toml.startsWith("s3://")) {
        """
        aws s3 cp ${toml} ./local_toml_file.toml
        toml="./local_toml_file.toml"
        """
    }
    if (vcf.startsWith("s3://")) {
        """
        aws s3 cp ${vcf} ./local_vcf_file.vcf.gz
        vcf="./local_vcf_file.vcf.gz"
        """
    }
    if (resources.startsWith("s3://")) {
        """
        aws s3 cp ${resources} ./local_resources_dir --recursive
        resources="./local_resources_dir"
        """
    }

    """
    vcfanno \\
        -p ${task.cpus} \\
        ${args} \\
        ${lua_cmd} \\
        ${toml} \\
        ${vcf} \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' '))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' '))
    END_VERSIONS
    """
}
