//
// Run SpaceRanger count
//

process SPACERANGER_COUNT {
    tag "${sample}"

    container "nfcore/spaceranger:1.3.0"

    input:
    tuple val(sample), path(fastq_dir), path(image), val(slide), val(area)
    path(reference)
    path(probeset)

    output:
    path "spaceranger-${sample}", type: "dir", emit: sr_dir
    path "versions.yml"                      , emit: versions

    script:
    """
    spaceranger count \
        --id=spaceranger-${sample} \
        --sample=${sample} \
        --fastqs=${fastq_dir} \
        --image=${image} \
        --slide=${slide} \
        --area=${area} \
        --transcriptome=${reference} \
        --probe-set=${probeset} \
        --localcores=${task.cpus}

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SpaceRanger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """
}
