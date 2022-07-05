//
// Run SpaceRanger count
//

process SPACERANGER_COUNT {
    tag "${meta.id}"

    container "nfcore/spaceranger:1.3.0"

    input:
    tuple val(meta), path(fastq_dir), path(image), val(slide), val(area)
    path(reference)
    path(probeset)

    output:
    path "spaceranger-${meta.id}", type: "dir", emit: sr_dir
    path "versions.yml"                       , emit: versions

    script:
    """
    spaceranger count \
        --id=spaceranger-${meta.id} \
        --sample=${meta.id} \
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
