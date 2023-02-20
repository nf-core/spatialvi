//
// Run SpaceRanger count
//

process SPACERANGER_COUNT {

    tag "${meta.id}"
    label "process_high"

    container "nfcore/spaceranger:1.3.0"

    input:
    tuple val(meta), path(fastq_dir), path(image), val(slide), val(area)
    path(reference)
    path(probeset)

    output:
    path "spaceranger-${meta.id}", type: "dir"               , emit: sr_dir
    tuple val  (meta),
        path ("*/outs/spatial/tissue_positions_list.csv"),
        path ("*/outs/spatial/tissue_lowres_image.png"),
        path ("*/outs/spatial/tissue_hires_image.png"),
        path ("*/outs/spatial/scalefactors_json.json"),
        path ("*/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),
        path ("*/outs/raw_feature_bc_matrix/features.tsv.gz"),
        path ("*/outs/raw_feature_bc_matrix/matrix.mtx.gz")  , emit: sr_out
    path "versions.yml"                                      , emit: versions

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
