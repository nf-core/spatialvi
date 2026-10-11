process STARDIST {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d964e0bef867bb2ff1a309c9c087d8d83ac734ce3aa315dd8311d4c1bfdafd8e/data' :
        'community.wave.seqera.io/library/python_pip_imagecodecs_nvidia-cublas-cu12_pruned:b668bcb6d531d350' }"

    input:
    tuple val(meta), path(image)
    tuple val(model_name), path(model_path)

    output:
    tuple val(meta), path("*.stardist.tif"), emit: mask
    tuple val("${task.process}"), val('stardist'), eval("pip show stardist | sed -n 's/^Version: //p'"), topic: versions, emit: versions_stardist
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tensorflow'), eval("pip show tensorflow | sed -n 's/^Version: //p'"), topic: versions, emit: versions_tensorflow
    tuple val("${task.process}"), val('tifffile'), eval("pip show tifffile | sed -n 's/^Version: //p'"), topic: versions, emit: versions_tifffile

    when:
    task.ext.when == null || task.ext.when

    // The `stardist-predict2d` script declares `--n_tiles` with `nargs=2`, so it
    // cannot tile RGB input, which needs one tile count per axis. Untiled, a
    // full-size H&E allocates its ray array whole and is OOM-killed.
    script:
    // Extra keyword arguments to `predict_instances`, e.g. 'prob_thresh=0.6, nms_thresh=0.3'
    def args      = task.ext.args ? ", ${task.ext.args}" : ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def model_ref = model_path ? "${model_path}" : "${model_name}"
    """
    python <<'PY'
    import logging
    import pathlib

    import numpy as np
    from csbdeep.utils import normalize
    from imageio import imread
    from stardist.models import StarDist2D
    from tifffile import imwrite

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logger = logging.getLogger(__name__)

    img = imread("${image}")
    if img.ndim not in (2, 3):
        raise ValueError("only 2d and 3d images are supported, got %dd" % img.ndim)

    # Roughly 2000 px per tile
    n_tiles = max(1, int(np.ceil(max(img.shape[:2]) / 2000)))
    # One tile count per axis; never tile the channel axis of YXC input.
    tiles = (n_tiles, n_tiles) if img.ndim == 2 else (n_tiles, n_tiles, 1)
    logger.info(f"Image {img.shape[1]} x {img.shape[0]}, n_tiles={tiles}")

    model_ref = "${model_ref}"
    if pathlib.Path(model_ref).is_dir():
        model = StarDist2D(None, name=model_ref)
    else:
        model = StarDist2D.from_pretrained(model_ref)
    if model is None:
        raise ValueError("unknown model: %s" % model_ref)

    img = normalize(img, 1, 99.8)
    labels, _ = model.predict_instances(img, n_tiles=tiles${args})
    imwrite("${prefix}.stardist.tif", labels, compression="zlib")
    PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stardist.tif
    """
}
