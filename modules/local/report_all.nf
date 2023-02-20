import groovy.json.JsonSlurper

//
// Report
//
process REPORT_ALL {

    // TODO: Create the final report script
    // TODO: Change this process to correspond to final report
    // TODO: Add Conda/container directive
    // TODO: Export versions

    label "process_low"

    input:
    val sample_state
    val outdir

    output:
    tuple env(sample_id), env(outpath)
    // path("versions.yml"), emit: versions

    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))

    """
    #!/bin/bash

    sample_id=${sample_id_gr}

    dname=${outdir}/\${sample_id}

    echo \${dname}/
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}
