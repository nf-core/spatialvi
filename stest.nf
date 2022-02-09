#!/usr/bin/env nextflow

nextflow.enable.dsl=2

WorkflowMain.initialise(workflow, params, log)

include { IN_TEST4   } from './modules/local/test_tasks'

Channel
.from(file(params.input))
.splitCsv(header: true)
.map{ row-> tuple(row.sample_id, row.species, row.st_data_dir, row.sc_data_dir) }
.set{input_samples}


process IN_TEST1 {

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    
    output:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir), env(outpath)
    
    """   
    sleep 1
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


process IN_TEST2 {

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir), path(p)
    
    output:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir), env(outpath)
    
    """   
    sleep 6
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


process IN_TEST3 {

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir), path(p)
    val(s)
    
    output:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir), env(outpath)
    
    """   
    sleep 2
    echo main task ${params.dryrun} 
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


process IN_TEST5 {

    echo true
   
    """
    echo IN_TEST5 main task ${params.emptyrun} 
    """
}


workflow {

    IN_TEST5()  
    IN_TEST4()
  
}






