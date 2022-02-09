#!/usr/bin/env nextflow

nextflow.enable.dsl=2

WorkflowMain.initialise(workflow, params, log)

Channel
.from(file(params.input))
.splitCsv(header: true)
.map{ row-> tuple(row.sample_id, row.species, row.st_data_dir, row.sc_data_dir) }
.set{input_samples}


process IN_TEST1 {

    echo false

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val(s)
    
    output:
    tuple val('state_IN1'), val(sample_id)
    
    """
    
    echo IN_TEST1 $sample_id $species $s
    
    """
}

process IN_TEST2 {

    echo false

    input: 
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val(s)
    
    output:
    tuple val('state_IN2'), val(sample_id)
    
    """
    
    echo IN_TEST2 $s
    
    """
}

process IN_TEST3 {

    echo false

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val(s)
    
    output:
    tuple val('state_IN3'), val(sample_id)
    
    """
    
    echo IN_TEST3 $sample_id $species $s
    
    """
}


workflow W1 {

    take:
      input_
    
    main:
      IN_TEST1( input_, 'A' )
      IN_TEST2( input_, IN_TEST1.out )
      IN_TEST3( input_, IN_TEST1.out )
      
    emit:
      IN_TEST2.out
      .combine(IN_TEST3.out)
      .combine(IN_TEST1.out)

}

process IN_TEST4 {

    echo false

    input: 
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val(s)
    
    output:
    tuple val('state_IN4'), val(sample_id)
    
    """
    
    echo IN_TEST4 $s
    
    """
}

process IN_TEST5 {

    echo false

    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val(s)
    
    output:
    tuple val('state_IN5'), val(sample_id)
    
    """
    
    echo IN_TEST5 $sample_id $species $s
    
    """
}


workflow W2 {

    take:
      input_
      out_
    
    main:
      IN_TEST4( input_, out_ )
      IN_TEST5( input_, IN_TEST4.out )
      
    emit:
      IN_TEST5.out

}


workflow {

    W1( input_samples )
    W2( input_samples, W1.out )
    
    W2.out.view()
    
}


