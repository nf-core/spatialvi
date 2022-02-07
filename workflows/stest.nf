#!/usr/bin/env nextflow

process TEST_REPORT1 {
     
    label "python_process"
    
    output:
    val 'state1' into test_channel_1
    
    script:  
    """
    #!/bin/bash
      
    python --version
    echo Start1
    sleep 5   
    echo End1
    echo 0
    """
}

process TEST_REPORT2 {
     
    label "r_process"
    
    input:
    val state from test_channel_1
    
    output:
    val 'state2' into test_channel_2
    
    script:  
    """
    #!/bin/bash
      
    R --version
    echo Start2
    sleep 7   
    echo End2
    echo 0
    """
}