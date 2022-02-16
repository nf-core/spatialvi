#!/usr/bin/env nextflow

nextflow.enable.dsl=2

WorkflowMain.initialise(workflow, params, log)

include { IN_TEST4   } from './modules/local/test_tasks'


Channel
.from(file(params.input))
.splitCsv(header:true, sep:',')
.map{ prep_input_csv_files(it) }
.set{sample_ids}

import org.apache.commons.csv.CSVPrinter
import org.apache.commons.csv.CSVFormat
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def prep_input_csv_files(LinkedHashMap row) {
    
    def fileName = String.format("%s/sample_%s", params.outdir, row.sample_id)    
    def FILE_HEADER = row.keySet() as String[];
    
    new File(fileName + ".csv").withWriter { fileWriter ->
        def csvFilePrinter = new CSVPrinter(fileWriter, CSVFormat.DEFAULT)
        csvFilePrinter.printRecord(FILE_HEADER)
        csvFilePrinter.printRecord(row.values())
    }
    
    File file = new File(fileName + ".json")
    file.write(JsonOutput.toJson(row))

    return row.sample_id
}

process IN_TEST0 {

    echo true

    input:
    val(sample_id)
    
    output:
    tuple val(sample_id), env(outpath) 
    
    script:
    def fileName = String.format("%s/sample_%s.json", params.outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    echo Species:${sample_info.species}
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


workflow {

    IN_TEST0( sample_ids )

    IN_TEST0.out.view()
  
}


def aaa(LinkedHashMap row) {

    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}






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









