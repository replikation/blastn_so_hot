#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Author: christian.jena@gmail.com
*/

if ( !nextflow.version.matches('20.+') ) {
    println "This workflow requires Nextflow version 20.X or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "WtP intended for Nextflow-version: 20.01.0"
println "Starting time: $nextflow.timestamp"
println "Workdir location [--workdir]:"
println "  $workflow.workDir"
println "Output location [--output]:"
println "  $params.output"
println "\033[2mDatabase location [--databases]:"
println "  $params.databases\u001B[0m"
if (workflow.profile.contains('singularity')) {
println "\033[2mSingularity cache location [--cachedir]:"
println "  $params.cachedir"
println "\u001B[33m  WARNING: Singularity image building sometimes fails!"
println "  Rerun WtP via -resume to retry the failed image build"
println "  Manually remove faulty images in $params.cachedir for a rebuild\u001B[0m"
}
println " "
println "\033[2mCPUs to use: $params.cores\033[0m"
println " "

/************* 
* ERROR HANDLING
*************/
// profiles
if ( workflow.profile == 'standard' ) { exit 1, "NO VALID EXECUTION PROFILE SELECTED, use e.g. [-profile local,docker]" }

if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker')
    ) { "engine selected" }
else { exit 1, "No engine selected:  -profile EXECUTER,ENGINE" }

if (
    workflow.profile.contains('local') 
    ) { "executer selected" }
else { exit 1, "No executer selected:  -profile EXECUTER,ENGINE" }

// params tests
if ( !params.fasta ) {
    exit 1, "input missing, use [--fasta]"}


/************* 
* INPUT HANDLING
*************/

// fasta input or via csv file
    if (params.fasta && params.list) { fasta_input_ch = Channel
            .fromPath( params.fasta, checkIfExists: true )
            .splitCsv()
            .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
                }
    else if (params.fasta) { fasta_input_ch = Channel
            .fromPath( params.fasta, checkIfExists: true)
            .map { file -> tuple(file.baseName, file) }
                }

/************* 
* MODULES
*************/



/************* 
* SUB WORKFLOWS
*************/  

workflow blast_against_NCBI {
    take:   fasta
    main:   
    emit:   
}

workflow blast_against_own_DB {
    take:   fasta
    main:   
    emit:   
}

workflow plot_blast_output {
    take:   fasta
    main:   
    emit:   
}

/************* 
* Main Workflow
*************/ 