#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Author: christian.jena@gmail.com
*/


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
    if (params.fasta) { 
        fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

// references
    if (params.references) { 
        references_input_ch = Channel
        .fromPath( params.references, checkIfExists: true)
    }

/************* 
* MODULES
*************/

include { blastn_NCBI                       } from './modules/blastn_NCBI' 
include { makeblastDB; blastn_local         } from './modules/blast' 
include { plot_xml                          } from './modules/plot_xml' 
include { split_multi_fasta                 } from './modules/split_multi_fasta'

/************* 
* DATABASES
*************/

workflow make_blast_DB {
    take:   references
    main:   makeblastDB(references)   
    emit:   makeblastDB.out
}

/************* 
* SUB WORKFLOWS
*************/  



workflow blast_against_NCBI_wf {
    take:   fasta
    main:   if (params.multifasta) {
                split_multi_fasta(fasta)
                mapped_channel=split_multi_fasta.out.flatten().map{it -> [it.baseName, it]}
                blastn_NCBI(mapped_channel)

                // report
                report_ch = blastn_NCBI.out.status.view { name, status -> "$name got NCBI response: $status" }

            }
            else {
                blastn_NCBI(fasta)

                report_ch = blastn_NCBI.out.status.view { name, status -> "$name got NCBI response: $status" }
            
            }
    
            
    emit:   blastn_NCBI.out.xml
}

workflow blast_against_own_DB_wf {
    take:   fasta
            database
    main:   if (params.multifasta) {
                split_multi_fasta(fasta)
                mapped_channel=split_multi_fasta.out.flatten().view() //.map { it -> tuple(it[1].baseName, it[1]) }
                blastn_local(mapped_channel, database)
            }
            else {
                blastn_local(fasta, database)
            }
    emit:   blastn_local.out
}

workflow plot_blast_output_wf {
    take:   xml
    main:   plot_xml(xml)   
    emit:   plot_xml.out
}

/************* 
* Main Workflow
*************/ 

workflow {
    if (params.fasta && !params.references) {
        plot_blast_output_wf(
            blast_against_NCBI_wf(fasta_input_ch))
    }

    if (params.fasta && params.references) {
        plot_blast_output_wf(blast_against_own_DB_wf(fasta_input_ch, make_blast_DB(references_input_ch)))
    }
}

/*************  
* --help
*************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    ${c_green}blastn so hot${c_reset} | A Nextflow blastn workflow for convinience
    
    ${c_yellow}Usage examples:${c_reset}
    nextflow run replikation/blastn_so_hot --fasta 'sample_01.fasta' --references references.fasta --cores 14 -profile local,docker
    ${c_dim}or${c_reset}
    nextflow run replikation/blastn_so_hot --fasta 'sample_01.fasta' --references references.fasta --cores 14 -profile local,singularity

    ${c_yellow}Inputs (choose one):${c_reset}

    --fasta         your query sequence

    --references    multiple references (will be converted to blast DB)


    ${c_yellow}Parameters${c_reset}
    --eValue        (not implemented)
    --hsplength     (not implemented)
    --multifasta    multifasta file, analysed each sequence separately 


    ${c_yellow}Options:${c_reset}
    --cores         max cores for local use [default: $params.cores]
    --memory        available memory [default: $params.memory]
    --output        name of the result folder [default: $params.output]
    --cachedir      defines the path where singularity images are cached
                    [default: $params.cachedir] 

    ${c_yellow}Execution/Engine profiles:${c_reset}
    This workflow supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.:
     -profile ${c_green}local${c_reset},${c_blue}docker${c_reset}

      ${c_green}Executer${c_reset} (choose one):
      local
      ${c_blue}Engines${c_reset} (choose one):
      docker
      singularity
    """.stripIndent()
}
