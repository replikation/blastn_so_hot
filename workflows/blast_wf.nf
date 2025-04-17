include { makeblastDB }         from './process/blast.nf'
include { blastn_local }        from './process/blast.nf'
include { blastn_NCBI }         from './process/blastn_NCBI.nf'
include { plot_xml }            from './process/plot_xml.nf'
include { split_multi_fasta }   from './process/split_multi_fasta.nf'


                     
//             // plasflow bekommt nur die fasta zuerst, dann auf den motif ch matchen
//             plasflow_in_ch=motif_ch.map( {it -> tuple(it[0], it[3])}) // find out what is fasta
//             plasflow()


/************* 
* DATABASES
*************/


workflow make_blast_DB_wf {
    take:   references_in_ch
    main:   makeblastDB(references_in_ch)   
    emit:   makeblastDB.out
}

/************* 
* WORKFLOWS
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
