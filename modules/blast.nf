process makeblastDB {
        label 'blast'
    input:
        path(references) 
    output:
	    path("blast_database", type: 'dir') 
    script:
    """
    makeblastdb -in ${references} -dbtype nucl -parse_seqids -out blast_database -title reference_db
    mkdir blast_database && mv blast_database.* blast_database/
    """
}

process blastn_local {
        label 'blast'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta)
        path(database)
    output:
	    tuple val(name), path("${name}.xml") 
    script:
    """
    blastn -query ${fasta} -db ${database}/${database} -out ${name}.xml -outfmt 5 -num_threads ${task.cpus} -evalue 10E-120 -qcov_hsp_perc 10 -max_hsps 10
    """
}