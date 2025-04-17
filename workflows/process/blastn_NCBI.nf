process blastn_NCBI {
        label 'biopython'
        publishDir "${params.output}/${name}/", mode: 'copy'
        maxForks 1
    input:
        tuple val(name), path(fasta)
  	output:
    	tuple val(name), path("${name}.xml"), emit: xml
        tuple val(name), env(STATUS), emit: status
  	script:
    """
    export PYTHONHTTPSVERIFY=0
	blast_ncbi.py --fasta ${fasta}
	mv my_blast.xml ${name}.xml

    STATUS=\$(grep "<Iteration_message>" ${name}.xml | cut -f2 -d ">" | cut -f1 -d "<")
    """
}
