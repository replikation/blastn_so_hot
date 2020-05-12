process blastn_NCBI {
        label 'biopython'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(fasta)
  	output:
    	tuple val(name), path("${name}.xml")
  	script:
    """
    export PYTHONHTTPSVERIFY=0
	blast_ncbi.py --fasta ${fasta}
	mv my_blast.xml ${name}.xml
    """
}
