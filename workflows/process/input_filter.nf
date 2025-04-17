process input_filter {
      label 'seqkit'
    input:
      tuple val(name), path(fasta) 
    output:
      path("${name}.fasta") 
    script:
      """

      
      mv ${fasta} "${name}".fasta
      """
}