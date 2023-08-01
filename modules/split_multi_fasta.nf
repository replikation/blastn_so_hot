process split_multi_fasta {
      label 'seqkit'
    input:
      tuple val(name), path(fasta) 
    output:
      path("${name}/*.fasta") 
    script:
      """
      seqkit split2 --by-size 1 ${fasta} -O ${name}
      """
}