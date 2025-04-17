process plot_xml {
        label 'altair'
        publishDir "${params.output}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(xml)
  	output:
    	tuple val(name), path("${name}.html"), path("chart.svg") 
  	script:
    """
	parse_and_plot.py --xml ${xml}
	mv chart.html ${name}.html
    """
}
