#!/usr/bin/env python
import altair as alt
import pandas as pd 
import xml.etree.ElementTree as et 
import argparse
from altair_saver import save

alt.renderers.enable('html')

parser = argparse.ArgumentParser(description='plotting of xml blast outputs')

parser.add_argument(
    '--xml', required=True,
    help='blastn xml file')

args = parser.parse_args()
xml_file = args.xml

## PARSE XML

xtree = et.parse(xml_file)
xroot = xtree.getroot() 

df_cols = ["Hit_def", "Hit_accession", "Hit_ID", "hsp_start", "hsp_end", "value", "Hsp_num"]
rows = []

for node in xroot.iter('Hit'):
    # collecting all <Hit> nodes with certain content:
    for hitnode in node: 
        if ( hitnode.tag == 'Hit_def' ): 
          s_def = hitnode.text
        elif ( hitnode.tag == 'Hit_accession' ): 
          s_accession = hitnode.text
        elif ( hitnode.tag == 'Hit_id' ): 
          s_id = hitnode.text
        elif ( hitnode.tag == 'Hit_hsps' ):
            # collecting all the HSPs
            for hsp in hitnode.iter():
                if ( hsp.tag == 'Hsp_query-from' ):
                    s_hsp_start = hsp.text
                elif ( hsp.tag == 'Hsp_query-to' ):
                    s_hsp_end = hsp.text
                elif ( hsp.tag == 'Hsp_num' ):
                    s_hsp_num = hsp.text
            # dropping the collected infos into a row after the Hsp_query-to and per hsp
                elif ( hsp.tag == 'Hsp_align-len' ):
                    value = ( int(s_hsp_end) - int(s_hsp_start) )
                    rows.append({"Hit_def": s_def, "Hit_accession": s_accession, 
                                    "Hit_ID": s_id, "hsp_start": s_hsp_start, 
                                    "hsp_end": s_hsp_end, #"value": value, 
                                    "Hsp_num": s_hsp_num})

out_df = pd.DataFrame(rows, columns = df_cols)

## PLOT

source = out_df

y_axis = alt.Axis(
    title='Accession',
    offset=100,
    ticks=False,
    minExtent=100,
    labelLimit=800,
    domain=False
)

#slider = alt.binding_range(min=1, max=4, step=1)
#select_year = alt.selection_single(name="hsp", fields=['Hsp_num'],
#                                   bind=slider, init={'Hsp_num': 4})

render = alt.Chart(source).mark_bar(opacity=0.8).encode(
    #x='hsp_start:Q',
    x= alt.X('hsp_start:Q', title='Length in bp'),
    x2='hsp_end:Q',
    y=alt.Y('Hit_ID:N', axis=y_axis),
   
    color=alt.Color(
        'Hsp_num:N',
        legend=alt.Legend( title='HSP-hits'),
        )
    ).properties(width=800) #    height=200#.add_selection(select_year)


render2 = alt.layer(render, data=source).facet(
    column='Hsp_num:N'
)




render.save('chart.html')

save(render, "chart.svg")