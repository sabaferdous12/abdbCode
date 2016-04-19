#!/bin/bash
#*********************************************************************        
# Description: This bash script creates a table for displaying statistics of 
#              unprocessed antibodies on an html page (.tt file). 
#             
#              It reads the data on command line to put in the table          
# Usage:      ./statsUnprocessed.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12    
#                                                                             
#**********************************************************************

args=("$@")

cat << __EOF__
<div class='abdb_stats2'>
<h2>Database Statistics - Unprocessed Data</h2>
<p>The following statistics show the type and number of unprocessed antibody related structures that are not present in AbDb</p>
<table class="mytable">

<tr>
<th valign="middle">Fc Fragments</th>
<th valign="middle">Numbering Failed</th>
<th valign="middle">Superseded</th>
<th valign="middle">Single Chain Antibody (scFVs)</th>
</tr>

<tr>
<td valign="middle">${args[0]}</td>
<td valign="middle">${args[1]}</td>
<td valign="middle">${args[2]}</td>
<td valign="middle">${args[3]}</td>
</tr>

</table>
<p><a href="http://www.bioinf.org.uk/abs/abdb/Data/Kabat_Failed.list">List of antibodies that failed the numbering</a></p>
<p><a href="http://www.bioinf.org.uk/abs/abdb/Data/Martin_logs/multi-chain.list">List of antibodies bound to multiple antigens</a></p>
<p><a href="http://www.bioinf.org.uk/abs/abdb/Data/Martin_logs/scFV.list">List of single chain antibodies (scFVs)</a></p>

</table>
</div>

__EOF__
