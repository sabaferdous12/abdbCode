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
<p>The following statistics show the type and number of antibody structures that are not present in AbDb at the moment</p>
<table class="mytable">

<tr>
<th valign="middle">Light Only</th>
<th valign="middle">Light-Antigen</th>
<th valign="middle">Heavy Only</th>
<th valign="middle">Heavy-Antigen</th>
<th valign="middle">Fc Fragments</th>
<th valign="middle">Numbering Failed</th>
<th valign="middle">Superseded</th>

</tr>

<tr>
<td valign="middle">${args[0]}</td>
<td valign="middle">${args[1]}</td>
<td valign="middle">${args[2]}</td>
<td valign="middle">${args[3]}</td>
<td valign="middle">${args[4]}</td>
<td valign="middle">${args[5]}</td>
<td valign="middle">${args[6]}</td>
</tr>

</table>
<p><a href="http://www.bioinf.org.uk/abs/abdb/Data/Kabat_Failed.list">List of antibodies that failed the numbering</a></p>
</div>

</table>

</div>

__EOF__
