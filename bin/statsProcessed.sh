#!/bin/bash
#*********************************************************************
# Description: This bash script creates a table for displaying statistics of 
#              processed antibodies on an html page (.tt file).
#              It reads the data on command line to put in the table 
#
# Usage:      ./statsProcessed.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12
#
#**********************************************************************
#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ' -> echo $1 $2 $3 $4 $5 $6 $7 $8 $9'

args=("$@")


cat << __EOF__
<div class='abdb_stats'>
<h2>Database Statistics - Processed Data</h2>

<table class="mytable">
<tr>
<th>Datasets</th>
<th>Complex Type</th>
<th>Processed PDB Files </th>
<th>Resultant Antibodies </th>
<th>Non-Redundant Antibodies</th>
</tr>


<tr class='firstrow'>
<td rowspan="4" valign="top">Complete Antibody</td>
<td>Protein</td>
<td>${args[0]}</td>
<td>${args[1]}</td>
<td>${args[2]}</td>
</tr>

<tr>                                          
<td>Non-protein</td>                                                   
<td>${args[3]}</td>
<td>${args[4]}</td>
<td>${args[5]}</td>
</tr>

<tr>
<td>Free Antibody</td>
<td>${args[6]}</td>
<td>${args[7]}</td>
<td>${args[8]}</td>
</tr>

<tr class='secondrow'>
<td>Complete Dataset</td>
<td>${args[9]}</td>
<td>${args[10]}</td>
<td>${args[11]}</td>
</tr>


<tr class='firstrow'>
 <td rowspan="4" valign="top">Light Chains</td>
 <td>Protein</td>
 <td>${args[12]}</td>
 <td>${args[13]}</td>
 <td>${args[14]}</td>
</tr>

<tr>
<td>Non-protein</td>
<td>${args[15]}</td>
<td>${args[16]}</td>
<td>${args[17]}</td>
</tr>

<tr>
<td>Light Only</td>
<td>${args[18]}</td>
<td>${args[19]}</td>
<td>${args[20]}</td>
</tr>

<tr class='secondrow'>
<td>Complete Dataset</td>
<td>${args[21]}</td>
<td>${args[22]}</td>
<td>${args[23]}</td>
</tr>

<tr class='firstrow'>
<td rowspan="4" valign="top">Heavy Chains</td>
<td>Protein</td>
<td>${args[24]}</td>
<td>${args[25]}</td>
<td>${args[26]}</td>
</tr>

<tr>
<td>Non-protein</td>
<td>${args[27]}</td>
<td>${args[28]}</td>
<td>${args[29]}</td>
</tr>

<tr>
<td>Heavy Only</td>
<td>${args[30]}</td>
<td>${args[31]}</td>
<td>${args[32]}</td>
</tr>

<tr class='secondrow'>
<td>Complete Dataset</td>
<td>${args[33]}</td>
<td>${args[34]}</td>
<td>${args[35]}</td>
</tr>


</table>

</div>

__EOF__
