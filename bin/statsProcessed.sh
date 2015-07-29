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
<p>The database only contains structures of free and complexed antibodies where both heavy and light chains are pres\
ent.</p>

<table class="mytable">
<tr>
<th>Antibody Structures</th>
<th>Total Processed PDBs </th>
<th>Total Resultant Complexes </th>
<th>Non-Redundant Antibodies</th>
</tr>


<tr>
<td>Complexes (Protein)</td>
<td>${args[0]}</td>
<td>${args[1]}</td>
<td>${args[2]}</td>
</tr>

<tr>                                          
<td>Complexes (Non-protein)</td>                                                   
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

<tr>
<td>Complete Dataset</td>
<td>${args[9]}</td>
<td>${args[10]}</td>
<td>${args[11]}</td>
</tr>

</table>

</div>

__EOF__
