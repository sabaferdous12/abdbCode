#!/bin/bash
#*************************************************************************
#
#   Program:    AbDBpipeline
#   File:       AbDBpipeline.sh
#
#   Version:    V1.1
#   Date:       22.04.14
#   Function:   This is an automatic pipeline to process the Data for AbDb and
#               It generates data just for all numbering schemes
#   Usage:      ./AbDBpipeline.sh <list of antibody PDB codes>
#               Note: abdb.cfg file should be in current working dir 
#
#
#   Copyright:  (c) UCL, Saba Ferdous, 2014
#   Author:     Miss Saba Ferdous
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      saba@bioinf.org.uk
#
#*************************************************************************
. ./abdb.cfg
lib=$abdb_HOME/lib
bin=$abdb_HOME/bin

function process
{
     # parameter
    free=$1
    proAntigen=$2
    npAntigen=$3
    Combined=$4
    scheme=$5
    label=$6

    data="Data"
    mkdir -p $data
        
    Redundant="Redundant_files"
    mkdir -p $Combined
    mkdir -p $Redundant

    # copy all the complexes in Combined directory
    cp ./$free/* ./$Combined
    cp ./$proAntigen/* ./$Combined
    cp ./$npAntigen/* ./$Combined

    # Preparing redundant clusters for all 4 datasets
    for f in `ls -d $label"_"*`;
    do
        cd $f
        echo "Calculating Redundant Clusters for $free";
        perl $bin/getRedundantClustersAntibody.pl
        mv *.txt ../$Redundant

        # Preparaing non-redundant data for each dataset
        mkdir -p "NR_"$f
        echo "Preparing Non-Redundant Data";
        perl $bin/prepareNRAbData.pl ../$Redundant/"Redundant_"$f".txt" "NR_"$f
        if [ $f = $Combined ] ; then
            `ls "NR_"$f/*.pdb | grep "_" | cut -f1 -d. | awk -F "/" '{print \$2}' >../"NR_"$label"_Combined.txt"`
        else    
        `ls "NR_"$f/*.pdb | grep "_" | cut -f1 -d. | awk -F "/" '{print \$2}' >>../"NR_"$label"_Merged.txt"`
        fi
        mv "NR_"$f ../
        cd ..
     done

    # This shell command finds diffrence between Combined.txt and Merged.txt 
    `comm -13 <(sort NR_"$label"_Combined.txt) <(sort NR_"$label"_Merged.txt) >$label"_difference.list"`

    mv $free ./$data
    mv $proAntigen ./$data
    mv $npAntigen ./$data
    mv $Combined ./$data
    mv "NR_"* ./$data
    cd $data

    combineData $scheme
    cd ..
    } # Function ends
##############################
# This function puts together all the data (Complete Antibody, Light and
# Heavy chains) for each numbering scheme
function combineData
{
    scheme=$1;
    
    mkdir -p "ALL_"$scheme
    cp ./"LH_Combined_"$scheme/* ./ALL_$scheme
    cp ./"L_Combined_"$scheme/* ./ALL_$scheme
    cp ./"H_Combined_"$scheme/* ./ALL_$scheme
    cd "ALL_"$scheme
    cd ..
}

function compress
{
    array=(*/)
    for directory in  "${array[@]}"
    do
        dirname=$(basename "$directory")
        tar -jcvf $dirname.tar.bz2 $dirname
    done
}


function runProg
{
    scheme=$1
    schemeFlag=$2
    # Running processAntibodyPDBs program for 3 numbering schemes
    
    free="LH_Free_"$scheme
    proAntigen="LH_Protein_"$scheme
    npAntigen="LH_NonProtein_"$scheme
    Combined="LH_Combined_"$scheme
    label="LH"

    echo "processAntibodyPDBs program is running for $scheme numbering";
    perl $bin/processAntibodyPDBs.pl $schemeFlag $3 # $3 is input file
    
    process $free $proAntigen $npAntigen $Combined $scheme $label
    perl $bin/FreeComplexedAntibody.pl -$label ./LH_difference.list ./Redundant_files/"Redundant_"$Combined".txt"
exit;    
    free="L_Free_"$scheme
    proAntigen="L_Protein_"$scheme
    npAntigen="L_NonProtein_"$scheme
    Combined="L_Combined_"$scheme
    label="L"
    process $free $proAntigen $npAntigen $Combined $scheme $label
    perl $bin/FreeComplexedAntibody.pl -$label ./L_difference.list ./Redundant_files/"Redundant_"$Combined".txt"

    free="H_Free_"$scheme
    proAntigen="H_Protein_"$scheme
    npAntigen="H_NonProtein_"$scheme
    Combined="H_Combined_"$scheme
    label="H"
    process $free $proAntigen $npAntigen $Combined $scheme $label
    perl $bin/FreeComplexedAntibody.pl -$label ./H_difference.list ./Redundant_files/"Redundant_"$Combined".txt"
    # Combine the LH,L and H combined redundant clusters in to Redundant_ALL file
cat ./Redundant_files/"Redundant_LH_Combined_"$scheme".txt" ./Redundant_files/"Redundant_L_Combined_"$scheme".txt" ./Redundant_files/"Redundant_H_Combined_"$scheme".txt" >./Redundant_files/"Redundant_ALL_"$scheme".txt"
    }






############################################
echo "Main program running";
############################################

scheme="Martin"
schemeFlag="-a"
runProg $scheme $schemeFlag $1
mkdir -p $scheme"_logs"
mv *.list *.dat ./$scheme"_logs"

exit; 
scheme="Kabat"
schemeFlag="-k"
runProg $scheme $schemeFlag $1
mkdir -p  $scheme"_logs"
mv *.list *.dat ./$scheme"_logs"

scheme="Chothia"
schemeFlag="-c"
runProg $scheme $schemeFlag $1
mkdir -p $scheme"_logs"
mv *.list *.dat ./$scheme"_logs"


mv Redundant_files Data
cd Data/Redundant_files
# To sort all the files in directory
shopt -s nullglob
filearray=( * ) # Reading directory files into an array
 for i in "${filearray[@]}"
 do
     sort $i -o $i
 done
cd ..
#mv *.txt ../$log

# Stats for processed data
perl $bin/getprocessedDataStats.pl
mv *.tt ../
compress
cd ..

# Stats for unprocessed data
cd Martin_logs

grep -r "multi-chain" ../Dataprep*Martin/ | awk -F "/" '{print $3}' | sort | uniq >multi-chain.list
grep -r "scFV" ../Dataprep*Martin/ | awk -F "/" '{print $3}' | sort | uniq >scFV.list

kabatError=`awk 'END {print NR}' Kabat_Error.list`;
Fc=`awk 'END {print NR}' FC.list`;
superceded=`awk 'END {print NR}' Superceded.list`;
scFV=`awk 'END {print NR}' scFV.list`;
multichains=`awk 'END {print NR}' multi-chain.list`;

realKabatError=$(($kabatError-$scFV));

bash $bin/statsUnprocessed.sh $Fc $realKabatError $superceded $scFV >../stats_unprocessed.tt

# To merge 2 consective lines into one
awk 'NR%2{printf $0" ";next;}1' header.dat >headerProcessed.dat

# grep -r "scFV" . | awk -F "/" '{print $2}' > scFv.list
# grep -r "multi-chain" . | awk -F "/" '{print $2}' | sort | uniq


#echo "Moving Data to destination"
#mv Data $data_dest
#mv *.tt $data_dest