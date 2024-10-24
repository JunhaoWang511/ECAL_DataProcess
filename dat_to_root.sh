#!/bin/bash
# execute $build/ECALdig2root ./matriceI-J/matriceI-J.txt matriceI-J.root
execute_program='/ustcfs/HICUser/jhwang/bt2024/ECAL/build/ECALdig2root'
appname='_mu_5GeV.root'
for((i=2;i<5;i++))
do
    for((j=2;j<5;j++))
    do
        dir_name="matrice$i-$j"
        txt_name="${dir_name}/${dir_name}.txt"
        rootname="${dir_name}${appname}"
        `$execute_program ${txt_name} ${rootname} >/dev/null 2>&1 &`
    done
done