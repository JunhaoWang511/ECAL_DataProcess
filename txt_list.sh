#!/bin/bash
# generate matriceI-J.txt in each matriceI-J/ folder
for((i=2;i<5;i++))
do
    for((j=2;j<5;j++))
    do
        dir_name="matrice$i-$j"
        txt_name="${dir_name}/${dir_name}.txt"
        > $txt_name
        for datfile in `ls $dir_name`
        do
            if [ ${datfile##*.} = "dat" ];then
                datpath=`realpath "$dir_name/$datfile"`
                echo $datpath >> $txt_name
            fi
        done
    done
done