#!/bin/bash
# make directories from matrice2-2/ to matrice4-4/
for((i=2;i<5;i++))
do
    for((j=2;j<5;j++))
    do
        dir_name="matrice$i-$j"
        mkdir -p $dir_name
    done
done