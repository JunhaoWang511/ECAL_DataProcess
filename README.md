# ECAL data process program

## What can this program do:

This program is to process beam test data of ECAL in CERN at 2024a, it has two parts: *calibration procedure* and *hardon data analysis*.

## Calibration procedure

1. Make a directory `mkdir data` and change to it `cd data/`, and run `../make_directory.sh` to make 'matrice folders'.
2. Copy .dat file from '/ustcfs3/stcf/BT2024' to corresponding 'matrice folder'.
3. Run `../txt_list.sh` in data directory to generate .txt file under each 'matrice folder'.
4. Run `../dat_to_root.sh` to convert binary file to root file, program will run in background and wait until finished. (you need to specify the path of executable *ECALdig2root* in shell script, and it is generated from [bt2024](https://git.ustc.edu.cn/stcf1/bt2024))
5. Compile 'process_all.cpp' in root environment and execute `process_all()` function, finally 'merged.root' is genenrated.

## Hadron data analysis