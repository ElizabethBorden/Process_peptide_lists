#!/bin/bash

for directory in HLAA2.1;
    do cd /scratch/eknodel/carreno_data/$directory;
    cat $directory"_netCTL.out" >> ../"All_netCTL.out";
    cat $directory"_netMHC.out" >> ../"All_netMHC.out";
    cat $directory"_netMHCstab.out" >> ../"All_netMHCstab.out";
    cat $directory".WTmatch_netMHC.out" >> ../"All.WTmatch_netMHC.out";
    cat $directory".epitopes.annotated.tsv" >> ../"All.epitopes.annotated.tsv";
    cd ~/Carreno_data;
done
