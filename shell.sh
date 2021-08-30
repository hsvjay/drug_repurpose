#!/bin/bash
sudo apt-get -f install datamash
for i in `ls *.txt`;do grep -v logFC $i|awk '{print $7"\t"$6}'|sed 's/"//g'|awk  -F'\t' 'x$1'|sort|datamash -g 1 mean 2 >> proc_"$i";done

