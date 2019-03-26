#!/bin/bash

for param in *.para;do
if [ -f $param ]; then
    new_dir="`basename $param .para`"
    mkdir $new_dir
    mv $param new_dir
    cd $new_dir
    echo "$new_dir"
    mcfost $param > $new_dir/$new_dir.txt
    mcfost $param -img 0.6 -only_scatt >> $new_dir/$new_dir.txt
    cd ..
fi
done
