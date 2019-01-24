#!/bin/bash
# extract file ending
EXT=$1
# make directory to store new files
mkdir $EXT
# rename new files with ending and move to new directory
mv data_th $EXT
mv data_0.5 $EXT
mv _dust_prop_th.tmp $EXT
mv ref3.0_$EXT.para $EXT
mv ref3.0_$EXT.para~ $EXT
# unzip SED + image files (don't need if using python)
# gunzip $EXT/data_th/sed_rt.fits.gz
# gunzip $EXT/data_0.5/RT.fits.gz
