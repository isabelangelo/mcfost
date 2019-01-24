"""
The purpose of this code is to reduce SED files taken from CDS Portal
by removing unnecessary data points

*code should be run in directory with all .sed files in text format, all files with the 
.sed ending will be reduced and saved as fname_reduced.sed in same directory

Written: Isabel Angelo (2018)
"""

import sys
import numpy as np
import glob

# retrieve all SED files to be reduced
filenames = glob.glob('*.sed')

# parse lines for data entries
def parseline(line):
    """
    parses a line and converts to a line of floats
    Args:
        line(str): string version of line from .txt file
    Returns:
        l2(list): list with elements corresponding to columns in original text file
        in format [freq, flux, % err, telescope, ref]
    """
    line = line.replace('\t',' ')
    l1 = line.split(' ')
    l2 = list(filter(None, l1))
    return l2

# reduce SED file and save as reduced file
for fname in filenames:
    # read lines of input sed file
    lines = open(fname).readlines()

    # remove non-data lines in file
    start=0
    for i in range(len(lines)):
        if lines[i][0:9]=='Frequency':
            start = i
            break
        
    toplines = lines[start:start+2]
    lines = lines[start+2:]

    # edit out lines with bad data    
    keep_lines = [] 
    for line in lines:
        l = parseline(line)
        # remove error > 100%
        if np.abs(float(l[2])) >= 100:
            pass
        # remove Gaia/Johnson/POSS
        elif 'Gaia' in l[3] or 'Johnson' in l[3] or 'POSS' in l[3]:
            pass
        # remove non-PAN-STARRS visible data
        elif (0.4 <= float(l[0]) <= 0.7) & ('PAN-STARRS' not in l[3]):
            pass
        # remove non-2MASS/UKIDSS IR data 
        elif (1 <= float(l[0]) <=2) & ('2MASS' not in l[3]) & ('UKIDSS' not in l[3]):
            pass
        else:
            keep_lines.append(line)

    # get Karl's data
    added_lines = []
    lines2 = open(fname).readlines()
    for i in range(start):
        line = lines2[i]
        if line.count('.')==3:
            added_lines.append(line)
    added_lines = added_lines[1:]

    # add Karl's data to line array
    final_lines = np.concatenate((toplines, keep_lines, added_lines))

    # write to text file
    with open(fname[:-4]+'_reduced.sed', 'a') as f:
        for line in final_lines:
            f.write(line)
    
    f.close()
    print(fname)