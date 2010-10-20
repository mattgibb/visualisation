#!/usr/bin/python

#
# Extracts the node indices of the PMJs.
#
# PMJs are defined as the end points of the Purkinje tree, so we look for nodes that are the
# end of a connection but not the start of one.
#

import sys
import os

if __name__ == "__main__":
    
    input_file_root = sys.argv[1];
    
    #Copy over the points file 
    cnnx_file = open(input_file_root + '.cnnx');
    

    line = cnnx_file.readline();
    
    cnnx_ends = []
    
    #Obtain the start of each connection and the ends (the potential pmjs)
    for line in cnnx_file:
        entries = line.split()
        if entries:
            cnnx_ends.append(entries[0])
            cnnx_ends.append(entries[1])

    #Count the frequency of each entry
    counted_ends = [(end, cnnx_ends.count(end)) for end in set(cnnx_ends)];
            
    #Write out the pmjs
    pmj_file = open(input_file_root + '.pmj_half', 'w');
    
    #pmj_file.write(str(len(counted_ends))  +" 1 1 1\n")
    
    for end in counted_ends:
        if(end[1] == 1):
            pmj_file.write(end[0] + "\n")

    #Close up files
    cnnx_file.close()
    pmj_file.close()
