#!/usr/bin/python

#
# Due to the way the skeletonizers works raw skeletonized files contain duplicates of nodes.
# This file removes the duplicates and adjusts the corresponding cnnx file
#
#

import sys
import os
from sets import Set

if __name__ == "__main__":
    
    input_file_root = sys.argv[1];
    output_file_root = sys.argv[2];
    
    original_pts_file = open(input_file_root + '.pts');

    line = original_pts_file.readline();
    entries = line.split()
    size = int(entries[0]);

    #we want node_index -> new_node_index

    #We'll create a coords -> node_index

    #Create a set of coordinates
    coords = []
    unique_coords = Set()

    line_number = 0
    for line in original_pts_file:
        entries = line.split()
        if entries:
            coords.append(line)
            unique_coords.add(line)

    unique_coords = list(unique_coords)

    print "Total Coordinates: "
    print len(coords)
    
    print "Unique Coordinates: "
    print len(unique_coords)
    
    #Output the unique coordinates file
    new_pts_file = open(output_file_root + '.pts', 'w');
    
    new_pts_file.write(str(len(unique_coords)) + "\n")
    
    i = 0
    for line in unique_coords:
        new_pts_file.write(unique_coords[i])
        i = i + 1
        
    original_pts_file.close()
    new_pts_file.close()
        
    #Now fix the cnnx file
    #first create a map coords_index -> unique_coords_index
    coordinates_dict = {}
    for i in range(0,len(coords)):
        coordinates_dict[i] = unique_coords.index(coords[i])
        
    #Read the existing cnnx file and output the new one
    original_cnnx_file = open(input_file_root + '.cnnx');
    new_cnnx_file = open(output_file_root + '.cnnx','w');
    
    
    line = original_cnnx_file.readline();
    new_cnnx_file.write(line);
    
    for line in original_cnnx_file:
        entries = line.split()
        numbers = [0,0]
        if entries:
            numbers[0] = coordinates_dict[int(entries[0])]
            numbers[1] = coordinates_dict[int(entries[1])]
            
            #It is possible that there are connections between identical nodes
            #we avoid writing those out
            if(numbers[0] != numbers[1]):
                new_cnnx_file.write(str(numbers[0]) + " " + str(numbers[1]) + "\n")
            
    original_cnnx_file.close()
    new_cnnx_file.close()
