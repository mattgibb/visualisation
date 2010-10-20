#!/usr/bin/python

#
# Translates the raw skeletonized Purkinje representation in voxel coordinates to 
# spatial coordinates compatible with Martin Bishop's rabbit heart meshes. All files
# are in Meshalyzer format.
#
# Note that the x & y axes are interchanged between the two representations
#

import sys
import os

# voxel resolution
#For full resolution
#xres = 0.0052; #cm
#yres = 0.0052; #cm
#zres = 0.0024; #cm

#Offsets for original prototype Pkj, note that x & y must be interchanged to match M.B. mesh
xres = 26.4; #cm
yres = 26.4; #cm
zres = 24.4; #cm

xoffset=90;
yoffset=110;
zoffset=0;

pkj_xoffset = 0;
pkj_yoffset = 0;
pkj_zoffset = 0;

# Offsets for full resolution lv purkinje, note that no interchange of x & y is required!
#xres = 1;
#yres = 1;
#zres = 1;

#xoffset = 90*26.4; 
#yoffset = 110*26.4;
#zoffset = 0#-50*24.4;

# Offsets for decimated 001 lv purkinje
#pkj_xoffset = 54*4*26.4; 
#pkj_yoffset = 84*4*26.4;
#pkj_zoffset = 270*12.2;#200*12.2;



if __name__ == "__main__":
    
    input_file_root = sys.argv[1];
    output_file_root = sys.argv[2];
    
    #Copy over the points file 
    original_file = open(input_file_root + '.pts');
    new_file = open(output_file_root + '.pts', 'w');

    line = original_file.readline();
    new_file.write(line);

    for line in original_file:
        entries = line.split();
        numbers = [];
        if entries:
            # Note that x & y axes are interchanged!
            numbers.append((float(entries[1]) - xoffset + pkj_xoffset)*yres);
            numbers.append((float(entries[0]) - yoffset + pkj_yoffset)*xres);
            numbers.append((float(entries[2]) - zoffset + pkj_zoffset)*zres);

            new_file.write(str(numbers[0]) + " " + str(numbers[1]) + " " + str(numbers[2]) + "\n")

    original_file.close();
    new_file.close();
    
    #Copy over the cnnx file, no changes are required
    os.system('cp ' + input_file_root + '.cnnx ' + output_file_root + '.cnnx');
