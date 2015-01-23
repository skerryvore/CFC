#!/usr/bin/bash

input=$1
output=$2

# Verify the type of input and number of values
# Display an error message if the username (input) is not correct
# Exit the shell script with a status of 1 using exit 1 command.
[ $# -ne 2 ] && { echo "Usage: $0 input output"; exit 1; }

gdal_translate -of GMT $input tempo.grd
grd2xyz tempo.grd -ZBLa >> $output
rm tempo.grd

