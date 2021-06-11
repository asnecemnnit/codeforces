#! /bin/bash
ABSOLUTE_PATH_OF_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$ABSOLUTE_PATH_OF_SCRIPT"/../IO/
rm input.txt 
touch input.txt
rm output.txt 
touch output.txt
rm log.txt 
touch log.txt