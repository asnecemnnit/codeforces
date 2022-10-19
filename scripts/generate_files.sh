#! /bin/bash
mkdir ../"$1"
echo "$1" "created"
python3 generate_files.py "$1" $2 $3
echo "Generated" $2 "files"
subl ../../codeforces
subl ../src/utility.cpp
subl -p ../"$1"/*.cpp
sh ./clear_io.sh
echo "Cleared IO files"