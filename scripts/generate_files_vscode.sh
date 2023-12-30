#! /bin/bash
mkdir ../"$1"
echo "$1" "created"
python3 generate_files.py "$1" $2 $3
echo "Generated" $2 "files"
code ../../codeforces
code ../src/utility.cpp
code -p ../"$1"/*.cpp
sh ./clear_binaries.sh
sh ./clear_io.sh