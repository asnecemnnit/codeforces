#! /bin/bash
mkdir ../CurrentContest
python3 generate_files.py $1
subl ../../codeforces
subl ../src/utility.cpp
subl -p ../CurrentContest/*.cpp
sh ./clear_io.sh