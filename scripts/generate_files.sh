#! /bin/bash
mkdir ../CurrentContest
python3 generate_files.py $1
subl ../../codeforces
subl -p ../CurrentContest/*.cpp