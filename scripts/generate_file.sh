#! /bin/bash
python3 generate_file.py $1
subl ../../codeforces
subl ../src/utility.cpp
subl ../ProblemSet/$1.cpp