#! /bin/bash
python3 generate_file.py $1
subl ../../codeforces
subl ../ProblemSet/$1.cpp
