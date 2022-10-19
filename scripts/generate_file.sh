#! /bin/bash
python3 generate_file.py $1 $2 $3
subl ../../codeforces
subl ../src/utility.cpp
subl ../ProblemSet/$2/$1.cpp
sh ./clear_io.sh