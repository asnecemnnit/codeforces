#! /bin/bash
py generate_file.py $1 $2 $3
subl ../../codeforces
subl ../src/utility.cpp
if [ "$#" -eq 1 ]; then
	subl ../ProblemSet/PS_Codeforces/$1.cpp
else
	subl ../ProblemSet/$2/$1.cpp
fi
sh ./clear_io.sh