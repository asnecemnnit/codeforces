#! /bin/bash
python3 generate_file.py $1 $2 $3
code ../../codeforces
code ../src/utility.cpp
if [ "$#" -eq 1 ]; then
	code ../ProblemSet/PS_Codeforces/$1.cpp
else
	code ../ProblemSet/$2/$1.cpp
fi
sh ./clear_binaries.sh
sh ./clear_io.sh