#! /bin/bash
g++ -std=c++11 -o $1 $1.cpp
time ./$1
