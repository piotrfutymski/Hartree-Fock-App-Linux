#!/bin/bash
mkdir int
./build/HartreeFockApp  $1 $2
python3 script.py
rm -r int