#!/bin/bash
mkdir int
./build/HartreeFockApp  $1 $2
python3 -W ignore script.py
rm -r int