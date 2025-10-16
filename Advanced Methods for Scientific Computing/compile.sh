#!/bin/bash

filename="$1"
suffix=".cpp"
g++ -Wall ${filename} -o ${filename%"$suffix"}