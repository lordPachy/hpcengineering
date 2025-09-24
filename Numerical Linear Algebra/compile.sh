#!/bin/bash

filename="$1"
suffix=".cpp"
g++ -I ${mkEigenInc} ${filename} -o ${filename%"$suffix"}
