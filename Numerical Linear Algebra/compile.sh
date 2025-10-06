#!/bin/bash

filename="$1"
suffix=".cpp"
g++ -I${mkLisInc} -L${mkLisLib} -llis -I${mkEigenInc} ${filename} -o ${filename%"$suffix"}