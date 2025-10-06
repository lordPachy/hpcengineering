#!/bin/bash

filename="$1"
suffix=".cpp"
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis -I${mkEigenInc} ${filename} -o ${filename%"$suffix"}