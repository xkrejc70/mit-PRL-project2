#!/bin/bash

if [ ${#1} -eq 1 ]; then
	np=1
else
	np=$(((2*${#1})-2))
fi

#preklad cpp zdrojaku
mpic++ --prefix /usr/local/share/OpenMPI -o pro pro.cpp

#pocet parametru
if [ $# -ne 1 ]; then
	echo "Program requires exactly one argument!"
	exit 0
fi

#spusteni
mpirun --prefix /usr/local/share/OpenMPI -oversubscribe -np $np pro $1

#uklid
rm -f pro
