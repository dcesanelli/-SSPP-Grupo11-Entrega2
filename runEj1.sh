#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output.txt
#SBATCH -e errores.txt
./ejercicio1 $1 $2 2>&1 > "output-ej1-$1-$2.txt"

