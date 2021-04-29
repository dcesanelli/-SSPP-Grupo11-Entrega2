#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output.txt
#SBATCH -e errores.txt
./ejercicio3 $1 $2 $3 2>&1 > "output-ej3-$1-$2-$3.txt"

