#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=32
#SBATCH --time=1:00:00
#SBATCH --job-name=finceng

# 5a)
mpirun ./5a 5a.dat 500 10000000 1000 1.0

# 5c)
for lambda in 0.25 0.50 0.90; do
  mpirun ./5c "5c-N500-l$lambda.dat" 500 10000000 1000 1.0 "$lambda"
done

# 5d)
for N in 500 1000; do
  for lambda in 0.00 0.25; do
    for alpha in 0.5 1.0 1.5 2.0; do
      mpirun ./5d "5d-N$N-l$lambda-a$alpha.dat" "$N" 10000000 1000 1.0 "$lambda" "$alpha" 100
    done
  done
done

# 5d)
for lambda in 0.00 0.25; do
  for alpha in 0.5 1.0 1.5 2.0; do
    for gamma in 0.0 1.0 2.0 3.0 4.0; do
      mpirun ./5d "5e-N1000-l$lambda-a$alpha-g$gamma.dat" "$N" 10000000 1000 1.0 "$lambda" "$alpha" "$gamma" 100
    done
  done
done
