#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=32
#SBATCH --time=1:00:00
#SBATCH --job-name=finceng

K=10000000
R=1000

# 5a)
mpirun ./5a 5a-N500.dat 500 "$K" "$R" 1.0

# 5c)
for lambda in 0.25 0.50 0.90; do
 mpirun ./5c "5c-N500-l$lambda.dat" 500 "$K" "$R" 1.0 "$lambda"
done

# 5d)
# N = 1000
# λ = 0.0
# α = 0.5
for alpha in 0.5 1.0 1.5 2.0; do
  for N in 500 1000; do
    for lambda in 0.00 0.25; do
      mpirun ./5d "5d-N$N-l$lambda-a$alpha.dat" "$N" "$K" "$R" 1.0 "$lambda" "$alpha" 10 2>&1
    done
  done
done

# 5e)
N=1000
for alpha in 0.5:10 1.0:25 1.5:50 2.0:50; do
  OLD_IFS="$IFS"
  IFS=":"
  set $alpha
  IFS="$OLD_IFS"

  alpha=$1
  S=$2

  for lambda in 0.00 0.25; do
    for gamma in 0.0 1.0 2.0 3.0 4.0; do
      mpirun ./5e "5e-N$N-l$lambda-a$alpha-g$gamma.dat" "$N" "$K" "$R" 1.0 "$lambda" "$alpha" "$gamma" "$S" 2>&1
    done
  done
done
