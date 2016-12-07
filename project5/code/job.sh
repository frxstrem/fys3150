#!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=32
#SBATCH --time=3:00:00
#SBATCH --job-name=finceng

MPIRUN=( mpirun )

K=10000000
R=1000

# 5a)
"${MPIRUN[@]}" ./5a 5a-N500.dat 500 "$K" "$R" 1.0

# 5c)
for lambda in 0.25 0.50 0.90; do
 "${MPIRUN[@]}" ./5c "5c-N500-l$lambda.dat" 500 "$K" "$R" 1.0 "$lambda"
done

# 5d)
for alpha in 0.5:2 1.0:25 1.5:500 2.0:2500; do
  OLD_IFS="$IFS"
  IFS=":"
  set $alpha
  IFS="$OLD_IFS"

  alpha=$1
  S=$2

  for N in 500 1000; do
    for lambda in 0.00 0.25; do
      "${MPIRUN[@]}" ./5d "5d-N$N-l$lambda-a$alpha.dat" "$N" "$K" "$R" 1.0 "$lambda" "$alpha" "$S" 2>&1
    done
  done
done

# 5e)
N=1000
for alpha in 1.0:25 2.0:2500; do
  OLD_IFS="$IFS"
  IFS=":"
  set $alpha
  IFS="$OLD_IFS"

  alpha=$1
  S=$2

  for lambda in 0.00 0.25; do
    for gamma in 0.0 1.0 2.0 3.0 4.0; do
      "${MPIRUN[@]}" ./5e "5e-N$N-l$lambda-a$alpha-g$gamma.dat" "$N" "$K" "$R" 1.0 "$lambda" "$alpha" "$gamma" "$S" 2>&1
    done
  done
done
