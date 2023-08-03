#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=lmp
#SBATCH --nodes=1
#SBATCH --nodelist=cn2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
##SBATCH --threads-per-core=1

module load intel/2020

#module load mkl/2022.0.2
#module load mpi/2021.5.1

  exe_dir="/data/home/hfhuang/software/Lammps_for_PWMLFF-master/lammps_neigh_mlff_20230508/src"
  input="lammps.in"
  fid="test"
 
  echo "SLURM_NPROCS ${SLURM_NPROCS}" 
  mpirun -np ${SLURM_NPROCS} ${exe_dir}/lmp_mpi -in ${input} > ./data_out/out.${fid}
  #mpirun -np 1 ${exe_dir}/lmp_mpi -in ${input} > ./data_out/out.${fid}


