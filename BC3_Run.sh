#!/bin/bash
#$ -cwd
# request resources:
#PBS -l nodes=3:ppn=16
#PBS -l walltime=0:10:00
# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR
# run your program, timing it for good measure:
#time ./my-multi-threaded-program

#Load Abaqus modules
#check modules using modules avail
module add apps/abaqus-6.14
#module add intel/compiler/64/12.1/2011_sp1.12.361
module add intel/compiler/64/12.1/2011_sp1.8.273

#Single core user sub
#Load the intel compiler stored in ABAQCOMPVER
module load $ABAQCOMPVER
abaqus job=RVE_0921_Edit input=/newhome/mm20962/RVE_0921/RVE_0921_Edit.inp user=/newhome/mm20962/RVE_0921/ABQ_Lengthscale.for cpus=16 scratch=$TMPDIR mp_mode=mpi memory=95% interactive

