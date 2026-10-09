#!/bin/bash 
#$ -cwd -pe nproc 1 -N NAME

source /usr/local/run/gromacs-dftbplus-machine-learning/bin/GMXRC

export OMP_NUM_THREADS=1
export GMX_QMMM_VARIANT=1 # defines type of electrostatics 0= ewald , 1= ewald mesh 2-4= different cut-offs
#export GMX_DFTB_ESP=100
#export GMX_DFTB_CHARGES=100
#export GMX_DFTB_QM_COORD=100


current=$PWD
# add gll topol to cwd
workdir=/scratch/$USER/GMX_$JOB_ID

mkdir -vp $workdir
cp dftb_in.hsd md.tpr $workdir
cd $workdir

gmx_d mdrun -ntomp 1 -deffnm md -plumed plumed.dat 
#gmx_d mdrun -ntomp 1 -deffnm md 
#gmx_d mdrun -ntomp 1 -deffnm metad_$i -cpi metad_$i.cpt -plumed plumed_$i.dat -noappend &>/dev/null & 

wait

cd $current

rsync -avruz $workdir/ .
