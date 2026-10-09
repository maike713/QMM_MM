#!/bin/bash 
#$ -cwd -pe nproc 24 -N NAME
n_walker=24

source /usr/local/run/gromacs-dftbplus-machine-learning/bin/GMXRC
declare -x LD_LIBRARY_PATH=/usr/local/run/plumed-2.9.0-openblas/lib
ulimit -c unlimited

export OMP_NUM_THREADS=1 # openMP threads/processes?
export GMX_QMMM_VARIANT=1 # defines type of electrostatics 0= ewald , 1= ewald mesh 2-4= different cut-offs
#export GMX_DFTB_ESP=1000
#export GMX_DFTB_CHARGES=1000
#export GMX_DFTB_QM_COORD=1000
#export GMX_DFTB_MM_COORD=100
#export GMX_DFTB_MM_COORD_FULL=100

current=$PWD
# add gll topol to cwd
workdir=/scratch/$USER/GMX_$JOB_ID
start_date=$(date)
path_top_='/data/jboeser/CONFIG/top'
path_joblog='/data/jboeser/CONFIG/joblog.txt'


mkdir -vp $workdir
#cp HILLS* $workdir
cp -r -v WALKER_* $workdir
cd $workdir

for i in $(seq 0 $((${n_walker}-1)))
do
  cd WALKER_$i
     #gmx_d mdrun -ntomp 1 -deffnm metad_$i -cpi metad_$i.cpt -plumed plumed_$i.dat -noappend &>/dev/null & 
     gmx_d mdrun -ntomp 1 -deffnm metad_$i -plumed wtmetad_$i.dat &>/dev/null &
     sleep 2
     rm dftb_in.out
  cd ..
done

wait

cd $current

rsync -avruz $workdir/ .
# to log file
echo $JOB_ID >> $path_joblog
echo ${start_date} >> $path_joblog
date >> $path_joblog
echo $current >> $path_joblog
echo $workdir >> $path_joblog
echo $HOSTNAME >> $path_joblog
echo 'DELETE IT!!! [ ], copy_analyse_data  [ ]' >> $path_joblog
echo -e "\n" >> $path_joblog
# to file
echo $JOB_ID >> job-ID.info
echo ${start_date} >> job-ID.info
date >> job-ID.info
echo $current >> job-ID.info
echo $workdir >> job-ID.info
echo $HOSTNAME >> job-ID.info
