#!/bin/bash

for i in {0..23}
do

mkdir WALKER_$i
mv metad_$i.tpr WALKER_$i/metad_$i.tpr
#if append traj
#mv metad_$i.cpt WALKER_$i/metad_$i.cpi

# cp reference_renum.pdb WALKER_$i/
cp dftb_in.hsd WALKER_$i/
cp wtmetad.dat WALKER_$i/wtmetad_$i.dat

sed -i -e 's/üü/'${i}'/g' WALKER_$i/wtmetad_$i.dat
done
