#!/bin/bash

for i in {0..23};do
	gmx_d grompp -f r02*.mdp -n index.ndx -c md.gro -o metad_${i}.tpr -p *.top -maxwarn 2
done
