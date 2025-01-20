
source /share/scripts/gromacs-XXX.env


# Based on the geometric shape of the PIEZO2 dome (PDB: 6kg7, from OPM database), we created a hemisphere structure comprising 2129 dummy particles (dome shape), and a torus structure with a minor radius of 2 nm and a major radius of 16 nm, made up of 1800 particles (doughnut shape). These two structures were positioned 16 nm apart along the Z-axis. 
python pyscript/genSSMB.py input_pdb/6kg7.pdb 600 600         # box X: 600, box Y: 600
# It will return the minimum box size for Z based on the hemisphere and doughnut generated from 6kg7.
# system.top : topology file.
# dome.pdb   : structure file.
# ssmb.itp   : hemisphere force field parameter file.


# Generate force field parameter file for donuts
grep "DON" dome.pdb > don.pdb
python pyscript/genTop.py don.pdb 0.5 3 > don.itp   # The elastic bond ranges from 0.5 nm to 3 nm. 

# Generate MARTINI bilayer membrane
python pyscript/insane.py  -o cg_lipid.gro -p cg_lipid.top -pbc square -box 60,60,44.3 -l POPC -u POPC -sol W  -dm 0  
gmx_mpi solvate -cp dome.pdb -cs cg_lipid.gro -o system.gro -p system.top -radius 0.47


## Minimazation 
gmx_mpi grompp -f mdp/min.mdp -c system.gro -p system.top -o min.tpr -maxwarn 2 -r system.gro
gmx_mpi mdrun -deffnm min -ntomp 16 

## Equilibration
gmx_mpi grompp -f mdp/eq.mdp  -c min.gro  -p system.top -o eq.tpr -maxwarn 2 -r system.gro
gmx_mpi mdrun -deffnm eq -ntomp 16 


## Kept the position of the torus fixed, designating the Z-axis as the reaction coordinate, and applied a force constant of 100,000 kJ·mol−1·nm−2 to pull the hemisphere downward by 15 nm. 
sed -i 's/1000  1000  1000/1000  1000  0/g' ssmb.itp
sed -i 's/1000  1000  1000/1000  1000  1000/g' don.itp
gmx_mpi grompp -f mdp/pull_1.mdp   -c eq.gro  -p system.top -o pull_1.tpr -maxwarn 2 -r system.gro
gmx_mpi mdrun -deffnm pull_1 -cpt 0.5 -ntomp 16  


## With the hemisphere fixed, the same force was applied to push the torus upward until the membrane structure resembles the PIEZO2 shape, and removed the structured particles
sed -i 's/1000  1000  0/1000  1000  1000/g' ssmb.itp
sed -i 's/1000  1000  1000/1000  1000  0/g' don.itp
gmx_mpi grompp -f mdp/pull_2.mdp   -c pull_1.gro  -p system.top -o pull_2.tpr -maxwarn 2 -r pull_1.gro
gmx_mpi mdrun -deffnm pull_2 -cpt 0.5 -ntomp 16  


