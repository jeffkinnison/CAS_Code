#!/bin/bash

###
# TODO: set main directory, walker directory, and molecular dynamics simulation program directory
# and number of nodes requested, number of cores per node, number of gpu's requested.
# num_sim_walkers should be equal to num_nodes * num_cpu, since one walker will run with one node.
###

export MAIN_DIRECTORY=/scratch/users/sahn1/Triazine
export WALKER_DIRECTORY=/scratch/users/sahn1/Triazine/CAS
export GROMACS=/home/sahn1/gromacs/4.6.4/bin

num_nodes=1
num_cpu=16
num_gpu=8
num_sim_walkers=16

# get sequence of walker indices from sh_input.txt
cd $MAIN_DIRECTORY
output=$(cat bash_script_input_file.txt)
first_walker=$(echo $output | awk '{$NF=""}1')
last_walker=$(echo $output | cut -d'_' -f 2)

# assume walkers run in separate directories, walker0, walker1, ...
cd $WALKER_DIRECTORY
counter=0
pin=0
for i in `seq $first_walker $last_walker`;
do
	if [ $counter -lt $num_sim_walkers ];
	then
		cd walker$i/
                let cpu_pin=($i%$num_cpu) #pin=$(($i/$num_nodes))
                let gpu_pin=($i%$num_gpu)
                # write hostfile for i-th job to use
                let lstart=($counter-1)*${num_nodes}+2
                let lend=${lstart}+${num_nodes}-1
                sed -n ${lstart},${lend}'p' < ../../nodefilelist.txt > nodefile$counter
                ssh $(cat nodefile$counter) bash -c "`
                cd ${WALKER_DIRECTORY}

                cd walker$i

                export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/

                ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol_ions.top -n ../../index.ndx -c minim.gro -o run.tpr -maxwarn 1

                ${GROMACS}/mdrun -gpu_id $gpu_pin -ntmpi 1 -ntomp 1 -pin on -pinstride 1 -pinoffset $cpu_pin -v -deffnm run
                `" &
                echo "running walker$i"
                cd ..
                let counter=counter+1
	fi
  	if [ $counter == $num_sim_walkers ];  # limit to num_sim_walkers concurrent walkers
  	then
      		let counter=0 # reset counter
      		wait
  	fi
done
wait

# post-process data after simulations are done. this needs to be modified depending on how the output files are named
# and what collective variables are being collected from the simulations.
for i in `seq 0 $last_walker`;
do
    cd walker$i
    echo "post-processing walker$i"
    while [ ! -f "run.gro" ];  # in case a simulation failed to run
    do
        sed -n 1,1'p' < ../../nodefilelist.txt > nodefile0
        ssh $(cat nodefile0) bash -c "`
        cd ${WALKER_DIRECTORY}

        cd walker$i

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/

        ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol_ions.top -n ../../index.ndx -c minim.gro -o run.tpr -maxwarn 1

        ${GROMACS}/mdrun -v -deffnm run
        `" &
    done
    wait
    rm -rf nodefile*
    rm -rf mdout.mdp
    mv run.gro minim.gro
    mv run.tpr minim.tpr
    echo 0 | ${GROMACS}/trjconv -s minim.tpr -f run.xtc -b 100.0 -o run.xtc
    echo 25 28 | ${GROMACS}/g_mindist -s minim.tpr -f run.xtc -n ../../index.ndx
    awk 'END {print $2*10}' mindist.xvg > coordinates.out
    #awk '{printf $0" ";next;}' coordinates.out > new_coordinates.out
    #mv new_coordinates.out coordinates.out
    rm -rf mindist*
    if [ -f "traj.xtc" ];  # concatenating trajectories after the first step
    then
             ${GROMACS}/trjcat -f traj.xtc run.xtc -o out.xtc -cat
    fi
    mv run.xtc traj.xtc
    mv out.xtc traj.xtc
    rm -rf run*
    rm -rf \#*
    cd ..
done
exit
