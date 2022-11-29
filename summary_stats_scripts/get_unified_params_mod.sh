#!/bin/bash

set -e

find_good_resudual_likelihood=/export-data/users/vsharov/ExpDyn/ExpDyn_PhyDyn/phy_script_dir/find_good_resudual_likelihood.py

CUR_DIR=`pwd`
OUT_DIR=`pwd`/unified_param
mkdir -p ${OUT_DIR}
rm -rf ${OUT_DIR}/*

for i in $(seq 1 11);
do
    j=$((i+1))
    
    WORK_DIR=${CUR_DIR}/p${i}_p${j}
    cd ${WORK_DIR}
    mkdir -p ./unified_params_exp
    rm -rf ./unified_params_exp/*
    mkdir -p ./unified_params_phy
    rm -rf ./unified_params_phy/*
    
    #cd ${WORK_DIR}/work
    cd ${WORK_DIR}/work/exp_dyn
    n=`ls -d exp_dyn*/ | wc -l`
    for d in $(seq 1 $n);
    do
    	cd ${WORK_DIR}/work/exp_dyn/exp_dyn_${d}
    	echo ${WORK_DIR}/work/exp_dyn/exp_dyn_${d}
	if [ -f ./unified_params_residuals*.txt ]; then
	    cp ./unified_params_residuals*.txt ${WORK_DIR}/unified_params_exp
	fi
    	#mv ${WORK_DIR}/unified_params_exp/unified_params_residuals.txt ${WORK_DIR}/unified_params_exp/unified_params_residuals_${n}.txt
    	#n=$((n+1))
    	cd ${WORK_DIR}/work/exp_dyn
    done
    
    cd ${WORK_DIR}/work/phy_dyn
    n=`ls -d beast*/ | wc -l`
    for d in $(seq 1 $n);
    do	
    	cd ${WORK_DIR}/work/phy_dyn/beast_${d}
	echo ${WORK_DIR}/work/phy_dyn/beast_${d}
	if [ -f ./unified_params_residuals*.txt ]; then
    	    awk '{print $6}' Posterior_init_parameters_variables_${d}.txt > treeLikelihood_col.txt && paste unified_params_residuals*.txt treeLikelihood_col.txt > ${WORK_DIR}/unified_params_phy/unified_params_residuals_${d}.txt
    	fi
    	#cp ./unified_params_residuals.txt ${WORK_DIR}/unified_params_phy
    	#mv ${WORK_DIR}/unified_params_phy/unified_params_residuals.txt ${WORK_DIR}/unified_params_phy/unified_params_residuals_${n}.txt
    	#n=$((n+1))
    	cd ${WORK_DIR}/work/phy_dyn
    done
    cp ${WORK_DIR}/target_traj_* ${WORK_DIR}/unified_params_phy
    cp ${WORK_DIR}/target_traj_* ${WORK_DIR}/unified_params_exp
    
    ##Find better residual and treeLikelihood
    cd ${WORK_DIR}/unified_params_phy/
    mkdir -p ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood
    rm -rf ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/*
    module load python3
    python3 $find_good_resudual_likelihood -p `ls unified_params_residuals* | paste -d, -s -` -n 10
    
    for mcmc_step_gl_step in `cat good_residual_likelihood_result.txt | sed '1d' | awk '{print $1"_"$23}'` 
    do
	mcmc_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $1}'`
	gl_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $2}'`
	cp ${WORK_DIR}/work/phy_dyn/beast_${gl_step}/solve/${mcmc_step}_phy_trajectory.txt ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt
	cp ${WORK_DIR}/work/phy_dyn/beast_${gl_step}/tree_pictures/STATE_${mcmc_step}.png ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_tree.png
    done
    
    cd ${WORK_DIR}
    
    mkdir -p ${OUT_DIR}/p${i}_p${j}
    rm -rf  ${OUT_DIR}/p${i}_p${j}/*
    cp -r  ./unified_params_exp  ${OUT_DIR}/p${i}_p${j}/
    
    cd ${WORK_DIR}
    cp -r  ./unified_params_phy  ${OUT_DIR}/p${i}_p${j}/
    cd ${WORK_DIR}
    
    #exit 1;
    
done
    