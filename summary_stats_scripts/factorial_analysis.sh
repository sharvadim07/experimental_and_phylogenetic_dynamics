#!/bin/bash

set -e

find_good_resudual_likelihood=/export-data/users/vsharov/ExpDyn/ExpDyn_PhyDyn/phy_script_dir/find_good_resudual_likelihood.py
fact_analysis_script_dir=/export-data/users/vsharov/FactorAnalysis/fact_analysis_scripts_new

module load python3

CUR_DIR=`pwd`
OUT_DIR=`pwd`/unified_param
mkdir -p ${OUT_DIR}

fact_an_all_pass_dir=${OUT_DIR}/fact_an_all_pass
mkdir -p ${fact_an_all_pass_dir}
rm -rf ${fact_an_all_pass_dir}/*

repl_dict=/export-data/users/vsharov/ExpDyn/ExpDyn_PhyDyn/tests/zikv_all_paired_passages_upd_eq_v3_traj_fix/repl_dict_ODE4R_zikv.txt

for i in $(seq 1 11);
do
    j=$((i+1))
    
    echo "Working with p${i}_p${j}!"
    WORK_DIR=${CUR_DIR}/p${i}_p${j}
    cd ${WORK_DIR}
    mkdir -p ./unified_params_phy
    
    cd ${WORK_DIR}/unified_params_phy/
    mkdir -p ${WORK_DIR}/unified_params_phy/factorial_analysis
    rm -rf ${WORK_DIR}/unified_params_phy/factorial_analysis/*
    #mkdir -p ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood
    #rm -rf ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/*
    #module load python3
    #python3 $find_good_resudual_likelihood -p `ls unified_params_residuals* | paste -d, -s -` -n 10
    
    k=1
    for mcmc_step_gl_step in `cat good_residual_likelihood_result.txt | sed '1d' | awk '{print $1"_"$23}'` 
    do
	mcmc_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $1}'`
	gl_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $2}'`
	#cp ${WORK_DIR}/work/phy_dyn/beast_${gl_step}/solve/${mcmc_step}_phy_trajectory.txt ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt
	#cp ${WORK_DIR}/work/phy_dyn/beast_${gl_step}/tree_pictures/STATE_${mcmc_step}.png ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_tree.png
	date
	echo "Run factorial for ${k} good residual!"
	###Factorial analysis
	mkdir ${WORK_DIR}/unified_params_phy/factorial_analysis/${k}_GL_${gl_step}_STATE_${mcmc_step}
	cd ${WORK_DIR}/unified_params_phy/factorial_analysis/${k}_GL_${gl_step}_STATE_${mcmc_step}
	#Get params
	cat ${WORK_DIR}/unified_params_phy/good_residual_likelihood_result.txt | grep -e "^${mcmc_step}" | grep -e "${gl_step}$" | cut -f 2-20 > params.txt
	#Get init values
	head -n 2 ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt | sed '1d' | cut -f 2- | sed 's/\t$//' > init_val.txt
	#Preparation
	echo "Preparation step..."
	$fact_analysis_script_dir/partial_design_prep.sh \
	${WORK_DIR}/unified_params_phy/factorial_analysis/${k}_GL_${gl_step}_STATE_${mcmc_step} \
	params.txt \
	init_val.txt > /dev/null
	#Consistency
	echo "Consistency step..."
	cp ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt ./
	$fact_analysis_script_dir/consistency.sh \
	${WORK_DIR}/unified_params_phy/factorial_analysis/${k}_GL_${gl_step}_STATE_${mcmc_step} \
	GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt > /dev/null
	date
	echo "End consistency."
	
	mv FactorialDesignMatr_Res.pdf p${i}_p${j}_${k}_GL_${gl_step}_STATE_${mcmc_step}_FactorialDesignMatr_Res.pdf
	mv FactorialDesignMatr_Res.txt p${i}_p${j}_${k}_GL_${gl_step}_STATE_${mcmc_step}_FactorialDesignMatr_Res.txt
	
	#Modify results
	python3 $fact_analysis_script_dir/mod_factorial_an_res.py \
	-f p${i}_p${j}_${k}_GL_${gl_step}_STATE_${mcmc_step}_FactorialDesignMatr_Res.txt \
	-r ${repl_dict}
	
	cp mod_p${i}_p${j}_${k}_GL_${gl_step}_STATE_${mcmc_step}_FactorialDesignMatr_Res.txt ${fact_an_all_pass_dir}/
	echo "p${i}_p${j}_${k}_GL_${gl_step}_STATE_${mcmc_step}_FactorialDesignMatr_Res.txt" >> ${fact_an_all_pass_dir}/${k}_all_pass.list
	
	###
	cd ${WORK_DIR}/unified_params_phy/
	k=$((k+1))
    done
    
    cd ${WORK_DIR}
    cp -rn  ./unified_params_phy  ${OUT_DIR}/p${i}_p${j}/
    cd ${WORK_DIR}
    
    exit 1;
    
done
    