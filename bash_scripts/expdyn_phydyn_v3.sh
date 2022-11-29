#!/bin/bash
#
#SETTINGS
#
INPUT_DIR=$1 # Absolute path input dir
OUT_DIR=$2 # Absolute path out dir
EQUATION_FILE=$3 # EQUATIONS and init min and max values for variables and parameters file absolute path
MSA_FASTA_OUT_FILE=$4 # MSA result in fasta format file absolute path
FASTA_INFO_FILE=$5 # File with information about all fasta sequences (sequence_name,time(as year.part_of_year),variable_name_in_ODE)
MCMC_ITERATIONS=$6 # Number of MCMC iterations
LOG_STEP=$7 # Size of MCMC logging step
MUTATION_RATE=$8 # Mutation rate for sequences in MCMC (default=1.1)
TARGET_FILE=$9 # Target trajectory file for Experimental dynamic
shift;
GLOB_ITERATIONS=$9 # Num of iterations of exp_dyn and phy_dyn (default=1)
shift;
FIRST_DYN=$9 # Set which dynamics will started first: Experimental dynamics (expdyn) or Phylodynamics (phydyn) (default=expdyn)
shift;
NUM_GOOD_RES_TR_L=$9 # The number of parameter sets that have a good residual and tree probability and have maximum variability
shift;
wR=$9 # Coeficient of Residual/min_residual ratio for weighted average sum (wR)


#shift;
#MULTIPLE_RANGE=$9 # Multiplier and divider value for min and max value for paramater range


WORK_DIR=${INPUT_DIR}/work
LOG=${WORK_DIR}/pipe.log
PHY_DYN_DIR=${WORK_DIR}/phy_dyn
MAKE_PHYTREE_DIR=${PHY_DYN_DIR}/make_phytree

EXP_DYN_DIR=${WORK_DIR}/exp_dyn


#
#SOFT
#
module unload tbrc-production-set
module load python3
module load R
#module load clustalw2
#PhyDyn
script_dir=$(dirname $(readlink -f $0))
phy_script_dir=${script_dir}/phy_script_dir
phy_dyn_gen_xml=${phy_script_dir}/phy_dyn_gen_xml.py
beast_tree_constr_gen_xml=${phy_script_dir}/beast_tree_constr_gen_xml.py
draw_phy_trees=${phy_script_dir}/draw_phy_trees.py
find_good_resudual_likelihood=${phy_script_dir}/find_good_resudual_likelihood.py
#beauti_command_line=${phy_script_dir}/BEAUTi_command_line

phydyn_res_to_equation=${phy_script_dir}/phydyn_res_to_equation.py
#ExpDyn
exp_script_dir=${script_dir}/exp_script_dir
param_parser=${exp_script_dir}/deBInfer_ParseParams.R
deBInfer_script=${exp_script_dir}/debinfer.R
expdyn_gen_debinfer_in=${exp_script_dir}/expdyn_gen_debinfer_in.py
expdyn_res_to_equation=${exp_script_dir}/expdyn_res_to_equation.py
solve=${exp_script_dir}/solve.py
find_residual=${exp_script_dir}/find_residual.py
get_two_cols=${exp_script_dir}/awk_by_col_name.awk


#beast=/export/dome/vsharov/SOFT/beast/bin/beast
#clustalw=/export/dome/vsharov/SOFT/clustalw-2.1/bin/clustalw2

set -e

make_phytree(){

set -e
#load SOFT
#module load BEAST/1.10.4
module load BEAST2

MAKE_PHYTREE_DIR=${PHY_DYN_DIR}/make_phytree
rm -rf ${MAKE_PHYTREE_DIR}
mkdir -p ${MAKE_PHYTREE_DIR}
cd ${MAKE_PHYTREE_DIR}

date >> ${LOG}
echo "Start constructing nexus tree by CLustalW for ${MSA_FASTA_OUT_FILE}" >> ${LOG}
cp ${MSA_FASTA_OUT_FILE} ${MAKE_PHYTREE_DIR}
MSA_FASTA_OUT_FILE_NAME=`basename ${MSA_FASTA_OUT_FILE}`

#$beauti_command_line -fa ${MAKE_PHYTREE_DIR}/${MSA_FASTA_OUT_FILE_NAME} -out gen_tree_for_beast -NA -Constant
local MCMC_ITER_FOR_INIT_TREE=10000000
python3 $beast_tree_constr_gen_xml -m ${MAKE_PHYTREE_DIR}/${MSA_FASTA_OUT_FILE_NAME} -i ${FASTA_INFO_FILE} -l ${MCMC_ITER_FOR_INIT_TREE} \
|| { echo "Error while generating XML file for BEAST tree construction for ${MSA_FASTA_OUT_FILE}!" >> ${LOG}; exit 1; }
beast ./gen_tree_coalescent_const.xml
mv ./*.trees ./beast_init.trees
grep -B 10000000 "tree STATE_0" ./beast_init.trees  | sed "$ d" > ./beast_init_last.tree
grep -A 1 "tree STATE_${MCMC_ITER_FOR_INIT_TREE}" ./beast_init.trees  >> ./beast_init_last.tree
sed 's/\[\&rate=[0-9]*\.[0-9]*\]//g' beast_init_last.tree | sed 's/\[\&lnP.*\]/=/g' > ./beast_init_last_corr.tree
#cp ${MAKE_PHYTREE_DIR}/beast_init_last_corr.tree ${OUT_DIR}
mkdir -p ./init_tree_picture
cd ./init_tree_picture
python3 $draw_phy_trees -t ${MAKE_PHYTREE_DIR}/beast_init_last_corr.tree \
|| { echo "Error in draw_phy_trees.py!" >> ${LOG}; exit 1; }
    
cd ${MAKE_PHYTREE_DIR}
#cp -r ./init_tree_picture ${OUT_DIR}
date >> ${LOG}

module unload BEAST2
#module unload BEAST/1.10.4
#clustalw2  -tree -infile=${MAKE_PHYTREE_DIR}/${MSA_FASTA_OUT_FILE_NAME} -outputtree=nexus \
#|| { echo "Error while contructing tree by ClustalW FOR ${MSA_FASTA_OUT_FILE}!" >> ${LOG}; exit 1; }
#cp ${MAKE_PHYTREE_DIR}/*.tre ${OUT_DIR}

}

gen_xml_for_beast(){

set -e

local CUR_ITER=$1
local IN_EQ_FILE=$2

CUR_GEN_XML_DIR=${PHY_DYN_DIR}/gen_xml_${CUR_ITER}
rm -rf ${CUR_GEN_XML_DIR}
mkdir ${CUR_GEN_XML_DIR}
cd ${CUR_GEN_XML_DIR}

date >> ${LOG}
echo "Start phy_dyn_gen_xml. Iter ${CUR_ITER}." >> ${LOG}
#OUT is pattern_model_seq_phy_dyn_v2_updated.xml
python3 $phy_dyn_gen_xml \
 -m ${MSA_FASTA_OUT_FILE} \
 -i ${FASTA_INFO_FILE} \
 -e ${IN_EQ_FILE} \
 -t ${MAKE_PHYTREE_DIR}/beast_init_last_corr.tree \
 -l ${MCMC_ITERATIONS} \
 -s ${LOG_STEP} \
 -r ${MUTATION_RATE} \
|| { echo "Error in phy_dyn_gen_xml!" >> ${LOG}; exit 1; }

mv ${CUR_GEN_XML_DIR}/pattern_model_seq_phy_dyn_v2_updated.xml ${CUR_GEN_XML_DIR}/pattern_model_seq_phy_dyn_v2_updated_${CUR_ITER}.xml


}

beast_start(){
set -e

local CUR_ITER=$1
local IN_EQ_FILE=$2
local ST_T=$3
local END_T=$4
    

local CUR_GEN_XML_DIR=${PHY_DYN_DIR}/gen_xml_${CUR_ITER}
local CUR_BEAST_DIR=${PHY_DYN_DIR}/beast_${CUR_ITER}
rm -rf ${CUR_BEAST_DIR}
mkdir ${CUR_BEAST_DIR}
cd ${CUR_BEAST_DIR}

cp ${IN_EQ_FILE} ${CUR_BEAST_DIR}/
mv ${CUR_BEAST_DIR}/${IN_EQ_FILE##*/} ${CUR_BEAST_DIR}/cur_initial_eq.txt
local IN_EQ_FILE=${CUR_BEAST_DIR}/cur_initial_eq.txt

module load BEAST2
date >> ${LOG}
echo "Start BEAST for cunstructed model. Iter ${CUR_ITER}." >> ${LOG}
beast ${CUR_GEN_XML_DIR}/pattern_model_seq_phy_dyn_v2_updated_${CUR_ITER}.xml > ${CUR_BEAST_DIR}/beast_${CUR_ITER}.log \
|| { echo "Error in BEAST!" >> ${LOG}; echo "BEAST failed"; return 1; }
module unload BEAST2

grep -A 1000000000 "Sample" ${CUR_BEAST_DIR}/MultiDemeBirthDeathPhyDyn_m1.log > ${CUR_BEAST_DIR}/Posterior_init_parameters_variables_${CUR_ITER}.txt
mv ${CUR_BEAST_DIR}/MDBDPhyDyn_m1.traj ${CUR_BEAST_DIR}/Variables_Trejectory_${CUR_ITER}.traj
mv ${CUR_BEAST_DIR}/MDBDPhyDyn_m1.trees ${CUR_BEAST_DIR}/All_Trees_${CUR_ITER}.tree

mkdir -p ${CUR_BEAST_DIR}/tree_pictures
cd ${CUR_BEAST_DIR}/tree_pictures
python3 $draw_phy_trees -t ${CUR_BEAST_DIR}/All_Trees_${CUR_ITER}.tree \
|| { echo "Error in draw_phy_trees.py!" >> ${LOG}; exit 1; }
cd ${CUR_BEAST_DIR}
#cp -r ./tree_pictures ${OUT_DIR}

date >> ${LOG}
#cp ${CUR_BEAST_DIR}/Variables_Trejectory.traj ${OUT_DIR}
#cp ${CUR_BEAST_DIR}/All_Trees.tree ${OUT_DIR}
#cp ${CUR_BEAST_DIR}/Posterior_init_parameters_variables.txt ${OUT_DIR}
#cp ${CUR_BEAST_DIR}/beast.log ${OUT_DIR}
rm -rf ${CUR_BEAST_DIR}/*.state
rm -rf ${CUR_BEAST_DIR}/MultiDemeBirthDeathPhyDyn*

echo "Running Solver for phydyn params. Iteration=${CUR_ITER}" >> ${LOG}
python3 $expdyn_gen_debinfer_in -e ${IN_EQ_FILE} \
|| { echo "Error in expdyn_gen_debinfer_in.py!" >> ${LOG}; exit 1; }
rm ${CUR_BEAST_DIR}/updated_equation.c ${CUR_BEAST_DIR}/parsed_userpriors.txt ${CUR_BEAST_DIR}/initial_values_deBInfer_formatted.txt ${CUR_BEAST_DIR}/repl_dict.txt
g++ ${CUR_BEAST_DIR}/updated_solve.cpp -o ${CUR_BEAST_DIR}/updated_solve

mkdir ${CUR_BEAST_DIR}/solve
cp ${IN_EQ_FILE} ${CUR_BEAST_DIR}/solve/updated_eq.txt
cat ${CUR_BEAST_DIR}/Posterior_init_parameters_variables_${CUR_ITER}.txt | sed '1d' | awk '{print $1}' > ${CUR_BEAST_DIR}/solve/step_list.txt
mv ${CUR_BEAST_DIR}/updated_solve ${CUR_BEAST_DIR}/solve/updated_solve
cd ${CUR_BEAST_DIR}/solve
local UPDATED_EQ_FILE=${CUR_BEAST_DIR}/solve/updated_eq.txt
python3 $phydyn_res_to_equation -e ${UPDATED_EQ_FILE} -s ${CUR_BEAST_DIR}/solve/step_list.txt -p ${CUR_BEAST_DIR}/Posterior_init_parameters_variables_${CUR_ITER}.txt \
|| { echo "Error in phydyn_res_to_equation.py!" >> ${LOG}; exit 1; }
rm -rf ${CUR_BEAST_DIR}/solve/updated_eq.txt
echo -e "step\tresidual" > ${CUR_BEAST_DIR}/solve/residuals.txt
for step in `cat ${CUR_BEAST_DIR}/solve/step_list.txt` ;
do
    python3 $solve -e ${CUR_BEAST_DIR}/solve/${step}_${UPDATED_EQ_FILE##*/} -f ${ST_T} -l ${END_T} -s 0.06 -a 0.001 > ${CUR_BEAST_DIR}/solve/${step}_phy_trajectory.txt \
    || { echo "Error in solve.py!" >> ${LOG}; exit 1; }
    echo -e "${step}\t`python3 $find_residual -t ${TARGET_FILE} -c ${step}_phy_trajectory.txt`" >> ${CUR_BEAST_DIR}/solve/residuals.txt
    local AFTER_PHY_DYN_EQ_FILE=`pwd`/${step}_${UPDATED_EQ_FILE##*/} # Get parameters from last step of phydyn mcmc
done
awk '{print $2}' ${CUR_BEAST_DIR}/solve/residuals.txt | paste ${CUR_BEAST_DIR}/solve/unified_params.txt - -d "\t" > ${CUR_BEAST_DIR}/unified_params_residuals_${CUR_ITER}.txt
rm ${CUR_BEAST_DIR}/solve/residuals.txt ${CUR_BEAST_DIR}/solve/unified_params.txt
rm ${CUR_BEAST_DIR}/solve/updated_solve
cd ${CUR_BEAST_DIR}

echo ${AFTER_PHY_DYN_EQ_FILE}

}

phydyn(){
    set -e

    local CUR_ITER=$1
    local IN_EQ_FILE=$2
    local ST_T=$3
    local END_T=$4
    
    mkdir -p ${PHY_DYN_DIR}
    cd ${PHY_DYN_DIR}
    echo "Start phydyn. Iteration=${CUR_ITER}" >> ${LOG}
    gen_xml_for_beast ${CUR_ITER} ${IN_EQ_FILE}
    local AFTER_PHY_DYN_EQ_FILE=$(beast_start ${CUR_ITER} ${IN_EQ_FILE} ${ST_T} ${END_T})
    echo "End phydyn. Iteration=${CUR_ITER}" >> ${LOG}
    
    echo ${AFTER_PHY_DYN_EQ_FILE}
}

expdyn(){

set -e

local CUR_ITER=$1
local IN_EQ_FILE=$2
local TARGET_FILE=$3
local ST_T=$4
local END_T=$5
local NUM_OF_POINTS_FROM_TARGET=`sed '1d' ${TARGET_FILE} | wc -l`
#Debug
echo ${TARGET_FILE} >> ${LOG}

#Parent dir
mkdir -p ${EXP_DYN_DIR}
cd ${EXP_DYN_DIR}

#Change to current dir
local CUR_EXPDYN_DIR=${EXP_DYN_DIR}/exp_dyn_${CUR_ITER}
rm -rf ${CUR_EXPDYN_DIR}
mkdir -p ${CUR_EXPDYN_DIR}
cd ${CUR_EXPDYN_DIR}

cp ${IN_EQ_FILE} ${CUR_EXPDYN_DIR}/
mv ${CUR_EXPDYN_DIR}/${IN_EQ_FILE##*/} ${CUR_EXPDYN_DIR}/cur_initial_eq.txt
local IN_EQ_FILE=${CUR_EXPDYN_DIR}/cur_initial_eq.txt


echo "Start expdyn. Iteration=${CUR_ITER}" >> ${LOG}
# Generate input files for deBInfer
date >> ${LOG}
echo "Generate input files for deBInfer ${IN_EQ_FILE}. Iteration=${CUR_ITER}." >> ${LOG}
python3 $expdyn_gen_debinfer_in -e ${IN_EQ_FILE} -t ${TARGET_FILE} -n ${NUM_OF_POINTS_FROM_TARGET} \
|| { echo "Error in expdyn_gen_debinfer_in.py!" >> ${LOG}; exit 1; }
g++ ${CUR_EXPDYN_DIR}/updated_solve.cpp -o ${CUR_EXPDYN_DIR}/updated_solve
# First convert parameter ranges to the required format
#echo "Formalizing a prior distribution based on the specified parameter ranges. Iteration=${CUR_ITER}" >> ${LOG}
#Rscript $param_parser ${CUR_EXPDYN_DIR}/ParamRanges.txt ${CUR_EXPDYN_DIR}/parsed_userpriors.txt > /dev/null

# Run deBInfer
echo "Running deBInder. Iteration=${CUR_ITER}" >> ${LOG}
R CMD SHLIB ${CUR_EXPDYN_DIR}/updated_equation.c > /dev/null
Rscript $deBInfer_script ${MCMC_ITERATIONS} ${CUR_EXPDYN_DIR}/updated_target_trajectory.txt > /dev/null

head -n 2 ${CUR_EXPDYN_DIR}/deBInfer_output.txt > ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output_Nth.txt
mv ${CUR_EXPDYN_DIR}/deBInfer_output.txt ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output.txt
sed '1d' ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output.txt | awk -v s=${LOG_STEP} 'NR % s == 0' >> ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output_Nth.txt
date >> ${LOG}

#tail -n 1 ./${CUR_ITER}_deBInfer_output_Nth.txt | awk '{print $(NF)}' > ./step_list.txt

echo "Running Solver for expdyn result params. Iteration=${CUR_ITER}" >> ${LOG}
mkdir ${CUR_EXPDYN_DIR}/solve
cp ${IN_EQ_FILE} ${CUR_EXPDYN_DIR}/solve/updated_eq.txt
cat ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output_Nth.txt | sed '1d' | awk '{print $(NF)}' > ${CUR_EXPDYN_DIR}/solve/step_list.txt
mv ${CUR_EXPDYN_DIR}/updated_solve ${CUR_EXPDYN_DIR}/solve/updated_solve
cd ${CUR_EXPDYN_DIR}/solve
local UPDATED_EQ_FILE=${CUR_EXPDYN_DIR}/solve/updated_eq.txt
python3 $expdyn_res_to_equation -e ${UPDATED_EQ_FILE} -r ${CUR_EXPDYN_DIR}/repl_dict.txt -s ${CUR_EXPDYN_DIR}/solve/step_list.txt -d ${CUR_EXPDYN_DIR}/${CUR_ITER}_deBInfer_output.txt \
|| { echo "Error in expdyn_res_to_equation.py!" >> ${LOG}; exit 1; }
rm -rf ${CUR_EXPDYN_DIR}/solve/updated_eq.txt
echo -e "step\tresidual" > ${CUR_EXPDYN_DIR}/solve/residuals.txt
for step in `cat step_list.txt` ;
do
    python3 $solve -e ${step}_${UPDATED_EQ_FILE##*/} -f ${ST_T} -l ${END_T} -s 0.06 -a 0.001 > ${CUR_EXPDYN_DIR}/solve/${step}_exp_trajectory.txt \
    || { echo "Error in solve.py!" >> ${LOG}; exit 1; }
    echo -e "${step}\t`python3 $find_residual -t ${TARGET_FILE} -c ${step}_exp_trajectory.txt`" >> ${CUR_EXPDYN_DIR}/solve/residuals.txt
    local AFTER_EXP_DYN_EQ_FILE=`pwd`/${step}_${UPDATED_EQ_FILE##*/} # Get parameters from last step of expdyn mcmc
done
awk '{print $2}' ${CUR_EXPDYN_DIR}/solve/residuals.txt | paste ${CUR_EXPDYN_DIR}/solve/unified_params.txt - -d "\t" > ${CUR_EXPDYN_DIR}/unified_params_residuals_${CUR_ITER}.txt
rm ${CUR_EXPDYN_DIR}/solve/residuals.txt ${CUR_EXPDYN_DIR}/solve/unified_params.txt
rm ${CUR_EXPDYN_DIR}/solve/updated_solve
cd ${CUR_EXPDYN_DIR}
rm -rf ${CUR_EXPDYN_DIR}/updated_equation.o ${CUR_EXPDYN_DIR}/updated_equation.so

echo "End expdyn. Iteration=${CUR_ITER}" >> ${LOG}
echo ${AFTER_EXP_DYN_EQ_FILE}
}

get_unified_params_mod(){
    cd ${WORK_DIR}
    mkdir -p ./unified_params_exp
    rm -rf ./unified_params_exp/*
    mkdir -p ./unified_params_phy
    rm -rf ./unified_params_phy/*
    
    #cd ${WORK_DIR}/work
    cd ${WORK_DIR}/exp_dyn
    n=`ls -d exp_dyn*/ | wc -l`
    for d in $(seq 1 $n);
    do
    	cd ${WORK_DIR}/exp_dyn/exp_dyn_${d}
    	echo ${WORK_DIR}/exp_dyn/exp_dyn_${d}
	if [ -f ./unified_params_residuals*.txt ]; then
	    cp ./unified_params_residuals*.txt ${WORK_DIR}/unified_params_exp
	fi
    	#mv ${WORK_DIR}/unified_params_exp/unified_params_residuals.txt ${WORK_DIR}/unified_params_exp/unified_params_residuals_${n}.txt
    	#n=$((n+1))
    	cd ${WORK_DIR}/exp_dyn
    done
    
    cd ${WORK_DIR}/phy_dyn
    n=`ls -d beast*/ | wc -l`
    for d in $(seq 1 $n);
    do	
    	cd ${WORK_DIR}/phy_dyn/beast_${d}
	echo ${WORK_DIR}/phy_dyn/beast_${d}
	if [ -f ./unified_params_residuals*.txt ]; then
    	    awk '{print $6}' Posterior_init_parameters_variables_${d}.txt > treeLikelihood_col.txt && paste unified_params_residuals*.txt treeLikelihood_col.txt > ${WORK_DIR}/unified_params_phy/unified_params_residuals_${d}.txt
    	fi
    	#cp ./unified_params_residuals.txt ${WORK_DIR}/unified_params_phy
    	#mv ${WORK_DIR}/unified_params_phy/unified_params_residuals.txt ${WORK_DIR}/unified_params_phy/unified_params_residuals_${n}.txt
    	#n=$((n+1))
    	cd ${WORK_DIR}/phy_dyn
    done
    #cp ${WORK_DIR}/target_traj_* ${WORK_DIR}/unified_params_phy
    #cp ${WORK_DIR}/target_traj_* ${WORK_DIR}/unified_params_exp
    
    ##Find better residual and treeLikelihood
    cd ${WORK_DIR}/unified_params_phy/
    #mkdir -p ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood
    #rm -rf ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/*
    module load python3
    #for g in $(seq 1 9);
    #do
    #if [ -f ./unified_params_residuals*.txt ]; then
	python3 $find_good_resudual_likelihood -r ${wR} -p `ls unified_params_residuals* | paste -d, -s -` -n ${NUM_GOOD_RES_TR_L}
	mkdir -p ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood
	rm -rf ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/*
	GOOD_PARAM_SETS=good_residual_likelihood_result.txt
	for mcmc_step_gl_step in `awk -f ${get_two_cols} c1=mcmc_step c2=glob_iter ${GOOD_PARAM_SETS}`
	do
	    mcmc_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $1}'`
	    gl_step=`echo ${mcmc_step_gl_step} | awk -F_ '{print $2}'`
	    cp ${WORK_DIR}/phy_dyn/beast_${gl_step}/solve/${mcmc_step}_phy_trajectory.txt ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_phy_traj.txt
	    cp ${WORK_DIR}/phy_dyn/beast_${gl_step}/tree_pictures/STATE_${mcmc_step}.png ${WORK_DIR}/unified_params_phy/traj_and_tree_good_residual_likelihood/GL_${gl_step}_STATE_${mcmc_step}_tree.png
	done
    #else
	#echo "No suitable parameters found to start phylodynamics! Please change initial parameters or initial variables values!" > ${LOG}
    #fi
    #done 
    
    cd ${WORK_DIR}

}


if  [[ $1 && $2 && $3 && $4 && $5 && $6 && $7 && $8 && $9 ]]
then
    mkdir -p ${OUT_DIR}
    rm -rf ${OUT_DIR}/*
    mkdir -p ${WORK_DIR}
    cd ${WORK_DIR}
    #rm -rf ${WORK_DIR}/*
    date >> ${LOG}
    echo "Start pipe with params $1 $2 $3 $4 $5 $6 $7 $8 $9" >> ${LOG}
    
    ##take start, end time
    ST_T=`sort -t, -nk2 ${FASTA_INFO_FILE} | head -n 1 | awk -F, '{print $2}' | bc`
    END_T=`sort -t, -nk2 ${FASTA_INFO_FILE} | tail -n 1 | awk -F, '{print $2}' | bc`
    #ST_T=`echo $ST_T-0.1 | bc -l`
    #END_T=`echo $END_T+0.1 | bc -l`
    # First generate phylogenetic tree
    make_phytree
    for i in $(seq 1 $GLOB_ITERATIONS); 
    do 
	if [[ ${FIRST_DYN} == "expdyn" ]]
	then
	EQUATION_FILE=$(expdyn $i ${EQUATION_FILE} ${TARGET_FILE} ${ST_T} ${END_T})
	EQUATION_FILE=$(phydyn $i ${EQUATION_FILE} ${ST_T} ${END_T})
	if [[ ${EQUATION_FILE} == "BEAST failed" ]]
	then
	    break;
	fi
	else
	EQUATION_FILE=$(phydyn $i ${EQUATION_FILE} ${ST_T} ${END_T})
	if [[ ${EQUATION_FILE} == "BEAST failed" ]]
	then
	    break;
	fi
	EQUATION_FILE=$(expdyn $i ${EQUATION_FILE} ${TARGET_FILE} ${ST_T} ${END_T})
	fi
	#get_unified_params_mod
    done
    get_unified_params_mod
    cd ${WORK_DIR}
    cp -r ./* ${OUT_DIR}
    #End
    rm -rf ${WORK_DIR}
else
    echo "Parameters is not specified!"; exit 1;
fi
