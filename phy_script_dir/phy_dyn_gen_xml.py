import logging
import argparse
import xml.etree.cElementTree as ET
import expdyn_gen_debinfer_in

#
#            "args": ["-mzikv_shorah2_group_1_hapl_MSA.fasta", "-tbeast_init_last_corr.tree", "-ifasta_info.txt", "-etest_full_4_equation.txt", "-l20000", "-s5000", "-r1.1"]
parser = argparse.ArgumentParser(description='Generating XML for PhyDyn BEAST')
parser.add_argument('-m','--msa', type=str,                   
                    help='Input MSA file.', required=True)
parser.add_argument('-t','--phytree', type=str,                   
                    help='Phylogenetic Tree in nexus format.', required=True)
parser.add_argument('-i','--fasta_info', type=str,                   
                    help='MSA sequences info. Demes(as in the equations) and date(year.part_of_year)', required=True)
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
parser.add_argument('-l','--mcmc_iter', type=int,                   
                    help='Number of mcmc iterations.', required=True)
parser.add_argument('-s','--log_step', type=int,                   
                    help='Size of logging step.', required=True)
parser.add_argument('-r','--init_mutation_rate', type=float,                   
                    help='Mutation rate for sequences in MCMC.', required=True)
# parser.add_argument('-v','--init_values', type=str,                   
#                     help='Initial values file for equations (ODE).', required=True)


# parser.add_argument('-g','--genbank', type=str,                    
#                     help='Genbank file for reference sequence.', required=True)

  

#Replace A -> A1 ...
#repl_dict = {}

def read_fasta_msa_file(fasta_msa_file_name):
    fasta_sequnces_dict = {}    
    with open(fasta_msa_file_name, 'r') as fasta_msa_file:
        cur_seq_name = ''
        for line in fasta_msa_file:
            if line[0] == '>':
                cur_seq_name = line.strip()[1:] 
                fasta_sequnces_dict[cur_seq_name] = ''
            else:
                fasta_sequnces_dict[cur_seq_name] += line.strip()
    return fasta_sequnces_dict

def read_fasta_names_info(fasta_names_info_file_name, repl_dict):
    fasta_names_info_dict = {}    
    with open(fasta_names_info_file_name, 'r') as fasta_names_info_file:
        for line in fasta_names_info_file:
            line_list = line.strip().split(',')
            seq_name = line_list[0]
            time = line_list[1]
            #Replace A -> A1 ...
            deme = repl_dict[line_list[2]]
            fasta_names_info_dict[seq_name] = [time,deme]
    return fasta_names_info_dict


# def read_equations_file(equations_file_name):
#     equations_dict = {}
#     repl_dict = {}
#     with open(equations_file_name, 'r') as equations_file:
#         import re
#         for line in equations_file:
#             if line[:2] == 'dx':
#                 eq_name_reres = re.search(r'dx\[\,\'(.+)\'\]=', line)
#                 eq_name = eq_name_reres.group(1) + '1'
#                 repl_dict[eq_name_reres.group(1)] = eq_name_reres.group(1) + '1'
#                 #eq_name = eq_name_reres.group(1).replace(eq_name_reres.group(1), eq_name_reres.group(1) + '1')
#                 equations_dict[eq_name] = []
#                 delimiters = "+(p", "+(-p"
#                 regexPattern = '(' + '|'.join(map(re.escape, delimiters)) + ')'
#                 line_str_terms = re.split(regexPattern, line.strip())[1:]
#                 line_str_terms_upd = []
#                 for i, str_term in enumerate(line_str_terms):
#                     if i % 2 == 0:
#                         line_str_terms_upd.append(line_str_terms[i] + line_str_terms[i + 1])
#                 line_str_terms = line_str_terms_upd
#                 for str_term in line_str_terms:
#                     #equations_dict[eq_name].append(re.sub('[\[\]\'p\(\)\,x ]', '', str_term))
#                     term = re.sub('[\[\]\'p\,x ]', '', str_term.strip())[2:-1]
#                     if term.count('(') != term.count(')'):
#                         raise ValueError('Count of left brackets does not correpond count of right brackets!')
#                     equations_dict[eq_name].append(term)  

#         for eq in equations_dict:
#             for ch in repl_dict:
#                 equations_dict[eq] = [ term.replace(ch, repl_dict[ch]) for term in equations_dict[eq]]
#     return equations_dict, repl_dict

# def read_init_values_file(init_values_file_name, repl_dict):
#     par_init_values_dict = {}
#     var_init_values_dict = {}
#     par_name_list = []
#     param_counter = 0
#     var_name_list = []
#     with open(init_values_file_name, 'r') as init_values_file:
#         StartValuesOfparameters = False
#         StartValuesOfequations = False        
#         #var_counter = 0
#         for line in init_values_file:
#             if line == '\n':
#                 StartValuesOfparameters = False
#                 StartValuesOfequations = False
#             if line.strip() == '[StartValuesOfparameters]':
#                 StartValuesOfparameters = True
#                 StartValuesOfequations = False
#                 continue
#             if line.strip() == '[StartValuesOfequations]':
#                 StartValuesOfparameters = False
#                 StartValuesOfequations = True                
#                 continue
#             if StartValuesOfparameters == True and StartValuesOfequations == False:                           
#                 raw_name, value = line.strip().split('=')
#                 par_spec = raw_name[0] + raw_name[2:] 
#                 par_name = par_spec.split('_')[0]               
#                 if par_name not in par_name_list:
#                      param_counter += 1   
#                      par_name_list.append(par_name)
#                      repl_dict[par_name] = 'a[' + str(param_counter) + ']'  
#                 try:
#                     par_init_values_dict[par_spec] = float(value)
#                 except ValueError:
#                     par_init_values_dict[par_spec] = value
#             if StartValuesOfparameters == False and StartValuesOfequations == True:
#                 #var_counter += 1
#                 raw_name, value = line.strip().split('=')
#                 var_spec = repl_dict[raw_name.split('_')[0]] + '_' + raw_name.split('_')[2]
#                 var_name = var_spec.split('_')[0]
#                 if var_name not in var_name_list:
#                     var_name_list.append(var_name)
#                 try:
#                     var_init_values_dict[var_spec] = float(value)
#                 except ValueError:
#                     var_init_values_dict[var_spec] = value    
                                            
#     # Add missing specs for parameters
#     for par_name in par_name_list:
#         if par_name + '_min' not in par_init_values_dict:
#             #init_values_dict[par_var_name + '_min'] = init_values_dict[par_var_name + '_initval']            
#             par_init_values_dict[par_name + '_min'] = par_init_values_dict[par_name + '_initval']/10.0
#         if par_name + '_max' not in par_init_values_dict:    
#             #init_values_dict[par_var_name + '_max'] = init_values_dict[par_var_name + '_initval']
#             par_init_values_dict[par_name + '_max'] = par_init_values_dict[par_name + '_initval']*10.0
#     # Add missing specs for initialvariabes
#     for var_name in var_name_list:
#         if var_name + '_min' not in var_init_values_dict:
#             #init_values_dict[par_var_name + '_min'] = init_values_dict[par_var_name + '_initval']            
#             var_init_values_dict[var_name + '_min'] = var_init_values_dict[var_name + '_initval']/10.0
#         if var_name + '_max' not in var_init_values_dict:    
#             #init_values_dict[par_var_name + '_max'] = init_values_dict[par_var_name + '_initval']
#             var_init_values_dict[var_name + '_max'] = var_init_values_dict[var_name + '_initval']*10.0
#     return par_init_values_dict, var_init_values_dict

def find_term_birth(equations_dict, variable_i, variable_j):
    #Need find term with max priority
    import re
    term_priority_list = []
    for k, term in enumerate(equations_dict[variable_i]):
        finded_vars_in_term = re.findall(r"[\w']+", term)
        finded_vars_in_term = [var for var in finded_vars_in_term if not var.isdigit()]
        if  len(term) != 0 and term[0] != '-' and variable_j in finded_vars_in_term:      
            term_priority_list.append([term, finded_vars_in_term.index(variable_j), k])                              
    if len(term_priority_list) == 0:
        return ''
    max_priority = min([term_prio[1] for term_prio in term_priority_list])
    #Find priority for other var
    for variable in equations_dict:
        for term_priority in term_priority_list.copy():
            finded_vars_in_term = re.findall(r"[\w']+", term_priority[0])
            if variable in finded_vars_in_term and finded_vars_in_term.index(variable) < max_priority:
                term_priority_list.remove(term_priority)
    if len(term_priority_list) == 0:
        return ''
    max_priority_terms = ''
    for term_priority in term_priority_list:
        if term_priority[1] == max_priority:
            max_priority_terms += '+' + term_priority[0]    
            equations_dict[variable_i][term_priority[2]] = ''                 
    return max_priority_terms    

def find_term_death(full_list, variable):    
    import re
    term_priority_list = []
    for k, term in enumerate(full_list):
        finded_vars_in_term = re.findall(r"[\w']+", term)
        finded_vars_in_term = [var for var in finded_vars_in_term if not var.isdigit()]
        if  len(term) != 0 and variable in finded_vars_in_term:      
            term_priority_list.append([term, finded_vars_in_term.index(variable), k])                              
    if len(term_priority_list) == 0:
        return ''
    #Find the term with max priority   
    max_priority = min([term_prio[1] for term_prio in term_priority_list])
    max_priority_terms = ''
    for term_priority in term_priority_list:
        if term_priority[1] == max_priority:
            max_priority_terms += '+' + term_priority[0]    
            full_list[term_priority[2]] = ''                 
    return max_priority_terms   

#TODO:Check 1
def make_birth_matrix(equations_dict):    
    birth_matrix = [['' for x in range(len(equations_dict))] for y in range(len(equations_dict))]
    for i, variable_i in enumerate(equations_dict):
        birth_matrix[i][i] = find_term_birth(equations_dict, variable_i, variable_i)
        for j, variable_j in enumerate(equations_dict):
            if i != j:
                birth_matrix[i][j] = find_term_birth(equations_dict, variable_i, variable_j)                          
    #Del first +
    for i in range(len(equations_dict)):
        for j in range(len(equations_dict)):
            if len(birth_matrix[i][j]) > 0 and birth_matrix[i][j][0] == '+':
                birth_matrix[i][j] = birth_matrix[i][j][1:]
            elif len(birth_matrix[i][j]) == 0:
                birth_matrix[i][j] = '0'
    return birth_matrix
#TODO:Check 2
def make_death_vector(equations_dict):    
    death_vector_dict = {}
    for variable_i in equations_dict:
        death_vector_dict[variable_i] = ''
        for k, term in enumerate(equations_dict[variable_i]):
            if  len(term) != 0 and term[0] == '-':
                death_vector_dict[variable_i] += term
                equations_dict[variable_i][k] = ''  
    #With priority
    # full_list = []
    # for eq in equations_dict:
    #     for i, term in enumerate(equations_dict[eq]):
    #         full_list.append(term)
    # for variable_i in equations_dict:
    #     death_vector_dict[variable_i] = ''    
    # for variable_i in equations_dict:
    #     # Find max proirity terms from all terms
    #     death_vector_dict[variable_i] = find_term_death(full_list, variable_i)
    #With priority

    #Change - by + and del first +
    for variable_i in death_vector_dict:    
        if len(death_vector_dict[variable_i]) == 0:
            death_vector_dict[variable_i] = '0'
        death_vector_dict[variable_i] = death_vector_dict[variable_i].replace('-','+')
        death_vector_dict[variable_i] = death_vector_dict[variable_i].replace('++','+')
        if death_vector_dict[variable_i][0] == '+':
            death_vector_dict[variable_i] = death_vector_dict[variable_i][1:]                     
    return death_vector_dict

def add_sequences_from_alignment(tree_root, fasta_names_info_dict, fasta_sequences_dict):
    data = tree_root.find('data')
    for sequence in data.findall('sequence'):
        data.remove(sequence)
    for sequence in fasta_sequences_dict:
        seq_xml = ET.Element('sequence')
        seq_xml.attrib['id'] = sequence + '_' + fasta_names_info_dict[sequence][1] #fasta_names_info_dict[1] it is the deme
        seq_xml.attrib['taxon'] = sequence + '_' + fasta_names_info_dict[sequence][1] #fasta_names_info_dict[1] it is the deme
        seq_xml.attrib['totalcount'] = '4'
        seq_xml.attrib['value'] = fasta_sequences_dict[sequence]
        seq_xml.tail = '\n'
        data.append(seq_xml)
    #ET.dump(data)
    trait = tree_root.find('run').find('state').find('tree').find('trait')
    trait.text = ''    
    for i, sequence in enumerate(fasta_names_info_dict):        
        trait.text += sequence + '_' + fasta_names_info_dict[sequence][1] + \
                      '=' + \
                      fasta_names_info_dict[sequence][0]  #fasta_names_info_dict[0] it is time
        if i < len(fasta_names_info_dict) - 1:
            trait.text += ',\n'

def add_ODE_matrices_to_model(tree_root, birth_matrix, death_vector_dict, equations_dict):
    model = tree_root.find('model')
    for definition in model.findall('definition'):
        model.remove(definition)
    for matrixeq in model.findall('matrixeq'):
        model.remove(matrixeq)
    variables_list = equations_dict.keys()
    #Birth
    for i, var_i in enumerate(variables_list):
        for j, var_j in enumerate(variables_list):
            matrixeq_xml = ET.Element('matrixeq')
            matrixeq_xml.attrib['destination'] = var_i
            matrixeq_xml.attrib['origin'] = var_j
            matrixeq_xml.attrib['spec'] = 'MatrixEquation'
            matrixeq_xml.attrib['type'] = 'birth'
            matrixeq_xml.text = birth_matrix[i][j]
            matrixeq_xml.tail = '\n'
            model.append(matrixeq_xml)
    #Death
    for var_i in death_vector_dict:
        matrixeq_xml = ET.Element('matrixeq')
        matrixeq_xml.attrib['origin'] = var_i
        matrixeq_xml.attrib['spec'] = 'MatrixEquation'
        matrixeq_xml.attrib['type'] = 'death'
        matrixeq_xml.text = death_vector_dict[var_i]
        matrixeq_xml.tail = '\n'
        model.append(matrixeq_xml)

def add_newick_tree(tree_root, newick_phy_tree_file_name, fasta_names_info_dict):
    newick_tree = ''
    with open(newick_phy_tree_file_name,'r') as newick_phy_tree_file:
        newick_tree = ''.join(newick_phy_tree_file.readlines()).strip()
    for seq_name in fasta_names_info_dict:
        newick_tree = newick_tree.replace(seq_name, seq_name + '_' + fasta_names_info_dict[seq_name][1]) #fasta_names_info_dict[1] it is the deme
    init = tree_root.find('run').find('init')
    init.attrib['newick'] = newick_tree

def change_mcmc_num_iterations(tree_root, mcmc_iter):
    run = tree_root.find('run')
    run.attrib['chainLength'] = str(mcmc_iter)


def add_init_values(tree_root, par_init_values_dict, var_init_values_dict, equations_dict, fasta_names_info_dict, log_step, init_mutation_rate, integration_steps = 1000):
    import math
    #rates
    rates = tree_root.find('rates')    
    for param in rates.findall('param'):
        rates.remove(param)
    # Parameters
    for init_name in par_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            param_xml = ET.Element('param')
            param_xml.attrib['spec'] = 'ParamValue'
            param_xml.attrib['names'] = init_name
            param_xml.attrib['values'] = '@' + init_name
            param_xml.tail = '\n'
            rates.append(param_xml)
    
    #trajparams
    trajparams = tree_root.find('trajparams')
    #Find min time 
    min_time = min([float(fasta_names_info_dict[sequence][0]) for sequence in fasta_names_info_dict])
    trajparams.attrib['t0'] = str(min_time)
    trajparams.attrib['integrationSteps'] = str(integration_steps)
    for initialValue in trajparams.findall('initialValue'):
        trajparams.remove(initialValue)
    # Variables
    for init_name in var_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            initialValue_xml = ET.Element('initialValue')
            initialValue_xml.attrib['spec'] = 'ParamValue'
            initialValue_xml.attrib['names'] = init_name
            initialValue_xml.attrib['values'] = '@init' + init_name + '0'
            initialValue_xml.tail = '\n'
            trajparams.append(initialValue_xml)
   
    #run state
    run_state = tree_root.find('run').find('state')
    run_state.attrib['storeEvery'] = str(log_step)
    standart_parameters_tuple = ('clockRate.c', 
                           'kappa.s', 
                           'mutationRate.s', 
                           'proportionInvariant.s',
                           'gammaShape.s',
                           'freqParameter.s')
    for parameter in run_state.findall('parameter'):
        if parameter.attrib['id'].split(':')[0] not in standart_parameters_tuple:
            run_state.remove(parameter)
        elif parameter.attrib['id'].split(':')[0] == 'mutationRate.s':
            parameter.text = str(round(float(init_mutation_rate),2))
    # Parameters
    for init_name in par_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            parameter_xml = ET.Element('parameter')
            parameter_xml.attrib['id'] = init_name
            min_val = float(par_init_values_dict[init_name + '_min'])
            max_val = float(par_init_values_dict[init_name + '_max'])
            if min_val != max_val:
                parameter_xml.attrib['lower'] = str(min_val)
                parameter_xml.attrib['upper'] = str(max_val)
            parameter_xml.attrib['name'] = 'stateNode'
            parameter_xml.text = str(float(par_init_values_dict[init_name + '_initval']))
            parameter_xml.tail = '\n'
            run_state.append(parameter_xml)
    # Variables
    for init_name in var_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            parameter_xml = ET.Element('parameter')
            parameter_xml.attrib['id'] = 'init' + init_name + '0'
            min_val = float(var_init_values_dict[init_name + '_min'])
            max_val = float(var_init_values_dict[init_name + '_max'])
            if min_val != max_val:
                parameter_xml.attrib['lower'] = str(min_val)
                parameter_xml.attrib['upper'] = str(max_val)
            parameter_xml.attrib['name'] = 'stateNode'
            parameter_xml.text = str(float(var_init_values_dict[init_name + '_initval']))
            parameter_xml.tail = '\n'
            run_state.append(parameter_xml)
   
    #run distribution distribution
    run_distribution_distribution = tree_root.find('run').find('distribution').find('distribution')
    standart_priors_tuple = ('ClockPrior.c', 
                            'GammaShapePrior.s', 
                            'KappaPrior.s',
                            'PropInvariantPrior.s')
    for prior in run_distribution_distribution.findall('prior'):
        if prior.attrib['id'].split(':')[0] not in standart_priors_tuple:
            run_distribution_distribution.remove(prior)
    # Variables
    for init_name in var_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            min_val = float(var_init_values_dict[init_name + '_min'])
            max_val = float(var_init_values_dict[init_name + '_max'])
            if min_val != max_val:
                init_name = init_name.split('_')[0]
                prior_xml = ET.Element('prior')
                prior_xml.attrib['id'] = 'init' + init_name + '0prior'
                prior_xml.attrib['x'] = '@init' + init_name + '0'
                prior_xml.attrib['name'] = 'distribution'
                if var_init_values_dict[init_name + '_distribution'] == 'log-normal':
                    LogNormal_xml = ET.Element('LogNormal')           
                    LogNormal_xml.attrib['id'] = 'LogNormal:' + 'init' + init_name + '0'  
                    #Calc sd and mean
                    gmean = (min_val * max_val) ** 0.5
                    span = math.log(max_val/gmean)
                    mean = math.log(gmean)
                    sd = span/3
                    LogNormal_xml.attrib['M'] = str(mean)
                    LogNormal_xml.attrib['S'] = str(sd)
                    LogNormal_xml.attrib['name'] = 'distr'
                    LogNormal_xml.tail = '\n'
                    prior_xml.append(LogNormal_xml)
                elif var_init_values_dict[init_name + '_distribution'] == 'normal':
                    pass            
                prior_xml.tail = '\n'
                run_distribution_distribution.append(prior_xml)
    # Parameters    
    for init_name in par_init_values_dict:
        if '_initval' in init_name:
            init_name = init_name.split('_')[0]
            min_val = float(par_init_values_dict[init_name + '_min'])
            max_val = float(par_init_values_dict[init_name + '_max'])
            if min_val != max_val:
                init_name = init_name.split('_')[0]
                prior_xml = ET.Element('prior')
                prior_xml.attrib['id'] = init_name + 'prior'
                prior_xml.attrib['x'] = '@' + init_name             
                prior_xml.attrib['name'] = 'distribution'
                if par_init_values_dict[init_name + '_distribution'] == 'log-normal':
                    LogNormal_xml = ET.Element('LogNormal')           
                    LogNormal_xml.attrib['id'] = 'LogNormal:' + init_name
                    #Calc sd and mean
                    gmean = (min_val * max_val) ** 0.5
                    span = math.log(max_val/gmean)
                    mean = math.log(gmean)
                    sd = span/3
                    LogNormal_xml.attrib['M'] = str(mean)
                    LogNormal_xml.attrib['S'] = str(sd)
                    LogNormal_xml.attrib['name'] = 'distr'
                    LogNormal_xml.tail = '\n'
                    prior_xml.append(LogNormal_xml)
                elif par_init_values_dict[init_name + '_distribution'] == 'normal':
                    pass            
                prior_xml.tail = '\n'
                run_distribution_distribution.append(prior_xml)

    #run 
    run = tree_root.find('run')
    for operator in run.findall('operator'):
        if 'rwoperator' in operator.attrib['id']:
            run.remove(operator)
    
    # Variables
    for init_name in var_init_values_dict:
        if '_initval' in init_name:            
            init_name = init_name.split('_')[0]
            min_val = float(var_init_values_dict[init_name + '_min'])
            max_val = float(var_init_values_dict[init_name + '_max'])
            operator_xml = ET.Element('operator')
            init_name = init_name.split('_')[0]
            operator_xml.attrib['id'] = 'rwoperator:init' + init_name + '0'
            operator_xml.attrib['parameter'] = '@init' + init_name + '0'
            operator_xml.attrib['spec'] = 'RealRandomWalkOperator'
            operator_xml.attrib['windowSize'] = str(abs(max_val-min_val))
            operator_xml.attrib['useGaussian'] = 'true'
            operator_xml.attrib['weight'] = '1'
            operator_xml.tail = '\n'
            run.append(operator_xml)
    # Parameters
    for init_name in par_init_values_dict:
        if '_initval' in init_name:            
            init_name = init_name.split('_')[0]
            min_val = float(par_init_values_dict[init_name + '_min'])
            max_val = float(par_init_values_dict[init_name + '_max'])
            operator_xml = ET.Element('operator')
            init_name = init_name.split('_')[0]
            operator_xml.attrib['id'] = 'rwoperator:' + init_name
            operator_xml.attrib['parameter'] = '@' + init_name
            operator_xml.attrib['spec'] = 'RealRandomWalkOperator'
            operator_xml.attrib['windowSize'] = str(abs(max_val-min_val))
            operator_xml.attrib['useGaussian'] = 'true'
            operator_xml.attrib['weight'] = '1'
            operator_xml.tail = '\n'
            run.append(operator_xml)
    
    #run loggers
    loggers = tree_root.find('run').findall('logger')
    for logger in loggers:
        if logger.attrib['id'] in ('tracelog', 'screenlog'):
            # Variables
            for init_name in var_init_values_dict:
                if '_initval' in init_name:
                    init_name = init_name.split('_')[0]
                    log_xml = ET.Element('log')
                    log_xml.attrib['idref'] = 'init' + init_name + '0'
                    log_xml.tail = '\n'
                    logger.append(log_xml)    
            # Parameters
            for init_name in par_init_values_dict:
                if '_initval' in init_name:
                    init_name = init_name.split('_')[0]
                    log_xml = ET.Element('log')
                    log_xml.attrib['idref'] = init_name
                    log_xml.tail = '\n'
                    logger.append(log_xml)  
    for logger in loggers:
        logger.attrib['logEvery'] = str(log_step)

def editXML(filename, 
            fasta_names_info_dict, 
            fasta_sequnces_dict, 
            birth_matrix, 
            death_vector_dict, 
            equations_dict, 
            par_init_values_dict, 
            var_init_values_dict,
            mcmc_iter, 
            log_step, 
            init_mutation_rate,
            newick_phy_tree_file_name):
    tree = ET.ElementTree(file=filename)
    tree_root = tree.getroot()
    add_sequences_from_alignment(tree_root, fasta_names_info_dict, fasta_sequnces_dict)
    add_ODE_matrices_to_model(tree_root, birth_matrix, death_vector_dict, equations_dict)
    change_mcmc_num_iterations(tree_root, mcmc_iter)
    add_newick_tree(tree_root, newick_phy_tree_file_name, fasta_names_info_dict)
    add_init_values(tree_root, par_init_values_dict, var_init_values_dict, equations_dict, fasta_names_info_dict, log_step, init_mutation_rate)
    new_tree = ET.ElementTree(tree_root)
    #with open("pattern_model_seq_phy_dyn.xml_updated.xml", "w") as f:
    new_tree.write("pattern_model_seq_phy_dyn_v2_updated.xml")

def convert_phy_tree(tree_file_name):
    from Bio import Phylo
    Phylo.convert(tree_file_name,'nexus',tree_file_name + '.newick','newick')
    return tree_file_name + '.newick'

def main():
    import os
    args = parser.parse_args()

    equations_dict, repl_dict = expdyn_gen_debinfer_in.read_equations_file(args.equations)
    par_init_values_dict, var_init_values_dict = expdyn_gen_debinfer_in.read_init_values_file(args.equations, repl_dict)
    #check_init_val(equations_dict, init_values_dict)
    birth_matrix = make_birth_matrix(equations_dict)
    death_vector_dict = make_death_vector(equations_dict)
    fasta_names_info_dict = read_fasta_names_info(args.fasta_info, repl_dict)
    fasta_sequnces_dict = read_fasta_msa_file(args.msa)
    newick_phy_tree_file_name = convert_phy_tree(args.phytree)

    editXML(os.path.dirname(os.path.realpath(__file__)) + "/pattern_model_seq_phy_dyn_v2.xml", 
                        fasta_names_info_dict, 
                        fasta_sequnces_dict, 
                        birth_matrix, 
                        death_vector_dict, 
                        equations_dict, 
                        par_init_values_dict, 
                        var_init_values_dict, 
                        args.mcmc_iter, 
                        args.log_step, 
                        args.init_mutation_rate,
                        newick_phy_tree_file_name)

if __name__ == "__main__":
    main()
