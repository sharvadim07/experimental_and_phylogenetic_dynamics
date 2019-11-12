import logging
import argparse
import os
import re
import xml.etree.cElementTree as ET

#
#            "args": ["-mzikv_shorah2_group_1_hapl_MSA.fasta", "-tbeast_init_last_corr.tree", "-ifasta_info.txt", "-etest_full_4_equation.txt", "-l20000", "-s5000", "-r1.1"]
parser = argparse.ArgumentParser(description='Generating XML for PhyDyn BEAST')
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
# parser.add_argument('-l','--mcmc_iter', type=int,                   
#                     help='Number of mcmc iterations.', required=True)

args = parser.parse_args()      

#Replace A -> x[1] ...
repl_dict = {}


def read_equations_file(equations_file_name):
    equations_dict = {}    
    with open(equations_file_name, 'r') as equations_file:        
        eq_counter = 0
        for line in equations_file:
            if line[:2] == 'dx':
                eq_counter += 1
                eq_name_reres = re.search(r'dx\[\,\'(.+)\'\]=', line)
                eq_name = eq_name_reres.group(1)               
                equations_dict[eq_name] = []
                repl_dict[eq_name] = 'x[' + str(eq_counter) + ']'
                line_str_terms = line.strip().split('+(')[1:]
                for str_term in line_str_terms:
                    #equations_dict[eq_name].append(re.sub('[\[\]\'p\(\)\,x ]', '', str_term))
                    equations_dict[eq_name].append(re.sub('[\[\]\'p\,x ]', '', str_term)[:-1])        
    return equations_dict

def read_init_values_file(init_values_file_name):
    par_init_values_dict = {}
    var_init_values_dict = {}
    par_name_list = []
    param_counter = 0
    #var_name_list = []
    with open(init_values_file_name, 'r') as init_values_file:
        StartValuesOfparameters = False
        StartValuesOfequations = False        
        #var_counter = 0
        for line in init_values_file:
            if line == '\n':
                StartValuesOfparameters = False
                StartValuesOfequations = False
            if line.strip() == '[StartValuesOfparameters]':
                StartValuesOfparameters = True
                StartValuesOfequations = False
                continue
            if line.strip() == '[StartValuesOfequations]':
                StartValuesOfparameters = False
                StartValuesOfequations = True                
                continue
            if StartValuesOfparameters == True and StartValuesOfequations == False:                           
                raw_name, value = line.strip().split('=')
                par_spec = raw_name[0] + raw_name[2:] 
                par_name = par_spec.split('_')[0]               
                if par_name not in par_name_list:
                     param_counter += 1   
                     par_name_list.append(par_name)
                     repl_dict[par_name] = 'a[' + str(param_counter) + ']'  
                try:
                    par_init_values_dict[par_spec] = float(value)
                except ValueError:
                    par_init_values_dict[par_spec] = value
            if StartValuesOfparameters == False and StartValuesOfequations == True:
                #var_counter += 1
                raw_name, value = line.strip().split('=')
                var_spec = raw_name.split('_')[0] + '_' + raw_name.split('_')[2]
                #var_name = var_spec.split('_')[0]
                #var_name_list.append(var_name)
                try:
                    var_init_values_dict[var_spec] = float(value)
                except ValueError:
                    var_init_values_dict[var_spec] = value    
                                            
    # Add missing specs for parameters
    for par_name in par_name_list:
        if par_name + '_min' not in par_init_values_dict:
            #init_values_dict[par_var_name + '_min'] = init_values_dict[par_var_name + '_initval']            
            par_init_values_dict[par_name + '_min'] = par_init_values_dict[par_name + '_initval']/10.0
        if par_name + '_max' not in par_init_values_dict:    
            #init_values_dict[par_var_name + '_max'] = init_values_dict[par_var_name + '_initval']
            par_init_values_dict[par_name + '_max'] = par_init_values_dict[par_name + '_initval']*10.0
    return par_init_values_dict, var_init_values_dict, par_name_list

def get_Jacobian(equations_dict):
    import symengine
    #Temporarly replace E variable to RE
    new_eq_dict = {}
    for eq in equations_dict:
        if eq == 'E':            
            new_eq_dict['RE'] =  [term.replace('E', 'RE') for term in equations_dict[eq]]
        else:
            new_eq_dict[eq] = [term.replace('E', 'RE') for term in equations_dict[eq]]

    variables = symengine.symbols(' '.join(new_eq_dict)) # Define variables
    f = symengine.sympify(['+'.join(new_eq_dict[eq]).replace('+-','-') for eq in new_eq_dict]) # Define function
    #J = symengine.zeros(len(f),len(variables)) # Initialise Jacobian matrix
    J_str = [['' for x in range(len(variables))] for y in range(len(f))] # Initialise Jacobian matrix with strings
    # Fill Jacobian matrix with entries
    for i, fi in enumerate(f):
        for j, s in enumerate(variables):
            J_str[i][j] = str(symengine.diff(fi, s)).replace('RE', 'E')
    return J_str

def replace_py_pow_to_c_pow(eq, re_res_positive_pow):
    for py_pow in re_res_positive_pow:
        index_py_pow = eq.index(py_pow)
        cur_py_pow = py_pow
        if py_pow[0] == ')':
            bracket_right_counter = 1
            bracket_left_counter = 0            
            for i, ch in enumerate(eq[index_py_pow - 1::-1]):
                if ch == ')':
                    bracket_right_counter += 1
                elif ch == '(':
                    bracket_left_counter += 1
                cur_py_pow = ch + cur_py_pow
                if bracket_right_counter == bracket_left_counter and \
                    eq[index_py_pow - i - 2] in '*/+-':
                    break                   
        else:
            for i, ch in enumerate(eq[index_py_pow - 1::-1]):                
                if eq[index_py_pow - i - 1] in '*/+-':
                    break 
                else:
                    cur_py_pow = ch + cur_py_pow
        exp_number = cur_py_pow.split('**')[0]
        exponent = cur_py_pow.split('**')[1]
        new_j_eq = eq.replace(cur_py_pow, 'pow(' + exp_number + ',' + exponent + ')')
    return new_j_eq

def edit_equation_C_code(c_code_pattern_file_name, equations_dict, params_counter):    
    with open(c_code_pattern_file_name, 'r') as c_code_pattern_file:
        with open('updated_' + os.path.basename(c_code_pattern_file_name), 'w') as updated_c_code_file:
            need_add_num_eq = False
            need_add_derivs_eq = False
            need_add_derivs_sum_vars = False
            need_add_jacobian = False
            for c_code_patt_line in c_code_pattern_file:
                if '// global num of parameters' in c_code_patt_line:
                    need_add_num_eq = True
                    updated_c_code_file.write(c_code_patt_line) 
                elif '// derivs Equations start' in c_code_patt_line:
                    need_add_derivs_eq = True
                    updated_c_code_file.write(c_code_patt_line) 
                elif '// derivs Sum of vars' in c_code_patt_line:
                    need_add_derivs_sum_vars = True
                    updated_c_code_file.write(c_code_patt_line)
                elif '//jac Jacobian start' in c_code_patt_line:
                    need_add_jacobian = True
                    updated_c_code_file.write(c_code_patt_line)  
                elif need_add_num_eq == True:
                    updated_c_code_file.write('#define NUM_OF_PAR ' + str(params_counter + 1) + '\n')
                    need_add_num_eq = False
                elif need_add_derivs_eq == True:
                    for eq in equations_dict:
                        new_eq = '+'.join(equations_dict[eq]).replace('+-','-')
                        for repl_ch in repl_dict:
                            new_eq = new_eq.replace(repl_ch, repl_dict[repl_ch])
                        new_eq = 'res' + repl_dict[eq][1:] + '=' + new_eq
                        updated_c_code_file.write(new_eq + ';\n')
                    need_add_derivs_eq = False
                elif need_add_derivs_sum_vars == True:
                    updated_c_code_file.write('yout[0]=' + '+'.join([repl_dict[eq] for eq in equations_dict]) + ';\n')
                    need_add_derivs_sum_vars = False
                elif need_add_jacobian == True:
                    jac_mat = get_Jacobian(equations_dict)
                    for i in range(len(equations_dict) + 1):
                        for j in range(len(equations_dict) + 1):
                            if i == 0 or j == 0:
                                updated_c_code_file.write('J('+ str(i) +', '+ str(j) +') = 0;\n')
                            else:
                                new_j_eq = str(jac_mat[i - 1][j - 1])
                                #re_res_positive_pow = re.findall(r'.*(\(.+\)\*\*[0-9]+).*', new_j_eq)                                
                                re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                while len(re_res_positive_pow) != 0:                                
                                    new_j_eq = replace_py_pow_to_c_pow(new_j_eq, re_res_positive_pow)
                                    re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                for repl_ch in repl_dict:
                                    new_j_eq = new_j_eq.replace(repl_ch, repl_dict[repl_ch])
                                updated_c_code_file.write('J('+ str(i) +', '+ str(j) +') = '+ new_j_eq +';\n')
                    need_add_jacobian = False
                else:
                    updated_c_code_file.write(c_code_patt_line) 

def generate_param_ranges_file(param_ranges_file_name, par_init_values_dict, par_name_list):
    with open(param_ranges_file_name, 'w') as param_ranges_file:
        param_ranges_file.write('Param\tMin\tMax\t# Comment\n')
        for par in par_name_list:
            param_ranges_file.write(repl_dict[par].replace('[','').replace(']','') + '\t' + \
                                    str(par_init_values_dict[par + '_min']) + '\t' + \
                                    str(par_init_values_dict[par + '_max']) + '\t' + \
                                    '# ' + par + '\n')
                        
def gen_initial_values_deBInfer(initial_values_deBInfer_file_name, var_init_values_dict):
    with open(initial_values_deBInfer_file_name, 'w') as initial_values_deBInfer_file:
        for var_init_val in var_init_values_dict:
            if '_initval' in var_init_val:
                initial_values_deBInfer_file.write(repl_dict[var_init_val.split('_')[0]].replace('[','').replace(']','') + \
                                                    '=' + str(var_init_values_dict[var_init_val]) + '\n')

def write_names_dict(names_dict_file_name):
    with open(names_dict_file_name, 'w') as names_dict_file:
        for repl in repl_dict:
            names_dict_file.write(repl + '\t' + repl_dict[repl].replace('[','').replace(']','') + '\n')

equations_dict = read_equations_file(args.equations)
par_init_values_dict, var_init_values_dict, par_name_list = read_init_values_file(args.equations)


edit_equation_C_code(os.path.dirname(os.path.realpath(__file__)) + "/equation.c", equations_dict, len(par_name_list))
generate_param_ranges_file('ParamRanges.txt', par_init_values_dict, par_name_list)
gen_initial_values_deBInfer('initial_values_deBInfer_formatted.txt', var_init_values_dict)
write_names_dict('repl_dict.txt')