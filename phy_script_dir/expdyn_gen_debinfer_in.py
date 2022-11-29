import logging
import argparse
import os
import re
import xml.etree.cElementTree as ET

#"args": ["-ezODE4r_equation_v2_fixed_init_var.txt", "-tMDBDPhyDyn_m1_first.traj", "-n10"]
parser = argparse.ArgumentParser(description='Generating files for ExpDyn and Solver')
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
parser.add_argument('-t','--target_traj', type=str,                   
                    help='Target trajectory in PhyDyn format.', required=False)
parser.add_argument('-n','--num_of_target_points', type=str,                   
                    help='Number of time points from trajectory.', required=False)
# parser.add_argument('-l','--mcmc_iter', type=int,                   
#                     help='Number of mcmc iterations.', required=True)



#Replace A -> x[1] ...
#repl_dict = {}


def read_equations_file(equations_file_name):
    equations_dict = {}    
    repl_dict = {}
    with open(equations_file_name, 'r') as equations_file:        
        eq_counter = 0
        for line in equations_file:
            if line[:2] == 'dx':
                eq_counter += 1
                eq_name_reres = re.search(r'dx\[\,\'(.+)\'\]=', line)
                eq_name = eq_name_reres.group(1)               
                equations_dict[eq_name] = []
                repl_dict[eq_name] = 'x[' + str(eq_counter) + ']'
                delimiters = "+(p", "+(-p"
                regexPattern = '(' + '|'.join(map(re.escape, delimiters)) + ')'
                line_str_terms = re.split(regexPattern, line.strip())[1:]
                line_str_terms_upd = []
                for i, str_term in enumerate(line_str_terms):
                    if i % 2 == 0:
                        line_str_terms_upd.append(line_str_terms[i] + line_str_terms[i + 1])
                line_str_terms = line_str_terms_upd
                for str_term in line_str_terms:
                    #equations_dict[eq_name].append(re.sub('[\[\]\'p\(\)\,x ]', '', str_term))
                    #term = re.sub('[\[\]\'p\,x ]', '', str_term.strip())[2:-1]
                    term = re.sub(r'(p\[\')|(x\[,\')|(\'\])|(^\+\()|(\)$)', '', str_term.strip())
                    term = re.sub(r'(\+-)|(-\+)', '-', term)
                    if term.count('(') != term.count(')'):
                        raise ValueError('Count of left brackets does not correpond count of right brackets!')
                    equations_dict[eq_name].append(term)        
    return equations_dict, repl_dict

def read_init_values_file(init_values_file_name, repl_dict):
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
                raw_name_splitted = raw_name.split('_')
                par_spec = raw_name
                par_name = '_'.join(raw_name_splitted[:-1])               
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
                var_spec = raw_name
                #var_spec = raw_name.split('_')[0] + '_' + raw_name.split('_')[2]
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
    eq_list = ['+'.join(new_eq_dict[eq]).replace('+-','-') for eq in new_eq_dict]
    f = symengine.sympify(eq_list) # Define function
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
                    eq[index_py_pow - i - 2] in '*/+-(':
                    break                   
        elif index_py_pow != 0:
            for i, ch in enumerate(eq[index_py_pow - 1::-1]):                
                if eq[index_py_pow - i - 1] in '*/+-(':
                    break 
                else:
                    cur_py_pow = ch + cur_py_pow
        exp_number = cur_py_pow.split('**')[0]
        exponent = cur_py_pow.split('**')[1]
        new_j_eq = eq.replace(cur_py_pow, 'pow(' + exp_number + ',' + exponent + ')')
    return new_j_eq

def edit_equation_C_code(c_code_pattern_file_name, equations_dict, params_counter, repl_dict):    
    with open(c_code_pattern_file_name, 'r') as c_code_pattern_file:
        with open('updated_' + os.path.basename(c_code_pattern_file_name), 'w') as updated_c_code_file:
            need_add_num_par = False
            need_add_derivs_eq = False
            need_add_derivs_sum_vars = False
            need_add_jacobian = False
            for c_code_patt_line in c_code_pattern_file:
                if '// global num of parameters' in c_code_patt_line:
                    need_add_num_par = True
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
                elif need_add_num_par == True:
                    updated_c_code_file.write('#define NUM_OF_PAR ' + str(params_counter + 1) + '\n')
                    need_add_num_par = False
                elif need_add_derivs_eq == True:
                    for eq in equations_dict:
                        new_eq = '+'.join(equations_dict[eq]).replace('+-','-')
                        for repl_ch in repl_dict:
                            #new_eq = new_eq.replace(repl_ch, repl_dict[repl_ch])
                            new_eq = re.sub(rf"(^|[\*\+\-\(\)\/\,]|\s)({repl_ch})($|[\*\+\-\(\)\/\,]|\s)", rf'\1{repl_dict[repl_ch]}\3', new_eq)
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
                                # W*k1 + G*M*a1/E**2 + G*(G + J)*r1/(E**2*(Lm + (G + J)/E)) - G*(G + J)**2*r1/(E**3*(Lm + (G + J)/E)**2)                               
                                #Old implementation
                                # re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                # while len(re_res_positive_pow) != 0:                                
                                #     new_j_eq = replace_py_pow_to_c_pow(new_j_eq, re_res_positive_pow)
                                #     re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 

                                re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                re_res_positive_pow += re.findall(r'.*(.{1,3}\*\*\(\-[0-9]+\)).*', new_j_eq)
                                while len(re_res_positive_pow) != 0:                                
                                    new_j_eq = replace_py_pow_to_c_pow(new_j_eq, re_res_positive_pow)
                                    re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq)
                                    re_res_positive_pow += re.findall(r'.*(.{1,3}\*\*\(\-[0-9]+\)).*', new_j_eq)
                                for repl_ch in repl_dict:
                                    #new_j_eq = new_j_eq.replace(repl_ch, repl_dict[repl_ch])
                                    new_j_eq = re.sub(rf"(^|[\*\+\-\(\)\/\,]|\s)({repl_ch})($|[\*\+\-\(\)\/\,]|\s)", rf'\1{repl_dict[repl_ch]}\3', new_j_eq)
                                updated_c_code_file.write('J('+ str(i) +', '+ str(j) +') = '+ new_j_eq +';\n')
                    need_add_jacobian = False
                else:
                    updated_c_code_file.write(c_code_patt_line) 

def edit_solve_C_code(c_code_solve_pattern_file_name, equations_dict, params_counter, var_init_values_dict, repl_dict): 
    with open(c_code_solve_pattern_file_name, 'r') as c_code_solve_pattern_file:
        with open('updated_' + os.path.basename(c_code_solve_pattern_file_name), 'w') as updated_c_code_solve_file:
            need_add_num_par = False
            need_add_num_eq = False
            need_add_derivs_eq = False
            need_add_derivs_sum_vars = False
            need_add_jacobian = False
            need_add_dfdt = False
            need_add_print_header = False
            need_add_init_var_val = False
            for c_code_patt_line in c_code_solve_pattern_file:
                if '// global num of parameters' in c_code_patt_line:
                    need_add_num_par = True
                    updated_c_code_solve_file.write(c_code_patt_line) 
                if '// global num of equations' in c_code_patt_line:
                    need_add_num_eq = True
                    updated_c_code_solve_file.write(c_code_patt_line) 
                elif '// derivs Equations start' in c_code_patt_line:
                    need_add_derivs_eq = True
                    updated_c_code_solve_file.write(c_code_patt_line) 
                elif '// derivs Sum of vars' in c_code_patt_line:
                    need_add_derivs_sum_vars = True
                    updated_c_code_solve_file.write(c_code_patt_line)
                elif '//jac Jacobian start' in c_code_patt_line:
                    need_add_jacobian = True
                    updated_c_code_solve_file.write(c_code_patt_line)  
                elif '//dfdt start' in c_code_patt_line:
                    need_add_dfdt = True
                    updated_c_code_solve_file.write(c_code_patt_line)  
                elif '// print header traj file' in c_code_patt_line:
                    need_add_print_header = True
                    updated_c_code_solve_file.write(c_code_patt_line) 
                elif '//initial variables values start' in c_code_patt_line:
                    need_add_init_var_val = True
                    updated_c_code_solve_file.write(c_code_patt_line) 
                elif need_add_num_par == True:
                    updated_c_code_solve_file.write('#define param_num  ' + str(params_counter) + '\n')
                    need_add_num_par = False
                elif need_add_num_eq == True:
                    updated_c_code_solve_file.write('#define dim ' + str(len(equations_dict)) + '\n')
                    need_add_num_eq = False
                elif need_add_derivs_eq == True:
                    for eq in equations_dict:
                        new_eq = '+'.join(equations_dict[eq]).replace('+-','-')
                        for repl_ch in repl_dict:
                            new_eq = re.sub(rf"(^|[\*\+\-\(\)\/\,]|\s)({repl_ch})($|[\*\+\-\(\)\/\,]|\s)", rf'\1{repl_dict[repl_ch]}\3', new_eq)
                        new_eq = 'res' + repl_dict[eq][1:] + '=' + new_eq
                        updated_c_code_solve_file.write(new_eq + ';\n')
                    need_add_derivs_eq = False
                elif need_add_derivs_sum_vars == True:
                    updated_c_code_solve_file.write('yout[0]=' + '+'.join([repl_dict[eq] for eq in equations_dict]) + ';\n')
                    need_add_derivs_sum_vars = False
                elif need_add_jacobian == True:
                    jac_mat = get_Jacobian(equations_dict)
                    for i in range(len(equations_dict) + 1):
                        for j in range(len(equations_dict) + 1):
                            if i == 0 or j == 0:
                                updated_c_code_solve_file.write('J('+ str(i) +', '+ str(j) +') = 0;\n')
                            else:
                                new_j_eq = str(jac_mat[i - 1][j - 1])
                                #re_res_positive_pow = re.findall(r'.*(\(.+\)\*\*[0-9]+).*', new_j_eq) 
                                # W*k1 + G*M*a1/E**2 + G*(G + J)*r1/(E**2*(Lm + (G + J)/E)) - G*(G + J)**2*r1/(E**3*(Lm + (G + J)/E)**2)                               
                                re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                re_res_positive_pow += re.findall(r'.*(.{1,3}\*\*\(\-[0-9]+\)).*', new_j_eq)
                                while len(re_res_positive_pow) != 0:                                
                                    new_j_eq = replace_py_pow_to_c_pow(new_j_eq, re_res_positive_pow)
                                    re_res_positive_pow = re.findall(r'.*(.{1,3}\*\*[0-9]+).*', new_j_eq) 
                                    re_res_positive_pow += re.findall(r'.*(.{1,3}\*\*\(\-[0-9]+\)).*', new_j_eq)
                                for repl_ch in repl_dict:
                                    #new_j_eq = new_j_eq.replace(repl_ch, repl_dict[repl_ch])
                                    new_j_eq = re.sub(rf"(^|[\*\+\-\(\)\/\,]|\s)({repl_ch})($|[\*\+\-\(\)\/\,]|\s)", rf'\1{repl_dict[repl_ch]}\3', new_j_eq)
                                updated_c_code_solve_file.write('J('+ str(i) +', '+ str(j) +') = '+ new_j_eq +';\n')
                    need_add_jacobian = False
                elif need_add_dfdt == True:
                    for i in range(len(equations_dict) + 1):
                        updated_c_code_solve_file.write('dfdt[' + str(i) + '] = 0;\n')
                    need_add_dfdt = False
                elif need_add_print_header == True:
                    updated_c_code_solve_file.write('cout << \"#time\\t' + \
                                                    '\\t'.join(equations_dict) + '\" << endl;')
                    need_add_print_header = False
                elif need_add_init_var_val == True:
                    for eq in equations_dict:                         
                        updated_c_code_solve_file.write(repl_dict[eq] + '=' + str(var_init_values_dict[eq + '_initval']) + ';\n' )
                    need_add_init_var_val = False
                else:
                    updated_c_code_solve_file.write(c_code_patt_line) 



def generate_param_ranges_file(param_ranges_file_name, par_init_values_dict, par_name_list, repl_dict):
    with open(param_ranges_file_name, 'w') as param_ranges_file:
        param_ranges_file.write('Param\tMin\tMax\t# Comment\n')
        for par in par_name_list:
            param_ranges_file.write(repl_dict[par].replace('[','').replace(']','') + '\t' + \
                                    str(par_init_values_dict[par + '_min']) + '\t' + \
                                    str(par_init_values_dict[par + '_max']) + '\t' + \
                                    '# ' + par + '\n')

def generate_user_priors_file(user_priors_file_name, par_init_values_dict, par_name_list, repl_dict):
    with open(user_priors_file_name, 'w') as user_priors_file:        
        import math
        for par in par_name_list:            
            min_val = float(par_init_values_dict[par + '_min'])
            max_val = float(par_init_values_dict[par + '_max'])
            gmean = 1.0
            span = 1.0
            distribution_type = 'norm'
            if (min_val - max_val) == 0:
                distribution_type = 'fixed'
                gmean = min_val
                span = 1.0
            elif min_val > max_val:
                raise ValueError ('ERROR: minimum value exceeds maximum value for a parameter!')
            elif min_val < 0:
                raise ValueError ('ERROR: negative value of a minimum value for a parameter!')
            else:
                gmean = (min_val * max_val) ** 0.5
                span = math.log10(max_val/gmean)            
            user_priors_file.write(repl_dict[par].replace('[','').replace(']','') + '\t' + \
                                    distribution_type + '\t' + \
                                    str(gmean) + '\t' + \
                                    str(span) + '\n')
                        
def gen_initial_values_deBInfer(initial_values_deBInfer_file_name, var_init_values_dict, repl_dict):
    with open(initial_values_deBInfer_file_name, 'w') as initial_values_deBInfer_file:
        for var_init_val in var_init_values_dict:
            if '_initval' in var_init_val:
                initial_values_deBInfer_file.write(repl_dict[var_init_val.split('_')[0]].replace('[','').replace(']','') + \
                                                    '=' + str(var_init_values_dict[var_init_val]) + '\n')

def write_names_dict(names_dict_file_name, repl_dict):
    with open(names_dict_file_name, 'w') as names_dict_file:
        for repl in repl_dict:
            names_dict_file.write(repl + '\t' + repl_dict[repl].replace('[','').replace(']','') + '\n')

def edit_target_traj_file(target_traj_file_name, num_of_target_points, repl_dict):
    num_lines_in_target_traj_file = sum(1 for line in open(target_traj_file_name, 'r'))
    with open(target_traj_file_name, 'r') as in_target_traj_file:        
        step = int((num_lines_in_target_traj_file-1)/int(num_of_target_points))
        with open('updated_target_trajectory.txt', 'w') as out_target_traj_file:
            line_in_list = []
            for i, line_in in enumerate(in_target_traj_file):
                line_in_splitted = line_in.strip().split('\t')
                if i == 0:
                    new_line_in_splitted = [repl_dict[re.sub('[0-9]', '', el_line_in_splitted)].replace('[','').replace(']','') 
                                        if el_line_in_splitted != 't' 
                                        else 'time'  
                                        for el_line_in_splitted in line_in_splitted[1:]]
                    out_target_traj_file.write('\t'.join(new_line_in_splitted) + '\n')
                elif i%step == 0 :
                    line_in_list.append(line_in_splitted[1:])
                    #out_target_traj_file.write('\t'.join(line_in_splitted[1:]) + '\n')
            line_in_list = sorted(line_in_list, key=lambda x: float(x[0]))
            for line_in_splitted in line_in_list:
                out_target_traj_file.write('\t'.join(line_in_splitted) + '\n')

def main():
    args = parser.parse_args()      
    equations_dict, repl_dict = read_equations_file(args.equations)
    par_init_values_dict, var_init_values_dict, par_name_list = read_init_values_file(args.equations, repl_dict)


    edit_equation_C_code(os.path.dirname(os.path.realpath(__file__)) + "/equation.c", equations_dict, len(par_name_list), repl_dict)
    edit_solve_C_code(os.path.dirname(os.path.realpath(__file__)) + "/solve.cpp", equations_dict, len(par_name_list), var_init_values_dict, repl_dict)
    #generate_param_ranges_file('ParamRanges.txt', par_init_values_dict, par_name_list)
    generate_user_priors_file('parsed_userpriors.txt', par_init_values_dict, par_name_list, repl_dict)
    gen_initial_values_deBInfer('initial_values_deBInfer_formatted.txt', var_init_values_dict, repl_dict)
    write_names_dict('repl_dict.txt', repl_dict)

    if args.target_traj != None and args.num_of_target_points != None:
        edit_target_traj_file(args.target_traj, args.num_of_target_points, repl_dict)

if __name__ == "__main__":
    main()

