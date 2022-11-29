import argparse

parser = argparse.ArgumentParser(description='debinfer output to unified equation file')
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
parser.add_argument('-s','--step_list', type=str,                   
                    help='List of number steps.', required=True)
parser.add_argument('-p','--phydyn_output', type=str,                   
                    help='Phydyn output with one set parameters for ODE.', required=True)
#parser.add_argument('-m','--multiple_range', type=str,                   
#                    help='Multiplier and divider value for min and max value for paramater range.', required=True)

args = parser.parse_args()      

def read_phydyn_output(phydyn_output_file_name, step_list_file_name):    
    step_list = []
    with open(step_list_file_name, 'r') as step_list_file:
        for step_line in step_list_file:
            step_list.append(int(step_line.strip()))
    with open(phydyn_output_file_name, 'r') as phydyn_output_file:
        par_init_values_dict = {}
        par_names_tuple = None
        par_val_tuple = None
        for i, phydyn_output_line in enumerate(phydyn_output_file):
            phydyn_output_line_list = phydyn_output_line.strip().split()
            if phydyn_output_line_list[0] == 'Sample':
                par_names_tuple = tuple(phydyn_output_line_list[1:])
            elif phydyn_output_line_list[0].isdigit() and int(phydyn_output_line_list[0]) in step_list:
                step_num = phydyn_output_line_list[0]
                par_val_tuple = (phydyn_output_line_list[1:])
                if len(par_names_tuple) != len(par_val_tuple):
                    raise ValueError('Phydyn output file format is bad!')
                for i, name in enumerate(par_names_tuple):
                    par_init_values_dict[str(step_num) + '_' + name] = float(par_val_tuple[i])
    return par_init_values_dict, step_list

def write_unified_params_file(unified_params_file_name, par_init_values_dict, step_list, par_names_list):
    with open(unified_params_file_name, 'w') as unified_params_file:
        unified_params_file.write('mcmc_step\t' + '\t'.join(par_names_list) + '\n')
        for step in step_list:
            unified_params_file.write(str(step) + '\t' + \
                                      '\t'.join([str(par_init_values_dict[str(step) + '_' + par_name]) for par_name in par_names_list]) + \
                                       '\n')

def find_min_max_value(init_values_line_list, raw_name):    
    max_val = ''
    min_val = ''
    for line in init_values_line_list:
        if raw_name + '_max' in line:
            max_val = line.strip().split('=')[1]
        elif raw_name + '_min' in line:
            min_val = line.strip().split('=')[1]
    return float(min_val), float(max_val)

def read_write_init_values_file(init_values_file_name, par_init_values_dict, step_list):
    init_values_line_list = []
    par_names_list = []
    with open(init_values_file_name, 'r') as init_values_file:
        for line in init_values_file:
            init_values_line_list.append(line)
    import os
    for step in step_list:
        with open(str(step) + '_' + os.path.basename(init_values_file_name), 'w') as out_eq_init_values_file:            
            StartValuesOfparameters = False
            StartValuesOfequations = False        
            #var_counter = 0
            for line in init_values_line_list:
                if line == '\n':
                    StartValuesOfparameters = False
                    StartValuesOfequations = False
                if line.strip() == '[StartValuesOfparameters]':
                    StartValuesOfparameters = True
                    StartValuesOfequations = False
                    out_eq_init_values_file.write(line) 
                    continue
                if line.strip() == '[StartValuesOfequations]':
                    StartValuesOfparameters = False
                    StartValuesOfequations = True  
                    out_eq_init_values_file.write(line)               
                    continue
                if StartValuesOfparameters == True and StartValuesOfequations == False:                           
                    # raw_name, value = line.strip().split('=')
                    # par_spec = raw_name[0] + raw_name[2:] 
                    # par_name = par_spec.split('_')[0] 
                    raw_name, value = line.strip().split('=')
                    raw_name_splitted = raw_name.split('_')
                    par_spec = raw_name
                    par_name = '_'.join(raw_name_splitted[:-1]) 
                    if '_initval' in raw_name and str(step) + '_' + par_name in par_init_values_dict:                                      
                        out_eq_init_values_file.write(raw_name + '=' + \
                                                      str(par_init_values_dict[str(step) + '_' + par_name]) + \
                                                      '\n')
                        if par_name not in par_names_list:
                            par_names_list.append(par_name)
                    elif '_max' in raw_name and str(step) + '_' + par_name in par_init_values_dict:
                        min_val, max_val = find_min_max_value(init_values_line_list, raw_name.replace('_max', ''))
                        mul_val = (max_val/min_val)**0.5
                        out_eq_init_values_file.write(raw_name + '=' + \
                                                      str(par_init_values_dict[str(step) + '_' + par_name] * mul_val) + \
                                                      '\n')
                    elif '_min' in raw_name and str(step) + '_' + par_name in par_init_values_dict:
                        min_val, max_val = find_min_max_value(init_values_line_list, raw_name.replace('_min', ''))
                        mul_val = (max_val/min_val)**0.5
                        out_eq_init_values_file.write(raw_name + '=' + \
                                                      str(par_init_values_dict[str(step) + '_' + par_name] / mul_val) + \
                                                      '\n')
                    else:
                        out_eq_init_values_file.write(line)                    
                else:
                    out_eq_init_values_file.write(line)
    return par_names_list

par_init_values_dict, step_list = read_phydyn_output(args.phydyn_output, args.step_list)
par_names_list = read_write_init_values_file(args.equations, par_init_values_dict, step_list)
write_unified_params_file('unified_params.txt', par_init_values_dict, step_list, par_names_list)
