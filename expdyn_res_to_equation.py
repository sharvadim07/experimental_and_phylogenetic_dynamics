import argparse

#"args": ["-ezODE4r_equation_v2_fixed_init_var.txt", "-rrepl_dict.txt", "-sstep_list.txt", "-ddeBInfer_output_Nth.txt"]
parser = argparse.ArgumentParser(description='debinfer output to unified equation file')
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
parser.add_argument('-r','--repl_dict', type=str,                   
                    help='Replacement dictionary file for parameters names.', required=True)
parser.add_argument('-s','--step_list', type=str,                   
                    help='List of number steps.', required=True)
parser.add_argument('-d','--debinfer_output', type=str,                   
                    help='Debinfer output with one set parameters for ODE.', required=True)

args = parser.parse_args()      

def read_debinfer_output(debinfer_output_file_name, repl_dict_file_name, step_list_file_name):
    repl_dict = {}
    step_list = []
    with open(step_list_file_name, 'r') as step_list_file:
        for step_line in step_list_file:
            step_list.append(int(step_line.strip()))
    with open(repl_dict_file_name, 'r') as repl_dict_file:
        for repl_line in repl_dict_file:
            repl_dict[repl_line.strip().split('\t')[1]] = repl_line.strip().split('\t')[0]
    with open(debinfer_output_file_name, 'r') as debinfer_output_file:
        par_init_values_dict = {}
        par_names_tuple = None
        par_val_tuple = None
        for i, debinfer_output_line in enumerate(debinfer_output_file):
            debinfer_output_line_list = debinfer_output_line.strip().split()
            if i == 0:
                par_names_tuple = tuple(repl_dict[deb_infer_par_name.replace('#','')] \
                                            for deb_infer_par_name in debinfer_output_line_list[:-1])
            elif int(debinfer_output_line_list[-1]) in step_list:
                step_num = debinfer_output_line_list[-1]
                par_val_tuple = (debinfer_output_line_list[:-1])
                if len(par_names_tuple) != len(par_val_tuple):
                    raise ValueError('debinfer output file format is bad!')
                for i, name in enumerate(par_names_tuple):
                    par_init_values_dict[str(step_num) + '_' + name] = 10.0**float(par_val_tuple[i]) # LOG10 to simple number
    return par_init_values_dict, step_list, par_names_tuple

def write_unified_params_file(unified_params_file_name, par_init_values_dict, step_list, par_names_tuple):
    with open(unified_params_file_name, 'w') as unified_params_file:
        unified_params_file.write('mcmc_step\t' + '\t'.join(par_names_tuple) + '\n')
        for step in step_list:
            unified_params_file.write(str(step) + '\t' + \
                                      '\t'.join([str(par_init_values_dict[str(step) + '_' + par_name]) for par_name in par_names_tuple]) + \
                                       '\n')

def read_write_init_values_file(init_values_file_name, par_init_values_dict, step_list):
    init_values_line_list = []
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
                    raw_name, value = line.strip().split('=')
                    par_spec = raw_name[0] + raw_name[2:] 
                    par_name = par_spec.split('_')[0] 
                    if '_initval' in raw_name and str(step) + '_' + par_name in par_init_values_dict:                                      
                        out_eq_init_values_file.write(raw_name + '=' + \
                                                      str(par_init_values_dict[str(step) + '_' + par_name]) + \
                                                      '\n')
                    else:
                        out_eq_init_values_file.write(line)                    
                else:
                    out_eq_init_values_file.write(line)

par_init_values_dict, step_list, par_names_tuple =  read_debinfer_output(args.debinfer_output, args.repl_dict, args.step_list)
write_unified_params_file('unified_params.txt', par_init_values_dict, step_list, par_names_tuple)
read_write_init_values_file(args.equations, par_init_values_dict, step_list)