import argparse

#"args": ["-ezODE4r_equation_v2_fixed_init_var.txt", "-f1", "-l60", "-s0.06", "0.000001"]
parser = argparse.ArgumentParser(description='Generating XML for PhyDyn BEAST')
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=True)
parser.add_argument('-f','--first_time', type=str,                   
                    help='First time point.', required=True)
parser.add_argument('-l','--last_time', type=str,                   
                    help='Last time point.', required=True)
parser.add_argument('-s','--step', type=str,                   
                    help='Step.', required=True)
parser.add_argument('-a','--accuracy', type=str,                   
                    help='Accuracy', required=True)

args = parser.parse_args()      

def read_init_values_file_run_solve(init_values_file_name):
    par_init_values_dict = {}    
    with open(init_values_file_name, 'r') as init_values_file:
        StartValuesOfparameters = False
        StartValuesOfequations = False        
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
                  
                try:
                    par_init_values_dict[par_spec] = float(value)
                except ValueError:
                    par_init_values_dict[par_spec] = value
        par_init_val_list = [par_init_values_dict[par] for par in par_init_values_dict if '_initval' in par]
        return par_init_val_list

def run_solve(t1, t2, step, accuracy, par_init_val_list):
    import subprocess
    param_str = [str(t1), str(t2), str(step), str(accuracy)]
    param_str = param_str + list(map(str, par_init_val_list))
    subprocess.call(['./updated_solve'] + param_str)

par_init_val_list = read_init_values_file_run_solve(args.equations)
run_solve(args.first_time, args.last_time, args.step, args.accuracy, par_init_val_list)
