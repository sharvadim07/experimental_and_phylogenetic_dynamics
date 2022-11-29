import argparse
import re

parser = argparse.ArgumentParser(description='Calculate residual for target and calculated trajectories.')
parser.add_argument('-t','--target_traj', type=str,                   
                    help='Target trajectory file.', required=True)
parser.add_argument('-c','--calc_traj', type=str,                   
                    help='Calculated trajectory file.', required=True)

args = parser.parse_args()  

def read_trajectory_file(traj_file):
    var_name_tuple = None
    time_values_list = []
    time_col_index = -1           
    sample_col_index = -1     
    for i, traj_line in enumerate(traj_file):
        traj_line_splitted = traj_line.strip().split('\t')
        if i == 0:   
            if 't' in traj_line_splitted:
                time_col_index = traj_line_splitted.index('t')
            elif 'time' in traj_line_splitted: 
                time_col_index = traj_line_splitted.index('time')
            elif '#time' in traj_line_splitted: 
                time_col_index = traj_line_splitted.index('#time')
            else:
                raise ValueError('Trajectory file not contains time column!')
            if 'Sample' in traj_line_splitted:
                sample_col_index = traj_line_splitted.index('Sample')
            var_name_tuple = tuple(re.sub('[0-9]', '', col_head) for col_head in traj_line_splitted if col_head != 'Sample' and \
                                                                                         col_head != 'time' and \
                                                                                         col_head != '#time' and \
                                                                                         col_head != 't')
        else:            
            time_values_list.append([float(traj_line_splitted[time_col_index]), i, \
                                     tuple(float(col_val) for i, col_val in enumerate(traj_line_splitted) if i != time_col_index and \
                                                                                        i != sample_col_index)])
    return var_name_tuple, time_values_list

def find_residual(target_traj_file_name, calc_traj_file_name):
    with open(target_traj_file_name, 'r') as target_traj_file, \
        open(calc_traj_file_name, 'r') as calc_traj_file:
        target_var_name_tuple, target_time_values_list = read_trajectory_file(target_traj_file)
        calc_var_name_tuple, calc_time_values_list = read_trajectory_file(calc_traj_file)
        target_calc_time_pairs_list = []     
        for k, target_time in enumerate([target_time_values[0] for target_time_values in target_time_values_list]):
            target_time = float(target_time)
            min_diff = -1
            min_diff_index = -1
            for i, calc_time in enumerate([calc_time_values[0] for calc_time_values in calc_time_values_list]):   
                calc_time = float(calc_time)
                if i == 0:
                    min_diff = abs(target_time - calc_time)    
                    min_diff_index = i        
                elif abs(target_time - calc_time) < min_diff:
                    min_diff = abs(target_time - calc_time)    
                    min_diff_index = i
            target_calc_time_pairs_list.append((k, min_diff_index))

        if len(target_calc_time_pairs_list) != len(target_time_values_list):
            raise ValueError('Times in target trajectory not matched times in calc trajectory!')
        
        #Calc residial for target columns
        sum_sq_diff = 0
        for target_var_name in target_var_name_tuple:
            col_index_target = target_var_name_tuple.index(target_var_name) 
            col_index_calc = calc_var_name_tuple.index(target_var_name)  
            for target_calc_time_pair in target_calc_time_pairs_list:
                sum_sq_diff += (target_time_values_list[target_calc_time_pair[0]][2][col_index_target] - \
                                calc_time_values_list[target_calc_time_pair[1]][2][col_index_calc])**2
        sum_sq_diff_rt = sum_sq_diff**0.5
        print (sum_sq_diff_rt)

find_residual(args.target_traj, args.calc_traj)

            

        