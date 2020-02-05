import argparse

#"args": ["-ezODE4r_equation_v2_fixed_init_var.txt", "-rrepl_dict.txt", "-sstep_list.txt", "-ddeBInfer_output_Nth.txt"]
parser = argparse.ArgumentParser(description='Generate trajectories for every replicate and passage')
parser.add_argument('-t','--table', type=str,                   
                    help='Table wiith cells numbers.', required=True)
parser.add_argument('-r','--replicate', type=str,                   
                    help='Replicate list separated by comma.', required=False)

args = parser.parse_args()   

pass_dict = {}
repl_list = []
with open(args.table, 'r') as table_file:
    for i, line in enumerate(table_file):
        if i == 0:
            repl_list = line.strip().split('\t')[1:]
        else:
            line_splitted = line.strip().split('\t')
            pass_dict[line_splitted[0]] = line_splitted[1:]

if args.replicate != None:   
    import math
    for passage in pass_dict:
        with open('target_traj/target_traj_' + str(passage) + '_W_D_A_E_B_var.txt', 'w') as out_for_passage_file:
            out_for_passage_file.write('Sample\tt\tW1\tD1\tA1\tE1\tB1\n')
            out_for_passage_file.write('0\t1.101\t83500\t1670\t50000\t1\t1\n')
            for replicate in args.replicate.split(','):
                W_val = int(pass_dict[passage][repl_list.index(replicate)])
                D_val = int(float(W_val)* 0.02)
                E_val = int(50000*math.log(W_val)/math.log(max(map(int, pass_dict[passage]))) - 5000)
                A_val = int((50000 - E_val)/7 + 1000)
                B_val = int(50000 - (A_val + E_val))
                out_for_passage_file.write('0\t59.999\t' + str(W_val) + \
                                                    '\t' + str(D_val) + \
                                                    '\t' + str(A_val) + \
                                                    '\t' + str(E_val) + \
                                                    '\t' + str(B_val) + '\n')
else:
    for passage in pass_dict:
        for replicate in repl_list:
            with open('target_traj/target_traj_' + str(passage) + '_' + str(replicate) + '_W_D_A_E_B_var.txt', 'w') as out_for_passage_file:
                out_for_passage_file.write('Sample\tt\tW1\tD1\tA1\tE1\tB1\n')
                out_for_passage_file.write('0\t1.101\t83500\t1670\t50000\t1\t1\n')
                W_val = int(pass_dict[passage][repl_list.index(replicate)])
                D_val = int(float(W_val)* 0.02)
                E_val = int(50000*math.log(W_val)/math.log(max(map(int, pass_dict[passage]))) - 5000)
                A_val = int((50000 - E_val)/7 + 1000)
                B_val = int(50000 - (A_val + E_val))
                out_for_passage_file.write('0\t59.999\t' + str(W_val) + \
                                                    '\t' + str(D_val) + \
                                                    '\t' + str(A_val) + \
                                                    '\t' + str(E_val) + \
                                                    '\t' + str(B_val) + '\n')
                