import logging
import argparse
import os
import re
import numpy
import phy_dyn_gen_xml

parser = argparse.ArgumentParser(
    description='Find good residual for PhyDyn results.')
parser.add_argument('-p', '--par_tables', type=str,
                    help='Parameters tables.', required=True)
parser.add_argument('-n', '--num_tuples', type=str,
                    help='Num tuples for search.', required=True)
parser.add_argument('--only_residual',
                    help='Use only residual column for search.', action='store_true', required=False)
parser.add_argument('--only_likelihood',
                    help='Use only tree likelihood column for search.', action='store_true', required=False)
parser.add_argument('-r', '--residual_coef', type=str,
                    help='residual_coef in weighted_average', required=False)
parser.add_argument('-e','--equations', type=str,                   
                    help='Equations (ODE).', required=False)
# parser.add_argument('-t', '--treeLikelihood_coef', type=str,
#                     help='treeLikelihood_coef in weighted_average', required=False)


class ParTable:
    def __init__(self, par_table_file_name, special_features_cols, add_features_name_val_list):
        self.header = None
        self.param_tuple_list = []
        self.special_features_cols_idx_dict = {}
        with open(par_table_file_name, 'r') as par_table_file:
            for i, line in enumerate(par_table_file):
                if i == 0:
                    self.header = tuple(line.strip().split('\t') + [add_feature_name_val[0]
                                                                    for add_feature_name_val in add_features_name_val_list])
                    for special_feature in special_features_cols:
                        self.special_features_cols_idx_dict[special_feature] = \
                            self.header.index(special_feature)
                else:
                    splitted_line_tuple = tuple(line.strip().split('\t') + [add_feature_name_val[1]
                                                                            for add_feature_name_val in add_features_name_val_list])
                    if len(splitted_line_tuple) != len(self.header):
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')
                    self.param_tuple_list.append(splitted_line_tuple)

    def append(self, par_table_file_name, add_features_name_val_list):
        with open(par_table_file_name, 'r') as par_table_file:
            for i, line in enumerate(par_table_file):
                if i != 0:
                    splitted_line_tuple = tuple(line.strip().split('\t') + [add_feature_name_val[1]
                                                                            for add_feature_name_val in add_features_name_val_list])
                    if len(splitted_line_tuple) != len(self.header):
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')
                    self.param_tuple_list.append(splitted_line_tuple)

    def __len__(self):
        return len(self.param_tuple_list)

    def __getitem__(self, position):
        return self.param_tuple_list[position]

# def find_good_residual_likelihood(par_table, init_residual_range = 0.01, treeLikelihood_range = 0.2):
#     residual_range = init_residual_range
#     param_tuple_sorted_list = sorted(par_table.param_tuple_list, key=lambda x: (float(x[par_table.header.index('residual')])), reverse=False)
#     max_treeLikelihood = max([float(cur_tuple[par_table.header.index('treeLikelihood')]) for cur_tuple in param_tuple_sorted_list])

#     param_tuple_sorted_in_residual_range_list = [cur_tuple for cur_tuple in param_tuple_sorted_list \
#         if float(cur_tuple[par_table.header.index('residual')]) <= \
#             float(param_tuple_sorted_list[0][par_table.header.index('residual')]) + \
#                 float(param_tuple_sorted_list[0][par_table.header.index('residual')])*residual_range]
#     treeLikelihood_in_residual_range_list = [float(cur_tuple[par_table.header.index('treeLikelihood')]) \
#         for cur_tuple in param_tuple_sorted_in_residual_range_list]
#     max_treeLikelihood_in_residual_range = max(treeLikelihood_in_residual_range_list)
#     while not (max_treeLikelihood - max_treeLikelihood*treeLikelihood_range) >= \
#         max_treeLikelihood_in_residual_range >= \
#         (max_treeLikelihood + max_treeLikelihood*treeLikelihood_range):
#         residual_range = residual_range + 0.01

#         param_tuple_sorted_in_residual_range_list = [cur_tuple for cur_tuple in param_tuple_sorted_list \
#         if float(cur_tuple[par_table.header.index('residual')]) <= \
#             float(param_tuple_sorted_list[0][par_table.header.index('residual')]) + \
#                 float(param_tuple_sorted_list[0][par_table.header.index('residual')])*residual_range]
#         treeLikelihood_in_residual_range_list = [float(cur_tuple[par_table.header.index('treeLikelihood')]) \
#             for cur_tuple in param_tuple_sorted_in_residual_range_list]
#         max_treeLikelihood_in_residual_range = max(treeLikelihood_in_residual_range_list)

#     max_treeLikelihood_in_residual_range_index = treeLikelihood_in_residual_range_list.index(max_treeLikelihood_in_residual_range)
#     return param_tuple_sorted_in_residual_range_list[max_treeLikelihood_in_residual_range_index]


def find_good_residual_likelihood(par_table, residual_coef=0.8):
    treeLikelihood_coef = 1 - residual_coef
    param_tuple_sorted_list = sorted(par_table.param_tuple_list, key=lambda x: (
        float(x[par_table.header.index('residual')])), reverse=False)
    max_treeLikelihood = max([float(cur_tuple[par_table.header.index(
        'treeLikelihood')]) for cur_tuple in param_tuple_sorted_list])
    min_residual = float(param_tuple_sorted_list[0][par_table.header.index(
        'residual')])

    residual_ratio_list = \
        [float(cur_tuple[par_table.header.index('residual')]) /
         min_residual for cur_tuple in param_tuple_sorted_list]
    treeLikelihood_ratio_list = \
        [float(cur_tuple[par_table.header.index('treeLikelihood')]) /
         max_treeLikelihood for cur_tuple in param_tuple_sorted_list]

    weighted_average_sum_list = \
        [residual_ratio*residual_coef + treeLikelihood_ratio_list[i]*treeLikelihood_coef
            for i, residual_ratio in enumerate(residual_ratio_list)]

    return param_tuple_sorted_list[weighted_average_sum_list.index(min(weighted_average_sum_list))]


def check_parameters_bounds(header, good_residual_likelihood):
    # return True if check not used
    return True
    b1_init_par = 0.036
    d1_init_par = 0.022
    r1_init_par = 2.095
    r2_init_par = 3.210
    g1_init_par = 0.121
    g2_init_par = 0.064
    g3_init_par = 0.052
    t1_init_par = 34.768

    b1_par = float(good_residual_likelihood[header.index('b1')])
    d1_par = float(good_residual_likelihood[header.index('d1')])
    r1_par = float(good_residual_likelihood[header.index('r1')])
    r2_par = float(good_residual_likelihood[header.index('r2')])
    g1_par = float(good_residual_likelihood[header.index('g1')])
    g2_par = float(good_residual_likelihood[header.index('g2')])
    g3_par = float(good_residual_likelihood[header.index('g3')])
    t1_par = float(good_residual_likelihood[header.index('t1')])

    if b1_init_par/4 <= b1_par <= b1_init_par*4 and \
            d1_init_par/15 <= d1_par <= d1_init_par*15 and \
            r1_init_par/40 <= r1_par <= r1_init_par*40 and \
            r2_init_par/40 <= r2_par <= r2_init_par*40 and \
            g1_init_par/40 <= g1_par <= g1_init_par*40 and \
            g2_init_par/40 <= g2_par <= g2_init_par*40 and \
            g3_init_par/40 <= g3_par <= g3_init_par*40 and \
            t1_init_par/40 <= t1_par <= t1_init_par*40:
        return True
    else:
        return False


def get_good_residual_likelihood_list(wanted_num_tuples, only_residual, only_likelihood, residual_coef, par_tables_list):
    par_table = None
    spec_features = ['mcmc_step', 'glob_iter', 'residual', 'treeLikelihood']
    for i, par_table_file_name in enumerate(par_tables_list):
        global_iteration_num = int(
            re.search(r'_([0-9]+)', os.path.basename(par_table_file_name)).groups()[0])
        add_features_name_val_list = [['glob_iter', str(global_iteration_num)]]
        if only_residual:
            add_features_name_val_list.append(['treeLikelihood', '1.0'])
        if only_likelihood:
            add_features_name_val_list.append(['residual', '1.0'])
        if i == 0:
            par_table = ParTable(par_table_file_name,
                                 spec_features, add_features_name_val_list)
        else:
            par_table.append(par_table_file_name, add_features_name_val_list)

    min_residual = min([float(cur_tuple[par_table.header.index('residual')])
                        for cur_tuple in par_table.param_tuple_list])
    max_treeLikelihood = max([float(cur_tuple[par_table.header.index(
        'treeLikelihood')]) for cur_tuple in par_table.param_tuple_list])

    wanted_num_tuples_multiplier = 10
    while len(par_table.param_tuple_list)/2 < int(wanted_num_tuples) * wanted_num_tuples_multiplier \
            and wanted_num_tuples_multiplier > 0:
        wanted_num_tuples_multiplier -= 1

    good_residual_likelihood_list = []
    for i in range(int(wanted_num_tuples) * wanted_num_tuples_multiplier):
        writed = False
        while writed == False:
            good_residual_likelihood = find_good_residual_likelihood(par_table, residual_coef)
            if check_parameters_bounds(par_table.header, good_residual_likelihood):
                good_residual_likelihood_list.append(good_residual_likelihood)
                writed = True
            par_table.param_tuple_list.remove(good_residual_likelihood)
            if len(par_table.param_tuple_list) <= 0:
                raise ValueError('Not found paramters which stay in bounds.')
    return par_table, good_residual_likelihood_list, min_residual, max_treeLikelihood


def calc_sum_variability_all_params(par_table, good_residual_likelihood_list):
    sum_variability_all_params = 0.0
    for par_name in par_table.header:
        if par_name not in par_table.special_features_cols_idx_dict:
            ind_cur_par = par_table.header.index(par_name)
            cur_par_values = list(map(float, [good_residual_likelihood[ind_cur_par]
                                              for good_residual_likelihood in good_residual_likelihood_list]))
            variability_cur_par = numpy.std(
                cur_par_values) / numpy.average(cur_par_values)
            sum_variability_all_params = sum_variability_all_params + variability_cur_par
    return sum_variability_all_params


def select_better_variability_param_set(par_table, good_residual_likelihood_list, wanted_num_tuples, only_likelihood):
    import random
    final_good_residual_likelihood_list = []
    n_list = []
    # Get random params sets
    while len(final_good_residual_likelihood_list) != int(wanted_num_tuples):
        rand_n = random.randint(0, len(good_residual_likelihood_list) - 1)
        if rand_n not in n_list:
            n_list.append(rand_n)
            final_good_residual_likelihood_list.append(
                good_residual_likelihood_list[n_list[-1]])

    # Select better variability for ea
    for tuple_i in range(int(wanted_num_tuples)):
        variabilities_list_for_tuple_i = []
        for j in range(len(good_residual_likelihood_list)):
            final_good_residual_likelihood_list = \
                final_good_residual_likelihood_list[:tuple_i] + \
                [good_residual_likelihood_list[j]] + \
                final_good_residual_likelihood_list[tuple_i + 1:]
            variabilities_list_for_tuple_i.append(
                calc_sum_variability_all_params(par_table, final_good_residual_likelihood_list))
        index_of_max_variability = variabilities_list_for_tuple_i.index(
            max(variabilities_list_for_tuple_i))
        final_good_residual_likelihood_list = \
            final_good_residual_likelihood_list[:tuple_i] + \
            [good_residual_likelihood_list[index_of_max_variability]] + \
            final_good_residual_likelihood_list[tuple_i + 1:]
        good_residual_likelihood_list.pop(index_of_max_variability)
    if not only_likelihood:
        final_good_residual_likelihood_list = \
            sorted(final_good_residual_likelihood_list, key=lambda x: (
                float(x[par_table.header.index('residual')])), reverse=False)
    else:
        final_good_residual_likelihood_list = \
            sorted(final_good_residual_likelihood_list, key=lambda x: (
                float(x[par_table.header.index('treeLikelihood')])), reverse=True)
    return final_good_residual_likelihood_list


def print_good_residual_likelihood(par_table,
                                   good_residual_likelihood_list,
                                   min_residual,
                                   max_treeLikelihood,
                                   only_residual,
                                   only_likelihood,
                                   par_init_values_dict = None):
    if par_init_values_dict != None:
        good_residual_likelihood_list = [list(good_residual_likelihood) for good_residual_likelihood in good_residual_likelihood_list]
        par_names_list = list(set([par_name.split('_')[0] for par_name in list(par_init_values_dict.keys())]))
        new_header = list(par_table.header)
        for par_name in par_names_list:
            if par_init_values_dict[par_name + '_min'] == par_init_values_dict[par_name + '_max']:
                new_header.insert(1, par_name)
                for good_residual_likelihood in good_residual_likelihood_list:
                    good_residual_likelihood.insert(1, str(par_init_values_dict[par_name + '_min']))
        par_table.header = tuple(new_header)
    if only_residual:
        with open('good_residual_result.txt', 'w') as out_file:
            new_header = par_table.header
            new_header = [
                col_name for col_name in par_table.header if col_name != 'treeLikelihood']
            out_file.write('\t'.join(new_header) + '\t' +
                           '\t'.join(['min_residual']) + '\n')
            for good_residual_likelihood in good_residual_likelihood_list:
                new_good_residual_likelihood = \
                    [col_val for i, col_val in enumerate(good_residual_likelihood)
                        if i != par_table.header.index('treeLikelihood')]
                out_file.write('\t'.join(new_good_residual_likelihood) + '\t' +
                               '\t'.join([str(min_residual)]) + '\n')
    elif only_likelihood:
        with open('good_likelihood_result.txt', 'w') as out_file:
            new_header = par_table.header
            new_header = [
                col_name for col_name in par_table.header if col_name != 'residual']
            out_file.write('\t'.join(new_header) + '\t' +
                           '\t'.join(['max_treeLikelihood']) + '\n')
            for good_residual_likelihood in good_residual_likelihood_list:
                new_good_residual_likelihood = \
                    [col_val for i, col_val in enumerate(good_residual_likelihood)
                        if i != par_table.header.index('residual')]
                out_file.write('\t'.join(new_good_residual_likelihood) + '\t' +
                               '\t'.join([str(min_residual)]) + '\n')
    else:
        with open('good_residual_likelihood_result.txt', 'w') as out_file:
            out_file.write('\t'.join(par_table.header) + '\t' +
                           '\t'.join(['min_residual', 'max_treeLikelihood']) + '\n')
            for good_residual_likelihood in good_residual_likelihood_list:
                out_file.write('\t'.join(good_residual_likelihood) + '\t' +
                               '\t'.join([str(min_residual), str(max_treeLikelihood)]) + '\n')


def main():
    args = parser.parse_args()
    only_residual = args.only_residual
    only_likelihood = args.only_likelihood

    if args.only_residual:
        only_likelihood = False
    if args.only_likelihood:
        only_residual = False

    if args.residual_coef:
        residual_coef = float(args.residual_coef)
    else:
        residual_coef = 0.8

    par_table, good_residual_likelihood_list, min_residual, max_treeLikelihood = \
        get_good_residual_likelihood_list(args.num_tuples, 
                                            only_residual, 
                                            only_likelihood, 
                                            residual_coef,
                                            args.par_tables.split(','))
    final_good_residual_likelihood_list = \
        select_better_variability_param_set(
            par_table, good_residual_likelihood_list, args.num_tuples, only_likelihood)

    if args.equations != None:
        equations_dict, repl_dict = phy_dyn_gen_xml.read_equations_file(args.equations)
        par_init_values_dict, var_init_values_dict = phy_dyn_gen_xml.read_init_values_file(args.equations, repl_dict)
        print_good_residual_likelihood(par_table,
                                    final_good_residual_likelihood_list,
                                    min_residual,
                                    max_treeLikelihood,
                                    only_residual,
                                    only_likelihood,
                                    par_init_values_dict = par_init_values_dict)
    else:
        print_good_residual_likelihood(par_table,
                                final_good_residual_likelihood_list,
                                min_residual,
                                max_treeLikelihood,
                                only_residual,
                                only_likelihood)

if __name__ == "__main__":
    main()
