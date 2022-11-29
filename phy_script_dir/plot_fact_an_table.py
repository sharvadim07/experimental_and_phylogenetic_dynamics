import matplotlib
matplotlib.use('Agg')
import argparse
import re
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Modify factorial analysis results.')
parser.add_argument('-l','--fact_an_res_list', type=str,                   
                    help='Factorial analysis results tables names list to merge in txt file.', required=True)

args = parser.parse_args()


def print_df(df, pdf, max_x_pars=16, font_size=12, xlabel = 'Parameter'):
    df = df.drop(columns=['DepVariable'])
    df = df.transpose()
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df_dict = {}
    for col in df.columns:
        df_dict[col] = list(map(float, list(df[col].str.split(' '))[0])) 
    
    val_list = []
    for par_name in df_dict:
        val_list = val_list + df_dict[par_name]
    min_v = min(val_list)
    max_v = max(val_list)
    bot_lim = -0.05
    top_lim = 0.05
    if min_v < bot_lim:
        bot_lim =  min_v
    if max_v > top_lim:
        top_lim = max_v
    df_dict_keys_list = list(df_dict.keys())
    for i in range(0, len(df_dict_keys_list), max_x_pars):
        if i + max_x_pars >= len(df_dict_keys_list):
            end = len(df_dict_keys_list)
        else:
            end = i + max_x_pars
        df_dict_keys_ranged_list = df_dict_keys_list[i:end]
        
        plt.rcParams.update({'font.size': font_size})
        plt.rcParams['figure.figsize'] = 16, 10
        plt.title('Sensitivity of parameters for the characteristic TotalDeviation')
        plt.xlabel(xlabel)
        plt.ylabel('Regression estimation')
        
        plt.ylim(bot_lim, top_lim)   
        plt.boxplot(x=[df_dict[par_name] for par_name in df_dict if par_name in df_dict_keys_ranged_list], \
         labels=[par_name for par_name in df_dict if par_name in df_dict_keys_ranged_list])
        pdf.savefig(dpi = 600)
        plt.close()


def print_fact_an_tables(fact_an_res_list_file_name):
    #tabl_name_num_good_res = ''
    with open(fact_an_res_list_file_name, 'r') as fact_an_res_list_file:
        for fact_an_res_list_file in fact_an_res_list_file:
            fact_an_file_name = fact_an_res_list_file.strip()
            with open(fact_an_file_name, 'r') as fact_an_file:
                #reres = re.search(r'(p[0-9]+_p[0-9]+)_([0-9]+)', fact_an_file_name)
                #tabl_name_pass = reres.groups()[0]
                #tabl_name_num_good_res = reres.groups()[1]
                fact_an_tabl = pd.read_csv(fact_an_file, sep='\t', header = 0)
                df1 = fact_an_tabl
                df1 = df1.drop(columns=['Avg', 'StdDev', 'z-score'])
                df1 = df1[df1['DepVariable'].str.contains("TotalDeviation")]
                
                df2 = df1[df1['RegressionTerm'].str.match('^(.(?!(:)))*$')]
                df2 = df2[df2['RegressionTerm'].str.match('^(.(?!(sq)))*$')]

                df3 = df1[df1['RegressionTerm'].str.match('.*:.*')]
                df3 = df3[df3['RegressionTerm'].str.match('^(.(?!(sq)))*$')]

                df4 = df1[df1['RegressionTerm'].str.match('.*sq.*')]
                df4 = df4[df4['RegressionTerm'].str.match('^(.(?!(:)))*$')]

                df5 = df1[df1['RegressionTerm'].str.match('.*sq.*')]
                df5 = df5[df5['RegressionTerm'].str.match('.*:.*')]
                with PdfPages(os.path.basename(fact_an_file_name).split('.')[0] + '.pdf', ) as pdf:
                    print_df(df2, pdf)
                    print_df(df4, pdf, xlabel='Square Parameter')
                    print_df(df3, pdf, xlabel='Parameter Intaractions')                    
                    print_df(df5, pdf, xlabel='Square Parameter Intaractions', font_size = 10)

print_fact_an_tables(args.fact_an_res_list)