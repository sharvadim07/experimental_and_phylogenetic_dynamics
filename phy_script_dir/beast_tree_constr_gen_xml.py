import logging
import argparse
import xml.etree.cElementTree as ET


parser = argparse.ArgumentParser(description='Generating XML for initial tree construction with BEAST.')
parser.add_argument('-m','--msa', type=str,                   
                    help='Input MSA file.', required=True)
parser.add_argument('-i','--fasta_info', type=str,                   
                    help='MSA sequences info. Demes(as in the equations) and date(year.part_of_year)', required=True)
parser.add_argument('-l','--mcmc_iter', type=int,                   
                    help='Number of mcmc iterations.', required=True)

args = parser.parse_args()

def read_fasta_names_info(fasta_names_info_file_name):
    fasta_names_info_dict = {}    
    with open(fasta_names_info_file_name, 'r') as fasta_names_info_file:
        for line in fasta_names_info_file:
            line_list = line.strip().split(',')
            seq_name = line_list[0]
            time = line_list[1]
            #Replace A -> A1 ...
            deme = line_list[2]
            fasta_names_info_dict[seq_name] = [time,deme]
    return fasta_names_info_dict

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

def add_sequences_from_alignment(tree_root, fasta_names_info_dict, fasta_sequences_dict):
    data = tree_root.find('data')
    for sequence in data.findall('sequence'):
        data.remove(sequence)
    for sequence in fasta_sequences_dict:
        seq_xml = ET.Element('sequence')
        seq_xml.attrib['id'] = 'seq_' + sequence
        seq_xml.attrib['taxon'] = sequence
        seq_xml.attrib['totalcount'] = '4'
        seq_xml.attrib['value'] = fasta_sequences_dict[sequence]
        seq_xml.tail = '\n'
        data.append(seq_xml)
    #ET.dump(data)
    trait = tree_root.find('run').find('state').find('tree').find('trait')
    trait.text = ''    
    for i, sequence in enumerate(fasta_names_info_dict):        
        trait.text += sequence + \
                      '=' + \
                      fasta_names_info_dict[sequence][0]  #fasta_names_info_dict[0] it is time
        if i < len(fasta_names_info_dict) - 1:
            trait.text += ',\n'

def change_mcmc_num_iterations(tree_root, mcmc_iter):
    run = tree_root.find('run')
    run.attrib['chainLength'] = str(mcmc_iter)


def editXML(filename, 
            fasta_names_info_dict, 
            fasta_sequnces_dict, 
            mcmc_iter):
    tree = ET.ElementTree(file=filename)
    tree_root = tree.getroot()
    add_sequences_from_alignment(tree_root, fasta_names_info_dict, fasta_sequnces_dict)    
    change_mcmc_num_iterations(tree_root, mcmc_iter)    
    new_tree = ET.ElementTree(tree_root)
    #with open("pattern_model_seq_phy_dyn.xml_updated.xml", "w") as f:
    new_tree.write("gen_tree_coalescent_const.xml")

fasta_names_info_dict = read_fasta_names_info(args.fasta_info)
fasta_sequnces_dict = read_fasta_msa_file(args.msa)

import os
editXML(os.path.dirname(os.path.realpath(__file__)) + "/pattern_gen_tree_coalescent_const.xml", 
                        fasta_names_info_dict, 
                        fasta_sequnces_dict,                         
                        args.mcmc_iter)