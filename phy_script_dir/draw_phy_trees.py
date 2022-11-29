import matplotlib
matplotlib.use('Agg')
import argparse
from Bio import Phylo
import pylab

parser = argparse.ArgumentParser(description='Drawing specified by user trees for NEXUS tree format file.')
parser.add_argument('-t','--phytree', type=str,                   
                    help='Phylogenetic Tree file in nexus format.', required=True)
parser.add_argument('-l','--index_list', type=str,                   
                    help='Comma separated list of indexes of trees for drawing.', required=False)

args = parser.parse_args()
if args.index_list != None:
        list_of_index = [ind.strip() for ind in args.index_list.strip().split(',')]

for tree in Phylo.parse(args.phytree,'nexus'):
    if args.index_list == None or tree.name.split('_')[1] in list_of_index:
        pylab.rcParams['figure.figsize'] = 20,15
        pylab.rcParams['figure.dpi'] = 200
        Phylo.draw(tree,do_show=False)
        pylab.savefig(tree.name)
        pylab.close()

