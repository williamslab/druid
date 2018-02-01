import subprocess
import sys
import Bio
import itertools
import networkx as nx
import multiprocessing
import random
import copy
import numpy as np
from joblib import Parallel,delayed
import multiprocessing
import glob, os
from DRUID_functions import *
from DRUID_graph_interaction import *
import argparse
import os.path


parser=argparse.ArgumentParser(
    description='''DRUID v0.9.1a -- A pairwise relatedness estimator. 1/3/2017.''',
    epilog="""Output has extension *.DRUID and contains columns [ind1, ind2, estimated shared IBD proportion (for debugging purposes -- will be removed), DRUID's inferred degree of relatedness, Refined IBD's inferred degree of relatedness (for debugging purposes -- will be removed """)
parser.add_argument('-o', type=str, nargs=1, default=['out'], help='Output file prefix', metavar='out')
parser.add_argument('-i', type=str, nargs=1, required=True, help='Pairwise IBD1 & IBD2 proportions file', metavar='file.ibd12')
parser.add_argument('-s', type=str, nargs=1, required=True, help='Pairwise IBD segments file', metavar='file.seg')
parser.add_argument('-m', type=str, nargs=1, required=True, help='Map file (PLINK format), non-zero cM column required', metavar='file.map')
parser.add_argument('-f', type=str, nargs=1, default=['NA'], help="Known first (FS/P/C) and second degree relationships (AU/NN/GP/GC/HS), columns: ind1/ind2/ind1's relation to ind2",metavar='file.faminfo')
#parser.add_argument('-F', type=str, nargs=1, default=['0'], help='Force DRUID to only use relationships input with -f', metavar='0/1')
parser.add_argument('-u', type=str, nargs=1, default=['NA'], help='File containing individuals to include', metavar='file.inds')
parser.add_argument('-C', type=int, nargs=1, default=[0], help='Whether to run DRUID in normal mode (0) or conservative mode (1); default is 0', metavar='0/1')

#args=parser.parse_args()
#args=parser.parse_args(['-i','../testfiles/safs_fam40.kin','-m','../testfiles/safs_filter2_geno0.02_mind0.1_ALL_fixed.map','-s','../testfiles/safs_fam40.seg','-u','../testfiles/safs_fam40.inds','-C','0'])#,'-f','../testfiles/fam2_SAFSex.rel2'])
args=parser.parse_args(['-i','../testfiles/safs_fam40.kin','-m','../testfiles/safs_filter2_geno0.02_mind0.1_ALL_fixed.map','-s','../testfiles/safs_fam40_small.seg','-u','../testfiles/safs_fam40_small.inds','-C','0', '-f', '../testfiles/safs_fam40_small.rel'])#,'-f','../testfiles/fam2_SAFSex.rel2'])

inds = []

print("Using IBD12 file: "+args.i[0])
print("Using map file: "+args.m[0])
print("Using output prefix: "+args.o[0])
if args.f[0] != 'NA':
    print("Including family info from "+args.f[0])
if args.u[0] != 'NA':
    print("Including inds from "+args.u[0])

# Get map info
global total_genome, chrom_starts, chrom_ends
[total_genome, chrom_starts, chrom_ends] = getChrInfo(args.m[0])
print("Genome length: " + str(total_genome))

# Get IBD1/2 info
[all_rel,inds,first,second] = getAllRel(args.i[0],args.u[0])
print("Total number of individuals: " + str(len(inds)))

#make graph
rel_graph = nx.DiGraph()
if args.f[0] != 'NA':
    print("Reading in family info")
    faminfo = getFamInfo(args.f[0], inds)
    forceFamInfo(rel_graph,inds, faminfo)

    print("Inferring first degree relatives")
    inferFirstFaminfo(rel_graph, all_rel, first, second, int(args.C[0]))

    # infer second degree & aunts/uncles of sibling sets
    print("Inferring second degree relatives")
    inferSecondPath(rel_graph, all_rel, second, args.s[0], args.o[0], int(args.C[0]))
else:
    # infer and add siblings, parents; other first degrees are labeled as '1'
    print("Inferring first degree relatives")
    inferFirst(rel_graph, all_rel, first, second, int(args.C[0]))

    # infer second degree & aunts/uncles of sibling sets
    print("Inferring second degree relatives")

    inferSecondPath(rel_graph, all_rel, second, args.s[0], args.o[0], int(args.C[0]))



all_results = runDRUID(rel_graph, all_rel, inds, args)



outfile_results = open(args.o[0]+'.DRUID','w')
outfile_results.write("ind1\tind2\tDRUID\tRefinedIBD\tMethod\n")
for res in all_results:
    if res[2] == '1U':
        res[2] = '1'
    outfile_results.write("\t".join(map(str,res))+'\n')


outfile_results.close()


