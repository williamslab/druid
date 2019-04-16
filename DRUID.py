#!/usr/bin/env python

import networkx as nx
from DRUID_functions import *
from DRUID_graph_interaction import *
import argparse

version='v1.02.6'
update='16 Apr 2019'

#import cProfile, pstats
#pr = cProfile.Profile()
#pr.enable()


parser=argparse.ArgumentParser(
    description="DRUID " + version + " -- A multiway relatedness estimator. Last update " + update,
    epilog="""Outputs file with extension *.DRUID and contains columns [ind1, ind2, DRUID's inferred degree of relatedness, input pairwise IBD1/2 proportions file's inferred degree of relatedness.
     With flag -F, can also output PLINK format .fam files containing first and second degree relationships types which were inferred/provided, outputting multiple .fam files in cases when a parent-child pair is found, but there is not
     enough information to determine who is the parent and who is the child""")
parser.add_argument('-o', type=str, nargs=1, default=['out'], help='Output file prefix', metavar='out')
parser.add_argument('-i', type=str, nargs=1, required=True, help='Pairwise IBD1 & IBD2 proportions file', metavar='file.ibd12')
parser.add_argument('-s', type=str, nargs=1, required=True, help='Pairwise IBD segments file', metavar='file.seg')
parser.add_argument('-m', type=str, nargs=1, required=True, help='Map file (PLINK format), non-zero cM column required', metavar='file.map')
parser.add_argument('-f', type=str, nargs=1, default=[''], help="Known first (FS/P/C) and second degree relationships (AU/NN/GP/GC/HS), columns: ind1/ind2/ind1's relation to ind2",metavar='file.faminfo')
parser.add_argument('-u', type=str, nargs=1, default=[''], help='File containing individuals to include', metavar='file.inds')
parser.add_argument('-C', type=int, nargs=1, default=[0], help='Whether to run DRUID in normal mode (0) or conservative mode (1); default is 0', metavar='0/1')
parser.add_argument('-F', type=int, nargs=1, default=[0], help='Whether to output fam files (PLINK format) containing inferred/provided relationship types; default is 0', metavar='0/1')

args=parser.parse_args()

inds = []

print("DRUID " + version + " -- A multiway relatedness estimator.")
print("Last update " + update + "\n")

print("Using IBD12 file: "+args.i[0])
print("Using map file: "+args.m[0])
print("Using output prefix: "+args.o[0])
if args.f[0] != '':
    print("Including family info from "+args.f[0])
if args.u[0] != '':
    print("Including inds from "+args.u[0])
if args.F[0]:
    print("Print fam files = TRUE")
else:
    print("Print fam files = FALSE")

# Get map info
global total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs
[total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs] = getChrInfo(args.m[0])
print("Genome length: " + str(total_genome)+'\n')

# Get IBD1/2 info
[all_rel,inds,first,second,third] = getAllRel(args.i[0],args.u[0])
print("Total number of individuals: " + str(len(inds)))

#make graph
rel_graph = nx.DiGraph()
rel_graph_tmp = nx.DiGraph()
if args.f[0] != '':
    print("Reading in family info")
    faminfo = getFamInfo(args.f[0], inds)
    forceFamInfo(rel_graph_tmp, faminfo) #store relationship information in rel_graph_tmp

print("\nInferring first degree relatives")
inferFirst(rel_graph, rel_graph_tmp, all_rel, first, second, int(args.C[0]))

all_segs = readSegments(args.s[0])

# infer second degree & aunts/uncles of sibling sets
print("\nInferring second degree relatives")
inferSecondPath(rel_graph, rel_graph_tmp, all_rel, second, all_segs, args.o[0], int(args.C[0]))

outfile_results = open(args.o[0]+'.DRUID','w')
outfile_results.write("#ind1\tind2\tDRUID\tRefinedIBD\tMethod\n")

print("\n\nRunning main DRUID inference")
runDRUID(rel_graph, all_rel, inds, all_segs, args, outfile_results)

outfile_results.close()


if args.F[0] == 1:
    print("\nPrinting .fam files")
    fillInGraph(rel_graph)

print("\ndone")

#pr.disable()
#ps = pstats.Stats(pr).sort_stats('cumulative')
#ps.print_stats()
