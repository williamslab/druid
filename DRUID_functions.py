import subprocess
import sys
import Bio
import itertools
import networkx as nx
import multiprocessing
import random
import copy
import numpy as np
import os
from DRUID_graph_interaction import *


global total_genome

degrees = {'MZ': 1/2.0**(3.0/2), 1: 1/2.0**(5.0/2), 2: 1/2.0**(7.0/2), 3: 1/2.0**(9.0/2), 4: 1/2.0**(11.0/2), 5: 1/2.0**(13.0/2), 6: 1/2.0**(15.0/2), 7: 1/2.0**(17.0/2), 8: 1/2.0**(19.0/2), 9: 1/2.0**(21.0/2), 10: 1/2.0**(23.0/2), 11: 1/2.0**(25.0/2), 12: 1/2.0**(27/2.0), 13: 1/2.0**(29.0/2)}  # threshold values for each degree of relatedness

def forceFamInfo(rel_graph,inds, faminfo):
    #force the provided faminfo file information into rel_graph
    for i1 in faminfo.keys():
        for i2 in faminfo[i1].keys():
            rel_graph.add_edge(i1,i2)
            rel_graph[i1][i2]['type'] = faminfo[i1][i2]


def inferFirst(rel_graph,all_rel, first, second, C):
    #build graphs using first degree relative inferences
    for [ind1,ind2] in first+second: #iterate through currently inferred first and second degree pairs
        if ind1 < ind2:
            i1 = ind1
            i2 = ind2
        else:
            i1 = ind2
            i2 = ind1
        if i1 in all_rel.keys() and i2 in all_rel[i1].keys():
            # all_rel: IBD1, IBD2, K, D
            if float(all_rel[i1][i2][1]) > 1/2.0**(5/2.0): #if IBD2 meets minimum threshold
                rel_graph.add_edge(i1,i2)
                rel_graph.add_edge(i2,i1)
                if float(all_rel[i1][i2][2]) < 1/2.0**(3/2.0):
                    rel_graph[i1][i2]['type'] = 'FS'
                    rel_graph[i2][i1]['type'] = 'FS'
                else:
                    rel_graph[i1][i2]['type'] = 'T' #twin
                    rel_graph[i2][i1]['type'] = 'T'

            elif float(all_rel[i1][i2][1]) > 1/2.0**(7/2.0):
                rel_graph.add_edge(i1,i2)
                rel_graph.add_edge(i2,i1)
                if not C and float(all_rel[i1][i2][1]) > 1/2.0**(11/4.0):
                    print("Warning: "+i1+' and '+i2+' have low levels of IBD2 for siblings, may be 3/4 sibs')
                    rel_graph[i1][i2]['type'] = 'FS'
                    rel_graph[i2][i1]['type'] = 'FS'
                elif C and float(all_rel[i1][i2][1]) > 1/2.0**(5/2.0):
                    rel_graph[i1][i2]['type'] = 'FS'
                    rel_graph[i2][i1]['type'] = 'FS'
                else:
                    rel_graph[i1][i2]['type'] = '1U'
                    rel_graph[i2][i1]['type'] = '1U'
            else:
                if float(all_rel[i1][i2][3]) == 1:
                    #possible parent-child pair
                    rel_graph.add_edge(i1,i2)
                    rel_graph.add_edge(i2,i1)
                    rel_graph[i1][i2]['type'] = '1U'
                    rel_graph[i2][i1]['type'] = '1U'

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            siblings = getSibsFromGraph(rel_graph,node)
            siblings.add(node)
            for sib in siblings:
                checked.add(sib)
            [add_sibs,remove] = checkSiblingSubgraph(rel_graph,siblings.copy()) #edit siblings list, return removed sibs

            for sib in add_sibs:
                checked.add(sib)

            for ind in remove:
                for sib in siblings:
                    if rel_graph.has_edge(ind,sib):
                        rel_graph[ind][sib]['type'] = '1U'


            for [ind1,ind2] in itertools.combinations(add_sibs, 2):
                if not rel_graph.has_edge(ind1,ind2):
                    rel_graph.add_edge(ind1, ind2)
                    rel_graph.add_edge(ind2, ind1)
                rel_graph[ind1][ind2]['type'] = 'FS'
                rel_graph[ind2][ind1]['type'] = 'FS'



            for ind in remove:
                if ind in siblings:
                    siblings.remove(ind)

            #now that sibling set is completed, look for parents
            #find neighbors of first sibling set labeled as '1'
            neighbors = rel_graph.neighbors(list(siblings)[0])
            inds_to_check = set()
            for n in neighbors:
                if rel_graph.get_edge_data(list(siblings)[0],n)['type'] == '1U':
                    inds_to_check.add(n)

            #check if other siblings also have this neighbor and are labeled as '1'
            par = set()
            if len(inds_to_check):
                for ind in inds_to_check:
                    if checkIfParent(rel_graph, all_rel, siblings, ind):
                        par.add(ind)

            for ind in par:
                for sib in siblings:
                    if not rel_graph.has_edge(ind,sib):
                        rel_graph.add_edge(ind, sib)
                        rel_graph.add_edge(sib, ind)
                    rel_graph[ind][sib]['type'] = 'P'
                    rel_graph[sib][ind]['type'] = 'C'
                    #edge_graph[sib][ind]['type'] = 'C'





def inferFirstFaminfo(rel_graph, all_rel, first, second, C):
    #build graphs using first degree relative inferences as well as provided faminfo file

    rel_graph_tmp = nx.DiGraph() #inferred, to be compared to input relationships
    for [ind1, ind2] in first + second:  # iterate through currently inferred or provided first and second degree pairs
        if ind1 < ind2:
            i1 = ind1
            i2 = ind2
        else:
            i1 = ind2
            i2 = ind1

        if i1 in all_rel.keys() and i2 in all_rel[i1].keys():
            if float(all_rel[i1][i2][1]) > 1 / 2.0 ** (5 / 2.0):  #check if IBD2 high enough to be siblings
                rel_graph_tmp.add_edge(i1, i2)
                rel_graph_tmp.add_edge(i2, i1)
                if float(all_rel[i1][i2][2]) < 1 / 2.0 ** (3 / 2.0):
                    rel_graph_tmp[i1][i2]['type'] = 'FS'
                    rel_graph_tmp[i2][i1]['type'] = 'FS'
                else:
                    rel_graph_tmp[i1][i2]['type'] = 'T'  # twin
                    rel_graph_tmp[i2][i1]['type'] = 'T'
            elif float(all_rel[i1][i2][1]) > 1 / 2.0 ** (7 / 2.0): #check if IBD2 still fairly high
                rel_graph_tmp.add_edge(i1, i2)
                rel_graph_tmp.add_edge(i2, i1)
                if not C and float(all_rel[i1][i2][1]) > 1 / 2.0 ** (11 / 4.0):
                    print("Warning: " + i1 + ' and ' + i2 + ' have low levels of IBD2 for siblings, may be 3/4 sibs')
                    rel_graph_tmp[i1][i2]['type'] = 'FS'
                    rel_graph_tmp[i2][i1]['type'] = 'FS'
                elif C and float(all_rel[i1][i2][1]) > 1 / 2.0 ** (5 / 2.0):
                    rel_graph_tmp[i1][i2]['type'] = 'FS'
                    rel_graph_tmp[i2][i1]['type'] = 'FS'
                else:
                    if all_rel[i1][i2][3] == 1:
                        rel_graph_tmp[i1][i2]['type'] = '1U'
                        rel_graph_tmp[i2][i1]['type'] = '1U'
            elif all_rel[i1][i2][3] == 1: #if IBD2 low but individual is inferred to be first degree, save that relationship
                rel_graph_tmp.add_edge(i1, i2)
                rel_graph_tmp.add_edge(i2, i1)
                rel_graph_tmp[i1][i2]['type'] = '1U'
                rel_graph_tmp[i2][i1]['type'] = '1U'


    #ensure our sibling subgraph is adequately connected
    checked = set()
    for node in rel_graph_tmp.nodes():
        if not node in checked:
            siblings = getSibsFromGraph(rel_graph_tmp,node)
            siblings.append(node)
            for sib in siblings:
                checked.add(sib)
            [add_sibs,remove] = checkSiblingSubgraph(rel_graph_tmp,siblings.copy()) #edit siblings list, return removed sibs

            for sib in add_sibs:
                checked.add(sib)

            for ind in remove:
                for sib in siblings:
                    if rel_graph_tmp.has_edge(ind,sib):
                        rel_graph_tmp[ind][sib]['type'] = '1U'
                        rel_graph_tmp[sib][ind]['type'] = '1U'


            for [ind1,ind2] in itertools.combinations(add_sibs, 2):
                if not rel_graph_tmp.has_edge(ind1,ind2):
                    rel_graph_tmp.add_edge(ind1,ind2)
                    rel_graph_tmp.add_edge(ind2, ind1)
                    rel_graph_tmp[ind1][ind2]['type'] = 'FS'
                    rel_graph_tmp[ind2][ind1]['type'] = 'FS'



            for ind in remove:
                if ind in siblings:
                    siblings.remove(ind)

            #now that sibling set is completed, look for parents
            #find neighbors of first sibling set labeled as '1'
            neighbors = rel_graph_tmp.neighbors(siblings[0])
            inds_to_check = set()
            for n in neighbors:
                if rel_graph_tmp.get_edge_data(siblings[0],n)['type'] == '1U':
                    inds_to_check.add(n)

            #check if other siblings also have this neighbor and are labeled as '1'
            par = set()
            if len(inds_to_check):
                for ind in inds_to_check:
                    if checkIfParent(rel_graph_tmp, all_rel, siblings, ind):
                        par.add(ind)

            for ind in par:
                for sib in siblings:
                    if not rel_graph_tmp.has_edge(ind,sib):
                        rel_graph_tmp.add_edge(ind, sib)
                        rel_graph_tmp.add_edge(sib, ind)
                    rel_graph_tmp[ind][sib]['type'] = 'P'
                    rel_graph_tmp[sib][ind]['type'] = 'C'


    #compare inferred graph to provided graph
    for edge in rel_graph.edges():
        if not edge in rel_graph_tmp.edges() or not rel_graph.get_edge_data(edge[0],edge[1])['type'] == rel_graph_tmp.get_edge_data(edge[0],edge[1])['type']:
            print("Warning: Unable to confirm "+edge[0]+" and "+edge[1]+" as "+str(rel_graph.get_edge_data(edge[0],edge[1])['type'])+" but including as such")

    for edge in rel_graph_tmp.edges():
        if not edge in rel_graph.edges():
            rel_graph.add_edge(edge[0],edge[1])
            rel_graph[edge[0]][edge[1]]['type'] = rel_graph_tmp.get_edge_data(edge[0],edge[1])['type']




def inferSecondPath(rel_graph, all_rel, second, file_for_segments, outfile, C):
    # infer and add 2nd degree relationships
    for [ind1, ind2] in second:
        #if not rel_graph.has_edge(ind1,ind2) or not rel_graph.get_edge_data(ind1,ind2)['type'] in ['NN','AU']:
        if not rel_graph.has_edge(ind1, ind2): #and not rel_graph.get_edge_data(ind1,ind2)['type'] in ['NN','AU']):
            if ind1 < ind2:
                i1 = ind1
                i2 = ind2
            else:
                i1 = ind2
                i2 = ind1
            if i1 in all_rel.keys() and i2 in all_rel[i1].keys():
                if float(all_rel[i1][i2][1]) < 0.01: #proportion IBD2 less than 0.01
                    rel_graph.add_edge(i1,i2)
                    rel_graph.add_edge(i2,i1)
                    rel_graph[i1][i2]['type'] = '2'
                    rel_graph[i2][i1]['type'] = '2'
                elif float(all_rel[i1][i2][1]) > 1/2.0**(7/2.0): #high IBD2 even though second degree
                    sib1 = getSibsFromGraph(rel_graph,i1)
                    sib2 = getSibsFromGraph(rel_graph,i2)
                    sib1.append(i1)
                    sib2.append(i2)
                    for s1 in sib1:
                        for s2 in sib2:
                            if not rel_graph.get_edge_data(s1,s2)['type'] == 'FS':
                                if not rel_graph.has_edge(s1,s2):
                                    rel_graph.add_edge(s1, s2)
                                    rel_graph.add_edge(s2, s1)
                                rel_graph[s1][s2]['type'] = '1U'
                                rel_graph[s2][s1]['type'] = '1U'
                else: #if IBD2 is too high, label as DC so we don't use
                    rel_graph.add_edge(i1,i2)
                    rel_graph.add_edge(i2,i1)
                    rel_graph[i1][i2]['type'] = 'DC'
                    rel_graph[i2][i1]['type'] = 'DC'

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [sibs, halfsibs] = getSibsAndHalfSibsFromGraph(rel_graph, node)
            sibs.add(node)

            if len(sibs) > 1:
                if not C:
                    second = getSecondDegreeRelativesFromAllRel_lessRestrictive(all_rel, inds, sibs)
                else:
                    second = getSecondDegreeRelativesFromAllRel(all_rel, inds, sibs)
                print("Checking for aunts/uncles of "+node+" and his/her siblings")
                [avunc, avunc_hs1, avunc_hs2] = getAuntsUncles_IBD011_nonoverlapping_pairs(all_rel, sibs, halfsibs, second, file_for_segments, rel_graph)
                print("Done")

                # add the inferred avuncular relationships to graph
                for av in avunc:
                    for sib in sibs:
                        if not rel_graph.has_edge(sib,av): #if it doesn't have this edge, add it and label as AV
                            rel_graph.add_edge(sib, av)
                            rel_graph.add_edge(av, sib)
                            rel_graph[sib][av]['type'] = 'NN'
                            rel_graph[av][sib]['type'] = 'AU'
                        elif not rel_graph[sib][av]['type'] == 'NN': #if it does have this edge but it isn't labeled as AV, label as AV
                            rel_graph[sib][av]['type'] = 'NN'
                            rel_graph[av][sib]['type'] = 'AU'

                if len(avunc_hs1):
                    for av in avunc_hs1:
                        for sib in sibs+halfsibs[0]:
                            if not rel_graph.has_edge(sib,av): #if it doesn't have this edge, add it and label as AV
                                rel_graph.add_edge(sib, av)
                                rel_graph.add_edge(av, sib)
                                rel_graph[sib][av]['type'] = 'NN'
                                rel_graph[av][sib]['type'] = 'AU'
                            elif not rel_graph[sib][av]['type'] == 'NN': #if it does have this edge but it isn't labeled as AV, label as AV
                                rel_graph[sib][av]['type'] = 'NN'
                                rel_graph[av][sib]['type'] = 'AU'

                if len(avunc_hs2):
                    for av in avunc_hs2:
                        for sib in sibs+halfsibs[1]:
                            if not rel_graph.has_edge(sib,av): #if it doesn't have this edge, add it and label as AV
                                rel_graph.add_edge(sib, av)
                                rel_graph.add_edge(av, sib)
                                rel_graph[sib][av]['type'] = 'NN'
                                rel_graph[av][sib]['type'] = 'AU'
                            elif not rel_graph[sib][av]['type'] == 'NN': #if it does have this edge but it isn't labeled as AV, label as AV
                                rel_graph[sib][av]['type'] = 'NN'
                                rel_graph[av][sib]['type'] = 'AU'


            for sib in sibs:
                checked.add(sib)

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [siblings, avunc_bothsides, nn, par, child, gp, halfsib_sets, twins] = pullFamily(rel_graph, node)
            siblings.add(node)
            checkAuntUncleGPRelationships(rel_graph, siblings, par)
            checked = checked.union(siblings)


def getFamInfo(famfile, inds):
    #read in faminfo file
    global faminfo
    faminfo = {}
    file = open(famfile,'r')
    for line in file:
        l = str.split(line.rstrip())
        if l[0] in inds and l[1] in inds:
            if not l[0] in faminfo.keys():
                faminfo[l[0]] = {}
            if not l[1] in faminfo.keys():
                faminfo[l[1]] = {}
            faminfo[l[0]][l[1]] = l[2]
            if l[2] == 'FS' or l[2] == 'HS':
                faminfo[l[1]][l[0]] = l[2]
            elif l[2] == 'P':
                faminfo[l[1]][l[0]] = 'C'
            elif l[2] == 'C':
                faminfo[l[1]][l[0]] = 'P'
            elif l[2] == 'AU':
                faminfo[l[1]][l[0]] = 'NN'
            elif l[2] == 'NN':
                faminfo[l[1]][l[0]] = 'AU'
            elif l[2] == 'GC':
                faminfo[l[1]][l[0]] = 'GP'
            elif l[2] == 'GP':
                faminfo[l[1]][l[0]] = 'GC'
            else:
                file.close()
                raise ValueError(str(l[2]) + ' is not an accepted relationship type (FS, P, C, AU, NN, GC, GP, HS)')

    file.close()

    return faminfo

def getChrInfo(mapfile):
    #read in information from .map file
    global chrom_starts
    global chrom_ends
    global total_genome
    chrom_starts = {}
    chrom_ends = {}
    for chr in range(1,23):
        chrom_starts[chr] = 9999999
        chrom_ends[chr] = 0

    file = open(mapfile,'r')
    for line in file:
        l = str.split(line.rstrip())
        chr = float(l[0])
        if chrom_starts[chr] > float(l[2]):
            chrom_starts[chr] = float(l[2])
        if chrom_ends[chr] < float(l[2]):
            chrom_ends[chr] = float(l[2])

    file.close()

    total_genome = 0
    for chr in range(1,23):
        total_genome = total_genome + chrom_ends[chr] - chrom_starts[chr]

    return [total_genome, chrom_starts, chrom_ends]


def getInferredFromK(K):
    # Return inferred degree of relatedness using kinship coefficient K
    if K >= degrees['MZ']:
        return -1
    if K >= degrees[1]:
        return 1
    elif K >= degrees[2]:
        return 2
    elif K >= degrees[3]:
        return 3
    elif K >= degrees[4]:
        return 4
    elif K >= degrees[5]:
        return 5
    elif K >= degrees[6]:
        return 6
    elif K >= degrees[7]:
        return 7
    elif K >= degrees[8]:
        return 8
    elif K >= degrees[9]:
        return 9
    elif K >= degrees[10]:
        return 10
    elif K >= degrees[11]:
        return 11
    else:
        return 0


def getIBDsegments(ind1, ind2, file_for_segments):
    # get IBD segments between ind1 and ind2, sorting segments by IBD2, IBD1, and IBD0
    rand = random.randint(1, 1e10)
    filename = ind1 + '_' + ind2 + '_' + str(rand) + '.txt'
    call1 = 'grep "' + ind1 + '" ' + file_for_segments + ' | grep ' + ind2 + ' > ' + filename
    subprocess.call(call1, shell=True)

    all = {}
    IBD_file = open(filename, 'r')
    for line in IBD_file:
        l = str.split(line.rstrip())
        chr = int(l[2])
        #chr = int(str.split(l[4],"chr")[1])
        if not chr in all.keys():
            all[chr] = []
        all[chr].append([float(l[4]), float(l[5]), int(str.split(l[3],'IBD')[1])])  # collect all children segments and which haplotype they arose from
        #all[chr].append([int(l[5]), int(l[6]), ])

    IBD_file.close()
    call1 = 'rm ' + filename
    subprocess.call(call1, shell=True)

    #prepare IBD1 and IBD2 dicts
    IBD1 = {}
    IBD2 = {}
    for chr in all.keys():
        all[chr].sort()
        for seg in all[chr]:
            if seg[2] == 1:
                if not chr in IBD1.keys():
                    IBD1[chr] = []
                IBD1[chr].append(seg[0:2]) #only append start and end
            elif seg[2] == 2:
                if not chr in IBD2.keys():
                    IBD2[chr] = []
                IBD2[chr].append(seg[0:2]) #only append start and end

    IBD12 = {} #IBD12 = regions that are IBD (IBD1 or IBD2)
    for chr in IBD1.keys():
        IBD12[chr] = []
        for k in IBD1[chr]:
            IBD12[chr].append(k)

    for chr in IBD2.keys():
        if not chr in IBD12.keys():
            IBD12[chr] = []
        for k in IBD2[chr]:
            IBD12[chr].append(k)
        IBD12[chr] = mergeIntervals(IBD12[chr][:])

    IBD0 = {}
    for chr in IBD12.keys():
        IBD0[chr] = []
        if len(IBD12[chr]) > 0:
            if IBD12[chr][0][0] > chrom_starts[chr]:
                IBD0[chr].append([chrom_starts[chr], IBD12[chr][0][0] - 1])
            if len(IBD12[chr]) > 1:
                for k in range(1, len(IBD12[chr])):
                    if IBD12[chr][k - 1][1] != IBD12[chr][k][0]:
                        IBD0[chr].append([IBD12[chr][k - 1][1], IBD12[chr][k][0]])
                if IBD12[chr][k][1] < chrom_ends[chr]:
                    IBD0[chr].append([IBD12[chr][k][1], chrom_ends[chr]])
        IBD0[chr] = mergeIntervals(IBD0[chr][:])

    return [IBD0, IBD1, IBD2]  # outputs IBD0, IBD1, IBD2



def mergeIntervals(intervals):
    #given a list of intervals, merge them where they overlap
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    return merged



def combineIBDsegments(sibset):
    # for a given sibset, collect all of their pairwise IBD0/1/2 regions
    # then, greedily collect IBD0 regions, IBD1 regions, and then IBD2 regions for comparison with avuncular segments
    IBD0 = {}
    IBD1 = {}
    IBD2 = {}
    all = {}
    if len(sibset) > 1:
        for [ind1, ind2] in itertools.combinations(sibset, 2):
            tmp = getIBDsegments(ind1, ind2)
            for chr in tmp[0].keys():
                if not chr in IBD0.keys():
                    IBD0[chr] = []
                IBD0[chr].append(tmp[0][chr])

            for chr in tmp[1].keys():
                if not chr in IBD1.keys():
                    IBD1[chr] = []
                IBD1[chr].append(tmp[1][chr])

            for chr in tmp[2].keys():
                if not chr in IBD2.keys():
                    IBD2[chr] = []
                IBD2[chr].append(tmp[2][chr])

        IBD0_all = {}
        IBD0_new = {}  # as many IBD0 ranges as possible, as large as possible
        IBD2_all = {}
        IBD2_new = {}
        IBD1_all = {}
        IBD1_new = {}
        for chr in IBD0.keys():
            IBD0_all[chr] = []
            tmp = [y for x in IBD0[chr] for y in x]
            for k in tmp:
                IBD0_all[chr].append(k)
            IBD0_all[chr].sort()
            IBD0_new[chr] = mergeIntervals(IBD0_all[chr])

        for chr in IBD1.keys():
            IBD1_all[chr] = []
            tmp = [y for x in IBD1[chr] for y in x]
            for k in tmp:
                IBD1_all[chr].append(k)
            IBD1_new[chr] = mergeIntervals(IBD1_all[chr])

        for chr in IBD2.keys():
            IBD2_all[chr] = []
            tmp = [y for x in IBD2[chr] for y in x]
            for k in tmp:
                IBD2_all[chr].append(k)
            IBD2_new[chr] = mergeIntervals(IBD2_all[chr])

        all_IBD = {}
        for chr in IBD0_new.keys():
            all_IBD[chr] = []
            for k in IBD0_new[chr]:
                all_IBD[chr].append([k[0], k[1], 0])  # start with all IBD0 segments
            all_IBD[chr].sort()

        ### WORKS!! ###
        for chr in IBD1_new.keys():
            IBD1_new[chr].sort()
            if not chr in all_IBD.keys():
                all_IBD[chr] = []
            k = 0
            k2 = 0
            newadd = []
            full_prev = 0
            if len(all_IBD[chr]) > 0:
                while k < len(all_IBD[chr]) - 1 and k2 < len(IBD1_new[chr]):
                    # print(str(k)+'\t'+str(k2)+'\n')
                    if IBD1_new[chr][k2][0] <= all_IBD[chr][k][0] and IBD1_new[chr][k2][1] > all_IBD[chr][k][1]:
                        if full_prev != 1:
                            newadd.append([IBD1_new[chr][k2][0], min(IBD1_new[chr][k2][1], all_IBD[chr][k][0])])
                            full_prev = 1
                        newadd.append([all_IBD[chr][k][1], min(IBD1_new[chr][k2][1], all_IBD[chr][k + 1][0])])
                        if all_IBD[chr][k + 1][0] <= IBD1_new[chr][k2][1]:
                            k = k + 1

                    elif IBD1_new[chr][k2][0] < all_IBD[chr][k][0]:
                        newadd.append([IBD1_new[chr][k2][0], min(IBD1_new[chr][k2][1], all_IBD[chr][k][0])])
                        full_prev = 0
                        k2 = k2 + 1

                    elif IBD1_new[chr][k2][1] > all_IBD[chr][k][1]:
                        newadd.append([all_IBD[chr][k][1], min(IBD1_new[chr][k2][1], all_IBD[chr][k + 1][0])])
                        k = k + 1
                        full_prev = 1
                        if all_IBD[chr][k][1] > IBD1_new[chr][k2][1]:
                            k2 = k2 + 1

                    elif IBD1_new[chr][k2][0] >= all_IBD[chr][k][0] and IBD1_new[chr][k2][1] <= all_IBD[chr][k][1]:
                        k2 = k2 + 1
                        full_prev = 0

                    if k < len(all_IBD[chr]) - 1 and k2 < len(IBD1_new[chr]):
                        while IBD1_new[chr][k2][1] < all_IBD[chr][k + 1][0]:
                            k2 = k2 + 1
                            if k2 >= len(IBD1_new[chr]):
                                break

                while k2 < len(IBD1_new[chr]):
                    if IBD1_new[chr][k2][0] > all_IBD[chr][len(all_IBD[chr]) - 1][1]:
                        newadd.append([IBD1_new[chr][k2][0], IBD1_new[chr][k2][1]])
                    elif IBD1_new[chr][k2][1] > all_IBD[chr][len(all_IBD[chr]) - 1][1]:
                        newadd.append(
                            [max(all_IBD[chr][len(all_IBD[chr]) - 1][1], IBD1_new[chr][k2][0]), IBD1_new[chr][k2][1]])
                    k2 = k2 + 1

                for seg in newadd:
                    all_IBD[chr].append([seg[0], seg[1], 1])

            else:
                for seg in IBD1_new[chr]:
                    all_IBD[chr].append([seg[0], seg[1], 1])

            all_IBD[chr].sort()

        # Add IBD2
        for chr in range(1, 23):
            newadd = []
            if chr in all_IBD.keys():
                for k in range(0, len(all_IBD[chr]) - 1):
                    if all_IBD[chr][k][1] < all_IBD[chr][k + 1][0] - 1:
                        newadd.append([all_IBD[chr][k][1] + 1, all_IBD[chr][k + 1][0] - 1])
            else:
                all_IBD[chr] = []
                newadd.append([chrom_starts[chr], chrom_ends[chr]])

            for seg in newadd:
                all_IBD[chr].append([seg[0], seg[1], 2])

            all_IBD[chr].sort()
    else:
        all_IBD = {}
        for chr in range(1,23):
            all_IBD[chr] = []

    return all_IBD


def collectIBDsegments(sibset,file_for_segments):
    # collect pairwise IBD0,1,2 regions between all pairs of siblings
    IBD_all = {}
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        IBD0 = {}
        IBD1 = {}
        IBD2 = {}
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}

        if not ind2 in IBD_all[ind1].keys():
            IBD_all[ind1][ind2] = []

        tmp = getIBDsegments(ind1, ind2,file_for_segments)
        for chr in tmp[0].keys():
            if not chr in IBD0.keys():
                IBD0[chr] = []
            for seg in tmp[0][chr]:
                IBD0[chr].append(seg)
            IBD0[chr].sort()

        for chr in tmp[1].keys():
            if not chr in IBD1.keys():
                IBD1[chr] = []
            for seg in tmp[1][chr]:
                IBD1[chr].append(seg)
            IBD1[chr].sort()

        for chr in tmp[2].keys():
            if not chr in IBD2.keys():
                IBD2[chr] = []
            for seg in tmp[2][chr]:
                IBD2[chr].append(seg)
            IBD2[chr].sort()

        IBD_all[ind1][ind2] = [IBD0, IBD1, IBD2]

    return IBD_all


any_in = lambda a, b: any(i in b for i in a)

def collectAllIBDsegments(sibset):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD0 = {}
    IBD1 = {}
    IBD2 = {}
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        tmp = getIBDsegments(ind1, ind2)
        for chr in tmp[0].keys():
            if not chr in IBD0.keys():
                IBD0[chr] = []
            for seg in tmp[0][chr]:
                IBD0[chr].append(seg)

        for chr in tmp[1].keys():
            if not chr in IBD1.keys():
                IBD1[chr] = []
            for seg in tmp[1][chr]:
                IBD1[chr].append(seg)

        for chr in tmp[2].keys():
            if not chr in IBD2.keys():
                IBD2[chr] = []
            for seg in tmp[2][chr]:
                IBD2[chr].append(seg)

    for chr in IBD0.keys():
        IBD0[chr] = mergeIntervals(IBD0[chr][:])
    for chr in IBD1.keys():
        IBD1[chr] = mergeIntervals(IBD1[chr][:])
    for chr in IBD2.keys():
        IBD2[chr] = mergeIntervals(IBD2[chr][:])

    return [IBD0,IBD1,IBD2]


def collectIBDsegmentsSibsAvuncular(sibset, avunc,file_for_segments):  # n is number of individuals
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}
        for ind2 in avunc:
            if not ind2 in IBD_all[ind1].keys():
                IBD_all[ind1][ind2] = []
            IBD0 = {}
            IBD1 = {}
            IBD2 = {}
            tmp = getIBDsegments(ind1, ind2,file_for_segments)
            for chr in tmp[0].keys():
                if not chr in IBD0.keys():
                    IBD0[chr] = []
                for seg in tmp[0][chr]:
                    IBD0[chr].append(seg)
                IBD0[chr].sort()

            for chr in tmp[1].keys():
                if not chr in IBD1.keys():
                    IBD1[chr] = []
                for seg in tmp[1][chr]:
                    IBD1[chr].append(seg)
                IBD1[chr].sort()

            for chr in tmp[2].keys():
                if not chr in IBD2.keys():
                    IBD2[chr] = []
                for seg in tmp[2][chr]:
                    IBD2[chr].append(seg)
                IBD2[chr].sort()

            IBD_all[ind1][ind2] = [IBD0, IBD1, IBD2]

    return IBD_all


def collectIBDsegmentsSibsAvuncularCombine(sibset, avunc, file_for_segments):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    # also merge IBD0/1/2 intervals
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        IBD_all[ind1] = {}
        IBD_all[ind1]['A'] = []
        IBD0 = {}
        IBD1 = {}
        IBD2 = {}
        for ind2 in avunc:
            tmp = getIBDsegments(ind1, ind2, file_for_segments)
            for chr in tmp[0].keys():
                if not chr in IBD0.keys():
                    IBD0[chr] = []
                for seg in tmp[0][chr]:
                    IBD0[chr].append(seg)
                #IBD0[chr].sort()

            for chr in tmp[1].keys():
                if not chr in IBD1.keys():
                    IBD1[chr] = []
                for seg in tmp[1][chr]:
                    IBD1[chr].append(seg)
                #IBD1[chr].sort()

            for chr in tmp[2].keys():
                if not chr in IBD2.keys():
                    IBD2[chr] = []
                for seg in tmp[2][chr]:
                    IBD2[chr].append(seg)
                #IBD2[chr].sort()

        for chr in IBD0.keys():
            IBD0[chr] = mergeIntervals(IBD0[chr])
        for chr in IBD1.keys():
            IBD1[chr] = mergeIntervals(IBD1[chr])
        for chr in IBD2.keys():
            IBD2[chr] = mergeIntervals(IBD2[chr])

        IBD_all[ind1]['A'] = [IBD0, IBD1, IBD2]

    return IBD_all



def checkOverlap(range1, range2):
    #check if two numerical ranges overlap
    if not range1[1] <= range2[0] and not range2[1] <= range1[0]:  # one range doesn't end before start of other range
        return 1
    else:
        return 0



def findOverlap(sibseg, avsib, ss1, sa1, sa2, ranges, Eval):
    # Find regions of the genome which have sibling and avuncular IBD states as defined by ss1, sa1, sa2
    # ranges = ranges we already have in place and therefore cannot overlap; we update and return ranges with added info
    # Eval = expected amount of parent genome we get with this ss1/sa1/sa2 combination
    # sibseg = pairwise IBD segments between siblings
    # avsib = pairwise IBD segments between siblings and avuncular
    # ss1 = IBD type (0/1/2) between siblings
    # sa1 = IBD type (0/1/2) between one of those siblings and avuncular
    # sa2 = IBD type (0/1/2) between the other sibling and the avuncular
    # Eval=1
    # ss1 = 0
    # sa1 = 0
    # sa2 = 0
    # chr = 1
    #
    # For IBD2 between cousins' parents (siblings):
    # sibseg = collectIBDsegments(sib1, file_for_segments)
    # avsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, file_for_segments)
    # IBD011 = findOverlap(sibseg, avsib, 0, 1, 1, {}, 0.5)
    all_seg = {}
    for chr in range(1, 23):
        if not chr in ranges.keys():
            ranges[chr] = []
        if not chr in all_seg.keys():
            all_seg[chr] = []
        for sib1 in sibseg.keys():
            for sib2 in sibseg[sib1].keys():
                for av in avsib[sib1].keys():  # avsib[sib1].keys() and avsib[sib2].keys() are the same
                    ranges_to_append = []
                    ksib = 0
                    kav1 = 0
                    kav2 = 0
                    krange_cont = 0
                    while chr in sibseg[sib1][sib2][ss1].keys() and chr in avsib[sib1][av][sa1].keys() and chr in \
                            avsib[sib2][av][sa2].keys() and ksib < len(sibseg[sib1][sib2][ss1][chr]) and kav1 < len(
                            avsib[sib1][av][sa1][chr]) and kav2 < len(avsib[sib2][av][sa2][chr]):

                        if checkOverlap(sibseg[sib1][sib2][ss1][chr][ksib],
                                        avsib[sib1][av][sa1][chr][kav1]) and checkOverlap(
                                sibseg[sib1][sib2][ss1][chr][ksib], avsib[sib2][av][sa2][chr][kav2]) and checkOverlap(
                                avsib[sib2][av][sa2][chr][kav2], avsib[sib1][av][sa1][chr][kav1]):
                            # if all three segments overlap
                            range_add = [max(sibseg[sib1][sib2][ss1][chr][ksib][0], avsib[sib1][av][sa1][chr][kav1][0],
                                             avsib[sib2][av][sa2][chr][kav2][0]),
                                         min(sibseg[sib1][sib2][ss1][chr][ksib][1], avsib[sib1][av][sa1][chr][kav1][1],
                                             avsib[sib2][av][sa2][chr][kav2][1]), Eval, sib1, sib2, av]
                            to_append = [range_add[0],range_add[1],sib1,sib2,av]
                            all_seg[chr].append(to_append)
                            if not krange_cont:
                                krange = 0
                                while krange < len(ranges[chr]) and ranges[chr][krange][1] <= range_add[0]:
                                    krange = krange + 1

                            if krange < len(ranges[chr]):
                                if range_add[0:2] != ranges[chr][krange][0:2]:
                                    if checkOverlap(range_add, ranges[chr][krange]):
                                        range_new = []
                                        if range_add[0] < ranges[chr][krange][0]:  # new range starts before ranges[krange]
                                            if krange > 0:
                                                range_new.append([max(range_add[0], ranges[chr][krange - 1][1]),
                                                                  ranges[chr][krange][0], Eval, sib1, sib2, av])
                                            else:
                                                range_new.append(
                                                    [range_add[0], ranges[chr][krange][0], Eval, sib1, sib2, av])
                                        if range_add[1] > ranges[chr][krange][1]:  # new range ends after krange
                                            if krange < len(ranges[chr]) - 1:
                                                new_range = [ranges[chr][krange][1],
                                                             min(range_add[1], ranges[chr][krange + 1][0]), Eval, sib1,
                                                             sib2, av]
                                                if new_range[0] != new_range[1]:
                                                    range_new.append(new_range)
                                                if new_range[0] > new_range[1]:
                                                    print('ERROR: '+sib1+'\t'+sib2+'\t'+av+'\t'+str(chr) + '\t' + str(ranges[chr][krange][1]) + '\t' + str(
                                                        ranges[chr][krange + 1][0]) + '\n')
                                            else:
                                                range_new.append(
                                                    [ranges[chr][krange][1], range_add[1], Eval, sib1, sib2, av])
                                        # krange = krange + 1
                                        for seg in range_new:
                                            if not seg in ranges_to_append and not seg[0] == seg[1]:
                                                ranges_to_append.append(seg)
                                    else:  # no overlap between range_add and ranges[chr][krange]
                                        if krange > 0:
                                            range_add = [max(ranges[chr][krange - 1][1], range_add[0]), range_add[1],
                                                         Eval, sib1, sib2, av]
                                            if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                                ranges_to_append.append(range_add)
                                        else:
                                            if not range_add in ranges_to_append:
                                                ranges_to_append.append(range_add)
                            else:
                                if krange > 0:
                                    range_add = [max(ranges[chr][krange - 1][1], range_add[0]), range_add[1], Eval,
                                                 sib1, sib2, av]
                                    if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                        ranges_to_append.append(range_add)
                                else:
                                    if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                        ranges_to_append.append(range_add)
                            if krange < len(ranges[chr]):

                                if sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= \
                                                avsib[sib2][av][sa2][chr][kav2][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= ranges[chr][krange][1]:
                                    ksib = ksib + 1
                                    krange_cont = 0

                                elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][
                                            1] and avsib[sib1][av][sa1][chr][kav1][1] <= ranges[chr][krange][1]:
                                    kav1 = kav1 + 1
                                    krange_cont = 0

                                elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][
                                            1] and avsib[sib2][av][sa2][chr][kav2][1] <= ranges[chr][krange][1]:
                                    kav2 = kav2 + 1
                                    krange_cont = 0

                                elif ranges[chr][krange][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                ranges[chr][krange][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                ranges[chr][krange][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                                    krange = krange + 1
                                    krange_cont = 1

                            else:
                                if sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= \
                                                avsib[sib2][av][sa2][chr][kav2][1]:
                                    ksib = ksib + 1
                                    krange_cont = 0
                                elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][
                                            1]:
                                    kav1 = kav1 + 1
                                    krange_cont = 0
                                elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][
                                            1]:
                                    kav2 = kav2 + 1
                                    krange_cont = 0

                        elif sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                        sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                            ksib = ksib + 1
                            krange_cont = 0
                        elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                        avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                            kav1 = kav1 + 1
                            krange_cont = 0
                        elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                        avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][1]:
                            kav2 = kav2 + 1
                            krange_cont = 0

                    for seg in ranges_to_append:
                        ranges[chr].append(seg)
                    ranges[chr].sort()

    return ranges



def getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, file_for_segments):
    #get total IBD length between two sets of relatives (sib1+avunc1 and sib2+avunc2)
    #return length and number of individuals in each set with IBD segments
    sibandav = sib1.copy()
    for avunc in avunc1:
        sibandav.add(avunc)

    sibandav_rel = sib2.copy()
    for avunc in avunc2:
        sibandav_rel.add(avunc)


    all_seg_IBD1 = {}
    all_seg_IBD2 = {}
    has_seg_sib1 = [0 for x in range(0,len(sib1))]
    has_seg_sib2 = [0 for x in range(0, len(sib2))]
    has_seg_avunc1 = [0 for x in range(0,len(avunc1))]
    has_seg_avunc2 = [0 for x in range(0, len(avunc2))]
    sib1 = list(sib1)
    sib2 = list(sib2)
    for ind1 in sibandav:
        for ind2 in sibandav_rel:
            tmp = getIBDsegments(ind1, ind2, file_for_segments)
            # mark if these individuals have segments that were used
            if tmp != [{},{},{}]:
                if ind1 in sib1:
                    has_seg_sib1[sib1.index(ind1)] = 1
                elif ind1 in avunc1:
                    has_seg_avunc1[avunc1.index(ind1)] = 1
                if ind2 in sib2:
                    has_seg_sib2[sib2.index(ind2)] = 1
                elif ind2 in avunc2:
                    has_seg_avunc2[avunc2.index(ind2)] = 1
            for chr in tmp[1].keys():  # add IBD1
                if not chr in all_seg_IBD1.keys():
                    all_seg_IBD1[chr] = []
                for seg in tmp[1][chr]:
                    all_seg_IBD1[chr].append(seg)
            for chr in tmp[2].keys():  # add IBD2
                if not chr in all_seg_IBD2.keys():
                    all_seg_IBD2[chr] = []
                for seg in tmp[2][chr]:
                    all_seg_IBD2[chr].append(seg)

    IBD_sum = 0
    for chr in all_seg_IBD1:
        all_seg_IBD1[chr] = mergeIntervals(all_seg_IBD1[chr][:])
        for seg in all_seg_IBD1[chr]:
            IBD_sum = IBD_sum + seg[1] - seg[0]
    for chr in all_seg_IBD2:
        all_seg_IBD2[chr] = mergeIntervals(all_seg_IBD2[chr][:])
        for seg in all_seg_IBD2[chr]:
            IBD_sum = IBD_sum + 2.0*(seg[1] - seg[0])


    return [IBD_sum, sum(has_seg_sib1), sum(has_seg_sib2), sum(has_seg_avunc1), sum(has_seg_avunc2)]


def getInferredWithRel(total_IBD, pct_par, pct_par_rel):
    # using total length of IBD (in cM) and expected percentage of parent genome present in sibling set or percentage of grandparent genome present in sib + aunt/uncle set, calculate estimated K
    if pct_par != 0 and pct_par_rel != 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par * 1 / pct_par_rel
    elif pct_par == 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par_rel
    elif pct_par_rel == 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par
    return K  # input getSiblingRelativeIBDLength



def getExpectedGP(num_sibs,num_avunc):
    exp = (1.0 - 1.0/2.0**(num_avunc)) + (1.0/2.0**(num_avunc+1))*(1.0-1.0/2.0**num_sibs)
    return exp

def getExpectedPar(num_sibs):
    exp = (1.0-1.0/2.0**num_sibs)
    return exp


def combineBothGPsKeepProportionOnlyExpectation_final(sib1, avunc1, sib2, avunc2, proportion_par_exp, proportion_gp_exp, file_for_segments, results_file, rel_graph):
    # perform ancestral genome reconstruction between two groups of related individuals (sib1+avunc1 and sib2+avunc2)
    # infers relatedness between all individuals within the two groups
    if len(sib1) == 1 and len(sib2) == 1 and len(avunc1) == 0 and len(avunc2) == 0:
        if sib1[0] < sib2[0]:
            #return [all_rel[sib1[0]][sib2[0]][3]]
            return [[sib1[0],sib2[0],all_rel[sib1[0]][sib2[0]][3],all_rel[sib1[0]][sib2[0]][3]]]
        else:
            return [[sib2[0],sib1[0],all_rel[sib2[0]][sib1[0]][3],all_rel[sib2[0]][sib1[0]][3]]]

    cont = 1
    for av1 in avunc1:
        if av1 in avunc2:
            cont = 0
    if not cont:
        deg = 'NA'
    else:
        # returns total length of genome IBD between sibandav and sibandav_rel, number of sibs in sib1 with IBD segments, number of sibs in sib2 with IBD segments
        [tmpsibav, sib1_len, sib2_len, av1_len, av2_len] = getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, file_for_segments)

        if proportion_gp_exp == -1 and proportion_par_exp == -1:
            if av1_len != 0:
                proportion_gp_exp = getExpectedGP(len(sib1),len(avunc1))
                proportion_par_exp = 0
            elif sib1_len > 1:
                proportion_gp_exp = 0
                proportion_par_exp = getExpectedPar(len(sib1))
            else:
                proportion_par_exp = 0
                proportion_gp_exp = 0
        if av2_len != 0:
            proportion_gp_rel_exp = getExpectedGP(len(sib2), len(avunc2))
            proportion_par_rel_exp = 0
        elif sib2_len > 1:
            proportion_gp_rel_exp = 0
            proportion_par_rel_exp = getExpectedPar(len(sib2))

        else:
            proportion_par_rel_exp = 0
            proportion_gp_rel_exp = 0


        both = 1
        if proportion_gp_exp != 0:
            if proportion_gp_rel_exp != 0: #both grandparents reconstructed
                K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, proportion_gp_rel_exp)
                base = 4
            elif proportion_par_rel_exp != 0: #gp1 reconstructed, par2 reconstructed
                K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, proportion_par_rel_exp)
                base = 3
            else:
                K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, 0)
                base = 2
                both = 0
        elif proportion_par_exp != 0:
            if proportion_gp_rel_exp != 0:
                K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, proportion_gp_rel_exp)
                base = 3
            elif proportion_par_rel_exp !=0:
                K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, proportion_par_rel_exp)
                base = 2
            else:
                K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, 0)
                base = 1
                both = 0
        else:
            both = 0
            if proportion_gp_rel_exp != 0:
                K_exp = getInferredWithRel(tmpsibav, 0, proportion_gp_rel_exp)
                base = 2
            elif proportion_par_rel_exp != 0:
                K_exp = getInferredWithRel(tmpsibav, 0, proportion_par_rel_exp)
                base = 1
            else:
                if sib1[0] < sib2[0]:
                    K_exp = all_rel[sib1[0]][sib2[0]][2]
                else:
                    K_exp = all_rel[sib2[0]][sib1[0]][2]
                base = 0

        estimated_exp = getInferredFromK(K_exp)

        IBD2=0
        if both and estimated_exp != 1 and K_exp > 1/2.0**(9.0/2):  #check for IBD2 if both sides are being reconstructed, might be reconstructing two sibs; K_exp is 3rd degree or closer
            if len(avunc1) and len(avunc2):
                sibseg = collectIBDsegments(avunc1, file_for_segments)
                sibsib = collectIBDsegmentsSibsAvuncularCombine(avunc1, avunc2, file_for_segments)
                IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, {}, 0.5)
                if len(sib2) > 1:  # could be the case that we have one sib and his/her aunts/uncles
                    sibseg2 = collectIBDsegments(avunc2, file_for_segments)
                    sibsib2 = collectIBDsegmentsSibsAvuncularCombine(avunc2, avunc1, file_for_segments)
                    IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, {}, 0.5)
                    for chr in IBD011_2.keys():
                        if not chr in IBD011.keys():
                            IBD011[chr] = []
                        IBD011[chr] = IBD011[chr] + IBD011_2[chr]
                        IBD011[chr] = mergeIntervals(IBD011[chr])
                IBD2 = getTotalLength(IBD011)
            elif not len(avunc1) and not len(avunc2):
                sibseg = collectIBDsegments(sib1, file_for_segments)
                sibsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, file_for_segments)
                IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, {}, 0.5)
                if len(sib2) > 1: #could be the case that we have one sib and his/her aunts/uncles
                    sibseg2 = collectIBDsegments(sib2, file_for_segments)
                    sibsib2 = collectIBDsegmentsSibsAvuncularCombine(sib2, sib1, file_for_segments)
                    IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, {}, 0.5)
                    for chr in IBD011_2.keys():
                        if not chr in IBD011.keys():
                            IBD011[chr] = []
                        IBD011[chr] = IBD011[chr] + IBD011_2[chr]
                        IBD011[chr] = mergeIntervals(IBD011[chr])
                IBD2 = getTotalLength(IBD011)

            IBD2_add = IBD2
            if proportion_par_exp != 0:
                IBD2_add = IBD2_add / proportion_par_exp
            elif proportion_gp_exp != 0:
                IBD2_add = IBD2_add / proportion_gp_exp
            if proportion_par_rel_exp != 0:
                IBD2_add = IBD2_add / proportion_par_rel_exp
            elif proportion_gp_rel_exp != 0:
                IBD2_add = IBD2_add / proportion_gp_rel_exp

            if IBD2_add != 0:
                estimated_exp = getInferredFromK(K_exp + (IBD2_add) / total_genome / 4.0)  # add in IBD2, remove IBD1

        result = []
        if estimated_exp != 0:
            if estimated_exp == -1:
                estimated_exp = base
            else:
                estimated_exp = estimated_exp + base

        #get all combinations of pairs for results
        for sib_rel in sib2:
            for sib in sib1:
                if rel_graph.has_edge(sib_rel,sib):
                    estimated_out_exp = rel_graph.get_edge_data(sib_rel,sib)['type']
                else:
                    estimated_out_exp = estimated_exp
                if sib_rel < sib:
                    refined = all_rel[sib_rel][sib][3]
                else:
                    refined = all_rel[sib][sib_rel][3]
                to_add = [sib_rel, sib, estimated_out_exp, refined]
                result.append(to_add)

            for avunc in avunc1:
                if rel_graph.has_edge(sib_rel, avunc):
                    estimated_out_exp = rel_graph.get_edge_data(sib_rel,avunc)['type']
                else:
                    if estimated_exp != 0:
                        estimated_out_exp = estimated_exp - 1
                    else:
                        estimated_out_exp = 0
                if sib_rel < avunc:
                    refined = all_rel[sib_rel][avunc][3]
                else:
                    refined = all_rel[avunc][sib_rel][3]
                to_add = [sib_rel, avunc, estimated_out_exp, refined]
                result.append(to_add)

        for avunc_rel in avunc2:
            for sib in sib1:
                if rel_graph.has_edge(avunc_rel, sib):
                    estimated_out_exp = rel_graph.get_edge_data(avunc_rel,sib)['type']
                else:
                    if estimated_exp != 0:
                        estimated_out_exp = estimated_exp - 1
                    else:
                        estimated_out_exp = 0
                if sib < avunc_rel:
                    refined = all_rel[sib][avunc_rel][3]
                else:
                    refined = all_rel[avunc_rel][sib][3]
                to_add = [sib, avunc_rel, estimated_out_exp, refined]
                result.append(to_add)

            for avunc in avunc1:
                if rel_graph.has_edge(avunc_rel, avunc):
                    estimated_out_exp = rel_graph.get_edge_data(avunc_rel,avunc)['type']
                else:
                    if estimated_exp != 0:
                        estimated_out_exp = estimated_exp - 2
                    else:
                        estimated_out_exp = 0
                if avunc < avunc_rel:
                    refined = all_rel[avunc][avunc_rel][3]
                else:
                    refined = all_rel[avunc_rel][avunc][3]
                to_add = [avunc, avunc_rel, estimated_out_exp, refined]
                result.append(to_add)


    return result

def checkRelevantAuntsUncles(sibset1, sibset2, avunc1_bothsides, avunc2_bothsides, par1, par2):
    # check whether aunt/uncle should be included in analysis
    avunc1 = []
    avunc2 = []
    if len(avunc1_bothsides):
        for i in range(0, len(avunc1_bothsides)):
            avunc1 = avunc1 + avunc1_bothsides[i]

    if len(avunc2_bothsides):
        for i in range(0, len(avunc2_bothsides)):
            avunc2 = avunc2 + avunc2_bothsides[i]


    for av1 in avunc1:
        if av1 in avunc2 or av1 in par2:
            return ['same',[], [], []]
        if av1 in sibset2:
            return ['av',[], [], []]

    for av2 in avunc2: #
        if av2 in par1:
            return ['same', [], [], []]
        if av2 in sibset1:
            return ['av',[], [], []]


    allsibs = []
    for s1 in sibset1:
        for s2 in sibset2:
            if s1 < s2:
                if s1 in all_rel.keys() and s2 in all_rel[s1].keys():
                    allsibs.append(float(all_rel[s1][s2][2])) #K between s1 and s2
            else:
                if s2 in all_rel.keys() and s1 in all_rel[s2].keys():
                    allsibs.append(float(all_rel[s2][s1][2]))


    for av1 in avunc1:
        for av2 in avunc2:
            if av1 < av2:
                if all_rel[av1][av2][3] == '1':
                    if float(all_rel[av1][av2][0]) > 0.9:
                        return ['avparent',[av1,av2], [], []]
                    elif float(all_rel[av1][av2][1] > 1/2.0**(3.0/2)):
                        return ['sib', [av1,av2], [], []]


    minsib = min(allsibs)
    relavunc1 = []
    relavunc2 = []
    for s1 in sibset1:
        for a2 in avunc2:
        #for a2 in avunc2_keep:
            if s1 < a2:
                if all_rel[s1][a2][3] == '1':
                    if float(all_rel[s1][a2][0]) + float(all_rel[s1][a2][1]) > 0.9: #a2 is parent of s1
                        return ['sibparent',[s1,a2], [], []]
                else:
                    if float(all_rel[s1][a2][2]) > minsib:
                        relavunc2.append(a2)
            else:
                if all_rel[a2][s1][3] == '1':
                    if float(all_rel[a2][s1][0]) + float(all_rel[a2][s1][1]) > 0.9:  # a2 is parent of s1
                        return ['sibparent',[s1,a2], [], []]
                else:
                    if float(all_rel[a2][s1][2]) > minsib:
                        relavunc2.append(a2)

    for s2 in sibset2:
        #for a1 in avunc1_keep:
        for a1 in avunc1:
            if s2 < a1:
                if all_rel[s2][a1][3] == '1':
                    if float(all_rel[s2][a1][0]) + float(all_rel[s2][a1][1]) > 0.9:  # a2 is parent of s1
                        return ['sibparent',[s2,a1], [], []]
                else:
                    if float(all_rel[s2][a1][2]) > minsib:
                        relavunc1.append(a1)
            else:
                if all_rel[a1][s2][3] == '1':
                    if float(all_rel[a1][s2][0]) + float(all_rel[a1][s2][1]) > 0.9:  # a2 is parent of s1
                        return ['sibparent',[s2,a1], [], []]
                else:
                    if float(all_rel[a1][s2][2]) > minsib:
                        relavunc1.append(a1)

    relavunc1 = list(set(relavunc1))
    relavunc2 = list(set(relavunc2))

    if len(avunc1_bothsides):
        for i1 in relavunc1:
            for i in range(0,len(avunc1_bothsides)):
                if i1 in avunc1_bothsides[i]:
                    relavunc1 = relavunc1 + avunc1_bothsides[i]

    if len(avunc2_bothsides):
        for i2 in relavunc2:
            for i in range(0,len(avunc2_bothsides)):
                if i2 in avunc2_bothsides[i]:
                    relavunc2 = relavunc2 + avunc2_bothsides[i]

    relavunc1 = list(set(relavunc1))
    relavunc2 = list(set(relavunc2))



    unused1 = []
    unused2 = []
    for ind1 in avunc1:
        if not ind1 in relavunc1:
            unused1.append(ind1)

    for ind2 in avunc2:
        if not ind2 in relavunc2:
            unused2.append(ind2)

    return [relavunc1, relavunc2, unused1, unused2]


def getTotalLength(IBD):
    # get length of IBD segments
    total = 0
    for chr in IBD.keys():
        for seg in IBD[chr]:
            total = total + seg[1] - seg[0]

    return total

def getAllRel(results_file, inds_file):
    # read in results file:
    # all_rel: dict of ind1, dict of ind2, list of [IBD1, IBD2, K, D
    # store pairwise relatedness information
    global all_rel
    global inds
    first = [] #list of first degree relative pairs, according to Refined IBD results
    second = [] #list of second degree relative pairs, according to Refined IBD results
    inds = []
    if inds_file != 'NA':
        file = open(inds_file,'r')
        for line in file:
            l = str.split(line.rstrip())
            inds.append(l[0])

        file.close()

    all_rel = {}
    file = open(results_file,'r')
    if inds_file == 'NA':
        for line in file:
            l = str.split(line.rstrip())
            if not l[0] in inds:
                inds.append(l[0])
            if not l[1] in inds:
                inds.append(l[1])
            K = float(l[2])/4.0 + float(l[3])/2.0
            if l[0] < l[1]:
                if not l[0] in all_rel.keys():
                    all_rel[l[0]] = {} #IBD1, IBD2, K, D
                all_rel[l[0]][l[1]] = [float(l[2]),float(l[3]), K, degree]
                if degree == 1:
                    first.append([l[0],l[1]])
                elif degree == 2:
                    second.append([l[0],l[1]])
            else:
                if not l[1] in all_rel.keys():
                    all_rel[l[1]] = {}
                all_rel[l[1]][l[0]] = [float(l[2]), float(l[3]), K, degree]
                if degree == 1:
                    first.append([l[1],l[0]])
                elif degree == 2:
                    second.append([l[1],l[0]])
    else:
        for line in file:
            l = str.split(line.rstrip())
            if l[0] in inds and l[1] in inds:
                K = float(l[2]) / 4.0 + float(l[3]) / 2.0
                degree = getInferredFromK(K)
                if l[0] < l[1]:
                    if not l[0] in all_rel.keys():
                        all_rel[l[0]] = {}  # IBD1, IBD2, K, D
                    all_rel[l[0]][l[1]] = [float(l[2]), float(l[3]), K, degree]
                    if degree == 1:
                        first.append([l[0], l[1]])
                    elif degree == 2:
                        second.append([l[0], l[1]])
                else:
                    if not l[1] in all_rel.keys():
                        all_rel[l[1]] = {}
                    all_rel[l[1]][l[0]] = [float(l[2]), float(l[3]), K, degree]
                    if degree == 1:
                        first.append([l[1], l[0]])
                    elif degree == 2:
                        second.append([l[1], l[0]])

    file.close()

    for [ind1,ind2] in itertools.combinations(inds, 2):
        if ind1 < ind2:
            if not ind1 in all_rel.keys():
                all_rel[ind1] = {}
            if not ind2 in all_rel[ind1].keys():
                all_rel[ind1][ind2] = [-1, -1, -1, -1]
        else:
            if not ind2 in all_rel.keys():
                all_rel[ind2] = {}
            if not ind1 in all_rel[ind2].keys():
                all_rel[ind2][ind1] = [-1,-1,-1,-1]


    return [all_rel,inds,first,second]


def getSecondDegreeRelativesFromAllRel(all_rel,inds,sibset):
    # collect all second degree relatives
    # for DRUID_C
    second = []
    for ind in inds:
        if not ind in sibset:
            k = 0
            deg = 2
            while k < len(sibset) and deg == 2:
                if ind < list(sibset)[k]:
                    deg = all_rel[ind][list(sibset)[k]][3]
                else:
                    deg = all_rel[list(sibset)[k]][ind][3]
                k = k + 1

            if k == len(sibset) and deg == 2:
                second.append(ind)

    return second

def getSecondDegreeRelativesFromAllRel_lessRestrictive(all_rel,inds,sibset):
    # collect all inferred second and third degree relatives of any sibling (to check for aunts/uncles later)
    second = []
    for ind in inds:
        if not ind in sibset:
            degs = []
            for k in range(0,len(sibset)):
                if ind < list(sibset)[k]:
                    degs.append(all_rel[ind][list(sibset)[k]][3])
                else:
                    degs.append(all_rel[list(sibset)[k]][ind][3])
            if 2 in degs or 3 in degs:
                second.append(ind)

    return second



def getAuntsUncles_IBD011_nonoverlapping_pairs(all_rel, sibset, halfsibs, second, file_for_segments, rel_graph):
    # check whether individuals in list 'second' are likely aunts/uncles of 'sibset' and possibly also 'halfsibs'
    avunc = []
    remove_second = []
    avunc_hs1 = []
    avunc_hs2 = []
    if len(second):
        for ind in second:
            for sib in sibset:
                if ind < sib:
                    if float(all_rel[ind][sib][1]) > 0.01 or (rel_graph.has_edge(ind,sib) and rel_graph.get_edge_data(ind,sib)['type'] in ['P','HS','GP']): #if proportion of genome shared IBD2 is > 1%
                        remove_second.append(ind)
                        break
                else:
                    if float(all_rel[sib][ind][1]) > 0.01 or (rel_graph.has_edge(sib,ind) and rel_graph.get_edge_data(sib,ind)['type'] in ['P','HS','GP']):
                        remove_second.append(ind)
                        break

        for ind in remove_second:
            second.remove(ind)

        sibset = list(sibset)
        if len(second):
            numpairs = int(len(sibset) / 2)
            for i in range(0, numpairs):
                sibseg = collectIBDsegments(sibset[(i*2):(i*2+2)],file_for_segments)
                for av in second:
                    avsib = collectIBDsegmentsSibsAvuncular(sibset[(i*2):(i*2+2)], [av],file_for_segments)
                    IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, {}, 0.5))
                    if IBD011 > 50:
                        avunc.append(av)
                        avsibs = getSibsFromGraph(rel_graph, av)
                        for avsib in avsibs:
                            avunc.append(av)
                    elif IBD011 < 20:
                        break

            #check with halfsibs
            if len(halfsibs):
                numpairs = int(len(sibset) * len(halfsibs[0]))
                allsibs = sibset + halfsibs[0]
                for i in range(0, numpairs):
                    sibseg = collectIBDsegments(allsibs[(i * 2):(i * 2 + 2)], file_for_segments)
                    for av in second:
                        avsib = collectIBDsegmentsSibsAvuncular(allsibs[(i * 2):(i * 2 + 2)], [av], file_for_segments)
                        IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, {}, 0.5))
                        if IBD011 > 50:
                            avunc_hs1.append(av)
                            avsibs = getSibsFromGraph(rel_graph, av)
                            for avsib in avsibs:
                                avunc_hs1.append(av)
                            is_av = 1
                        elif IBD011 < 20:
                            break

                if len(halfsibs) > 1 and len(halfsibs[1]):
                    numpairs = int(len(sibset) * len(halfsibs[1]))
                    allsibs = sibset + halfsibs[1]
                    for i in range(0, numpairs):
                        sibseg = collectIBDsegments(allsibs[(i * 2):(i * 2 + 2)], file_for_segments)
                        is_av = 0
                        for av in second:
                            avsib = collectIBDsegmentsSibsAvuncular(allsibs[(i * 2):(i * 2 + 2)], [av], file_for_segments)
                            IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, {}, 0.5))
                            if IBD011 > 50:
                                avunc_hs2.append(av)
                                avsibs = getSibsFromGraph(rel_graph, av)
                                for avsib in avsibs:
                                    avunc_hs1.append(av)
                                is_av = 1
                            elif IBD011 < 20:
                                break



    return [avunc, avunc_hs1, avunc_hs2]


def checkUseHalfsibs(sibs,halfsib_sets,rel,all_rel):
    # check whether halfsibs should be included in analysis between 'sibs' and 'rel'
    # sibs = initial set of sibs
    # halfsib_sets = sibs' halfsibs
    # rel = distant relative

    if len(halfsib_sets):
        sibmin = 0
        for sib in sibs:
            if sib < rel:
                if all_rel[sib][rel][2] < sibmin or sib == sibs[0]:
                    sibmin = all_rel[sib][rel][2]
            else:
                if all_rel[rel][sib][2] < sibmin or sib == sibs[0]:
                    sibmin = all_rel[rel][sib][2]

        hsk = []
        hsk_all = []
        hsk_all_mean = []
        for hsset in halfsib_sets:
            hsmax = 0
            hsk_all_set = []
            for hs in hsset:
                if hs < rel:
                    if all_rel[hs][rel][2] > hsmax:
                        hsmax = all_rel[hs][rel][2]
                    hsk_all_set.append(all_rel[hs][rel][2])
                else:
                    if all_rel[rel][hs][2] > hsmax:
                        hsmax = all_rel[rel][hs][2]
                    hsk_all_set.append(all_rel[rel][hs][2])
            hsk.append(hsmax)
            hsk_all.append(hsk_all_set)
            hsk_all_mean.append(sum(hsk_all_set)/len(hsk_all_set))

        use_hs = -1
        use_hs_val = 0
        for i in range(0,len(hsk)):
            if hsk[i] > 3/4*sibmin and hsk_all_mean[i] > use_hs_val:
                use_hs = i
                use_hs_val = hsk_all_mean

        if use_hs != -1:
            return halfsib_sets[use_hs]
        else:
            return []
    else:
        return []






def runDRUID(rel_graph, all_rel, inds, args):
    # a chunky monkey
    all_results = []
    checked = []
    for [ind1,ind2] in itertools.combinations(inds,2): #test each pair of individuals
        if not [ind1,ind2] in checked and not [ind2,ind1] in checked: #if pair not yet tested
            print("Comparing "+ind1+" and "+ind2+"\n")
            results = []
            if rel_graph.has_edge(ind1, ind2):
                checked.append([ind1,ind2])
                if ind1 < ind2:
                    refined = all_rel[ind1][ind2][3]
                else:
                    refined = all_rel[ind2][ind1][3]
                type = rel_graph.get_edge_data(ind1, ind2)['type']
                if type == '1U':
                    type = '1'
                results = results + [[ind1, ind2, type, refined, 'graph']]
            else:
                #print(ind1+" and "+ind2+'\n')
                [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                sib1.add(ind1)
                sib2.add(ind2)
                hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                sib1 = sib1.union(set(hs1))
                sib2 = sib2.union(set(hs2))
                reltype = getRelationship(rel_graph, ind1, ind2)
                if reltype != -1:
                    if ind1 < ind2:
                        refined = all_rel[ind1][ind2][3]
                    else:
                        refined = all_rel[ind2][ind1][3]
                    results = results + [[ind1,ind2,reltype,refined, 'graph']]
                else: #no path between individuals
                    #check if ind1 has parents/grandparents more closely related to other set of individuals
                    ind1_original = ind1
                    ind1_new = checkForMoveUp(all_rel,ind1, sib1, par1.union(gp1), sib2)
                    moves1 = []
                    moves_inds1 = []
                    if ind1_new == 'same':
                        moves1.append(getRelationship(rel_graph, ind1_new, ind1))
                        moves_inds1.append(ind1)
                    while ind1 != ind1_new and ind1_new != 'same':
                        moves1.append(getRelationship(rel_graph,ind1_new, ind1))
                        moves_inds1.append(ind1)
                        ind1 = ind1_new
                        [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                        sib1.add(ind1)
                        ind1_new = checkForMoveUp(all_rel,ind1, sib1, gp1+par1, sib2)

                    if ind1_new == 'same':
                        ind1 = anyIn(gp1 + par1, sib2)[0]
                        [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                    else:
                        # check if ind2 has parents/grandparentsmore closely related to other set of individuals
                        ind2_original = ind2
                        ind2_new = checkForMoveUp(all_rel,ind2,sib2,gp2.union(par2),sib1)
                        moves2 = []
                        moves_inds2 = []
                        if ind2_new == 'same':
                            moves2.append(getRelationship(rel_graph, ind2_new, ind2))
                            moves_inds2.append(ind2)
                        while ind2 != ind2_new and ind2_new != 'same':
                            moves2.append(getRelationship(rel_graph, ind2_new, ind2))
                            moves_inds2.append(ind2)
                            ind2 = ind2_new
                            [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                            sib2.add(ind2)
                            ind2_new = checkForMoveUp(all_rel,ind2, sib2, gp2+par2, sib1)

                        if ind2_new == 'same':
                            ind2 = anyIn(gp2 + par2, sib1)[0]
                            [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)

                    # ind1 and ind2 can't be related via graph, otherwise they'd be considered above
                    # continue onto composite relatedness method

                    #switch focus to youngest generation if available
                    if ind1_new != 'same' and ind2_new != 'same':
                        if len(avunc1_bothsides) or len(avunc2_bothsides):
                            [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2)
                            if not len(relavunc1) and len(nn1):
                                # tmp = sib1[:]
                                # sib1 = nn1[:]
                                # avunc1_bothsides = [tmp[:]]
                                # nn1 = []  # doesn't matter
                                sib1 = getLargestSibsets(rel_graph, nn1)
                                ind1 = nn1[0]
                                [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                                sib1.add(ind1)
                                [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2)
                            if not len(relavunc2) and len(nn2):
                                # tmp = sib2[:]
                                # sib2 = nn2[:]
                                # avunc2_bothsides = [tmp[:]]
                                # nn2 = []  # doesn't matter
                                sib2 = getLargestSibsets(rel_graph,nn2)
                                ind2 = sib2[0]
                                [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                                sib2.add(ind2)
                                [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2)
                            # add these sibsets to checked
                            for i1 in sib1:
                                for i2 in sib2:
                                    checked.append([i1, i2])
                                    checked.append([i2, i1])
                            if len(unused1):
                                unused_check = []
                                for i1 in unused1:
                                    if not i1 in unused_check:
                                        [sib1u, avunc1u_bothsides, nn1u, par1u, child1u, gp1u, halfsib1u_sets, twins1u] = pullFamily(rel_graph, i1)
                                        sib1u.add(i1)
                                        unused_check = unused_check + sib1u
                                        results = results + combineBothGPsKeepProportionOnlyExpectation_final(sib1u, [], sib2, [], -1, -1, args.s[0], args.i[0], rel_graph)
                            if len(unused2):
                                unused_check = []
                                for i2 in unused2:
                                    if not i2 in unused_check:
                                        [sib2u, avunc2u_bothsides, nn2u, par2u, child2u, gp2u, halfsib2u_sets, twins2u] = pullFamily(rel_graph, i2)
                                        sib2u.add(ind2)
                                        unused_check = unused_check + sib2u
                                        results = results + combineBothGPsKeepProportionOnlyExpectation_final(sib1, [], sib2u, [], -1, -1, args.s[0], args.i[0], rel_graph)

                        else:
                            relavunc1 = []
                            relavunc2 = []
                        for av1 in relavunc1:
                            for av2 in relavunc2:
                                checked.append([av1,av2])
                                checked.append([av2,av1])
                            for s2 in sib2:
                                checked.append([av1,s2])
                                checked.append([s2,av1])
                        for av2 in relavunc2:
                            for s1 in sib1:
                                checked.append([av2,s1])
                                checked.append([s1,av2])
                        results = results + combineBothGPsKeepProportionOnlyExpectation_final(sib1, relavunc1, sib2, relavunc2, -1, -1, args.s[0], args.i[0], rel_graph)
                        if ind1_original != ind1 or ind2_original != ind2:
                            for res in results:
                                if (res[0] == ind1 and res[1] == ind2) or (res[0] == ind2 and res[1] == ind1):
                                    closest_result = res
                                    break
                            if closest_result[2] == '1U':
                                closest_result[2] = 1
                            total = int(closest_result[2])
                    if (ind1_new == 'same' or ind2_new == 'same'):
                        if ind1 == ind2:
                            closest_result = [ind1,ind2,0]
                        else: #siblings
                            closest_result = [ind1,ind2,1]

                    if ind1_original != ind1 or ind2_original != ind2:
                        for ii in range(len(moves1)-1,-1,-1):
                            if moves1[ii] in ['P','C']:
                                total = total + 1
                            else: #gp or gc
                                total = total + 2
                            if moves_inds1[ii] < ind2:
                                refined = all_rel[moves_inds1[ii]][ind2][3]
                            else:
                                refined = all_rel[ind2][moves_inds1[ii]][3]
                            results.append([moves_inds1[ii],ind2,total,refined,ind1_original,ind2_original, 'graph+inferred'])
                            #check for close relatives of moves_inds[ii]
                            [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[ii])
                            for s1 in sib1:
                                results.append([s1,ind2,total,refined,ind1_original,ind2_original, 'graph+inferred'])
                            sib1.append(moves_inds1[ii])
                            hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                            for h1 in hs1:
                                results.append([h1,ind2,total,refined,ind1_original,ind2_original, 'graph+inferred'])

                        total = int(closest_result[2])
                        for ii in range(len(moves2)-1,-1,-1):
                            if moves2[ii] in ['P','C']:
                                total = total + 1
                            else: #gp or gc
                                total = total + 2
                            if moves_inds2[ii] < ind1:
                                refined = all_rel[moves_inds2[ii]][ind1][3]
                            else:
                                refined = all_rel[ind1][moves_inds2[ii]][3]
                            results.append([ind1,moves_inds2[ii],total,refined,ind1_original,ind2_original, 'graph+inferred'])
                            #check for close relatives of moves_inds[ii]
                            [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[ii])
                            for s2 in sib2:
                                results.append([ind1,s2,total,refined,ind1_original,ind2_original, 'graph+inferred'])
                            sib2.append(moves_inds2[ii])
                            hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                            for h2 in hs2:
                                results.append([ind1,h2,total,refined,ind1_original,ind2_original, 'graph+inferred'])

                        if len(moves1) and len(moves2):
                            total = int(closest_result[2])
                            for i1 in range(len(moves1)-1,-1,-1):
                                if moves1[i1] in ['P', 'C']:
                                    total = total + 1
                                else:
                                    total = total + 2
                                for i2 in range(len(moves2)-1,-1,-1):
                                    if moves2[i2] in ['P','C']:
                                        total = total + 1
                                    else:
                                        total = total + 2
                                    if moves_inds1[i1] < moves_inds2[i2]:
                                        refined = all_rel[moves_inds1[i1]][moves_inds2[i2]][3]
                                    else:
                                        refined = all_rel[moves_inds2[i2]][moves_inds1[i1]][3]
                                    results.append([moves_inds1[i1],moves_inds2[i2],total,refined, 'graph+inferred'])

                                    # check for close relatives of moves_inds[ii]
                                    [sib1, avunc1_bothsides, nn1, par1, child1, gp1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[i1])
                                    [sib2, avunc2_bothsides, nn2, par2, child2, gp2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[i2])
                                    for s1 in sib1:
                                        results.append([s1, moves_inds2[i2], total, refined, 'graph+inferred'])
                                    for s2 in sib2:
                                        results.append([moves_inds1[i1], s2, total, refined, 'graph+inferred'])
                                        for s1 in sib1:
                                            results.append([s1, s2, total, refined, 'graph+inferred'])
                                    sib1 = list(sib1)
                                    sib2 = list(sib2)
                                    sib1.append(moves_inds1[i1])
                                    sib2.append(moves_inds2[i2])
                                    hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                                    hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                                    for h1 in hs1:
                                        results.append([h1, ind2, total, refined, 'graph+inferred'])
                                    for h2 in hs2:
                                        results.append([ind1, h2, total, refined, 'graph+inferred'])
                                        for h1 in hs1:
                                            results.append([h1,h2,total,refined, 'graph+inferred'])

                if len(twins1):
                    for res in results:
                        if sibs1[0] == res[0]:
                            for t1 in twins1:
                                results.append([t1,res[1],res[2],res[3],'graph'])
                        elif sibs1[0] == res[1]:
                            for t1 in twins1:
                                results.append([res[0],t1,res[2],res[3],'graph'])

                if len(twins2):
                    for res in results:
                        if sibs2[0] == res[0]:
                            for t2 in twins2:
                                results.append([t2, res[1], res[2], res[3],'graph'])
                        elif sibs2[0] == res[1]:
                            for t2 in twins2:
                                results.append([res[0], t2, res[2], res[3],'graph'])

            all_results = all_results + results

    return all_results
