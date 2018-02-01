import subprocess
import sys
import Bio
import itertools
import networkx as nx
import multiprocessing
import random
import copy
import numpy as np

def checkChangeLineage(tmp_graph,path):
    if len(path) > 2:
        for i in range(0,len(path)-2):
            indfirst = path[i]
            indsecond = path[i+1]
            indthird = path[i+2]
            type1 = tmp_graph.get_edge_data(indfirst, indsecond)['type']
            type2 = tmp_graph.get_edge_data(indsecond, indthird)['type']
            if (type1 == 'P' and type2 == 'C') or (type1 == 'GP' and type2 == 'GC') or (type1 == 'P' and type2 == 'GC') or (type1 == 'GP' and type2 == 'C') or (type1 == 'C' and type2 == 'P') or (type1 == 'GC' and type2 == 'GP') or (type1 == 'C' and type2 == 'GP') or (type1 == 'GC' and type2 == 'P') or (type1 == 'C' and type2 == 'AU') or (type1 == 'GC' and type2 == 'AU') or (type1 == 'P' and type2 == 'NN') or (type1 == 'GP' and type2 == 'NN'):
                return 1
        return 0
    else:
        return 0

def getRelationship(tmp_graph,ind1,ind2):
    #if ind1 anad ind2 have a path between them, we find their degree of relatedness/relationship type
    if ind1 in tmp_graph.nodes() and ind2 in tmp_graph.nodes() and nx.has_path(tmp_graph, ind1, ind2):
        #check all paths because some shortest paths may try to travel through the other lineage
        paths = nx.all_shortest_paths(tmp_graph,ind1,ind2)
        for path in paths:
            if len(path) == 2: #path only involves those two individuals
                return tmp_graph.get_edge_data(ind1,ind2)['type']
            elif len(path) > 2 and not checkChangeLineage(tmp_graph,path):
                total = 0 #degree of relatedness
                i = 2
                while i <= len(path):
                    indfirst = path[i-2]
                    indsecond = path[i-1]
                    type1 = tmp_graph.get_edge_data(indfirst, indsecond)['type']
                    if i != len(path):
                        indthird = path[i]
                        type2 = tmp_graph.get_edge_data(indsecond,indthird)['type']
                    else:
                        type2 = -1
                    if (type1 == 'P' and type2 == 'C') or (type1 == 'GP' and type2 == 'GC') or (type1 == 'P' and type2 == 'GC') or (type1 == 'GP' and type2 == 'C') or (type1 == 'C' and type2 == 'P') or (type1 == 'GC' and type2 == 'GP') or (type1 == 'C' and type2 == 'GP') or (type1 == 'GC' and type2 == 'P') or (type1 == 'C' and type2 == 'AU') or (type1 == 'GC' and type2 == 'AU') or (type1 == 'P' and type2 == 'NN') or (type1 == 'GP' and type2 == 'NN'): #traveling to other lineage, stop
                        total = -1
                        i = len(path) + 1
                    elif (type1 == 'AU' and type2 == 'C'):
                        total = -1
                        i = len(path) + 1
                    elif type1 in ['2','1U','HS'] or type2 in ['2','1U']:
                        total = -1
                        i = len(path) + 1
                    elif type1 in ['FS','P','C']:
                        total = total + 1
                        i = i + 1
                    elif type1 == 'NN':
                        if i == len(path) - 1:
                            if type2 == 'GP':
                                total = total + 4
                                i = i + 2
                        if type2 == 'P' or type2 == 'GP': #traveling to other lineage, stop
                            total = -1
                            i = len(path) + 1
                        elif type2 == 'AU' or type2 == 'NN':
                            total = total + 3
                            i = i + 2
                        else:
                            total = -1
                            i = len(path) + 1
                    elif type1 == 'AU':
                        if type2 == 'P': #grandparent
                            total = total + 1
                        elif type2 in ['C','GP']: #ind1 and ind2 are cousins or great-grandparent
                            total = total + 2
                        elif type2 in ['GC','HS','DC']: #cousins once removed, half-aunt/uncle, or complex
                            total = total + 3
                        elif type2 == 'NN': #great-aunt/uncle
                            total = total + 3
                        elif type2 == '-1':
                            total = total + 3
                        elif type2 == 'AU': #cousin
                            total = total + 3
                        i = i + 2
                    else:
                        total = -1
                        i = len(path) + 1

                if total != 0 and total != -1:
                    return total
        return -1

    else:
        return -1

def getLargestSibsets(tmp_graph,all_inds):
    sibsets = []
    checked = set()
    for ind in all_inds:
        if not ind in checked:
            [sib, avunc_bothsides, nn, par, child, gp, halfsib_sets, twins] = pullFamily(tmp_graph, ind)
            sib.add(ind)
            checked = checked.union(sib)
            sibsets.append(sib)

    sizes = [len(x) for x in sibsets]
    return sibsets[sizes.index(max(sizes))]



def checkIfSib(tmp_graph,ind1,ind2):
    #determine whether there's an edge labeled 'FS' between ind1 and ind2
    type = tmp_graph.get_edge_data(ind1,ind2)['type']
    if type == 'FS':
        return 1
    else:
        return 0


def checkIfParent(tmp_graph,all_rel,sibset,ind):
    #determine if ind is the parent of sibset
    if len(sibset) == 1:
        return 0
    else:
        par = 1
        for sib in sibset:
            if sib < ind:
                if not (sib,ind) in tmp_graph.edges(sib) or not tmp_graph.get_edge_data(sib,ind)['type'] == '1U' or not float(all_rel[sib][ind][1]) < 1/2.0**(7/2.0):
                    par = 0
                    break
            else:
                if not (sib,ind) in tmp_graph.edges(sib) or not tmp_graph.get_edge_data(sib,ind)['type'] == '1U' or not float(all_rel[ind][sib][1]) < 1/2.0**(7/2.0):
                    par = 0
                    break

    return par


def getSibsFromGraph(tmp_graph,ind):
    #return the siblings of ind
    sibs = set()
    if ind in tmp_graph:
        neighbors = tmp_graph.neighbors(ind)
        for x in neighbors:
            if checkIfSib(tmp_graph,ind,x):
                sibs.add(x)
    return sibs



def getSibsAndHalfSibsFromGraph(tmp_graph,ind):
    #return the siblings and half-sibs of ind
    hs = set()
    sibs = set()
    neighbors = tmp_graph.neighbors(ind)
    for x in neighbors:
        if tmp_graph.get_edge_data(ind,x)['type'] == 'HS':
            hs.add(x)
        elif tmp_graph.get_edge_data(ind,x)['type'] == 'FS':
            sibs.add(x)

    halfsib_sets = []
    while len(hs):
        hssibs = getSibsFromGraph(tmp_graph, hs[0])
        halfsib_set = [hs[0]]
        for hs in hssibs:
            halfsib_set.append(hs)
            hs.remove(hs)
        hs.remove(hs[0])
        halfsib_sets.append(halfsib_set)

    return [sibs, hs]



def checkSiblingSubgraph(tmp_graph,siblings):
    ###only used for file preparation###
    # check whether all nodes in a sibling subgraph have a direct edge between one-another;
    # if not, greedily remove the sibling node with the least direct edges connecting to other siblings until they're all connected
    remove = []
    all_sibs = []
    for ind in siblings: #for each sibling, collect his/her siblings currently in graph
        sibs_ind = getSibsFromGraph(tmp_graph,ind)
        all_sibs = all_sibs + list(sibs_ind)
    collected_sibs = list(set(all_sibs))
    sibs_add = []
    for x in collected_sibs:
        if not x in siblings:
            sibs_add.append(x)
    siblings = list(siblings) + sibs_add
    for ind in sibs_add:
        sibs_ind = getSibsFromGraph(tmp_graph,ind)
        all_sibs = all_sibs + list(sibs_ind)
    counts = [all_sibs.count(x) for x in all_sibs] #count number of times each reported sibling appears

    while not all(x >= round(len(list(set(all_sibs)))/2.0) for x in counts): #if any individual in all_sibs is found to be a sibling with less than half of the other siblings
        to_remove = getElementThatAppearsLeast(all_sibs)
        remove.append(to_remove)
        if to_remove in siblings:
            siblings.remove(to_remove)
        all_sibs = []
        for ind in siblings:
            sibs_ind = getSibsFromGraph(tmp_graph,ind)
            for x in sibs_ind:
                if not x in remove:
                    all_sibs.append(x)
        counts = [all_sibs.count(x) for x in all_sibs]
    return [list(set(all_sibs)), remove]

def anyIn(list1,list2):
    return [x for x in list1 if x in list2]

def checkForMoveUp(all_rel, ind, sibset, older_gen, third_party):
    #check if parent/grandparent is in dataset and is more related to third_party
    if len(older_gen):
        if anyIn(older_gen,third_party):
            return 'same'
        else:
            all_sib = set()
            for sib in sibset:
                for tp in third_party:
                    if sib < tp:
                        all_sib.add(float(all_rel[sib][tp][2]))
                    else:
                        all_sib.add(float(all_rel[tp][sib][2]))
            maxsib = max(all_sib)

            maxpar = set()
            for par in older_gen:
                all_par = set()
                for tp in third_party:
                    if par < tp:
                        all_par.add(float(all_rel[par][tp][2]))
                    else:
                        all_par.add(float(all_rel[tp][par][2]))
                maxpar.add(max(all_par))



        if max(maxpar) > maxsib:
            par_use = older_gen[maxpar.index(max(maxpar))]
            return par_use
        else:
            return ind
    else:
        return ind



def checkAuntUncleGPRelationships(tmp_graph,siblings,par):
    # ensure the siblings of 'par' are listed as aunts/uncles of 'siblings' (par = parent of siblings)
    if par!= []:
        for p in par:
            [sibpar, avunc_bothsides, nn, parpar, childpar, gppar, halfsib_sets, twins] = pullFamily(tmp_graph, p)
            for sib in siblings:
                for sp in sibpar:
                    if not tmp_graph.has_edge(sib, sp):
                        tmp_graph.add_edge(sib, sp)
                        tmp_graph.add_edge(sp, sib)
                    tmp_graph.get_edge_data(sib,sp)['type'] = 'NN'
                    tmp_graph.get_edge_data(sp,sib)['type'] = 'AU'
            if parpar != []:
                for sib in siblings:
                    for pp in parpar:
                        if not tmp_graph.has_edge(sib,pp):
                            tmp_graph.add_edge(sib,pp)
                            tmp_graph.add_edge(pp,sib)
                        tmp_graph.get_edge_data(sib, pp)['type'] = 'GC'
                        tmp_graph.get_edge_data(pp, sib)['type'] = 'GP'


def getAuntsUnclesFromGraph(tmp_graph,ind):
    #return the aunts/uncles of ind
    au = []
    neighbors = tmp_graph.neighbors(ind)
    for x in neighbors:
        if checkIfAuntUncle(tmp_graph,ind,x):
            au.append(x)
    return au



def getElementThatAppearsLeast(tmp_list):
    #return the element of a list that appears the least number of times in the list
    return min(set(tmp_list), key=tmp_list.count)

def checkAllNeighborsForSibs(tmp_graph,ind):
    neighbors = tmp_graph.neighbors(ind)
    sibs = []
    for n in neighbors:
        if tmp_graph.get_edge_data(ind,n)['type'] == 'FS':
            sibs.append(n)
    return sibs


def pullFamily(tmp_graph,ind):
    # get all possible connections in graph: siblings, aunts/uncles, parents, children, grandparents, half-siblings, and twins of ind
    edges = tmp_graph.edges(ind)
    parents = set()
    children = set()
    avunc = set()
    sib = set()
    grandparents = set()
    halfsibs = set()
    nn = set()
    twins = set()
    for edge in edges:
        edge_info = tmp_graph.get_edge_data(edge[0],edge[1])['type']
        if edge[0] == ind:
            if edge_info == 'P':
                children.add(edge[1])
            elif edge_info == 'C':
                parents.add(edge[1])
            elif edge_info == 'FS':
                if not edge[1] in sib:
                    sib.add(edge[1])
            elif edge_info == 'NN':
                avunc.add(edge[1])
            elif edge_info == 'GC':
                grandparents.add(edge[1])
            elif edge_info == 'HS':
                halfsibs.add(edge[1])
            elif edge_info == 'AU':
                nn.add(edge[1])
            elif edge_info == 'T':
                twins.add(edge[1])
    avunc_sets = []
    tmp=tmp_graph.subgraph(avunc)
    for x in nx.strongly_connected_components(tmp):
        avunc_sets.append(list(x))

    halfsib_sets = []
    while len(halfsibs):
        hssibs = getSibsFromGraph(tmp_graph, halfsibs[0])
        halfsib_set = [halfsibs[0]]
        for hs in hssibs:
            halfsib_set.append(hs)
            halfsibs.remove(hs)
        halfsibs.remove(halfsibs[0])
        halfsib_sets.append(halfsib_set)

    return [sib, avunc_sets, nn, parents, children, grandparents, halfsib_sets, twins]