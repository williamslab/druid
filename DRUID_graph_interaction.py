import subprocess
import sys
import Bio
import itertools
import networkx as nx
import multiprocessing
import random
import copy
import numpy as np



def getRelationship(tmp_graph,ind1,ind2):
    #if ind1 anad ind2 have a path between them, we find their degree of relatedness/relationship type
    if nx.has_path(tmp_graph, ind1, ind2):
        #check all paths because some shortest paths may try to travel through the other lineage
        paths = nx.all_shortest_paths(tmp_graph,ind1,ind2)
        for path in paths:
            if len(path) == 2: #path only involves those two individuals
                return tmp_graph.get_edge_data(ind1,ind2)['type']
            elif len(path) > 2:
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
                    elif type1 in ['2','1U','HS'] or type2 in ['2','1U']:
                        total = -1
                        i = len(path) + 1
                    elif type1 in ['FS','P','C']:
                        total = total + 1
                        i = i + 1
                    elif type1 == 'NN':
                        if type2 == 'P' or type2 == 'GP': #traveling to other lineage, stop
                            total = -1
                            i = len(path) + 1


                    elif type1 == 'AU':
                        if type2 == 'P': #grandparent
                            total = total + 1
                        elif type2 in ['C','GP']: #ind1 and ind2 are cousins or great-grandparent
                            total = total + 2
                        elif type2 in ['GC','HS','DC']: #cousins once removed, half-aunt/uncle, or complex
                            total = total + 3
                        elif type2 == 'NN': #niece/nephew of aunt/uncle = cousin
                            total = total + 3
                        elif type2 == '-1':
                            total = total + 2
                        i = i + 2
                    else:
                        total = -1
                        i = len(path) + 1

                if total != 0 and total != -1:
                    return total

        return -1

    else:
        return -1


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
    sibs = []
    if ind in tmp_graph:
        neighbors = tmp_graph.neighbors(ind)
        for x in neighbors:
            if checkIfSib(tmp_graph,ind,x):
                sibs.append(x)
    return sibs



def getSibsAndHalfSibsFromGraph(tmp_graph,ind):
    #return the siblings and half-sibs of ind
    hs = []
    sibs = []
    neighbors = tmp_graph.neighbors(ind)
    for x in neighbors:
        if tmp_graph.get_edge_data(ind,x)['type'] == 'HS':
            hs.append(x)
        elif tmp_graph.get_edge_data(ind,x)['type'] == 'FS':
            sibs.append(x)

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
        all_sibs = all_sibs + sibs_ind
    collected_sibs = list(set(all_sibs))
    sibs_add = []
    for x in collected_sibs:
        if not x in siblings:
            sibs_add.append(x)
    siblings = siblings + sibs_add
    for ind in sibs_add:
        sibs_ind = getSibsFromGraph(tmp_graph,ind)
        all_sibs = all_sibs + sibs_ind
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




def pullFamily(tmp_graph,ind):
    # get all possible connections in graph: siblings, aunts/uncles, parents, children, grandparents, half-siblings, and twins of ind
    edges = tmp_graph.edges(ind)
    parents = []
    children = []
    avunc = []
    sib = []
    grandparents = []
    halfsibs = []
    nn = []
    twins = []
    for edge in edges:
        edge_info = tmp_graph.get_edge_data(edge[0],edge[1])['type']
        if edge[0] == ind:
            if edge_info == 'P':
                children.append(edge[1])
            elif edge_info == 'C':
                parents.append(edge[1])
            elif edge_info == 'FS':
                if not edge[1] in sib:
                    sib.append(edge[1])
            elif edge_info == 'NN':
                avunc.append(edge[1])
            elif edge_info == 'GC':
                grandparents.append(edge[1])
            elif edge_info == 'HS':
                halfsibs.append(edge[1])
            elif edge_info == 'AU':
                nn.append(edge[1])
            elif edge_info == 'T':
                twins.append(edge[1])
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