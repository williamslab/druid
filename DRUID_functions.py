import itertools
import networkx as nx
import copy
from DRUID_graph_interaction import *
from DRUID_all_rel import *

global total_genome, chrom_name_to_idx, chrom_idx_to_name, num_chrs

degrees = {'MZ': 1/2.0**(3.0/2), 1: 1/2.0**(5.0/2), 2: 1/2.0**(7.0/2), 3: 1/2.0**(9.0/2), 4: 1/2.0**(11.0/2), 5: 1/2.0**(13.0/2), 6: 1/2.0**(15.0/2), 7: 1/2.0**(17.0/2), 8: 1/2.0**(19.0/2), 9: 1/2.0**(21.0/2), 10: 1/2.0**(23.0/2), 11: 1/2.0**(25.0/2), 12: 1/2.0**(27/2.0), 13: 1/2.0**(29.0/2)}  # threshold values for each degree of relatedness


def forceFamInfo(rel_graph, faminfo):
    #force the provided faminfo file information into rel_graph
    for i1 in faminfo.keys():
        for i2 in faminfo[i1].keys():
            rel_graph.add_edge(i1,i2)
            rel_graph[i1][i2]['type'] = faminfo[i1][i2]


def inferFirst(rel_graph, rel_graph_tmp, all_rel, first, second, C):
    #build graphs using first degree relative inferences
    for [i1,i2] in first+second: #iterate through currently inferred first and second degree pairs
        if C:
            fs_IBD2 = 1/2.0**(5/2.0)
            fs_kin = 1/2.0**(5/2.0)
        else:
            fs_IBD2 = 1/2.0**(11/4.0)
            fs_kin = 1/2.0**(7/2.0)
        ibd2 = getIBD2(i1, i2, all_rel)
        kinship = getPairwiseK(i1, i2, all_rel)
        if ibd2 >= fs_IBD2 and kinship >= fs_kin: #if IBD2 meets minimum threshold
            if kinship < 1/2.0**(3/2.0): #not twin
                if ibd2 < 1/2.0**(5/2.0): #lower IBD2 than expected
                    print("Warning: " + i1 + ' and ' + i2 + ' have low levels of IBD2 for siblings, may be 3/4 sibs')
                addEdgeType(i1, i2, 'FS', 'FS', rel_graph)
            else: #twin
                addEdgeType(i1, i2, 'T', 'T', rel_graph)
        else:
            if getPairwiseD(i1, i2, all_rel) == 1:
                addEdgeType(i1, i2, '1U', '1U', rel_graph)


    # ensure subgraphs of siblings are connected
    checked = set()
    sibsets = [] #collect list of sets of siblings
    for node in rel_graph.nodes():
        if not node in checked:
            #print(node+'\n')
            siblings = getSibsFromGraph(rel_graph,node)
            siblings.add(node)

            # get sibs to add to subgraph, and individuals to remove from sibling subgraph
            [add_sibs,remove] = checkSiblingSubgraph(rel_graph,siblings.copy(),C) #edit siblings list, return removed sibs


            # remove the individuals no longer believed to be siblings
            for ind in remove:
                for sib in siblings:
                    if rel_graph.has_edge(ind,sib):
                        rel_graph[ind][sib]['type'] = '1U'
                        rel_graph[sib][ind]['type'] = '1U'
                if ind in siblings:
                    siblings.remove(ind)

            # add in missing siblings
            for [ind1,ind2] in itertools.combinations(add_sibs, 2):
                addEdgeType(ind1, ind2, 'FS', 'FS', rel_graph)

            #get updated set of siblings, add to checked list
            siblings = getSibsFromGraph(rel_graph,node)
            siblings.add(node)

            for sib in siblings:
                checked.add(sib)

            sibsets.append(siblings)


            #now that sibling set is completed, look for parents
            #find neighbors of first sibling set labeled as '1'
            inds_to_check = set()
            for sib in siblings:
                neighbors = rel_graph.neighbors(sib)
                for sib_neighbor in neighbors:
                    if rel_graph.get_edge_data(sib,sib_neighbor)['type'] == '1U':
                        inds_to_check.add(sib_neighbor)


            #check if other siblings also have this neighbor and are labeled as '1'
            if len(inds_to_check):
                for ind in inds_to_check:
                    if checkIfParent(rel_graph, all_rel, siblings, ind, C):
                        if len(siblings) == 1: #only 1 parent or child, give the pair generic PC label
                            siblings_pc = getSibsFromGraph(rel_graph, ind)
                            if len(siblings_pc) and not any([x in inds_to_check for x in siblings]): #if the other individual has siblings, then the individual in "siblings" must be the child of "ind"
                                for s in siblings_pc: #add ind and his/her siblings as child of siblings[0]
                                    addEdgeType(list(siblings)[0], s, 'P', 'C', rel_graph)
                            else:
                                addEdgeType(ind, list(siblings)[0], 'PC', 'PC', rel_graph)
                        else:
                            for sib in siblings: #add ind as parent for each sibling
                                addEdgeType(ind, sib, 'P', 'C', rel_graph)

    for sibset in sibsets:
        sibset = list(sibset)
        pars = getParent(rel_graph,sibset[0]) #parents of sibset
        for par in pars:
            [sib_par,par_par] = getSibsParentsFromGraph(rel_graph,par) #parents of parents of sibset (gp)
            for sib in sibset:
                for sp in sib_par: #for each sibling of the parent
                    if not rel_graph.has_edge(sib,sp):
                        addEdgeType(sib,sp,'NN','AU',rel_graph)
                for pp in par_par: #for each parent of the parent
                    if not rel_graph.has_edge(sib, pp):
                        addEdgeType(sib,pp,'GC','GP',rel_graph)




    # compare inferred graph to provided graph
    for edge in rel_graph_tmp.edges():
        if not edge in rel_graph.edges():
            print("Warning: Unable to confirm " + edge[0] + " and " + edge[1] + " as " + str(rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']) + " but including as such")
            rel_graph.add_edge(edge[0],edge[1])
            rel_graph[edge[0]][edge[1]]['type'] = rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']
        elif rel_graph_tmp.get_edge_data(edge[0], edge[1])['type'] != rel_graph.get_edge_data(edge[0], edge[1])['type']:
            print("Warning: Unable to confirm " + edge[0] + " and " + edge[1] + " as " + str(rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']) + " but including as such")
            rel_graph[edge[0]][edge[1]]['type'] = rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']

    # #ensure sibsets have same relatives
    # for sibset in sibsets:
    #     #collect neighbors of the sibs
    #     neighbor_set = set()
    #     for ind in sibset:
    #         nei = rel_graph.neighbors(ind)
    #         for n in nei:
    #             neighbor_set = neighbor_set.union((set(n,rel_graph.get_edge_data(ind,n)['type'])))
    #     for n in neighbor_set:
    #         for ind in sibset:
    #             if not rel_graph.has_edge(ind,n):
    #                 addEdgeType


def readSegments(file_for_segments):
    all_segs = {}

    IBD_file = open(file_for_segments, 'r')
    for line in IBD_file:
        l = str.split(line.rstrip())
        if l[0] < l[1]:
            ids = [ l[0], l[1] ]
        else:
            ids = [ l[1], l[0] ]

        if not ids[0] in all_segs:
            all_segs[ ids[0] ] = \
                    { ids[1]: [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ] }
        elif not ids[1] in all_segs[ ids[0] ]:
            all_segs[ ids[0] ][ ids[1] ] = \
                                [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
        chrom_name = l[2]
        chr = chrom_name_to_idx[chrom_name]
        ibd_type = int(l[3][3]) # chop "IBD" off, get integer type IBD_1_ or 2
        all_segs[ ids[0] ][ ids[1] ][ibd_type - 1][chr].append([float(l[4]), float(l[5])])

    IBD_file.close()

    return all_segs


def inferSecondPath(rel_graph, rel_graph_tmp, all_rel, second, all_segs, outfile, C):
    # infer and add 2nd degree relationships
    if C:
        dc_lower = 1/2.0**(9/2.0) #minimum IBD2 for DC classification
    else:
        dc_lower = 1/2.0**(9/2.0)
    for [i1, i2] in second:
        if not rel_graph.has_edge(i1, i2):
            if getIBD2(i1, i2, all_rel) < dc_lower: #proportion IBD2 less than requirement for DC classification
                addEdgeType(i1, i2, '2', '2', rel_graph)
            else: #proportion IBD2 within requirement for DC classification
                sib1 = getSibsFromGraph(rel_graph,i1)
                sib2 = getSibsFromGraph(rel_graph,i2)
                sib1.add(i1)
                sib2.add(i2)
                #if one i1 is a DC of i2, then siblings of i1 are DC of siblings of i2 (and i2)
                for s1 in sib1:
                    for s2 in sib2:
                        addEdgeType(s1, s2, 'DC', 'DC', rel_graph)

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [sibs, halfsibs, par] = getSibsHalfSibsParentsFromGraph(rel_graph, node)
            sibs.add(node)

            if len(sibs) > 1:
                #print('TESTING '+" ".join(sibs)+'\n')
                second_av = getSecondDegreeRelatives(rel_graph, second, sibs, par, all_rel)
                [avunc, avunc_hs_all] = getAuntsUncles_IBD011_nonoverlapping_pairs(sibs, halfsibs, second_av, all_rel, all_segs, rel_graph)

                # add the inferred avuncular relationships to graph
                for av in avunc:
                    for sib in sibs:
                        if not rel_graph_tmp.has_edge(av,sib):
                            addEdgeType(av,sib,'AU','NN',rel_graph) # if the provided family information doesn't contain this relationship, add it
                        else:
                            print(av+" inferred as aunt/uncle of "+sib+", but will continue using provided relationship type "+rel_graph_tmp[av][sib]['type']+'\n')


                if len(avunc_hs_all):
                    for hs in range(0,len(avunc_hs_all)):
                        for av in avunc_hs_all[hs]:
                            for sib in sibs.union(set(halfsibs[hs])):
                                if not rel_graph_tmp.has_edge(av, sib):
                                    addEdgeType(av, sib, 'AU', 'NN', rel_graph)  # if the provided family information doesn't contain this relationship, add it
                                else:
                                    print(av + " inferred as aunt/uncle of " + sib + ", but will continue using provided relationship type " + rel_graph_tmp[av][sib]['type'] + '\n')

            for sib in sibs:
                checked.add(sib)

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [siblings, avunc_bothsides, nn, par, child, pc, gp, gc, halfsib_sets, twins] = pullFamily(rel_graph, node)
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
        if l != []:
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
            else:
                if not l[0] in inds:
                    print("Warning: "+l[0]+" not included in .inds file, not including "+l[2]+" relationship with "+l[1])
                if not l[1] in inds:
                    print("Warning: "+l[1]+" not included in .inds file, not including "+l[2]+" relationship with "+l[0])

    file.close()

    return faminfo

def getChrInfo(mapfile):
    #read in information from .map file
    global total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs
    chrom_name_to_idx = {}
    chrom_idx_to_name = []
    chrom_starts = []
    chrom_ends = []
    num_chrs = 0

    file = open(mapfile,'r')
    for line in file:
        l = str.split(line.rstrip())
        chr_name = l[0]
        if not chr_name in chrom_name_to_idx.keys():
            chrom_name_to_idx[chr_name] = chr = num_chrs
            chrom_idx_to_name.append(chr_name)
            num_chrs += 1
            chrom_starts.append(99999999)
            chrom_ends.append(0)
        else:
            chr = chrom_name_to_idx[chr_name]

        pos = float(l[2])
        if chrom_starts[chr] > pos:
            chrom_starts[chr] = pos
        if chrom_ends[chr] < pos:
            chrom_ends[chr] = pos

    file.close()

    total_genome = 0
    for chr in range(num_chrs):
        total_genome += chrom_ends[chr] - chrom_starts[chr]

    return [total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs]


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


def getIBDsegments(ind1, ind2, all_segs):
    # get IBD segments between ind1 and ind2, sorting segments by IBD2, IBD1, and IBD0
    if ind1 < ind2:
        ids = [ ind1, ind2 ]
    else:
        ids = [ ind2, ind1 ]
    if not ids[0] in all_segs.keys() or not ids[1] in all_segs[ ids[0] ].keys():
        return [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]

    return all_segs[ ids[0] ][ ids[1] ]


def getIBD0(IBD1,IBD2):
    #IBD12 = regions that are IBD (IBD1 or IBD2)
    IBD12 = { chr : mergeIntervals(IBD1[chr] + IBD2[chr])
              for chr in range(num_chrs) }

    IBD0 = { chr : [] for chr in range(num_chrs) }
    for chr in range(num_chrs):
        if len(IBD12[chr]) > 0:
            if IBD12[chr][0][0] > chrom_starts[chr]:
                IBD0[chr].append([chrom_starts[chr], IBD12[chr][0][0]])
            if len(IBD12[chr]) > 1:
                for k in range(1, len(IBD12[chr])):
                    if IBD12[chr][k - 1][1] != IBD12[chr][k][0]:
                        IBD0[chr].append([IBD12[chr][k - 1][1], IBD12[chr][k][0]])
                if IBD12[chr][k][1] < chrom_ends[chr]:
                    IBD0[chr].append([IBD12[chr][k][1], chrom_ends[chr]])
            else:
                IBD0[chr].append([IBD12[chr][0][1], chrom_ends[chr]])

    return IBD0



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


def collectIBDsegments(sibset, all_segs):
    # collect pairwise IBD0,1,2 regions between all pairs of siblings
    IBD_all = {}
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}

        IBD_all[ind1][ind2] = None

        tmp = getIBDsegments(ind1, ind2, all_segs)
        tmp0 = getIBD0(tmp[0],tmp[1])
        for chr in range(num_chrs):
            tmp0[chr].sort()
            tmp[0][chr].sort()
            tmp[1][chr].sort()

        IBD_all[ind1][ind2] = [tmp0, tmp[0], tmp[1]]

    return IBD_all


any_in = lambda a, b: any(i in b for i in a)

def collectAllIBDsegments(sibset):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD0 = { chr : [] for chr in range(num_chrs) }
    IBD1 = { chr : [] for chr in range(num_chrs) }
    IBD2 = { chr : [] for chr in range(num_chrs) }
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        tmp = getIBDsegments(ind1, ind2)
        tmp0 = getIBD0(tmp[0],tmp[1])
        for chr in range(num_chrs):
            IBD0[chr] += tmp0[chr]
            IBD1[chr] += tmp[0][chr]
            IBD2[chr] += tmp[1][chr]

    for chr in range(num_chrs):
        IBD0[chr] = mergeIntervals(IBD0[chr][:])
        IBD1[chr] = mergeIntervals(IBD1[chr][:])
        IBD2[chr] = mergeIntervals(IBD2[chr][:])

    return [IBD0,IBD1,IBD2]


def collectIBDsegmentsSibsAvuncular(sibset, avunc, all_segs):  # n is number of individuals
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}
        for ind2 in avunc:
            if not ind2 in IBD_all[ind1].keys():
                IBD_all[ind1][ind2] = None

            tmp = getIBDsegments(ind1, ind2, all_segs)
            #tmp0 = getIBD0(tmp[0],tmp[1])
            for chr in range(num_chrs):
                #tmp0[chr].sort()
                tmp[0][chr].sort()
                tmp[1][chr].sort()

            IBD_all[ind1][ind2] = [[],tmp[0],tmp[1]]

    return IBD_all


def collectIBDsegmentsSibsAvuncularCombine(sibset, avunc, all_segs):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    # also merge IBD0/1/2 intervals
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        IBD_all[ind1] = {}
        IBD_all[ind1]['A'] = []
        tmp_ind1 = [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
        for ind2 in avunc:
            tmp = getIBDsegments(ind1, ind2, all_segs) #[IBD1, IBD2]
            # for chr in range(num_chrs):
            #     tmp[0][chr].sort()
            # for chr in range(num_chrs):
            #     tmp[1][chr].sort()
            for chr in range(num_chrs):
                tmp_ind1[0][chr] += tmp[0][chr]
                tmp_ind1[1][chr] += tmp[1][chr]

        for chr in range(num_chrs):
            tmp_ind1[0][chr] = mergeIntervals(tmp_ind1[0][chr][:])
            tmp_ind1[1][chr] = mergeIntervals(tmp_ind1[1][chr][:])

        IBD_all[ind1]['A'] = [{},tmp_ind1[0],tmp_ind1[1]] #return IBD1, IBD2

    return IBD_all



def checkOverlap(range1, range2):
    #check if two numerical ranges overlap
    if not range1[1] <= range2[0] and not range2[1] <= range1[0]:  # one range doesn't end before start of other range
        return 1
    else:
        return 0



def findOverlap(sibseg, avsib, ss1, sa1, sa2, Eval):
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
    #
    # For IBD2 between cousins' parents (siblings):
    # sibseg = collectIBDsegments(sib1, all_segs)
    # avsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, all_segs)
    # IBD011 = findOverlap(sibseg, avsib, 0, 1, 1, 0.5)
    all_seg = { chr : [] for chr in range(num_chrs) }
    ranges = { chr : [] for chr in range(num_chrs) }
    for chr in range(num_chrs):
        for sib1 in sibseg.keys():
            for sib2 in sibseg[sib1].keys():
                for av in avsib[sib1].keys():  # avsib[sib1].keys() and avsib[sib2].keys() are the same
                    ranges_to_append = []
                    ksib = 0
                    kav1 = 0
                    kav2 = 0
                    krange_cont = 0
                    while ksib < len(sibseg[sib1][sib2][ss1][chr]) and kav1 < len(avsib[sib1][av][sa1][chr]) and kav2 < len(avsib[sib2][av][sa2][chr]):

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
                                                    chr_name = chrom_idx_to_name[chr]
                                                    print('ERROR: '+sib1+'\t'+sib2+'\t'+av+'\t'+ chr_name + '\t' + str(ranges[chr][krange][1]) + '\t' + str(
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

                    ranges[chr] += ranges_to_append
                    ranges[chr].sort()

    return ranges



def getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs):
    #get total IBD length between two sets of relatives (sib1+avunc1 and sib2+avunc2)
    #return length and number of individuals in each set with IBD segments
    sibandav = sib1.copy()
    for avunc in avunc1:
        sibandav.add(avunc)

    sibandav_rel = sib2.copy()
    for avunc in avunc2:
        sibandav_rel.add(avunc)


    all_seg_IBD1 = { chr : [] for chr in range(num_chrs) }
    all_seg_IBD2 = { chr : [] for chr in range(num_chrs) }
    has_seg_sib1 = [0 for x in range(len(sib1))]
    has_seg_sib2 = [0 for x in range(len(sib2))]
    has_seg_avunc1 = [0 for x in range(len(avunc1))]
    has_seg_avunc2 = [0 for x in range(len(avunc2))]
    sib1 = list(sib1)
    sib2 = list(sib2)
    avunc1 = list(avunc1)
    avunc2 = list(avunc2)
    for ind1 in sibandav:
        for ind2 in sibandav_rel:
            tmp = getIBDsegments(ind1, ind2, all_segs)
            for chr in range(num_chrs):  # add IBD1
                if len(tmp[0][chr]) > 0 or len(tmp[1][chr]) > 0:
                    # mark if these individuals have segments that were used
                    if ind1 in sib1:
                        has_seg_sib1[sib1.index(ind1)] = 1
                    elif ind1 in avunc1:
                        has_seg_avunc1[avunc1.index(ind1)] = 1
                    if ind2 in sib2:
                        has_seg_sib2[sib2.index(ind2)] = 1
                    elif ind2 in avunc2:
                        has_seg_avunc2[avunc2.index(ind2)] = 1
                    all_seg_IBD1[chr] += tmp[0][chr]
                    all_seg_IBD2[chr] += tmp[1][chr]

    IBD_sum = 0
    for chr in range(num_chrs):
        all_seg_IBD1[chr] = mergeIntervals(all_seg_IBD1[chr][:])
        for seg in all_seg_IBD1[chr]:
            IBD_sum += seg[1] - seg[0]
        all_seg_IBD2[chr] = mergeIntervals(all_seg_IBD2[chr][:])
        for seg in all_seg_IBD2[chr]:
            IBD_sum += 2.0*(seg[1] - seg[0])


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


def combineBothGPsKeepProportionOnlyExpectation(sib1, avunc1, pc1, sib2, avunc2, pc2, all_rel, all_segs, results_file, rel_graph):
# perform ancestral genome reconstruction between two groups of related individuals (sib1+avunc1 and sib2+avunc2)
# infers relatedness between all individuals within the two groups
    if len(sib1) == 1 and len(sib2) == 1 and len(avunc1) == 0 and len(avunc2) == 0:
        i1 = list(sib1)[0]
        i2 = list(sib2)[0]
        degree = getPairwiseD(i1, i2, all_rel)
        if i1 < i2:
            return [[i1,i2, degree, degree]]
        else:
            return [[i2,i1, degree, degree]]

    cont = 1
    for av1 in avunc1:
        if av1 in avunc2:
            cont = 0
    if not cont:
        deg = 'NA'
    else:
        # returns total length of genome IBD between sibandav and sibandav_rel, number of sibs in sib1 with IBD segments, number of sibs in sib2 with IBD segments
        [tmpsibav, sib1_len, sib2_len, av1_len, av2_len] = getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs)

        #get proportion of ancestor genome information expected on side 1
        if av1_len != 0:
            proportion_gp_exp = getExpectedGP(len(sib1),len(avunc1))
            proportion_par_exp = 0
        elif sib1_len > 1:
            proportion_gp_exp = 0
            proportion_par_exp = getExpectedPar(len(sib1))
        else:
            proportion_par_exp = 0
            proportion_gp_exp = 0

        #get proportion of ancestor genome information expectedo n side2
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
                i1 = list(sib1)[0]
                i2 = list(sib2)[0]
                K_exp = getPairwiseK(i1, i2, all_rel)
                base = 0

        estimated_exp = getInferredFromK(K_exp)

        IBD2=0
        if both and estimated_exp != 1 and K_exp > 1/2.0**(9.0/2):  #check for IBD2 if both sides are being reconstructed, might be reconstructing two sibs; K_exp is 3rd degree or closer
            if len(avunc1) and len(avunc2):
                sibseg = collectIBDsegments(avunc1, all_segs)
                sibsib = collectIBDsegmentsSibsAvuncularCombine(avunc1, avunc2, all_segs)
                IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, 0.5)
                if len(sib2) > 1:  # could be the case that we have one sib and his/her aunts/uncles
                    sibseg2 = collectIBDsegments(avunc2, all_segs)
                    sibsib2 = collectIBDsegmentsSibsAvuncularCombine(avunc2, avunc1, all_segs)
                    IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, 0.5)
                    for chr in range(num_chrs):
                        IBD011[chr] += IBD011_2[chr]
                        IBD011[chr] = mergeIntervals(IBD011[chr])
                IBD2 = getTotalLength(IBD011)
            elif not len(avunc1) and not len(avunc2):
                sibseg = collectIBDsegments(sib1, all_segs)
                sibsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, all_segs)
                IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, 0.5)
                if len(sib2) > 1: #could be the case that we have one sib and his/her aunts/uncles
                    sibseg2 = collectIBDsegments(sib2, all_segs)
                    sibsib2 = collectIBDsegmentsSibsAvuncularCombine(sib2, sib1, all_segs)
                    IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, 0.5)
                    for chr in range(num_chrs):
                        IBD011[chr] += IBD011_2[chr]
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
                estimated_exp = getInferredFromK(K_exp + (IBD2_add) / total_genome / 2.0)  # add in IBD2, remove IBD1

        result = []
        if estimated_exp != 0:
            if estimated_exp == -1:
                estimated_exp = base
            else:
                estimated_exp = estimated_exp + base

        #get all combinations of pairs for results
        for sib_rel in sib2:
            for sib in sib1:
                # if rel_graph.has_edge(sib_rel,sib):
                #     estimated_out_exp = rel_graph.get_edge_data(sib_rel,sib)['type']
                # else:
                estimated_out_exp = estimated_exp
                refined = getPairwiseD(sib_rel, sib, all_rel)
                to_add = [sib_rel, sib, estimated_out_exp, refined, 'combine1']
                result.append(to_add)

            for avunc in avunc1:
                # if rel_graph.has_edge(sib_rel, avunc):
                #     estimated_out_exp = rel_graph.get_edge_data(sib_rel,avunc)['type']
                # else:
                if estimated_exp != 0:
                    estimated_out_exp = estimated_exp - 1
                else:
                    estimated_out_exp = 0
                refined = getPairwiseD(sib_rel, avunc, all_rel)
                to_add = [sib_rel, avunc, estimated_out_exp, refined, 'combine2']
                result.append(to_add)

            # for p1 in pc1:
            #     if rel_graph.has_edge(sib_rel, p1):
            #         estimated_out_exp = rel_graph.get_edge_data(sib_rel,p1)['type']
            #     else:
            #         if estimated_exp != 0:
            #             estimated_out_exp = estimated_exp + 1
            #         else:
            #             estimated_out_exp = 0
            #     refined = getPairwiseD(sib_rel, p1, all_rel)
            #     to_add = [sib_rel, p1, estimated_out_exp, refined]
            #     result.append(to_add)

        for avunc_rel in avunc2:
            for sib in sib1:
                # if rel_graph.has_edge(avunc_rel, sib):
                #     estimated_out_exp = rel_graph.get_edge_data(avunc_rel,sib)['type']
                # else:
                if estimated_exp != 0:
                    estimated_out_exp = estimated_exp - 1
                else:
                    estimated_out_exp = 0
                refined = getPairwiseD(sib, avunc_rel, all_rel)
                to_add = [sib, avunc_rel, estimated_out_exp, refined, 'combine3']
                result.append(to_add)

            for avunc in avunc1:
                # if rel_graph.has_edge(avunc_rel, avunc):
                #     estimated_out_exp = rel_graph.get_edge_data(avunc_rel,avunc)['type']
                # else:
                if estimated_exp != 0:
                    estimated_out_exp = estimated_exp - 2
                else:
                    estimated_out_exp = 0
                refined = getPairwiseD(avunc, avunc_rel, all_rel)
                to_add = [avunc, avunc_rel, estimated_out_exp, refined, 'combine4']
                result.append(to_add)

            # for p1 in pc1:
            #     if rel_graph.has_edge(sib_rel, p1):
            #         estimated_out_exp = rel_graph.get_edge_data(sib_rel,p1)['type']
            #     else:
            #         if estimated_exp != 0:
            #             estimated_out_exp = estimated_exp + 1
            #         else:
            #             estimated_out_exp = 0
            #     refined = getPairwiseD(avunc_rel, p1, all_rel)
            #     to_add = [avunc_rel, p1, estimated_out_exp, refined]
            #     result.append(to_add)

        # for p1 in pc1:
        #     for p2 in pc2:
        #         if rel_graph.has_edge(p1, p2):
        #             estimated_out_exp = rel_graph.get_edge_data(p1, p2)['type']
        #         else:
        #             if estimated_exp != 0:
        #                 estimated_out_exp = estimated_exp + 2
        #             else:
        #                 estimated_out_exp = 0
        #         refined = getPairwiseD(p1, p2, all_rel)
        #         to_add = [p1, p2, estimated_out_exp, refined]
        #         result.append(to_add)

    return result

def checkRelevantAuntsUncles(sibset1, sibset2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel):
    # check whether aunt/uncle should be included in analysis
    avunc1 = []
    avunc2 = []
    if len(avunc1_bothsides):
        for i in range(0, len(avunc1_bothsides)):
            avunc1 += list(avunc1_bothsides[i])

    if len(avunc2_bothsides):
        for i in range(0, len(avunc2_bothsides)):
            avunc2 += list(avunc2_bothsides[i])


    for av1 in avunc1:
        if av1 in avunc2 or av1 in par2:
            return ['same',[], [], []] #same set of individuals
        if av1 in sibset2:
            return ['av',[], [], []] #sibset1's aunt/uncle is in sibset2 (sibset2 = aunts/uncles of sibset1)

    for av2 in avunc2: #
        if av2 in par1:
            return [[], 'same', [], []]
        if av2 in sibset1:
            return [[],'av', [], []]


    minsib = 50 # impossible value will be updated
    for s1 in sibset1:
        for s2 in sibset2:
            kinship = getPairwiseK(s1, s2, all_rel)
            if kinship < minsib:
                minsib = kinship


    # for av1 in avunc1:
    #     for av2 in avunc2:
    #         if getPairwiseD(av1, av2, all_rel) == 1:
    #             if getIBD1(av1, av2, all_rel) > 0.9:
    #                 return ['avparent',[av1,av2], [], []]
    #             elif getIBD2(av1, av2, all_rel) > 1/2.0**(3.0/2):
    #                 return ['sib', [av1,av2], [], []]


    relavunc1 = set()
    relavunc2 = set()
    for s1 in sibset1:
        for a2 in avunc2:
            # if getPairwiseD(s1, a2, all_rel) == 1:
            #     if getIBD1(s1, a2, all_rel) + getIBD2(s1, a2, all_rel) > 0.9: # a2 is parent of s1
            #         return ['sibparent',[s1,a2], [], []]
            kinship = getPairwiseK(s1, a2, all_rel)
            if kinship > minsib:
                relavunc2.add(a2)

    for s2 in sibset2:
        for a1 in avunc1:
            # if getPairwiseD(s2, a1, all_rel) == 1:
            #     if getIBD1(s2, a1, all_rel) + getIBD2(s2, a1, all_rel) > 0.9: # a1 is parent of s2
            #         return [[s2,a1],'sibparent', [], []]
            kinship = getPairwiseK(s2, a1, all_rel)
            if kinship > minsib:
                relavunc1.add(a1)

    if len(avunc1_bothsides):
        for i1 in relavunc1:
            for i in range(0,len(avunc1_bothsides)):
                if i1 in avunc1_bothsides[i]:
                    relavunc1 = relavunc1.union(avunc1_bothsides[i])

    if len(avunc2_bothsides):
        for i2 in relavunc2:
            for i in range(0,len(avunc2_bothsides)):
                if i2 in avunc2_bothsides[i]:
                    relavunc2 = relavunc2.union(avunc2_bothsides[i])



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
    for chr in range(num_chrs):
        for seg in IBD[chr]:
            total += seg[1] - seg[0]

    return total

def getAllRel(results_file, inds_file):
    # read in results file:
    # all_rel: dict of ind1, dict of ind2, list of [IBD1, IBD2, K, D
    # store pairwise relatedness information
    global inds
    first = [] #list of first degree relative pairs according to Refined IBD results
    second = [] #list of second degree relative pairs according to Refined IBD results
    third = [] #list of third degree relative pairs according to Refined IBD results
    inds = set()
    if inds_file != '':
        file = open(inds_file,'r')
        for line in file:
            l = str.split(line.rstrip())
            if len(l):
                inds.add(l[0])

        file.close()

    all_rel = {}
    file = open(results_file,'r')
    for line in file:
        l = str.split(line.rstrip())
        if inds_file == '':
            inds.add(l[0])
            inds.add(l[1])
        if l[0] in inds and l[1] in inds:
            ibd1 = float(l[2])
            ibd2 = float(l[3])
            K = ibd1/4.0 + ibd2/2.0
            degree = getInferredFromK(K)
            if l[0] < l[1]:
                ind1 = l[0]
                ind2 = l[1]
            else:
                ind1 = l[1]
                ind2 = l[0]
            if not ind1 in all_rel.keys():
                all_rel[ind1] = {} #IBD1, IBD2, K, D
            all_rel[ind1][ind2] = [ibd1,ibd2, K, degree]
            if degree == 1:
                first.append([ind1,ind2])
            elif degree == 2:
                second.append([ind1,ind2])
            elif degree == 3:
                third.append([ind1, ind2])


    file.close()


    for [ind1,ind2] in itertools.combinations(inds, 2):
        if ind1 < ind2:
            if not ind1 in all_rel.keys():
                all_rel[ind1] = {}
            if not ind2 in all_rel[ind1].keys():
                all_rel[ind1][ind2] = [0,0,0,0]
        else:
            if not ind2 in all_rel.keys():
                all_rel[ind2] = {}
            if not ind1 in all_rel[ind2].keys():
                all_rel[ind2][ind1] = [0,0,0,0]


    return [all_rel,inds,first,second,third]


def getSecondDegreeRelatives(rel_graph, second, sibset, par, all_rel):
    # collect all individuals we should check for being aunts/uncles of the sibset
    par = list(par)
    check_for_au = set()
    siblist = list(sibset)
    check_inds = set()
    for [ind1, ind2] in second:
        if not rel_graph.has_edge(ind1,ind2) or rel_graph.get_edge_data(ind1,ind2)['type'] in ['2','3']:
            if ind1 in sibset and not ind2 in par:
                check_inds.add(ind2)
            elif ind2 in sibset and not ind1 in par:
                check_inds.add(ind1)

    for ind in check_inds:
        if not ind in sibset:
            degs = set()
            for k in range(0,len(sibset)):
                degs.add( getPairwiseD(ind, siblist[k], all_rel) )
            if 2 in degs or 3 in degs:
                check_for_au.add(ind)

    return check_for_au




def getAuntsUncles_IBD011_nonoverlapping_pairs(sibset, halfsibs, second, all_rel, all_segs, rel_graph):
    # check whether individuals in list 'second' are likely aunts/uncles of 'sibset' and possibly also 'halfsibs'
    avunc = set()
    remove_second = set()
    avunc_hs_all = []
    if len(second):
        for ind in second:
            for sib in sibset:
                if getIBD2(ind, sib, all_rel) > 0.01 or (rel_graph.has_edge(ind,sib) and rel_graph.get_edge_data(ind,sib)['type'] in ['P','HS','GP']): #if proportion of genome shared IBD2 is > 1%
                    remove_second.add(ind)
                    break

        for ind in remove_second:
            second.remove(ind)

        second_original = list(second.copy()) #get copy for use with halfsibs later; we'll edit 'second' below
        second = list(second)
        sibset = list(sibset)
        if len(second):
            for [sib1, sib2] in itertools.combinations(sibset,2):
                sibseg = collectIBDsegments([sib1,sib2],all_segs)
                k = 0
                while k < len(second):
                    av = second[k]
                    avsib = collectIBDsegmentsSibsAvuncular([sib1,sib2], [av],all_segs)
                    IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, 0.5))
                    if IBD011 > 50:
                        avunc.add(av)
                        k = k + 1
                        # avsibs = getSibsFromGraph(rel_graph, av)
                        # second.remove(av)
                        # for avsib in avsibs:
                        #     avunc.add(avsib)
                        #     checkAndRemove(avsib, second)
                    elif IBD011 < 20:
                        second.remove(av)
                        checkAndRemove(av,avunc)
                        break
                    else:
                        k = k + 1

            #check with halfsibs
            second = second_original
            if len(halfsibs):
                avunc_hs = []
                for hs in range(0,len(halfsibs)): #hs = index of halfsib set
                    for [sib1,sib2] in itertools.product(sibset,halfsibs[hs]): #all pairs of [sib, halfsib]
                        sibseg = collectIBDsegments([sib1,sib2], all_segs)
                        k = 0
                        while k < len(second):
                            av = second[k]
                            avsib = collectIBDsegmentsSibsAvuncular([sib1,sib2], [av], all_segs)
                            IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, 0.5))
                            if IBD011 > 50:
                                avunc_hs.add(av)
                                # avsibs = getSibsFromGraph(rel_graph, av)
                                # second.remove(av)
                                # for avsib in avsibs:
                                #     possible_avunc_hs.add(avsib)
                                #     checkAndRemove(avsib, second)
                            elif IBD011 < 20:
                                second.remove(av)
                                checkAndRemove(av,avunc_hs)
                                break
                            else:
                                k = k + 1
                    if len(avunc_hs):
                        avunc_hs_all.append(avunc_hs)




    return [avunc, avunc_hs_all]


def checkUseHalfsibs(sibs, halfsib_sets, rel, all_rel):
    # check whether halfsibs should be included in analysis between 'sibs' and 'rel'
    # sibs = initial set of sibs
    # halfsib_sets = sibs' halfsibs
    # rel = distant relative

    sibs = list(sibs)

    if len(halfsib_sets):
        sibmin = 0
        for sib in sibs:
            kinship = getPairwiseK(sib, rel, all_rel)
            if kinship < sibmin or sib == sibs[0]:
                sibmin = kinship

        hsk = []
        hsk_all = []
        hsk_all_mean = []
        for hsset in halfsib_sets:
            hsmax = 0
            hsk_all_set = []
            for hs in hsset:
                kinship = getPairwiseK(hs, rel, all_rel)
                if kinship > hsmax:
                    hsmax = kinship
                hsk_all_set.append(kinship)

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


def checkSibs(sibs, rel, all_rel):
    #check if all siblings have same relatedness to rel
    all_pw = set()
    for sib in sibs:
        all_pw.add( getPairwiseD(sib, rel, all_rel) )
    if len(all_pw) == 1:
        return 1
    else:
        return 0

def checkAndRemove(x,setOrList):
    if x in setOrList:
        setOrList.remove(x)

def runDRUID(rel_graph, all_rel, inds, all_segs, args):
    # a chunky monkey
    all_results = []
    checked = set()
    for [ind1,ind2] in itertools.combinations(inds,2): #test each pair of individuals
        pair_name = getPairName(ind1, ind2)
        if not pair_name in checked: #if pair not yet tested
            #print("Comparing "+ind1+" and "+ind2)
            results = []
            [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
            [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
            sib1.add(ind1)
            sib2.add(ind2)
            if rel_graph.has_edge(ind1, ind2) and (rel_graph.get_edge_data(ind1,ind2)['type'] in ['FS','PC'] or ((len(sib1) > 1 and checkSibs(sib1,ind2, all_rel)) or len(sib1) == 1) and ((len(sib2) > 1 and checkSibs(sib2,ind1, all_rel)) or len(sib2) == 1)):
                #if the pair is already connected via graph, output that relationship
                checked.add(pair_name)
                refined = getPairwiseD(ind1, ind2, all_rel)
                type = rel_graph.get_edge_data(ind1, ind2)['type']
                if type == '1U':
                    type = '1'
                results.append([ind1, ind2, type, refined, 'graph1'])
            else:
                hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                sib1 = sib1.union(set(hs1))
                sib2 = sib2.union(set(hs2))

                reltype = getRelationship(rel_graph, ind1, ind2)
                if reltype != -1:
                    refined = getPairwiseD(ind1, ind2, all_rel)
                    results.append([ind1,ind2,reltype,refined, 'graph2'])
                else: #no path between individuals
                    #check if ind1 has parents/grandparents more closely related to other set of individuals
                    ind1_original = ind1
                    ind1_new = checkForMoveUp(all_rel, ind1, sib1, par1.union(gp1), pc1, sib2)
                    moves1 = []
                    moves_inds1 = []
                    if ind1_new == 'same': # same as ind2: shouldn't happen, but just in case
                        moves1.append(getRelationship(rel_graph, ind1_new, ind1))
                        moves_inds1.append(ind1)
                    while ind1 != ind1_new and ind1_new != 'same':
                        moves1.append(getRelationship(rel_graph,ind1_new, ind1))
                        moves_inds1.append(ind1)
                        ind1 = ind1_new
                        [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                        sib1.add(ind1)
                        ind1_new = checkForMoveUp(all_rel, ind1, sib1, gp1.union(par1), pc1, sib2)

                    if ind1_new == 'same': # ind1_new same as ind2?
                        same = anyIn(gp1.union(par1), sib2)
                        if len(same):
                            ind1 = same[0]
                        else:
                            ind1 = anyIn(pc1,sib2)[0]
                        [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                    else:
                        # check if ind2 has parents/grandparentsmore closely related to other set of individuals
                        ind2_original = ind2
                        ind2_new = checkForMoveUp(all_rel, ind2,sib2,gp2.union(par2), pc2, sib1)
                        moves2 = []
                        moves_inds2 = []
                        if ind2_new == 'same':
                            moves2.append(getRelationship(rel_graph, ind2_new, ind2))
                            moves_inds2.append(ind2)
                        while ind2 != ind2_new and ind2_new != 'same':
                            moves2.append(getRelationship(rel_graph, ind2_new, ind2))
                            moves_inds2.append(ind2)
                            ind2 = ind2_new
                            [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                            sib2.add(ind2)
                            ind2_new = checkForMoveUp(all_rel, ind2, sib2, gp2.union(par2), pc2, sib1)

                        if ind2_new == 'same':
                            same = anyIn(gp2.union(par2), sib1)
                            if len(same):
                                ind2 = same[0]
                            else:
                                ind2 = anyIn(pc2,sib1)[0]
                            [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)

                    # ind1 and ind2 can't be related via graph, otherwise they'd be considered above
                    # continue onto composite relatedness method

                    #switch focus to youngest generation if available
                    if ind1_new != 'same' and ind2_new != 'same':
                        if len(avunc1_bothsides) or len(avunc2_bothsides):
                            [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel)
                            # if relavunc1 == 'av':
                            #     results_tmp = []
                            #     for s1 in sib1:
                            #         for s2 in sib2:
                            #             refined = getPairwiseD(s1, s2, all_rel)
                            #             results_tmp.append([s1,s2,'A',refined,"graph"])
                            #     closest_result = [ind1,ind2,2,refined,"graph"] #refined may not be true Refined IBD inference for exact pair, but doesn't matter
                            # elif relavunc2 == 'av':
                            #     results_tmp = []
                            #     for s1 in sib1:
                            #         for s2 in sib2:
                            #             refined = getPairwiseD(s1, s2, all_rel)
                            #             results_tmp.append([s2,s1,'A',refined,"graph"])
                            #     closest_result = [ind1,ind2,2,refined,'graph']
                            # elif relavunc1 == 'sibparent':
                            #     results_tmp = []
                            #     for s1 in sib1:
                            #         for s2 in sib2:
                            #             refined = getPairwiseD(s1, s2, all_rel)
                            #             results_tmp.append([s1,s2,''])
                            #
                            # elif relavunc2 == 'sibparent':
                            #
                            # else:
                            if relavunc1 != 'av':
                                if not len(relavunc1) and len(nn1):
                                    # tmp = sib1[:]
                                    # sib1 = nn1[:]
                                    # avunc1_bothsides = [tmp[:]]
                                    # nn1 = []  # doesn't matter
                                    old_ind1 = ind1
                                    sib1 = getLargestSibsets(rel_graph, nn1)
                                    ind1 = list(sib1)[0]
                                    [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
                                    sib1.add(ind1)
                                    [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel)
                                    if relavunc1 != 'av' and not old_ind1 in relavunc1:
                                        relavunc1.add(old_ind1)
                                        checkAndRemove(old_ind1,unused1)
                                if not len(relavunc2) and len(nn2):
                                    # tmp = sib2[:]
                                    # sib2 = nn2[:]
                                    # avunc2_bothsides = [tmp[:]]
                                    # nn2 = []  # doesn't matter
                                    old_ind2 = ind2
                                    sib2 = getLargestSibsets(rel_graph,nn2)
                                    ind2 = list(sib2)[0]
                                    [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                                    sib2.add(ind2)
                                    [relavunc1, relavunc2, unused1, unused2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel)
                                    if relavunc2 != 'av' and not old_ind2 in relavunc2:
                                        relavunc2.add(old_ind2)
                                        checkAndRemove(old_ind2, unused2)

                                if len(unused1):
                                    unused_check = set()
                                    for i1 in unused1:
                                        if not i1 in unused_check:
                                            [sib1u, avunc1u_bothsides, nn1u, par1u, child1u, pc1u, gp1u, gc1u, halfsib1u_sets, twins1u] = pullFamily(rel_graph, i1)
                                            sib1u.add(i1)
                                            unused_check = unused_check.union(sib1u)
                                            results_to_add = combineBothGPsKeepProportionOnlyExpectation(sib1u, [], pc1, sib2, [], pc2, all_rel, all_segs, args.i[0], rel_graph)
                                            for item in results_to_add:
                                                this_pair = getPairName(item[0], item[1])
                                                if not this_pair in checked:
                                                    checked.add(this_pair)
                                                    item.append('inferred1')
                                                    results.append(item)
                                if len(unused2):
                                    unused_check = set()
                                    for i2 in unused2:
                                        if not i2 in unused_check:
                                            [sib2u, avunc2u_bothsides, nn2u, par2u, child2u, pc2u, gp2u, gc2u, halfsib2u_sets, twins2u] = pullFamily(rel_graph, i2)
                                            sib2u.add(ind2)
                                            unused_check = unused_check.union(sib2u)
                                            results_to_add = combineBothGPsKeepProportionOnlyExpectation(sib1, [], pc1, sib2u, [], pc2, all_rel, all_segs, args.i[0], rel_graph)
                                            for item in results_to_add:
                                                this_pair = getPairName(item[0], item[1])
                                                if not this_pair in checked:
                                                    checked.add(this_pair)
                                                    item.append('inferred2')
                                                    results.append(item)

                        else:
                            relavunc1 = []
                            relavunc2 = []

                        if relavunc1 != 'av' and relavunc2 != 'av':
                            results_tmp = combineBothGPsKeepProportionOnlyExpectation(sib1, relavunc1, pc1, sib2, relavunc2, pc2, all_rel, all_segs, args.i[0], rel_graph)
                        else:
                            # sibset1's aunt/uncle is in sibset2 (sibset2 = aunts/uncles of sibset1)
                            results_tmp = []
                        for resu in results_tmp:
                            this_pair = getPairName(resu[0], resu[1])
                            if not this_pair in checked:
                                resu.append('inferred3')
                                results.append(resu)
                                checked.add(this_pair)
                        if ind1_original != ind1 or ind2_original != ind2:
                            for res in results_tmp:
                                if (res[0] == ind1 and res[1] == ind2) or (res[0] == ind2 and res[1] == ind1):
                                    closest_result = res
                                    break
                            if closest_result[2] == '1U':
                                closest_result[2] = 1
                            elif closest_result[2] == 'A':
                                closest_result[2] = 2
                            total = int(closest_result[2])
                    if (ind1_new == 'same' or ind2_new == 'same'):
                        if ind1 == ind2: #catch bugs
                            closest_result = [ind1,ind2,0]
                        else: #siblings
                            closest_result = [ind1,ind2,1]

                    elif ind1_original != ind1 or ind2_original != ind2: #if we've travelled through the graph
                        for ii in range(len(moves1)-1,-1,-1):
                            #go through each move in moves1, add move length to total
                            #only add if total != 0 (i.e., there is a relationship)
                            if total > 0:
                                if moves1[ii] in ['P','C','PC']:
                                    total = total + 1
                                else: #gp or gc
                                    total = total + 2
                            (refined, this_pair) = getPairD_w_Name(moves_inds1[ii], ind2, all_rel)
                            if not this_pair in checked:
                                results.append([moves_inds1[ii],ind2,total,refined, 'graph+inferred1'])
                                checked.add(this_pair)
                            #check for close relatives of moves_inds[ii]
                            [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[ii])
                            for s1 in sib1:
                                (refined, this_pair) = getPairD_w_Name(s1, ind2, all_rel)
                                if not this_pair in checked:
                                    results.append([s1, ind2, total, refined, 'graph+inferred2'])
                                    checked.add(this_pair)
                            sib1.add(moves_inds1[ii])
                            hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                            for h1 in hs1:
                                (refined, this_pair) = getPairD_w_Name(h1, ind2, all_rel)
                                if not this_pair in checked:
                                    results.append([h1,ind2,total,refined, 'graph+inferred3'])
                                    checked.add(this_pair)
                            for p in pc1:
                                if not p in moves_inds1:  # if we didn't travel through this relationship already
                                    (refined, this_pair) = getPairD_w_Name(p, ind2, all_rel)
                                    if not this_pair in checked:
                                        if total > 0:
                                            results.append([p,ind2,total+1,refined,'graph+inferred4'])
                                        else:
                                            results.append([p,ind2,total,refined,'graph+inferred5'])
                                        checked.add(this_pair)

                        total = int(closest_result[2])
                        [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1) #get set of close relatives of ind1
                        sib1.add(ind1)
                        for ii in range(len(moves2)-1,-1,-1):
                            if total > 0:
                                if moves2[ii] in ['P','C','PC']:
                                    total = total + 1
                                else: #gp or gc
                                    total = total + 2
                            for s1 in sib1:
                                (refined, this_pair) = getPairD_w_Name(moves_inds2[ii], s1, all_rel)
                                if not this_pair in checked:
                                    results.append([s1,moves_inds2[ii],total,refined, 'graph+inferred6'])
                                    checked.add(this_pair)
                                #check for close relatives of moves_inds[ii]
                                [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[ii])
                                for s2 in sib2:
                                    (refined, this_pair) = getPairD_w_Name(s1, s2, all_rel)
                                    if not this_pair in checked:
                                        results.append([s1,s2,total,refined, 'graph+inferred7'])
                                        checked.add(this_pair)
                                sib2.add(moves_inds2[ii])
                                hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                                for h2 in hs2:
                                    (refined, this_pair) = getPairD_w_Name(h2, s1, all_rel)
                                    if this_pair in checked:
                                        results.append([s1,h2,total,refined,'graph+inferred8'])
                                        checked.add(this_pair)
                                for p in pc2:
                                    if not p in moves_inds2 and not p == ind2: #if we didn't travel through this relationship already
                                        (refined, this_pair) = getPairD_w_Name(p, s1, all_rel)
                                        if not this_pair in checked:
                                            if total > 0:
                                                results.append([s1, p, total+1, refined, 'graph+inferred9'])
                                            else:
                                                results.append([s1, p, total, refined, 'graph+inferred9'])
                                            checked.add(this_pair)

                        if len(moves1) and len(moves2):
                            total = int(closest_result[2])
                            for i1 in range(len(moves1)-1,-1,-1):
                                if total > 0:
                                    if moves1[i1] in ['P', 'C', 'PC']:
                                        total = total + 1
                                    else:
                                        total = total + 2
                                for i2 in range(len(moves2)-1,-1,-1):
                                    if total > 0:
                                        if moves2[i2] in ['P', 'C', 'PC']:
                                            total = total + 1
                                        else:
                                            total = total + 2
                                    (refined, this_pair) = getPairD_w_Name(moves_inds1[i1], moves_inds2[i2], all_rel)
                                    if not this_pair in checked:
                                        results.append([moves_inds1[i1],moves_inds2[i2],total,refined, 'graph+inferred10'])
                                        checked.add(this_pair)

                                    # check for close relatives of moves_inds[ii]
                                    [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[i1])
                                    [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[i2])
                                    for s1 in sib1:
                                        (refined, this_pair) = getPairD_w_Name(s1, moves_inds2[i2], all_rel)
                                        if not this_pair in checked:
                                            results.append([s1, moves_inds2[i2], total, refined, 'graph+inferred11'])
                                            checked.add(this_pair)
                                    for s2 in sib2:
                                        (refined, this_pair) = getPairD_w_Name(s2, moves_inds1[i1], all_rel)
                                        if not this_pair in checked:
                                            results.append([moves_inds1[i1], s2, total, refined, 'graph+inferred12'])
                                            checked.add(this_pair)
                                        for s1 in sib1:
                                            (refined, this_pair) = getPairD_w_Name(s1, s2, all_rel)
                                            if not this_pair in checked:
                                                results.append([s1, s2, total, refined, 'graph+inferred13'])
                                                checked.add(this_pair)
                                    sib1 = list(sib1)
                                    sib2 = list(sib2)
                                    sib1.append(moves_inds1[i1])
                                    sib2.append(moves_inds2[i2])
                                    hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                                    hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                                    for h1 in hs1:
                                        this_pair = getPairName(h1, ind2)
                                        if not this_pair in checked:
                                            results.append([h1, ind2, total, refined, 'graph+inferred14'])
                                            checked.add(this_pair)
                                    for h2 in hs2:
                                        this_pair = getPairName(h2, ind1)
                                        if not this_pair in checked:
                                            results.append([ind1, h2, total, refined, 'graph+inferred15'])
                                            checked.add(this_pair)
                                        for h1 in hs1:
                                            this_pair = getPairName(h1, h2)
                                            if not this_pair in checked:
                                                results.append([h1,h2,total,refined, 'graph+inferred16'])
                                                checked.add(this_pair)


                if len(twins1):
                    for res in results:
                        if sibs1[0] == res[0]:
                            for t1 in twins1:
                                results.append([t1,res[1],res[2],res[3],'graph3'])
                                this_pair = getPairName(t1, res[1])
                                checked.add(this_pair)
                        elif sibs1[0] == res[1]:
                            for t1 in twins1:
                                results.append([res[0],t1,res[2],res[3],'graph4'])
                                this_pair = getPairName(t1, res[0])
                                checked.add(this_pair)

                if len(twins2):
                    for res in results:
                        if sibs2[0] == res[0]:
                            for t2 in twins2:
                                results.append([t2, res[1], res[2], res[3],'graph5'])
                                this_pair = getPairName(t2, res[1])
                                checked.add(this_pair)
                        elif sibs2[0] == res[1]:
                            for t2 in twins2:
                                results.append([res[0], t2, res[2], res[3],'graph6'])
                                this_pair = getPairName(t2, res[0])
                                checked.add(this_pair)

            all_results += results

    return all_results
