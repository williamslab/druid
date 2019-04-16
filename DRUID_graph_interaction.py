import networkx as nx
from DRUID_functions import *
from DRUID_all_rel import *

def checkChangeLineage(tmp_graph,path):
    if len(path) > 2:
        for i in range(0,len(path)-2):
            indfirst = path[i]
            indsecond = path[i+1]
            indthird = path[i+2]
            type1 = tmp_graph.get_edge_data(indfirst, indsecond)['type']
            type2 = tmp_graph.get_edge_data(indsecond, indthird)['type']
            if (type1 == 'PC' and type2 =='PC') or (type1 == 'P' and type2 == 'PC') or (type1 == 'P' and type2 == 'C') or (type1 == 'GP' and type2 == 'GC') or (type1 == 'P' and type2 == 'GC') or (type1 == 'GP' and type2 == 'C') or (type1 == 'PC' and type2 == 'P') or (type1 == 'C' and type2 == 'P') or (type1 == 'GC' and type2 == 'GP') or (type1 == 'C' and type2 == 'GP') or (type1 == 'GC' and type2 == 'P') or (type1 == 'C' and type2 == 'AU') or (type1 == 'GC' and type2 == 'AU') or (type1 == 'P' and type2 == 'NN') or (type1 == 'GP' and type2 == 'NN'):
                return True
        return False
    else:
        return False

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
                    # check if we are (possibly) traveling to other lineage, and if so, stop
                    if (type1 == 'P' and type2 == 'C') or (type1 == 'GP' and type2 == 'GC') or (type1 == 'P' and type2 == 'GC') or (type1 == 'GP' and type2 == 'C') or (type1 == 'C' and type2 == 'P') or (type1 == 'GC' and type2 == 'GP') or (type1 == 'C' and type2 == 'GP') or (type1 == 'GC' and type2 == 'P') or (type1 == 'C' and type2 == 'AU') or (type1 == 'GC' and type2 == 'AU') or (type1 == 'P' and type2 == 'NN') or (type1 == 'GP' and type2 == 'NN') or (type1 == 'NN' and type2 == 'AU') or (type1 == 'AU' and type2 == 'NN'):
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
                    elif type1 == 'PC': #added 2/16/18
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
            [sib, avunc_bothsides, nn, par, child, pc, gp, gc, halfsib_sets, twins] = pullFamily(tmp_graph, ind)
            sib.add(ind)
            checked = checked.union(sib)
            sibsets.append(sib)

    sizes = [len(x) for x in sibsets]
    return sibsets[sizes.index(max(sizes))]



def checkIfSib(tmp_graph,ind1,ind2):
    #determine whether there's an edge labeled 'FS' between ind1 and ind2
    return tmp_graph.get_edge_data(ind1,ind2)['type'] == 'FS'


def checkIfParentInGraph(tmp_graph,ind1,ind2):
    #determine whether ind1 is the child of ind2
    return tmp_graph.get_edge_data(ind1,ind2)['type'] == 'C'


def checkIfParent(tmp_graph, all_rel, sibset, ind, C):
    #determine if ind is the parent of sibset
    if len(sibset) == 1:
        if C: #DRUID_C
            return False
        else:
            sib = list(sibset)[0]
            if getIBD2(sib, ind, all_rel) < 0.05: #very little IBD2
                #return True to give generic PC categorization to pair later
                return True
            else:
                return False
    else:
        for sib in sibset:
            if not (sib,ind) in tmp_graph.edges(sib) or not tmp_graph.get_edge_data(sib,ind)['type'] == '1U' or not float(getIBD2(sib, ind, all_rel)) < 1/2.0**(7/2.0):
                return False

    return True


def getSibsFromGraph(tmp_graph,ind):
    #return the siblings of ind
    sibs = set()
    if ind in tmp_graph:
        neighbors = tmp_graph.neighbors(ind)
        for x in neighbors:
            if checkIfSib(tmp_graph,ind,x):
                sibs.add(x)
    return sibs

def getParent(tmp_graph,ind):
    par = set()
    if ind in tmp_graph:
        neighbors = tmp_graph.neighbors(ind)
        for x in neighbors:
            if checkIfParentInGraph(tmp_graph,ind,x):
                par.add(x)

    return par


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
    hs = list(hs)
    while len(hs):
        hssibs = getSibsFromGraph(tmp_graph, hs[0])
        halfsib_set = [hs[0]]
        for ind_sib in hssibs:
            halfsib_set.append(ind_sib)
            hs.remove(ind_sib)
        hs.remove(hs[0])
        halfsib_sets.append(halfsib_set)

    return [sibs, halfsib_sets]


def getSibsParentsFromGraph(tmp_graph,ind):
    #return the siblings and half-sibs of ind
    par = set()
    sibs = set()
    neighbors = tmp_graph.neighbors(ind)
    for x in neighbors:
        if tmp_graph.get_edge_data(ind,x)['type'] == 'FS':
            sibs.add(x)
        elif tmp_graph.get_edge_data(ind,x)['type'] == 'C':
            par.add(x)

    return [sibs, par]

def getSibsHalfSibsParentsFromGraph(tmp_graph,ind):
    #return the siblings and half-sibs of ind
    hs = set()
    sibs = set()
    par = set()
    neighbors = tmp_graph.neighbors(ind)
    for x in neighbors:
        if tmp_graph.get_edge_data(ind,x)['type'] == 'HS':
            hs.add(x)
        elif tmp_graph.get_edge_data(ind,x)['type'] == 'FS':
            sibs.add(x)
        elif tmp_graph.get_edge_data(ind,x)['type'] == 'C':
            par.add(x)

    halfsib_sets = []
    hs = list(hs)
    while len(hs):
        hssibs = getSibsFromGraph(tmp_graph, hs[0])
        halfsib_set = [hs[0]]
        for ind_sib in hssibs:
            halfsib_set.append(ind_sib)
            hs.remove(ind_sib)
        hs.remove(hs[0])
        halfsib_sets.append(halfsib_set)

    return [sibs, halfsib_sets, par]



def checkSiblingSubgraph(tmp_graph,siblings,C):
    # check whether all nodes in a sibling subgraph have a direct edge between one-another;
    # if not, greedily remove the sibling node with the least direct edges connecting to other siblings until they're all connected

    #first, get all sibs (check neighbors) of given set 'siblings' and add to 'siblings'
    remove = []
    all_sibs = []
    for ind in siblings: #for each sibling, collect his/her siblings currently in graph
        sibs_ind = getSibsFromGraph(tmp_graph,ind)
        all_sibs = all_sibs + list(sibs_ind)
    #get sibs that weren't included in given set 'siblings'
    sibs_add = set(all_sibs).difference(siblings)
    #update 'siblings' to include those sibs
    siblings = siblings.union(sibs_add)

    #second, count the number of times each sib missing from the initial 'siblings' set is considered a sibling of another sib
    for ind in sibs_add:
        sibs_ind = getSibsFromGraph(tmp_graph,ind)
        all_sibs = all_sibs + list(sibs_ind)
    counts = [all_sibs.count(x) for x in all_sibs] #count number of times each reported sibling appears

    if not C: #if not DRUID_C
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
    else: #if DRUID_C
        while not all(x == len(set(all_sibs)) for x in counts):  # if any individual in all_sibs is found to be a sibling with less than half of the other siblings
            to_remove = getElementThatAppearsLeast(all_sibs)
            remove.append(to_remove)
            if to_remove in siblings:
                siblings.remove(to_remove)
            all_sibs = []
            for ind in siblings:
                sibs_ind = getSibsFromGraph(tmp_graph, ind)
                for x in sibs_ind:
                    if not x in remove:
                        all_sibs.append(x)
            counts = [all_sibs.count(x) for x in all_sibs]

    return [list(set(all_sibs)), remove]

def addEdgeType(ind1,ind2,type1,type2,rel_graph):
    if not rel_graph.has_edge(ind1, ind2):
        rel_graph.add_edge(ind1,ind2)
        rel_graph.add_edge(ind2,ind1)
    rel_graph[ind1][ind2]['type'] = type1
    rel_graph[ind2][ind1]['type'] = type2


def mean(nums):
    return sum(nums)/len(nums)

def thresholdK(pcD):
    return 0.75-0.025*pcD



def checkForMoveUp(all_rel, ind, sibset, older_gen, possible_par, third_party):
    # check if parent/grandparent is more related to third_party
    if len(older_gen):
        if len(older_gen.intersection(third_party)):
            return 'same'

        maxsibK = -1 # impossible value will be updated
        for sib in sibset:
            for tp in third_party:
                maxsibK = max(maxsibK, getPairwiseK(sib, tp, all_rel) )

        maxparK = -1
        the_max_par = -1
        for par in older_gen:
            for tp in third_party:
                parK = getPairwiseK(par, tp, all_rel)
                if parK > maxparK:
                    maxparK = parK
                    the_max_par = par

        if maxparK > maxsibK:
            # there's a parent/grandparent more closely related
            return the_max_par

    if len(possible_par):
        if len(possible_par.intersection(third_party)):
            return 'same'

        # check if any possible parents are more closely related
        max_pc_K = -1
        the_max_pc = -1
        found_max = False
        for pc in possible_par:
            should_consider_pc = True
            cur_max_pc_K = -1
            for tp in third_party:
                pcD = getPairwiseD(tp, pc, all_rel)
                pcK = getPairwiseK(tp, pc, all_rel)
                indK =  getPairwiseK(tp, ind, all_rel)

                if pcD < 0 or (thresholdK(pcD) * pcK <= indK or pcK / (indK + 1e-6) >= 20):
                    # Ensure that for each third party individual, pc's K is
                    # sufficiently larger than current individual's K.
                    # Because this person could be a child of ind, who is
                    # potentially related to tp through a different lineage,
                    # we impose an upper bound on the ratio of relatedness
                    # between pcK and indK: it must be < 20x higher
                    # (expectation is 2x)
                    should_consider_pc = False
                    break
                elif pcK > cur_max_pc_K:
                    cur_max_pc_K = pcK

            if should_consider_pc and cur_max_pc_K > max_pc_K:
                max_pc_K = cur_max_pc_K
                the_max_pc = pc
                found_max = True


        if found_max:
            return the_max_pc # use PC with largest K

    return ind


def checkAuntUncleGPRelationships(tmp_graph,siblings,par):
    # ensure the siblings of 'par' are listed as aunts/uncles of 'siblings' (par = parent of siblings)
    for p in par:
        [sibpar, avunc_bothsides, nn, parpar, childpar, pc, gppar, gcpar, halfsib_sets, twins] = pullFamily(tmp_graph, p)
        for sib in siblings:
            for sp in sibpar:
                addEdgeType(sib,sp,'NN','AU',tmp_graph)
        for pp in parpar:
            for sib in siblings:
                addEdgeType(sib,pp,'GC','GP',tmp_graph)


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
    grandchildren = set()
    halfsibs = set()
    nn = set()
    twins = set()
    pc = set()
    for edge in edges:
        edge_info = tmp_graph.get_edge_data(edge[0],edge[1])['type']
        if edge[0] == ind:
            if edge_info == 'P':
                children.add(edge[1])
            elif edge_info == 'C':
                parents.add(edge[1])
            elif edge_info == 'FS':
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
            elif edge_info == 'PC':
                pc.add(edge[1])
            elif edge_info == 'GP':
                grandchildren.add(edge[1])
    avunc_sets = []
    checked = set()
    for ind in avunc:
        if not ind in checked:
            av_sibs = getSibsFromGraph(tmp_graph, ind)
            av_sibs.add(ind)
            avunc_sets.append(av_sibs)
            checked = checked.union(av_sibs)



    halfsib_sets = []
    halfsibs = list(halfsibs)
    while len(halfsibs):
        hssibs = getSibsFromGraph(tmp_graph, halfsibs[0])
        halfsib_set = [halfsibs[0]]
        for hs in hssibs:
            halfsib_set.append(hs)
            halfsibs.remove(hs)
        halfsibs.remove(halfsibs[0])
        halfsib_sets.append(halfsib_set)

    return [sib, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins]


def getAllCloseRelationships(tmp_graph):
    fam = 1
    checked = set()
    for node in tmp_graph:
        if not node in checked:
            switches = []
            missing = 1
            # get close relatives
            [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph,node)

            # 'climb' pedigree by switching to grandparents or parents if available
            current_node = node
            while len(grandparents+parents+avunc_sets):
                if len(grandparents):
                    i = 0
                    while i < len(grandparents):
                        if not list(grandparents)[i] in checked:
                            switches.append(current_node)
                            current_node = list(grandparents)[i]
                            [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                            break
                        i += 1

                elif len(parents): #switch to parent
                    i = 0
                    while i < len(parents):
                        if not parents[i] in checked:
                            switches.append(current_node)
                            current_node = parents[i]
                            [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                            break
                        i += 1
                elif len(avunc_sets):
                    i = 0
                    while i < len(avunc_sets):
                        ii = 0
                        while ii < len(avunc_sets[i]):
                            if not avunc_sets[i][ii] in checked:
                                switches.append(current_node)
                                current_node = avunc_sets[i][ii]
                                [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                                break
                            ii += 1
                        i += 1
                else:
                    break #we've already considered all parents/grandparents

            # add all close relatives to list of checked individuals
            checked = checked.union(sibs,nn,parents,children,pc,grandparents,halfsib_sets,twins)
            for av in range(0,len(avunc_sets)):
                checked = checked.union(av)
            for hs in range(0,len(halfsib_sets)):
                checked = checked.union(hs)

            # add in missing parents
            if len(parents) == 1:
                parents.add(str(fam)+'_missing'+str(missing)) # e.g. '1_missing1' for missing individual 1 in fam 1
            if not len(parents):
                parents.add(str(fam)+'_missing'+str(missing))
                parents.add(str(fam)+'_missing'+str((missing+1)))
                missing += 2
            outfile.write(str(fam)+'\t'+current_node+'\t'+parents[0]+'\t'+parents[1]+'\n') #this individual's fam info

            for sib in sibs: #siblings
                outfile.write(str(fam)+'\t'+sib+'\t'+parents[0]+'\t'+parents[1]+'\n')

            # add nieces and nephews to outfile
            checked_nn = set()
            for n in nn: #nieces and nephews
                if not n in checked_nn:
                    [n_sibs,n_par] = getSibsParentsFromGraph(tmp_graph,n)
                    n_sibs.add(n)
                    if len(n_par) == 2: #sibling of current_node already in graph
                        for ns in n_sibs:
                            outfile.write(str(fam)+'\t'+ns+'\t'+n_par[0]+'\t'+n_par[1]+'\n')
                    elif len(n_par) == 1:
                        if n_par[0] in sibs: #sibling of current node already in graph
                            for ns in n_sibs:
                                outfile.write(str(fam)+'\t'+ns+'\t'+n_par[0]+'\t'+str(fam)+'_missing'+str(missing)+'\n')
                            missing += 1
                        else: #sibling of current node not already in graph
                            for ns in n_sibs:
                                outfile.write(str(fam) + '\t' + ns + '\t' + str(fam) + '_missing' + str(missing) + '\t' + str(fam) + '_missing' + str(missing+1) + '\n')
                            missing += 2
                    else:
                        for ns in n_sibs:
                            outfile.write(str(fam) + '\t' + ns + '\t' + str(fam) + '_missing' + str(missing) + '\t' + str(fam) + '_missing' + str(missing + 1) + '\n')
                        missing += 2
                    checked_nn = checked_nn.union(n_sibs)

            # add children to outfile
            for ch in children:
                outfile.write()

            if len(avunc_sets):
                if len(grandparents) == 1:
                    parent_gp
            for av_set in avunc_sets:
                for av in av_set:
                    outfile.write(str(fam))



def moveUpForFillIn(node,checked,tmp_graph):
    switches = []
    [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, node)
    current_node = node
    past_node = ''
    while (len(grandparents) or len(parents) or len(avunc_sets)) and not past_node == current_node:
        past_node = current_node
        if len(grandparents):
            grandparents = list(grandparents)
            i = 0
            while i < len(grandparents):
                print('grandparents')
                if not grandparents[i] in checked:
                    switches.append(current_node)
                    current_node = grandparents[i]
                    [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                    break
                i += 1

        elif len(parents):  # switch to parent
            print('parents')
            i = 0
            parents = list(parents)
            while i < len(parents):
                if not parents[i] in checked:
                    switches.append(current_node)
                    current_node = parents[i]
                    [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                    break
                i += 1
        elif len(avunc_sets):
            i = 0
            while i < len(avunc_sets):
                print('avunc')
                ii = 0
                avset = list(avunc_sets[i])
                while ii < len(avset):
                    if not avset[ii] in checked:
                        switches.append(current_node)
                        current_node = avset[ii]
                        [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                        break
                    ii += 1
                i += 1
        else:
            break  # we've already considered all parents/grandparents/aunts/uncles

    return [current_node,switches, sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins]

def fillInMissingParentGPFromAv(this_fam, grandparents, parents, avunc_set, av_par, missing):
    #compare grandparents in graph with parents of avuncular in graph, determine if the two lineages are the same
    #fillInMissingParentGPFromAv(this_fam, grandparents, avunc_sets[0], av1_par, missing)
    if len(grandparents): #if grandparents are available, check whether they are the parents of the av set
        common_inds = list(set(grandparents) & set(av_par))
        if not len(common_inds):  # none of the grandparents of current_node are the parents of the aunt/uncle set
            for avs in avunc_set:  # all aunts/uncles in this avunc set
                #addEdgeType(avs, 'missing_' + str(missing), 'C', 'P', this_fam)
                #addEdgeType(avs, 'missing_' + str(missing + 1), 'C', 'P', this_fam)
                addEdgeType(avs, 'missing_' + str(missing - 2), 'FS', 'FS', this_fam)  # add edge between avunc_set and the missing parent (their sib)
            addEdgeType(av_par[0], 'missing_' + str(missing - 2), 'P', 'C', this_fam)  # add edge between missing parent and his/her parent
            addEdgeType(av_par[1], 'missing_' + str(missing - 2), 'P', 'C', this_fam)  # add edge between missing parent and his/her parent
            missing += 2
        else:  # a parent/parents of the aunt/uncle set are in the grandparent set
            if len(common_inds) == 2:  # both parents of missing parent are in dataset
                addEdgeType('missing_' + str(missing - 2), common_inds[0], 'C', 'P', this_fam)  # add edge between missing parent and his/her parent
                addEdgeType('missing_' + str(missing - 2), common_inds[1], 'C', 'P', this_fam)  # add edge between missing parent and his/her parent
                for avs in avunc_sets:
                    addEdgeType('missing_' + str(missing - 2), avs, 'FS', 'FS', this_fam)  # add edge between missing parent and his/her parent
            elif len(common_inds) == 1:  # only one parent of missing parent is in dataset
                addEdgeType('missing_' + str(missing - 2), common_inds[0], 'C', 'P', this_fam)
                addEdgeType('missing_' + str(missing - 2), 'missing_' + str(missing), 'C', 'P', this_fam)  # create node for missing grandparent
                for avs in avunc_set:
                    addEdgeType(avs, 'missing_' + str(missing), 'C', 'P', this_fam)
                missing += 1
    #else, avset could be either lineage
    else:
        for avs in avunc_set:  # all aunts/uncles in this avunc set
            # addEdgeType(avs, 'missing_' + str(missing), 'C', 'P', this_fam)
            # addEdgeType(avs, 'missing_' + str(missing + 1), 'C', 'P', this_fam)
            addEdgeType(avs, list(parents)[0], 'FS', 'FS', this_fam)  # add edge between avunc_set and the missing parent (their sib)
        addEdgeType(list(av_par)[0], list(parents)[0], 'P', 'C', this_fam)  # add edge between missing parent and his/her parent
        addEdgeType(list(av_par)[1], list(parents)[0], 'P', 'C', this_fam)  # add edge between missing parent and his/her parent



def fillInGraph(tmp_graph):
    # add missing individuals as nodes to graph
    checked = set()
    fam = 1
    for node in tmp_graph.nodes():
        if not node in checked:
            this_fam = nx.ego_graph(tmp_graph, node)
            #keep a list of all tree possibilities
            all_trees_this_fam = []

            #iterate through network of close relatives, filling in missing individuals
            missing = 1  # keep track of what missing individual # we're on
            for fam_node in this_fam.nodes():
                if not fam_node in checked and not 'missing' in fam_node:
                    switches = [] #keep track of what switches we make as we 'climb' the pedigree

                    # 'climb' pedigree by switching to grandparents or parents if available; get relatives of 'current_node'
                    [current_node, switches, sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = moveUpForFillIn(fam_node,checked,this_fam)
                    sibs.add(current_node)
                    switches.append(current_node)

                    while len(switches):
                        current_node = switches.pop()
                        [sibs, avunc_sets, nn, parents, children, pc, grandparents, grandchildren, halfsib_sets, twins] = pullFamily(tmp_graph, current_node)
                        sibs.add(current_node)
                        checked = checked.union(sibs)

                        ## add missing parents ##
                        #if len(parents) == 2, then the aunts/uncles and grandparents should be taken care of already
                        if len(parents) == 1:
                            for hset in halfsib_sets:
                                hs = hset[0]
                                [hs_sibs, hs_par] = getSibsParentsFromGraph(this_fam, hs)
                                if any([x in parents for x in hs_par]): #if the parent is shared between sibset and halfsib set
                                    for h in hs:
                                        addEdgeType(h, 'missing_' + str(missing), 'C', 'P', this_fam)
                                        missing += 1
                            for s in sibs.union(twins):
                                addEdgeType(s, 'missing_'+str(missing), 'C', 'P', this_fam)
                            missing += 1

                            # connect missing parent to relevant aunts/uncles, if any
                            if len(avunc_sets):
                                if len(avunc_sets) == 1:
                                    if not list(parents)[0] in avunc_sets[0]: #the available parent isn't in the avunc_set --> the missing parent must be avunc_set's sib
                                        for avs in avunc_sets[0]:
                                            addEdgeType(avs, 'missing_' + str(missing-1), 'FS', 'FS', this_fam)
                                else:
                                    if not list(parents)[0] in avunc_sets[0]: #the available parent must be sib of second avunc set
                                        for avs in avunc_sets[0]:
                                            addEdgeType(avs, 'missing_' + str(missing-1), 'C', 'P', this_fam)
                                    else: #the available parent must be sib offirst avunc set
                                        for avs in avunc_sets[1]:
                                            addEdgeType(avs, 'missing_' + str(missing-1), 'C', 'P', this_fam)



                        elif len(parents) == 0:
                            for s in sibs.union(twins):
                                addEdgeType(s, 'missing_' + str(missing), 'C', 'P', this_fam)
                                addEdgeType(s, 'missing_' + str(missing+1), 'C', 'P', this_fam)
                            parents.add('missing_'+str(missing))
                            parents.add('missing_'+str(missing+1))
                            missing += 2

                            # if there are any avunc_sets, ensure they're connected to their respective siblings (missing parents)
                            # also add in any missing grandparents of current_node
                            if len(avunc_sets) == 2: #we have avunculars through both sides of family
                                [av1_sibs, av1_par] = getSibsParentsFromGraph(this_fam, list(avunc_sets[0])[0])
                                [av2_sibs, av2_par] = getSibsParentsFromGraph(this_fam, list(avunc_sets[1])[0])
                                fillInMissingParentGPFromAv(this_fam, grandparents, parents, avunc_sets[0], av1_par, missing)
                                fillInMissingParentGPFromAv(this_fam, grandparents, parents, avunc_sets[1], av2_par, missing)
                                missing += 2
                                if len(halfsib_sets):
                                    for hsset in halfsib_sets:
                                        hs = list(hsset)
                                        [hs_sibs, hs_avunc_sets, hs_nn, hs_parents, hs_children, hs_pc, hs_grandparents, hs_grandchildren, hs_halfsib_sets, hs_twins] = pullFamily(this_fam, hs[0])
                                        if len(hs_avunc_sets):
                                            if list(hs_avunc_sets[0])[0] in avunc_sets[0]: #if the halfsib avunc set is the same as the sib avunc set, the missing parent is also parent of hs set
                                                for halfsib in hs:
                                                    addEdgeType(halfsib, 'missing_' + str(missing-2), 'C', 'P', this_fam)
                                            elif len(hs_avunc_sets) == 2 and list(hs_avunc_sets[0])[0] in avunc_sets[1]:
                                                for halfsib in hs:
                                                    addEdgeType(halfsib, 'missing_' + str(missing-1), 'C', 'P', this_fam)
                            elif len(avunc_sets) == 1:
                                [av1_sibs, av1_par] = getSibsParentsFromGraph(this_fam, list(avunc_sets[0])[0])
                                fillInMissingParentGPFromAv(this_fam, grandparents, parents, avunc_sets[0], av1_par, missing) #add the grandparents of the sibling sets if missing
                                missing += 1
                                if len(halfsib_sets):
                                    for hsset in halfsib_sets:
                                        hs = list(hsset)
                                        [hs_sibs, hs_avunc_sets, hs_nn, hs_parents, hs_children, hs_pc, hs_grandparents, hs_grandchildren, hs_halfsib_sets, hs_twins] = pullFamily(this_fam, hs[0])
                                        if len(hs_avunc_sets):
                                            if any(x in avunc_sets[0] for x in hs_avunc_sets[0]):  # if the halfsib avunc set is the same as the sib avunc set, the missing parent is also parent of hs set
                                                for halfsib in hs+list(sibs):
                                                    addEdgeType(halfsib, 'missing_' + str(missing - 2), 'C', 'P', this_fam)
                                            elif len(hs_avunc_sets) == 2 and any(x in avunc_sets[0] for x in hs_avunc_sets[1]):
                                                for halfsib in hs:
                                                    addEdgeType(halfsib, 'missing_' + str(missing - 1), 'C', 'P', this_fam)

                            else: #no aunts/uncles; for halfsibs, just share a single parent with full sibs
                                if len(halfsib_sets):
                                    for hsset in halfsib_sets:
                                        hs = list(hsset)
                                        for halfsib in hsset:
                                            addEdgeType(halfsib, 'missing_'+str(missing-2),'C','P',this_fam)



            printFam(this_fam,fam,'fam_'+str(fam)+'.fam')
            fam = fam + 1





def printFam(this_fam,fam,outfile):
    checkedfam = set()
    output = open(outfile,'w')
    for node in this_fam.nodes():
        if not 'missing' in node and not node in checkedfam:
            [sibs, parents] = getSibsParentsFromGraph(this_fam,node)
            output.write(fam+'\t'+node+'\t'+list(parents)[0]+'\t'+list(parents)[1]+'\n')
            for sib in sibs:
                output.write(fam+'\t'+sib + '\t' + list(parents)[0] + '\t' + list(parents)[1] + '\n')
            checkedfam = checkedfam.union(sibs)
            checkedfam.add(node)

    output.close()



            #
            #
            #
            # #add missing siblings according to nieces/nephews
            # checked_nn = set()
            # for n in nn: #nieces and nephews
            #     if not n in checked_nn:
            #         [n_sibs,n_par] = getSibsParentsFromGraph(this_fam,n)
            #         n_sibs.add(n)
            #         if len(n_par) == 1:
            #             if not n_par[0] in sibs: #sibling of current node already in graph
            #                 for ns in n_sibs: #add parent-child relationships between these nn and their missing parent
            #                     addEdgeType(ns, 'missing_' + str(missing), 'C', 'P', this_fam)
            #                 for s in sibs: #add sibling relationship between the missing parent and their sibs
            #                     addEdgeType(s, 'missing_' + str(missing), 'FS', 'FS', this_fam)
            #                 missing += 1
            #         elif len(n_par) == 0:
            #             for ns in n_sibs:  # add parent-child relationships between these nn and their missing parent
            #                 addEdgeType(ns, 'missing_' + str(missing), 'C', 'P', this_fam)
            #                 addEdgeType(ns, 'missing_' + str(missing+1), 'C', 'P', this_fam)
            #             for s in sibs: #add sibling relationship between one of the missing parents and their sibs
            #                 addEdgeType(s, 'missing_' + str(missing), 'FS', 'FS', this_fam)
            #             missing += 2
            #         checked_nn = checked_nn.union(n_sibs)
            #
            # all_trees_this_fam.append(this_fam) #add this_fam now; relationships added beyond this point can have multiple possibilities
            # #add missing children according to grandchildren
            # checked_gc = set()
            # this_fam_tmp = this_fam.copy()
            # for gc in grandchildren:
            #     if not gc in checked_gc:
            #         [gc_sibs, gc_par] = getSibsParentsFromGraph(this_fam_tmp,gc)
            #         gc_sibs.add(gc)
            #         if len(gc_par): #if the grandchild has parents in the dataset
            #             if not any([x in children for x in gc_par]): #the grandchild's parent that is the grandparent's child isn't in the dataset
            #                 for gcs in gc_sibs:
            #                     addEdgeType(gcs, 'missing_' + str(missing), 'C', 'P', this_fam_tmp) #
            #        # else:
            #             #current_node and any other grandparents of the grandchildren in the dataset could be pairs
            #
