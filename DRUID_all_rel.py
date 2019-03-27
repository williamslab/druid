def getIBD1(ind1, ind2, all_rel):
    if ind1 < ind2:
        return all_rel[ind1][ind2][0]
    else:
        return all_rel[ind2][ind1][0]

def getIBD2(ind1, ind2, all_rel):
    if ind1 < ind2:
        return all_rel[ind1][ind2][1]
    else:
        return all_rel[ind2][ind1][1]

def getPairwiseK(ind1, ind2, all_rel):
    if ind1 < ind2:
        return all_rel[ind1][ind2][2]
    else:
        return all_rel[ind2][ind1][2]

def getPairwiseD(ind1, ind2, all_rel):
    if ind1 < ind2:
        return all_rel[ind1][ind2][3]
    else:
        return all_rel[ind2][ind1][3]

def getPairD_w_Name(ind1, ind2, all_rel):
    if ind1 < ind2:
        pair_name = ind1 + "$" + ind2
        return (all_rel[ind1][ind2][3], pair_name)
    else:
        pair_name = ind2 + "$" + ind1
        return (all_rel[ind2][ind1][3], pair_name)

def getPairName(ind1, ind2):
    if ind1 < ind2:
        pair_name = ind1 + "$" + ind2
        return pair_name
    else:
        pair_name = ind2 + "$" + ind1
        return pair_name
