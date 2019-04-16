def getIBD1(ind1, ind2, all_rel):
    if ind1 < ind2:
        i1 = ind1
        i2 = ind2
    else:
        i1 = ind2
        i2 = ind1

    if not i1 in all_rel.keys():
        return 0
    elif not i2 in all_rel[i1].keys():
        return 0
    else:
        return all_rel[i1][i2][0]


def getIBD2(ind1, ind2, all_rel):
    if ind1 < ind2:
        i1 = ind1
        i2 = ind2
    else:
        i1 = ind2
        i2 = ind1

    if not i1 in all_rel.keys():
        return 0
    elif not i2 in all_rel[i1].keys():
        return 0
    else:
        return all_rel[i1][i2][1]


def getPairwiseK(ind1, ind2, all_rel):
    if ind1 < ind2:
        i1 = ind1
        i2 = ind2
    else:
        i1 = ind2
        i2 = ind1

    if not i1 in all_rel.keys():
        return 0
    elif not i2 in all_rel[i1].keys():
        return 0
    else:
        return all_rel[i1][i2][2]


def getPairwiseD(ind1, ind2, all_rel):
    if ind1 < ind2:
        i1 = ind1
        i2 = ind2
    else:
        i1 = ind2
        i2 = ind1

    if not i1 in all_rel.keys():
        return -1
    elif not i2 in all_rel[i1].keys():
        return -1
    else:
        return all_rel[i1][i2][3]


def getPairD_w_Name(ind1, ind2, all_rel):
    if ind1 < ind2:
        pair_name = ind1 + "$" + ind2
        i1 = ind1
        i2 = ind2
    else:
        pair_name = ind2 + "$" + ind1
        i1 = ind2
        i2 = ind1

    if not i1 in all_rel.keys():
        return (-1, pair_name)
    elif not i2 in all_rel[i1].keys():
        return (-1, pair_name)
    else:
        return (all_rel[i1][i2][3], pair_name)


def getPairName(ind1, ind2):
    if ind1 < ind2:
        pair_name = ind1 + "$" + ind2
        return pair_name
    else:
        pair_name = ind2 + "$" + ind1
        return pair_name
