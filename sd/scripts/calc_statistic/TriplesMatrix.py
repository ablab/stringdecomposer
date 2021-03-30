import numpy as np

def BuildTriplesM(trp_cnt, mnid, mnlist):
    res = [[0] * len(mnlist) for i in range(len(mnlist))]
    for mns, cnts in trp_cnt.items():
        if mns[1] == mnid:
            if mns[0] in mnlist and mns[2] in mnlist:
                res[mnlist.index(mns[0])][mnlist.index(mns[2])] += cnts
    return res


def BuildPrefixTriplesM(trp_cnt, mnid, mnlist):
    res = [[0] * len(mnlist) for i in range(len(mnlist))]
    for mns, cnts in trp_cnt.items():
        if mns[2] == mnid:
            if mns[0] in mnlist and mns[1] in mnlist:
                res[mnlist.index(mns[0])][mnlist.index(mns[1])] += cnts
    return res


def BuildSuffixTriplesM(trp_cnt, mnid, mnlist):
    res = [[0] * len(mnlist) for i in range(len(mnlist))]
    for mns, cnts in trp_cnt.items():
        if mns[0] == mnid:
            if mns[1] in mnlist and mns[2] in mnlist:
                res[mnlist.index(mns[1])][mnlist.index(mns[2])] += cnts
    return res


def Norm(Trp):
    sqr_sum = np.sum(np.array(Trp)**2)**(1/2)
    res = (np.array(Trp)/sqr_sum).tolist()
    return res


def BuildNormTriplesM(trp_cnt, mnid, mnlist):
    Triples = BuildTriplesM(trp_cnt, mnid, mnlist)
    return Norm(Triples)


def BuildNormPreifxTriplesM(trp_cnt, mnid, mnlist):
    Triples = BuildPrefixTriplesM(trp_cnt, mnid, mnlist)
    return Norm(Triples)


def BuildNormSuffixTriplesM(trp_cnt, mnid, mnlist):
    Triples = BuildSuffixTriplesM(trp_cnt, mnid, mnlist)
    return Norm(Triples)


def getMnScore(trp1, trp2):
    res = np.sum((np.array(trp1) * np.array(trp2)))
    return res


def GetCenVec(mnlist, trp_cnt):
    trps = {}

    for mn in mnlist:
        Trp = BuildNormTriplesM(trp_cnt, mn, mnlist)
        trps[mn] = Trp
    return trps


def handleAllMn(trp_cnt, db_cnt):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > 100]

    trps = []

    for mn in mnlist:
        Trp = BuildNormTriplesM(trp_cnt, mn, mnlist)
        trps.append(Trp)

    PositionScore = {}
    for i in range(len(mnlist)):
        for j in range(len(mnlist)):
            PositionScore[(mnlist[i], mnlist[j])] = getMnScore(trps[i], trps[j])

    #print(PositionScore)
    return PositionScore, GetCenVec(mnlist, trp_cnt)


def PrefixPosScore(trp_cnt, db_cnt):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > 100]

    trps = []

    for mn in mnlist:
        Trp =  BuildNormPreifxTriplesM(trp_cnt, mn, mnlist)
        trps.append(Trp)

    PositionScore = {}
    for i in range(len(mnlist)):
        for j in range(len(mnlist)):
            PositionScore[(mnlist[i], mnlist[j])] = getMnScore(trps[i], trps[j])
    return PositionScore, GetCenVec(mnlist, trp_cnt)


def SuffixPosScore(trp_cnt, db_cnt):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > 100]

    trps = []

    for mn in mnlist:
        Trp = BuildNormSuffixTriplesM(trp_cnt, mn, mnlist)
        trps.append(Trp)

    PositionScore = {}
    for i in range(len(mnlist)):
        for j in range(len(mnlist)):
            PositionScore[(mnlist[i], mnlist[j])] = getMnScore(trps[i], trps[j])
    return PositionScore, GetCenVec(mnlist, trp_cnt)


def SplitAllMn(trp_cnt, db_cnt):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > 100]
    splTrp = {}

    for mn in mnlist:
        Trp = BuildNormTriplesM(trp_cnt, mn, mnlist)
        bst_lst = []
        for i in range(len(mnlist)):
            for j in range(len(mnlist)):
                bst_lst.append((i, j))
        bst_lst.sort(key=lambda x: -Trp[x[0]][x[1]])

        if mn == "mn_3":
            print("3: ", mnlist[bst[0]], mnlist[bst[1]], mnlist[scb[0]], mnlist[scb[1]])
        for i in range(1, len(bst_lst)):
            isIndependent = True
            scb = bst_lst[i]
            bst = bst_lst[0]
            if not Trp[bst[0]][bst[1]] < 10 * Trp[scb[0]][scb[1]]:
                continue
            for j in range(0, len(bst_lst)):
                    if i == j:
                        continue
                    bst = bst_lst[j]
                    if ((j < i and Trp[bst[0]][bst[1]] < 10 * Trp[scb[0]][scb[1]]) or
                       (i < j and Trp[bst[0]][bst[1]] * 10 > Trp[scb[0]][scb[1]])):
                        if bst[0] == scb[0] or bst[1] == scb[1]:
                            print("DEPEND", mn, mnlist[bst[0]], mnlist[bst[1]], mnlist[scb[0]], mnlist[scb[1]])
                            isIndependent = False

            if isIndependent:
                splTrp[(mnlist[scb[0]], mn, mnlist[scb[1]])] = mn + "." + str(i)

    return splTrp
