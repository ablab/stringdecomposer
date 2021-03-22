import numpy as np

def BuildTriplesM(trp_cnt, mnid, mnlist):
    res = [[0] * len(mnlist) for i in range(len(mnlist))]
    for mns, cnts in trp_cnt.items():
        if mns[1] == mnid:
            if mns[0] in mnlist and mns[2] in mnlist:
                res[mnlist.index(mns[0])][mnlist.index(mns[2])] += cnts
    return res


def Norm(Trp):
    sqr_sum = np.sum(np.array(Trp)**2)**(1/2)
    res = (np.array(Trp)/sqr_sum).tolist()
    return res

def BuildNormTriplesM(trp_cnt, mnid, mnlist):
    Triples = BuildTriplesM(trp_cnt, mnid, mnlist)
    return Norm(Triples)


def getMnScore(trp1, trp2):
    res = np.sum((np.array(trp1) * np.array(trp2)))
    return res


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
    return PositionScore

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
                bst = bst_lst[0]
                scb = bst_lst[i]
                if mn == "mn_43":
                    print(mnlist[bst[0]], mnlist[bst[1]], mnlist[scb[0]], mnlist[scb[1]], Trp[bst_lst[i][0]][bst_lst[i][1]])
                if Trp[bst_lst[0][0]][bst_lst[0][1]] < 5 * Trp[bst_lst[i][0]][bst_lst[i][1]]:
                    if bst[0] != scb[0] and bst[1] != scb[1]:
                        splTrp[(mnlist[scb[0]], mn, mnlist[scb[1]])] = mn + "." + str(i)

    return splTrp
