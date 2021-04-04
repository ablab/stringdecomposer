import numpy as np
import csv

def calc_mn_order_stat(sdtsv, cenid, maxk = 1, exchange=None, exchTrp=None):
    k_cnt = [{} for k in range(maxk)]
    rows = []
    with open(sdtsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if cenid not in row[0]:
                continue
            if row[2] == "start":
                continue
            if cenid == "cen1_" and row[1][-1] == "'":
                continue
            #print(row[1], end=" ")
            row[1] = row[1].rstrip("'")
            if exchange is not None and row[1] in exchange:
                row[1] = exchange[row[1]]
            rows.append(row)

    if exchTrp is not None:
        for i in range(1, len(rows) - 1):
            #if rows[i][1] == "mn_43":
            #    print((rows[i+1][1],rows[i][1],rows[i-1][1]))

            if (rows[i+1][1],rows[i][1],rows[i-1][1]) in exchTrp:
                #print("exchange")
                rows[i][1] = exchTrp[(rows[i+1][1],rows[i][1],rows[i-1][1])]

    for i, row in enumerate(rows):
        identity = float(row[4])
        mon = row[1]
        if row[-1] == '?':
            continue

        cur_mons = (mon,)

        for k in range(1, maxk + 1):
            if i - k >= 0 and rows[i - k] != []:
                pident = float(rows[i - k][4])
                pmon = rows[i - k][1]
                if pmon[-1] == "'":
                    pmon = pmon[:-1]

                if rows[i - k][-1] == '?':
                    break

                cur_mons = (*cur_mons, pmon)
                if cur_mons not in k_cnt[k - 1]:
                    k_cnt[k - 1][cur_mons] = 0
                k_cnt[k - 1][cur_mons] += 1
            else:
                break

    return k_cnt

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


def handleAllMn(trp_cnt, db_cnt, thr=100):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > thr]

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


def PrefixPosScore(trp_cnt, db_cnt, thr=100):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > thr]

    trps = []

    for mn in mnlist:
        Trp =  BuildNormPreifxTriplesM(trp_cnt, mn, mnlist)
        trps.append(Trp)

    PositionScore = {}
    for i in range(len(mnlist)):
        for j in range(len(mnlist)):
            PositionScore[(mnlist[i], mnlist[j])] = getMnScore(trps[i], trps[j])
    return PositionScore, GetCenVec(mnlist, trp_cnt)


def SuffixPosScore(trp_cnt, db_cnt, thr = 100):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > thr]

    trps = []

    for mn in mnlist:
        Trp = BuildNormSuffixTriplesM(trp_cnt, mn, mnlist)
        trps.append(Trp)

    PositionScore = {}
    for i in range(len(mnlist)):
        for j in range(len(mnlist)):
            PositionScore[(mnlist[i], mnlist[j])] = getMnScore(trps[i], trps[j])
    return PositionScore, GetCenVec(mnlist, trp_cnt)


def SplitAllMn(trp_cnt, db_cnt, thr=100):
    mncnt = {x[0]: 0 for x in db_cnt.keys()}
    for x, y in db_cnt.items():
        mncnt[x[0]] += y

    mnlist = [x for x, y in mncnt.items() if y > thr]
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
            if not Trp[bst[0]][bst[1]] < 5 * Trp[scb[0]][scb[1]]:
                continue
            for j in range(0, len(bst_lst)):
                    if i == j:
                        continue
                    bst = bst_lst[j]
                    if ((j < i and Trp[bst[0]][bst[1]] < 5 * Trp[scb[0]][scb[1]]) or
                       (i < j and Trp[bst[0]][bst[1]] * 5 > Trp[scb[0]][scb[1]])):
                        if bst[0] == scb[0] or bst[1] == scb[1]:
                            print("DEPEND", mn, mnlist[bst[0]], mnlist[bst[1]], mnlist[scb[0]], mnlist[scb[1]])
                            isIndependent = False

            if isIndependent:
                splTrp[(mnlist[scb[0]], mn, mnlist[scb[1]])] = mn + "." + str(i)

    return splTrp
