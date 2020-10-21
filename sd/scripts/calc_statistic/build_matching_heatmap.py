import sys

log_path = "data/monomers_mutations_with_isolates.log"
prefix_monomer_name = "S3CXH1L."
output_file = "output/heatmap_"

#MONOMER_NAME: [(MonomerId, mutations)...]
monomers_list = {}
with open(log_path) as fr:
    for line in fr:
        if "MonomerId" in line:
            monomerName = line.split()[0]
            monomerId = line.split("MonomerId= ")[-1].split()[0]
            mutationList = line.split('[ ')[-1].split(' ]')[0].split(', ')
            if mutationList[0] == '':
                mutationList = []
            
            for i in range(len(mutationList)):
                mutationList[i] = mutationList[i].split('(')[0]
            
            if monomerName not in monomers_list:
                monomers_list[monomerName] = []
            
            monomers_list[monomerName].append((monomerId, mutationList))


def get_mutation_graph(MonomerName, Pos, Len, Delta, MxVertCnt=2000):
    global monomers_list
    G = [[] for i in range(min(MxVertCnt, Len))]
    for i in range(min(Len, MxVertCnt)):
        for j in range(-Delta, Delta + 1):
            real_pos = Pos + i
            real_pos2 = Pos + i + Len + j
            
            if real_pos < 0 or real_pos >= len(monomers_list[MonomerName]):
                continue
                
            if real_pos2 < 0 or real_pos2 >= len(monomers_list[MonomerName]):
                continue
                
            
            if len(set(monomers_list[MonomerName][real_pos][1]) & set(monomers_list[MonomerName][real_pos2][1])) > 0:
                G[i].append((real_pos2 - Pos - Len + Delta, 1))
    return G
    
    
def get_max_matching_w(G):
    n = len(G)
    m = n;
    for i in range(n):
        for j in G[i]:
            m = max(m, j[0] + 1)
            
    a = [[1]*(m + 1) for i in range(n + 1)]
    for i in range(n):
        for j in range(len(G[i])):
            a[i + 1][G[i][j][0] + 1] = -G[i][j][1]
    u = [0]*(n + 1)
    v = [0]*(m + 1)
    p = [0]*(m + 1)
    way = [0] * (m + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [10**9]*(m + 1)
        used = [0] * (m + 1)
        while (p[j0] != 0):
            used[j0] = 1;
            i0 = p[j0]
            delta = 10**9
            j1 = 0
            for j in range(1, m + 1):
                if used[j] == 0:
                    cur = a[i0][j]-u[i0]-v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(0, m + 1):
                if used[j] == 1:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    
    res = 0
    for i in range(1, m + 1):
        if (p[i] != 0 and a[p[i]][i] < 0):
            res -= a[p[i]][i]
    return res

def build_pick_matrix(MonomerName, mx_dist, mxVert=2000):
    pick_h = [[0]*mx_dist for i in range(len(monomers_list[MonomerName]))]
    for i in range(len(monomers_list[MonomerName])):
        print(i)
        for Len in range(1, mx_dist):
            delta = min(int(0.2*Len), int(0.2*mxVert))
            pick_h[i][Len] = get_max_matching_w(get_mutation_graph(MonomerName, i, Len, delta, mxVert))
    return pick_h
 
monoName = prefix_monomer_name + sys.argv[1]
pick_h = build_pick_matrix(monoName, 1505, 20)
output_file += monoName + ".txt"

original_stdout = sys.stdout
with open(output_file, "w") as f:
    sys.stdout = f
    print(pick_h)
    sys.stdout = original_stdout
