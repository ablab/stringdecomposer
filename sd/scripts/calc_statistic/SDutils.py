import edlib

def get_dist(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    return result["editDistance"]


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def get_closest_mn(mn, mn_list, distf=get_dist):
    dists = [(distf(str(mn.seq), str(mnx.seq)), mnx.id) for mnx in mn_list]
    return min(dists)


def map_mn(mn_list1, mn_list2, distf=get_dist):
    return {mn1.id: get_closest_mn(mn1, mn_list2, distf) for mn1 in mn_list1}