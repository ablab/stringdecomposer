#!/usr/bin/env python3

import sys

"""
    Edmonds branching algorithm implementation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Computes a spanning arborescence of minimum weight.

    Input file format::

        {# of vertices} {# of edges} {root vertex}
        {edge source} {edge dest} {edge weight}
        {edge source} {edge dest} {edge weight}
        {edge source} {edge dest} {edge weight}
        ...

    Output file format::

        {# of vertices} {# of edges} {root vertex} {total weight of branching}
        {edge source} {edge dest} {edge weight}
        {edge source} {edge dest} {edge weight}
        {edge source} {edge dest} {edge weight}
        ...

    (C) 2016, CC-0, Lukas Prokop
"""


def read_input_graph(filepath):
    """Given a `filepath`, read a digraph file.
    The format is specified in this file's documentation.

    :param filepath:    filepath to a digraph file
    :type filepath:     str
    :return:            a tuple of vertices, edges and a root
    :rtype:             ([int], [(int, int, float)], int)
    """
    vertices = []
    edges = []
    root = None

    first = True
    with open(filepath, encoding='utf-8') as fd:
        for lineno, line in enumerate(fd):
            if any(line.startswith(c) for c in 'cb#'):
                continue
            if first:
                vals = tuple(map(int, line.split()))
                assert len(vals) >= 3, "first line must contain 3 integers"
                assert vals[0] > 0, "number of vertices must be positive"
                assert vals[1] >= 0, "number of edges must be non-negative"
                assert vals[2] > 0, "root must be an existing vertex"
                vertices = list(range(1, vals[0] + 1))
                num_edges = vals[1]
                root = vals[2]
                first = False
            else:
                vals = line.split()
                assert len(vals) == 3, "every edge line must contain 3 values"
                assert int(vals[0]) > 0 and int(vals[1]) > 0, \
                    "vertices must be 1-enumerated (1..n)"
                edges.append((int(vals[0]), int(vals[1]), float(vals[2])))

    assert not first, "file must not be empty"
    assert len(edges) == num_edges, "Actual # of edges differs from specified"
    assert root in vertices, "root id exceeds vertex enumeration"
    assert all(s in vertices and d in vertices for (s, d, w) in edges)

    return (vertices, edges, root)


def remove_multiedges(E):
    """Returns ``(s, d, w)`` with unique ``(s, d)`` values and `w` minimized.

    :param E:       a set of edges
    :type E:        [(int, int, float)]
    :return:        a subset of edges `E`
    :rtype:         [(int, int, float), ...]
    """
    result = []
    exclusion = set()
    for i, (si, di, wi) in enumerate(E):
        if i in exclusion:
            continue
        minimum = 0
        for j in range(i + 1, len(E)):
            if j in exclusion:
                continue
            sj, dj, wj = E[j]
            if si == sj and di == dj:
                if wi > wj:
                    exclusion.add(i)
                elif wi < wj:
                    exclusion.add(j)
        if i in exclusion:
            continue
        result.append(E[i])
    return result


def traverse(start, E):
    """Given a set of edges, find all DFS paths starting at vertex `start`.

    :param start:   root of paths to be traversed
    :type start:    int
    :param E:       set of edges to traverse
    :type E:        [(int, int, float)]
    :return:        generator for a sequence of edges
    :rtype:         [((int, int, float), ...)]
    """
    path = []
    current = start

    while True:
        options = list(filter(lambda e: e[0] == current, E))
        if len(options) >= 1:
            current = options[0][1]
            path.append(options)
        if not options or current in set(node[0][0] for node in path):
            if path:
                yield tuple(map(lambda node: node[0], path))
            while path and len(path[-1]) == 1:
                path.pop()
            if path:
                path[-1] = tuple(path[-1][1:])
            if not path:
                break


def find_cycle(E):
    """Does the given set of edges contain a cycle?
    If so, return a sequence of edges describing the cycle.
    If not, return False.

    :param E:       a set of edges
    :param E:       [(int, int, float), ...]
    :return:        a sequence of edges describing a cycle
    :rtype:         ((int, int, float), ...) or False
    """
    V = set(map(lambda e: e[0], E)).union(set(map(lambda e: e[1], E)))
    visited = set()
    while V.difference(visited):
        v = V.difference(visited).pop()
        for path in traverse(v, E):
            sp = set(node[0] for node in path)
            is_cycle = path[-1][1] in sp
            visited = visited.union(sp)
            if is_cycle:
                return path
        visited.add(v)
    return False


def cheapest_edges(root, E):
    """Given a set of edges (s, d, w), make d unique and minimize w.

    :param root:    a root vertex to start search from
    :type root:     int
    :param E:       a set of edges
    :type E:        [(int, int, float), ...]
    :return:        subset of edges with the cheapest edges entering any vertex
    :rtype:         [(int, int, float), ...]
    """
    result = {}
    for (s, d, w) in E:
        if d == root:
            continue
        if d in result and result[d][1] < w:
            src = result[d][0]
            weight = result[d][1]
        else:
            src = s
            weight = w

        result[d] = (src, weight)
    return [(s, d, w) for d, (s, w) in result.items()]


def pi(dest, E):
    """Return ``s`` in ``(s, d, w)`` with ``d`` == `dest` and `w` minimized.

    :param dest:    destination vertex
    :type dest:     int
    :param E:       a set of edges
    :type E:        [(int, int, float), ...]
    :return:        vertex with cheapest edge connected to `dest`
    :rtype:         int or None
    """
    src, weight = None, 0
    for (s, d, w) in E:
        if d == dest:
            if src is None or w < weight:
                src = s
                weight = w
    return src


def unique_edge(dest, E):
    """Return the unique edge pointing to `dest`.

    :param dest:    destination vertex
    :type dest:     int
    :param E:       a set of edges
    :type E:        [(int, int, float), ...]
    :return:        first edge connected to vertex `dest`
    :rtype:         (int, int, float) or None
    """
    for (s, d, w) in E:
        if d == dest:
            return (s, d, w)


def edmonds(V, E, root):
    """Recursive application of Edmonds' algorithm according
    to Wikipedia's description [0].

    [0] https://en.wikipedia.org/wiki/Edmonds'_algorithm#Description

    :param V:       set of vertices
    :type V:        [int, ...]
    :param E:       a set of edges
    :type E:        [(int, int, float), ...]
    :param root:    root vertex
    :type root:     int
    :return:        a subset of `E` representing a min-weight arborescence
    :rtype:         [(int, int, float), ...]
    """
    print(("c computing spanning arborescence of minimum weight "
           " for {} with root={}").format(E, root))

    E = list(filter(lambda e: e[1] != root, E))
    E = remove_multiedges(E)

    P = cheapest_edges(root, E)
    print("c P = {}".format(P))
    C = find_cycle(P)
    if not C:
        print("c found no cycle, returning {}".format(P))
        return P
    else:
        print("c found a cycle: {}".format(C))
    C_V = set(e[0] for e in C).union(set(e[1] for e in C))

    v_c = max(V) + 1
    E_prime = []
    correspondence = {}
    for (s, d, w) in E:
        if s not in C_V and d in C_V:
            fe = filter(lambda e: e[0] == pi(d, E) and e[1] == d, E)
            incoming_weight = list(map(lambda e: e[2], fe))[0]
            correspondence[s, v_c, w - incoming_weight] = (s, d, w)
            E_prime.append((s, v_c, w - incoming_weight))
        elif s in C_V and d not in C_V:
            correspondence[v_c, d, w] = (s, d, w)
            E_prime.append((v_c, d, w))
        elif s not in C_V and d not in C_V:
            correspondence[s, d, w] = (s, d, w)
            E_prime.append((s, d, w))

    D_prime = (V + [v_c], E_prime)
    A_prime = edmonds(D_prime[0], D_prime[1], root)

    u, v, w = correspondence[unique_edge(v_c, A_prime)]
    assert v in C_V
    A_prime_corr = list(map(lambda e: correspondence[e], A_prime))
    C_wo_pi_v = list(filter(lambda e: e[0] != pi(v, E) and e[1] != v, C))
    print("c returning arborescence {}".format(A_prime_corr + C_wo_pi_v))
    return A_prime_corr + C_wo_pi_v


def main(filepath):
    """Main routine.

    :param filepath:        A filepath to a digraph file
    :type filepath:         str
    """
    V, E, root = read_input_graph(filepath)
    max_branching = edmonds(V, E, root)

    vertices = max(max(e[0] for e in max_branching),
                   max(e[1] for e in max_branching))
    total_weight = sum(e[2] for e in max_branching)

    print(vertices, len(max_branching), root, total_weight)
    for (s, d, w) in E:
        print('b', s, d, int(w) if w % 1 == 0.0 else w)
    for (s, d, w) in max_branching:
        print(s, d, int(w) if w % 1 == 0.0 else w)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('./edmonds.py <filepath>')
        sys.exit(1)
    main(sys.argv[1])
