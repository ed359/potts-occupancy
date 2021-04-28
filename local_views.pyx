from cysignals.signals cimport sig_check
from cysignals.memory cimport sig_free, sig_malloc
from cytoolz.itertoolz cimport *

from numpy.random import binomial as sample_binomial

from sage.graphs.base.c_graph cimport CGraph
from sage.graphs.base.dense_graph cimport DenseGraph

from tqdm import tqdm_notebook

# Local views are built in stages
# 1) Choose G_Nv, the graph on the neighborhood of v.
# 2) Add boundary vertices B such that every vertex in Nv has degree d, 
#    and isolated color vertices C. This forms an uncolored G_L
# 3) Color the G_L by adding edges between B and C. If b is connected 
#    to c we say that b is colored with color c.
# 
# Stages 1 and 3 use a canonical augmentation algorithm from Brendan D
# McKay, Isomorph-Free Exhaustive Generation, Journal of Algorithms, 
# Volume 26, Issue 2, 1998, Pages 306-324.
#
# One can also estimate the number of local views with this process.
# The algortithm generates a tree of partially colored local views and if 
# every generated vertex at depth i is accepted with probability p_i
# then a generated local view (at depth b) appears with probability given
# by the product of the p_i. This yields an estimate of the total number.


def save_Ls(int d, max_q=None, int verbosity=0):
    cdef DenseGraph cG, sG
    cdef dict c_dict
    cdef int b, i, n, q, u, x, y
    cdef int i_GNv, i_L, size_Bu, next_b, cmax_q
    cdef list aut_gens, B, Bu, Bus, C, edges_cG, ps
    cdef str filename_cG
    cdef tuple cNu
    
    # The maximum boundary size is d*(d-1), hence we need no more colors
    cmax_q = d*(d-1) if max_q is None else min(d*(d-1), max_q)
    n  = d*d + cmax_q
    
    # Allocate space for boundary vertices and color vertices
    sG = DenseGraph(nverts=d, extra_vertices=d*(d-1)+cmax_q)
    # Compute aut gens and begin canonical augmentation algorithm. Here
    # aut_gens is list of generators for the symmetric group on range(d)
    aut_gens = search_tree(sG, [list(range(d))], False, False)

    # Stage 1) Generate the GNvs: all graphs on d vertices
    for i_GNv, cG in enumerate(tqdm_notebook(augment_edges(d, sG, aut_gens),
                                             desc='GNv loop',
                                             total=int(num_GNvs[d]))):

        # Stage 2) Add boundary vertices and edges to them
        next_b = d
        Bus = []
        B   = []
        for u in xrange(d):
            size_Bu = d - 1 - cG.out_degrees[u]
            Bu = range(next_b, next_b + size_Bu)
            
            # Optimised version of b in Bu
            for b in xrange(next_b, next_b + size_Bu):
                cG.add_vertex_unsafe(b)
                cG.add_arc_unsafe(u, b)
                cG.add_arc_unsafe(b, u)
            
            Bus.append(Bu)
            B.extend(Bu)
            next_b += size_Bu

        # Add color vertices
        q = min(cmax_q, len(B))
        C = range(d + len(B), d + len(B) + q)
        # optimised version of c in C
        for c in xrange(d + len(B), d + len(B) + q):
            cG.add_vertex_unsafe(c)

        # cG numbers Nv from 0 to d-1, but in the returned local views we shift indices to 1..d
        edges_cG = [(min(x+1, y+1), max(x+1, y+1)) for x, y in gen_all_edges(d) if cG.has_arc(x, y)]
        filename_cG = 'data/Ls_d{}_{}.ls'.format(d, '_'.join('{}{}'.format(x, y) for x, y in edges_cG))

        if verbosity > 0: 
            print('GNv %2d with edges %s.' % (i_GNv+1, edges_cG))
            print('Saving to {}'.format(filename_cG))
        if verbosity > 1:
            print 'B = {}.\nC = {}.\nStarting augmentation algorithm'.format(B, C)

        with open(filename_cG, 'w') as f:
            ps = [1] * (len(B)+1)    
            # Stage 3) Use canonical augemntation to color the boundary
            aut_gens = search_tree(cG, [range(d), B, C], False, False)
            for colored_cG in tqdm_notebook(augment_L(d, ps, B, C, cG, aut_gens, verbosity), 
                                            desc='Color loop {}'.format(i_GNv)):
                # Create dictionary of colors for the boundary, shifted to use colors 
                # starting from 0
                c_dict = {b: max(colored_cG.out_neighbors(b)) - d - len(B) for b in B}
                
                # Yield the local view in a concise form
                cNu = tuple(tuple(c_dict[b] for b in Bus[u]) for u in xrange(d))
                f.write('({}, {})\n'.format(cNu, edges_cG))

        # Clear existing boundary and color vertices ready for next iteration
        for v in xrange(d, n):
            cG.del_vertex(v)

        if verbosity > 0: print


def gen_Ls(int d, max_q=None, int verbosity=0):
    cdef DenseGraph cG, sG
    cdef dict c_dict
    cdef int b, i, n, q, u, x, y
    cdef int i_GNv, i_L, size_Bu, next_b, cmax_q
    cdef list aut_gens, B, Bu, Bus, C, edges_cG, ps
    cdef tuple cNu
    
    # The maximum boundary size is d*(d-1), hence we need no more colors
    cmax_q = d*(d-1) if max_q is None else min(d*(d-1), max_q)
    n  = d*d + cmax_q
    
    # Allocate space for boundary vertices and color vertices
    sG = DenseGraph(nverts=d, extra_vertices=d*(d-1)+cmax_q)
    # Compute aut gens and begin canonical augmentation algorithm. Here
    # aut_gens is list of generators for the symmetric group on range(d)
    aut_gens = search_tree(sG, [list(range(d))], False, False)

    # Stage 1) Generate the GNvs: all graphs on d vertices
    for i_GNv, cG in enumerate(tqdm_notebook(augment_edges(d, sG, aut_gens),
                                             desc='GNv loop',
                                             total=int(num_GNvs[d]))):

        # Stage 2) Add boundary vertices and edges to them
        next_b = d
        Bus = []
        B   = []
        for u in xrange(d):
            size_Bu = d - 1 - cG.out_degrees[u]
            Bu = range(next_b, next_b + size_Bu)
            
            # Optimised version of b in Bu
            for b in xrange(next_b, next_b + size_Bu):
                cG.add_vertex_unsafe(b)
                cG.add_arc_unsafe(u, b)
                cG.add_arc_unsafe(b, u)
            
            Bus.append(Bu)
            B.extend(Bu)
            next_b += size_Bu

        # Add color vertices
        q = min(cmax_q, len(B))
        C = range(d + len(B), d + len(B) + q)
        # Optimised version of c in C
        for c in xrange(d + len(B), d + len(B) + q):
            cG.add_vertex_unsafe(c)

        # cG numbers Nv from 0 to d-1, but in the returned local views we shift indices to 1..d
        edges_cG = [(min(x+1, y+1), max(x+1, y+1)) for (x, y) in gen_all_edges(d) if cG.has_arc(x, y)]

        if verbosity > 0: 
            print('GNv %2d with edges %s.' % (i_GNv+1, edges_cG))
        if verbosity > 1:
            print 'B = {}.\nC = {}.\nStarting augmentation algorithm'.format(B, C)
        
        ps = [1] * (len(B)+1)
        # Stage 3) Use canonical augemntation to color the boundary
        aut_gens = search_tree(cG, [range(d), B, C], False, False)
        for colored_cG in augment_L(d, ps, B, C, cG, aut_gens, verbosity):
            # Create dictionary of colors for the boundary, shifted to use colors 
            # starting from 0
            c_dict = {b: max(colored_cG.out_neighbors(b)) - d - len(B) for b in B}
            
            # Yield the local view in a concise form
            cNu = tuple(tuple(c_dict[b] for b in Bus[u]) for u in xrange(d))
            yield (cNu, edges_cG)

        # Clear existing boundary and color vertices ready for next iteration
        for v in xrange(d, n):
            cG.del_vertex(v)

        if verbosity > 0: print


def estimate_Ls(int d, list ps, max_q=None, int verbosity=0):
    cdef DenseGraph cG, sG
    cdef dict c_dict
    cdef float estimate, p
    cdef int b, i, n, q, u, x, y
    cdef int i_GNv, i_L, size_Bu, next_b, cmax_q
    cdef list aut_gens, B, Bu, Bus, C, edges_cG, estimates, ps_slice
    cdef tuple cNu
    
    # The maximum boundary size is d*(d-1), hence we need no more colors
    cmax_q = d*(d-1) if max_q is None else min(d*(d-1), max_q)
    n  = d*d + cmax_q
    estimates = []
    
    # Allocate space for boundary vertices and color vertices
    sG = DenseGraph(nverts=d, extra_vertices=d*(d-1)+cmax_q)
    # Compute aut gens and begin canonical augmentation algorithm. Here
    # aut_gens is list of generators for the symmetric group on range(d)
    aut_gens = search_tree(sG, [list(range(d))], False, False)

    # Stage 1) Generate the GNvs: all graphs on d vertices
    for i_GNv, cG in enumerate(tqdm_notebook(augment_edges(d, sG, aut_gens),
                                             desc='GNv loop',
                                             total=int(num_GNvs[d]))):

        # Stage 2) Add boundary vertices and edges to them
        next_b = d
        Bus = []
        B   = []
        for u in xrange(d):
            size_Bu = d - 1 - cG.out_degrees[u]
            Bu = range(next_b, next_b + size_Bu)
            
            # Optimised version of b in Bu
            for b in xrange(next_b, next_b + size_Bu):
                cG.add_vertex_unsafe(b)
                cG.add_arc_unsafe(u, b)
                cG.add_arc_unsafe(b, u)
            
            Bus.append(Bu)
            B.extend(Bu)
            next_b += size_Bu

        # Add color vertices
        q = min(cmax_q, len(B))
        C = range(d + len(B), d + len(B) + q)
        # Optimised version of c in C
        for c in xrange(d + len(B), d + len(B) + q):
            cG.add_vertex_unsafe(c)

        # cG numbers Nv from 0 to d-1, but in the returned local views we shift indices to 1..d
        edges_cG = [(min(x+1, y+1), max(x+1, y+1)) for (x, y) in gen_all_edges(d) if cG.has_arc(x, y)]

        if verbosity > 0: 
            print('GNv %2d with edges %s.' % (i_GNv+1, edges_cG))
        if verbosity > 1:
            print 'B = {}.\nC = {}.\nStarting augmentation algorithm'.format(B, C)
        
        ps_slice = ps[:len(B)+1]
        # Stage 3) Use canonical augemntation to color the boundary
        aut_gens = search_tree(cG, [range(d), B, C], False, False)
        count = sum(1 for _ in augment_L(d, ps_slice, B, C, cG, aut_gens, 0))
        estimate = count
        for p in ps_slice:
            estimate /= p
        estimates.append(estimate)

        # Clear existing boundary and color vertices ready for next iteration
        for v in xrange(d, n):
            cG.del_vertex(v)

        if verbosity > 0: 
            print 'count: %d, estimate: %.1f' % (count, estimate) 
            print
    if verbosity > 0: 
        print 'total estimate: %.1f' % sum(estimates)

    return estimates


# Both the augment_ functions in this file implement the same algorithm.
# Edges are added to a DenseGraph one at a time. The method is:
# 1) Yield the current graph X and terminate if it is complete.
# 2) Consider the possible 'uppers': ways of extending X and find
#    representations of each isomorphism class of uppers under the action
#    of aut(X). Here an upper is always a non-edge of X.
# 3) For each upper_rep e, form Y by adding e to X.
# 4) Let ystar_orig be a canonical 'lower' in Y, which is an edge
#    of Y. The choice is made via picking the 'last' lower ystar in a
#    canonical relabeling Ycan of Y (see search_tree).
# 5) If the upper_rep e is in the orbit of ystar_orig under the action of
#    aut(Y), continue the method with Y (recurse to step 1), else consider 
#    the next upper_rep
def augment_edges(int n, DenseGraph X, list aut_gens):
    """
    Generate a single representative of every unlableled graph on n vertices.

    INPUT:
    -  ``n`` - number of vertices.
    -  ``X`` - a DenseGraph object with room for at least n vertices.
    -  ``aut_gens`` - a list of generators for aut(X) in permutation list
       format, i.e. the generator g maps i to g[i].

    OUTPUT:
    -  A python generator which yields DenseGraph objects.

    See also implementations in the sage standard libraray (v8.0 2017-07-21):
        sage.groups.perm_gps.partn_ref.refinement_graphs
        sage.graphs.graph_generators
    """
    sig_check() # Respond to interrupts

    cdef DenseGraph Y, Ycan
    cdef list mY, upper_reps
    cdef set uppers
    cdef tuple e, ystar, ystar_orig
    
    yield X
    if X.num_arcs == n*(n-1):
        return
        
    uppers = {f for f in gen_all_edges(n) if not X.has_arc(f[0], f[1])}
    upper_reps = find_upper_reps(aut_gens, uppers, False)
    for e in upper_reps:

        Y = copy_DenseGraph(X)
        Y.add_arc(e[0], e[1])
        Y.add_arc(e[1], e[0])
        
        Y_aut_gens, Ycan, cr = search_tree(Y, [range(n)], True, True)
        icr = {v: k for k, v in cr.iteritems()} # inverse of relableling
        ystar = max(f for f in gen_all_edges(n) if Ycan.has_arc(f[0], f[1]))
        ystar_orig = tuple(sorted([icr[ystar[0]], icr[ystar[1]]]))

        mY = find_orbit(ystar_orig, Y_aut_gens)
        if e in mY:
            for Z in augment_edges(n, Y, Y_aut_gens):
                yield Z


# Minor modifications of the above algorithm can be used to color graphs in the
# the generation of local views. To color (a subset of) the vertices of a graph
# with q colors, we augment the graph with a set C of q vertices and represent
# coloring a vertex v with color c in C by an edge from c to v. Automorphisms 
# of the colored graph correspond to automorphisms of the augmented graph that
# fix C. The modifications to steps 1-5 above are as follows:
# 1) Only yield X when every boundary vertex is colored.
# 2) Uppers are edges from uncolored boundary vertices (B_unc) to allowed 
#    colors (C_allowed). As a minor optimisation we only generate colorings
#    that use an initial segment of the colors C.
# 3) No changes
# 4) Inform the search_tree method of the the partition [Nv, B, C] so that it
#    computes the automorphism group fixing these parts.
# 5) Chosing a canonical lower is similar, find Ycan and pick the 'last' edge
#    from a boundary vertex to a color vertex.
def augment_L(int d, list ps, list B, list C, 
              DenseGraph X, list aut_gens, int verbosity=0):
    """
    Generate local views. Add edges representing colors to the graph X.

    INPUT:
    -  ``d`` - the size of Nv. Nv is represented as [0, ..., d-1].
    -  ``ps`` - a list of probabilities with which each generation is accepted
    -  ``B`` - the list of boundary vertices.
    -  ``C`` - the list of color vertices.
    NOTE: Nv, B, and C must intervals of integers with no gaps between them.
    -  ``X`` - a DenseGraph object with vertices the union of Nv, B, and C.
    -  ``aut_gens`` - a list of generators for aut(X) in permutation list
       format, i.e. the generator g maps i to g[i].
    NOTE: Here aut(X) is a permutation of the union of Nv, B, and C that fixes
    the parts Nv, B, and C.

    OUTPUT:
    -  A python generator which yields DenseGraph objects.
    """
    sig_check() # Respond to interrupts

    cdef DenseGraph Y, Ycan
    cdef list B_unc, C_allowed, C_used, mY, upper_reps
    cdef set uppers
    cdef tuple e, ystar, ystar_orig
    cdef int depth
    
    depth = 0
    for b in B:
        if <int>X.out_degrees[b] > 1:
            depth += 1
    
    if not sample_bernoulli(ps[depth]):
        return

    B_unc  = [b for b in B if <int>X.out_degrees[b] == 1]

    if B_unc == []:
        if verbosity > 5: print "        yield"
        yield X
        return

    # Any used color is allowed, and if there are unused colors
    # the first of those is allowed
    C_allowed = [c for c in C if <int>X.out_degrees[c] > 0]
    if len(C) > len(C_allowed):
        C_allowed.append(C[len(C_allowed)])

    uppers = {(b, c) for b in B_unc for c in C_allowed}
    upper_reps = find_upper_reps(aut_gens, uppers, True)
   
    for e in upper_reps:

        if verbosity > 1:
            color_edges = [(b, c) for c in C for b in X.in_neighbors(c)]
            print "  partially colored cG: %s" % color_edges
        if verbosity > 4: print "    B_unc      = %s" % B_unc
        if verbosity > 4: print "    C_allowed  = %s" % C_allowed
        if verbosity > 3: print "    upper_reps = %s" % upper_reps
        if verbosity > 3: print "    e          = %s" % str(e)

        # Select Xupper = (X, e) where e is the first member of the orbit. 
        # Let Y be X with edge e added, so that Ylower = (Y, x) is in f'(Xupper)
        Y = copy_DenseGraph(X)
        Y.add_arc(e[0], e[1])
        Y.add_arc(e[1], e[0])
        
        Y_aut_gens, Ycan, cr = search_tree(Y, [range(d), B, C], True, True)
        icr = {v: k for k, v in cr.iteritems()}  # inverse of relableling
        Ccan = sorted(cr[c] for c in C)
        ystar = max((b, c) for c in Ccan for b in Ycan.in_neighbors(c))
        ystar_orig = tuple(sorted([icr[ystar[0]], icr[ystar[1]]]))

        mY = find_orbit(ystar_orig, Y_aut_gens)
        if verbosity > 3: 
            print "    ystar_orig = %s" % str(ystar_orig)
            print "    mY         = %s" % mY
        if e in mY:
            if verbosity > 4: print "      e in mY, recursing"
            for Z in augment_L(d, ps, B, C, Y, Y_aut_gens, verbosity):
                yield Z
        
        if verbosity > 4: print "      going to next upper_rep"
    if verbosity > 4: print "      upper_reps exhausted"


cpdef bint sample_bernoulli(p):
    return sample_binomial(1, p)


def list_GNv_edges(int d):
    cdef DenseGraph cG, sG
    cdef int x, y
    cdef list aut_gens
    
    sG = DenseGraph(nverts=d)
    # Compute aut gens and begin canonical augmentation algorithm. Here
    # aut_gens is list of generators for the symmetric group on range(d)
    aut_gens = search_tree(sG, [list(range(d))], False, False)

    return [[(min(x+1, y+1), max(x+1, y+1)) for x, y in gen_all_edges(d) if cG.has_arc(x, y)]
            for cG in augment_edges(d, sG, aut_gens)]

def list_GNv_filenames(int d):
    cdef int x, y
    cdef list edges_cG
    return ['data/Ls_d{}_{}.ls'.format(d, '_'.join('{}{}'.format(x, y) for x, y in edges_cG))
            for edges_cG in list_GNv_edges(d)]


# https://oeis.org/A000088
cpdef list num_GNvs = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864]


cdef inline gen_all_edges(int n):
    return ((i, j) for j in xrange(n) for i in xrange(j))


cpdef DenseGraph copy_DenseGraph(DenseGraph G):
    cdef int i, j, n_total, n_active
    cdef DenseGraph H
    n_total = <int>G.active_vertices.size
    n_active = G.num_verts
    H = DenseGraph(nverts=n_active, extra_vertices=n_total-n_active)
    for i in xrange(n_active):
        for j in G.out_neighbors(i):
            H.add_arc_unsafe(i, j)
    return H


cdef inline tuple permute_edge(tuple e, list gen):
    cdef tuple f = (gen[e[0]], gen[e[1]])
    return f if f[0] < f[1] else (f[1], f[0])


cpdef tuple reverse(tuple t):
  return (t[1], t[0])


cpdef list find_upper_reps(list aut_gens, set uppers, bint end_sort):
    """
    Compute a list of representatives of each isomorphism class of uppers under
    the action of the permutation group given by aut_gens

    INPUT:
    -  ``aut_gens`` - a list of generators for the group in permutation list.
       format, i.e. the generator g maps i to g[i].
    -  ``uppers`` - a set of edges (tuples of length 2) to permute
    -  ``end_sort`` - a flag indicating whether to sort the reperesentatives 
       by reverse lexicographic order.
    
    OUTPUT:
    -  A list of representatives
    """
    cdef list orbit, orbits
    cdef set sorbit
    
    orbits = []
    # WARNING the set that is iterated over gets modified in the loop
    while uppers:
        orbit = find_orbit(next(iter(uppers)), aut_gens)
        orbits.append(orbit)
        for e in orbit:
            uppers.discard(e)

    if end_sort:
        return sorted([orbit[0] for orbit in orbits], key=reverse)
    else:
        return sorted([orbit[0] for orbit in orbits])


cpdef list find_orbit(tuple e, list aut_gens):
    """
    Compute the orbit of an edge under the action of a permutation group

    INPUT:
    -  ``e`` - an edge (tuple of length 2).
    -  ``aut_gens`` - a list of generators for the group in permutation list 
       format, i.e. the generator g maps i to g[i].
    
    OUTPUT:
    -  The orbit of e as a list
    """
    cdef list lorbit = [e]
    cdef set  sorbit = {e}
    cdef int i = 0
    while i < len(lorbit):
        f = lorbit[i]
        for gen in aut_gens:
            g = permute_edge(f, gen)
            if g not in sorbit:
                sorbit.add(g)
                lorbit.append(g)
        i += 1
    lorbit.sort()
    return lorbit


# Import cython implementations of partition refinement algorithms from the
# sage standard library (v8.0 2017-07-21)
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport (
    aut_gp_and_can_lab, deallocate_agcl_output, get_aut_gp_and_can_lab)
from sage.groups.perm_gps.partn_ref.data_structures cimport (
    PartitionStack, PS_dealloc, PS_from_list, PS_num_cells)
from sage.groups.perm_gps.partn_ref.refinement_graphs cimport (
    all_children_are_equivalent, compare_graphs, GraphStruct, refine_by_degree)


cpdef search_tree(DenseGraph G, list partition, bint lab, bint certificate):
    """
    Compute automorphism groups and canonical labels of graphs. This function 
    is copied from sage.groups.perm_groups.partn_ref.refinement_graphs (v8.0 
    2017-07-21) but ancilliary code has been removed.

    INPUT:
    -  ``G`` - a DenseGraph object.
    -  ``partitions`` - a list of lists representing a partition of V(G). The
       regurned group fixes parts of this partition.
    -  ``lab`` - a flab indicating whether the canonically relabeled G should
       be returned as a DenseGraph.
    -  ``certificate`` - a flag indicating whether a dictionary of the 
       canonical relabeling should be returned.

    OUTPUT (as a tuple if multiple items returned):
    -  A list of generators for aut(G) (repsecting the partition)
    -  A DenseGraph which is the canoncial labeling of G (if lab)
    -  A dictionary cr such that cr[v] is the canonical label of v (if 
       certificate)
    """
    cdef int i, j, n
    cdef aut_gp_and_can_lab *output
    cdef PartitionStack *part

    n = G.num_verts

    cdef GraphStruct GS = GraphStruct()
    GS.G = <CGraph>G
    GS.directed = 0
    GS.loops = 0
    GS.use_indicator = 1

    if n == 0:
        return_tuple = [[]] # no aut gens
        if lab:
            G_C = DenseGraph(n)
            return_tuple.append(G_C)
        if certificate:
            return_tuple.append({})
        if len(return_tuple) == 1:
            return return_tuple[0]
        else:
            return tuple(return_tuple)

    GS.scratch = <int *> sig_malloc( (3*G.num_verts + 1) * sizeof(int) )
    part = PS_from_list(partition)
    if GS.scratch is NULL or part is NULL:
        PS_dealloc(part)
        sig_free(GS.scratch)
        raise MemoryError

    output = get_aut_gp_and_can_lab(<void *>GS, part, G.num_verts, 
                                    all_children_are_equivalent, 
                                    refine_by_degree, compare_graphs, 
                                    lab, NULL, NULL, NULL)
    sig_free(GS.scratch)

    # prepare output
    list_of_gens = []
    for i in xrange(output.num_gens):
        list_of_gens.append([output.generators[j+i*G.num_verts] 
                            for j in xrange(G.num_verts)])
    return_tuple = [list_of_gens]
    if lab:
        G_C = DenseGraph(n)
        for i in xrange(n):
            for j in G.out_neighbors(i):
                G_C.add_arc(output.relabeling[i],output.relabeling[j])
        return_tuple.append(G_C)
    if certificate:
        cr = {}
        for i in xrange(G.num_verts):
            cr[i] = output.relabeling[i]
        return_tuple.append(cr)
    PS_dealloc(part)
    deallocate_agcl_output(output)
    if len(return_tuple) == 1:
        return return_tuple[0]
    else:
        return tuple(return_tuple)
