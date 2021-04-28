from cython.view cimport array

from sage.arith.misc import factor
from sage.arith.misc import falling_factorial as sfallfact
from sage.calculus.var import var
from sage.combinat.partition import Partitions
from sage.combinat.tuple import Tuples
from sage.functions.other import binomial as sbinomial
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense as Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.vector_rational_dense cimport Vector_rational_dense as Vector
from sage.rings.integer_ring import Z as ZZ
from sage.rings.rational cimport Rational
from sage.rings.rational_field import Q as QQ
from sage.structure.sage_object import dumps, loads
from sage.symbolic.expression import Expression
from sage.symbolic.relation import solve
from sage.symbolic.ring import SR

from ast import literal_eval
from collections import Counter
from cytoolz.itertoolz import *
import itertools as its
from joblib import Parallel, delayed
import struct
from tqdm import tqdm_notebook


###########
# Data IO #
###########


# Reading local views
def read_Ls(filename, start=0, count=None):    
    with open(filename, 'r') as f:
        gen = drop(start, f) if count is None else take(count, drop(start, f))  
        for line in gen:
            yield literal_eval(line[:-1])

def read_all_Ls(d):
    for filename in L_filenames[d]:
        for L in read_Ls(filename):
            yield L

# IO for length-padded data in binary form
cpdef read_padded_raw(f):
    length_dump = f.read(4)
    if len(length_dump) == 0:
        return None # graceful failure, reached EOF
    elif len(length_dump) < 4:
        print('error reading file: %s' % str(f.name))
        return None

    (length,) = struct.unpack('I', length_dump)
    prop_dump = f.read(length)
    if len(prop_dump) < length:
        print('error reading file: %s' % str(f.name))
        return None

    return prop_dump

cpdef write_padded_raw(f, payload):
    length_dump = struct.pack('I', int(len(payload)))
    f.write(length_dump)
    f.write(payload)   

# Reading properties

def read_props(filename, start=0, count=None):
    with open(filename, 'rb') as f:
        iterator = its.count(0) if count is None else xrange(count+start)
        for i in iterator:
            props_dump = read_padded_raw(f)
            if props_dump is None:
                break
            if i >= start:
                yield loads(props_dump)

def read_all_props(int d):
    for filename in P_filenames[d]:
        for p in read_props(filename):
            yield p

# Writing properties

cpdef build_prop_dump(tuple L, q, lam):
    return dumps(build_props(L, q, lam, factor=True), compress=True)

def TqdmParallelExecutor(**joblib_args):
    def aprun(**tqdm_args):
        def tmp(op_iter):
            return Parallel(**joblib_args)(tqdm_notebook(op_iter, **tqdm_args))
        return tmp
    return aprun

def save_all_props(int d, int n_jobs=-1, leave_pbar=False):
    cdef int i, total
    cdef str desc, L_filename, P_filename
    q = var('q')
    lam = var('lam')

    print('Parallel property generation for d={}'.format(d))
    for i in range(num_GNvs[d]):
        L_filename = L_filenames[d][i]
        P_filename = P_filenames[d][i]
        total      = num_Ls_GNv[d][i]
        desc       = 'GNv_{}'.format(i)
        
        aprun = TqdmParallelExecutor(n_jobs=-1)
        pds = aprun(total=total, desc=desc, leave=leave_pbar)(delayed(build_prop_dump)(L, q, lam)
                                            for L in read_Ls(L_filename))

        print 'Writing to %s' % P_filename
        with open(P_filename, 'wb') as f:
            for pd in pds:
                write_padded_raw(f, pd)
        print


#######################################
# Partition functions for Kdd and Kd1 #
#######################################

def build_ZKdd(int d, q, lam):
    cdef int num_colorings, num_colors
    Z_acc = 0
    for X_color_freqs in Partitions(d):
        num_colors = len(X_color_freqs)
        color_choices = sfallfact(q, num_colors)
        num_colorings = num_set_partitions(X_color_freqs.to_list())
        Z_acc += color_choices*num_colorings*(sum(lam**f for f in X_color_freqs) + q - num_colors)**d
    return Z_acc


def build_UKdd(int d, q, lam):
    sq   = var('sq')
    slam = var('slam')
    sym_Z = build_ZKdd(d, sq, slam)
    return lam*sym_Z.derivative(slam).subs(sq=q, slam=lam)/sym_Z.subs(sq=q, slam=lam)/(2*d)
    

def build_ZKd1(int d, q, lam):
    cdef int num_colorings
    Z_acc = 0
    for color_freqs in Partitions(d+1):
        color_choices = sfallfact(q, len(color_freqs))
        num_colorings = num_set_partitions(color_freqs.to_list())
        Z_acc += color_choices*num_colorings*lam**(sum(cbinomial(f,2) for f in color_freqs))
    return Z_acc


def build_UKd1(int d, q, lam):
    sq   = var('sq')
    slam = var('slam')
    sym_Z = build_ZKd1(d, sq, slam)
    return lam*sym_Z.derivative(slam).subs(sq=q, slam=lam)/sym_Z.subs(sq=q, slam=lam)/(d+1)

####################################
# Helper functions for local views #
####################################

# Sage can't factor 0
def factor0(e):
    return 0 if e == 0 else factor(e)

# Return true if L can occur in Kdd
cpdef bint is_Kdd(L):
    cdef tuple Bu0 = L[0][0]
    if len(Bu0) != len(L[0]) - 1: return False
    for Bu in L[0][1:]:
        if Bu != Bu0: return False
    return True

# Return true if L is Kd1
cdef bint is_Kd1(L):
    return get_colors(L) == []

# Get multiset of colors used in L
cdef list get_colors(tuple L):
    cdef list ret = []
    for colors in L[0]:
        ret.extend(colors)
    return ret

# Find the number of colors used in a local view. Since the colors start
# at 0 and are contiguous this is 1 more than the maximum color
cpdef int find_qL(tuple L):
    cdef list colors = get_colors(L)
    return 1 + max(colors) if colors else 0

# valid_chi returns True if chi uses an initial seqment of the extra 
# colors {qL,...,max(chi)}
cdef bint valid_chi(chi, qL):
    extra_colors = set(xrange(qL, max(chi) + 1))
    return set(chi).intersection(extra_colors) == extra_colors

# Count the monochromatic edges in local view L with local coloring chi
# There are three terms, monochromatic edges incident to v, those from
# u in N(v) to an external neighbor, and those inside N(v)
cdef int mChi(tuple L, list chi):
    cdef int d = len(L[0])
    cdef int ret = 0
    cv = first(chi)
    cN = rest(chi)
    ret = cN.count(cv)
    ret += sum(L[0][j].count(cN[j]) for j in xrange(0, d))
    ret += sum(1 for (x, y) in L[1] if cN[x-1] == cN[y-1])
    return ret

# Count the monochromatic edges incident to v
cdef int mvChi(L, chi):
    cv = first(chi)
    cN = rest(chi)
    return cN.count(cv)

# Count the monochromatic edges incident to neighbors of v
cdef int mNChi(L, chi):
    cdef int d = len(L[0])
    cdef int ret = 0
    cv = first(chi)
    cN = rest(chi)
    ret = cN.count(cv)
    ret += sum(L[0][j].count(cN[j]) for j in xrange(0, d))
    ret += sum(2 for (x, y) in L[1] if cN[x-1] == cN[y-1])
    return ret

# zChi returns the total weight of (a class of) local colourings
cdef zChi(tuple L, int qL, q, lam, list chi):
    color_weight = sbinomial(q-qL, max(max(chi)+1-qL, 0))
    return color_weight * lam ** mChi(L, chi)

# Return the colors used on N(u) where u in {1,2,...,d} is a neighbor of v
cdef list cNu(L, chi, u): 
    cdef int d = len(L[0])
    neighbors = [w for w in xrange(1,d+1) if tuple(sorted([u,w])) in L[1]]
    return [first(chi)] + list(L[0][u-1]) + [chi[w] for w in neighbors]

# Return the frequencies of elements in N in descending order
cdef list H(N):
    return sorted(Counter(N).values(), reverse=True)

#############################################
# Computations for the minimization problem #
#############################################

# Generate the local views found in Kdd from scratch
def gen_Ls_Kdd(d, max_q=None):
    if max_q is None: max_q = d-1
    
    for color_freqs in Partitions(d-1, max_length=max_q):
        B = reduce(lambda a, t: a + [t[0]]*t[1], zip(its.count(0), color_freqs), [])
        yield ((tuple(B),)*d, ())

# Compute properties of local views in Kdd
def gen_props_Kdd(d, q, lam):
    gen = gen_Ls_Kdd(d, max_q=q) if q.is_integer() else gen_Ls_Kdd(d)
    for L in gen:
        yield build_props(L, q, lam)

# Compute properties of Kd1
def props_Kd1(d, q, lam):
    L = (tuple((() for i in xrange(d))), tuple((i, j) for i in xrange(1,d) for j in xrange(i+1,d+1)))
    return build_props(L, q, lam)

# # Solve for dual variables in minimization program. Returns a dictionary of the dual
# # variables Delta[S] in terms of q, lam for S in Ss
# def build_ZDelta(d, ZKdd, ZUKdd, Ss):
#     if d < 3 or len(Ss) < 1: # TODO: figure out correct len(Ss) constraint
#         return {}

#     Delta = list(var('Delta_%d' % i) for i in xrange(len(Ss)))

#     equations = [p['ZUv'] - ZUKdd*p['Z']/ZKdd
#                  - sum(Delta[i]*p['Zgamd', S] for i, S in enumerate(Ss)) == 0
#                  for p in gen_props_Kdd(d)]

#     solutions = solve(equations, *Delta)
#     return {S: factor(ZKdd*solutions[0][i].rhs()) for i, S in enumerate(Ss)}

# Solve for dual variables in minimization program. Returns a dictionary of the dual
# variables Delta[S] in terms of q, lam for S in Ss
def build_ZDelta(d, ZKdd, ZUKdd, Ss, props_Kdd):
    if d < 3 or len(Ss) < 1: # TODO: figure out correct len(Ss) constraint
        return {}

    Delta = list(var('Delta_%d' % i) for i in xrange(len(Ss)))

    equations = [p['ZUv'] - ZUKdd*p['Z']/ZKdd
                 - sum(Delta[i]*p['Zgamd', S] for i, S in enumerate(Ss)) == 0
                 for p in props_Kdd]

    solutions = solve(equations, *Delta)
    return {S: factor(ZKdd*solutions[0][i].rhs()) for i, S in enumerate(Ss)}

##############################
# Properties of a local view #
##############################


def subs_props(p, q, lam, fields):
    for field in fields:
        if p[field] != 0:
            p[field] = p[field].subs(q=q, lam=lam)


# The calculations of properties of a local view are fused into one loop.
cpdef dict build_props(tuple L, q, lam, minslack=None, maxslack=None, factor=False):
    cdef int d = len(L[0])
    cdef int qL = find_qL(L)
    cdef int i, u
    Ss = Partitions(d)

    acc_Z = 0
    acc_ZUv = 0
    acc_ZUN = 0
    acc_gamds = [0] * len(Ss)

    for chi in Tuples(range(qL + d + 1), d + 1):
        if not valid_chi(chi, qL):
            continue

        curr_zChi = zChi(L, qL, q, lam, chi)
        curr_HcNv = H(rest(chi))
        curr_HcNus = [H(cNu(L, chi, u)) for u in xrange(1, d+1)]

        acc_Z  += curr_zChi
        acc_ZUv += mvChi(L, chi) * curr_zChi/2
        acc_ZUN += mNChi(L, chi) * curr_zChi/2/d

        for i, S in enumerate(Ss):
            if curr_HcNv == S:
                acc_gamds[i] += curr_zChi
            for u in xrange(d):
                if curr_HcNus[u] == S:
                    acc_gamds[i] -= curr_zChi/d

    if factor:
        acc_ZUv = factor0(acc_ZUv)
        acc_ZUN = factor0(acc_ZUN)
        acc_Z = factor0(acc_Z)
        for i in xrange(len(Ss)):
            acc_gamds[i] = factor0(acc_gamds[i])

    p = {
        'L':     L,
        'isKdd': is_Kdd(L),
        'isKd1': is_Kd1(L),
        'ZUv':   acc_ZUv,
        'ZUN':   acc_ZUN,
        'Z':     acc_Z
    }
    for i, S in enumerate(Ss):
        p['Zgamd', S] = acc_gamds[i]
    
    p['minslack'] = None if minslack is None else minslack(p) 
    p['maxslack'] = None if maxslack is None else maxslack(p) 

    return p

#############
# Constants #
#############

# https://oeis.org/A000088
num_GNvs    = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864]
# My computations
num_Ls      = [0, 0, 3, 35, 3529, 7017118]
num_Ls_GNv  = [[], [], 
                 [2, 1], 
                 [23, 9, 2, 1],
                 [1636, 1288, 328, 23, 31, 148, 56, 9, 7, 2, 1], 
                 [2142277, 2692431, 957099, 69899, 1636, 44947, 578264, 245465, 26392, 9842, 2033, 109, 155091, 29448, 1288, 36764, 3323, 4557, 328, 190, 23, 1229, 588, 56, 66, 9, 7, 2, 1, 8458, 4557, 148, 560, 31]
                ]
L_filenames = [[], [], 
                 ['data/Ls_d2_.ls', 'data/Ls_d2_12.ls'], 
                 ['data/Ls_d3_.ls', 'data/Ls_d3_12.ls', 'data/Ls_d3_12_13.ls', 'data/Ls_d3_12_13_23.ls'],
                 ['data/Ls_d4_.ls', 'data/Ls_d4_12.ls', 'data/Ls_d4_12_13.ls', 'data/Ls_d4_12_13_14.ls', 'data/Ls_d4_12_13_23.ls', 'data/Ls_d4_12_34.ls', 'data/Ls_d4_12_13_34.ls', 'data/Ls_d4_12_13_14_34.ls', 'data/Ls_d4_12_13_24_34.ls', 'data/Ls_d4_12_13_14_24_34.ls', 'data/Ls_d4_12_13_23_14_24_34.ls'], 
                 ['data/Ls_d5_.ls', 'data/Ls_d5_12.ls', 'data/Ls_d5_12_13.ls', 'data/Ls_d5_12_13_14.ls', 'data/Ls_d5_12_13_14_15.ls', 'data/Ls_d5_12_13_23.ls', 'data/Ls_d5_12_34.ls', 'data/Ls_d5_12_13_34.ls', 'data/Ls_d5_12_13_14_34.ls', 'data/Ls_d5_12_13_24_34.ls', 'data/Ls_d5_12_13_14_24_34.ls', 'data/Ls_d5_12_13_23_14_24_34.ls', 'data/Ls_d5_12_34_15.ls', 'data/Ls_d5_12_13_34_15.ls', 'data/Ls_d5_12_13_14_34_15.ls', 'data/Ls_d5_12_23_34_15.ls', 'data/Ls_d5_12_13_23_34_15.ls', 'data/Ls_d5_12_23_14_34_15.ls', 'data/Ls_d5_12_13_23_14_34_15.ls', 'data/Ls_d5_12_23_14_34_15_35.ls', 'data/Ls_d5_12_13_23_14_34_15_35.ls', 'data/Ls_d5_12_23_34_15_45.ls', 'data/Ls_d5_12_13_23_34_15_45.ls', 'data/Ls_d5_12_13_23_14_34_15_45.ls', 'data/Ls_d5_12_13_23_24_34_15_45.ls', 'data/Ls_d5_12_13_23_14_24_34_15_45.ls', 'data/Ls_d5_12_13_23_24_34_15_25_45.ls', 'data/Ls_d5_12_13_23_14_24_34_15_25_45.ls', 'data/Ls_d5_12_13_23_14_24_34_15_25_35_45.ls', 'data/Ls_d5_12_34_15_25.ls', 'data/Ls_d5_12_13_34_15_25.ls', 'data/Ls_d5_12_13_14_34_15_25.ls', 'data/Ls_d5_12_13_23_34_15_25.ls', 'data/Ls_d5_12_13_23_34_15_25_35.ls']
                ]
P_filenames = [[], [], 
                 ['data/Ps_d2_.dat', 'data/Ps_d2_12.dat'], 
                 ['data/Ps_d3_.dat', 'data/Ps_d3_12.dat', 'data/Ps_d3_12_13.dat', 'data/Ps_d3_12_13_23.dat'],
                 ['data/Ps_d4_.dat', 'data/Ps_d4_12.dat', 'data/Ps_d4_12_13.dat', 'data/Ps_d4_12_13_14.dat', 'data/Ps_d4_12_13_23.dat', 'data/Ps_d4_12_34.dat', 'data/Ps_d4_12_13_34.dat', 'data/Ps_d4_12_13_14_34.dat', 'data/Ps_d4_12_13_24_34.dat', 'data/Ps_d4_12_13_14_24_34.dat', 'data/Ps_d4_12_13_23_14_24_34.dat'], 
                 ['data/Ps_d5_.dat', 'data/Ps_d5_12.dat', 'data/Ps_d5_12_13.dat', 'data/Ps_d5_12_13_14.dat', 'data/Ps_d5_12_13_14_15.dat', 'data/Ps_d5_12_13_23.dat', 'data/Ps_d5_12_34.dat', 'data/Ps_d5_12_13_34.dat', 'data/Ps_d5_12_13_14_34.dat', 'data/Ps_d5_12_13_24_34.dat', 'data/Ps_d5_12_13_14_24_34.dat', 'data/Ps_d5_12_13_23_14_24_34.dat', 'data/Ps_d5_12_34_15.dat', 'data/Ps_d5_12_13_34_15.dat', 'data/Ps_d5_12_13_14_34_15.dat', 'data/Ps_d5_12_23_34_15.dat', 'data/Ps_d5_12_13_23_34_15.dat', 'data/Ps_d5_12_23_14_34_15.dat', 'data/Ps_d5_12_13_23_14_34_15.dat', 'data/Ps_d5_12_23_14_34_15_35.dat', 'data/Ps_d5_12_13_23_14_34_15_35.dat', 'data/Ps_d5_12_23_34_15_45.dat', 'data/Ps_d5_12_13_23_34_15_45.dat', 'data/Ps_d5_12_13_23_14_34_15_45.dat', 'data/Ps_d5_12_13_23_24_34_15_45.dat', 'data/Ps_d5_12_13_23_14_24_34_15_45.dat', 'data/Ps_d5_12_13_23_24_34_15_25_45.dat', 'data/Ps_d5_12_13_23_14_24_34_15_25_45.dat', 'data/Ps_d5_12_13_23_14_24_34_15_25_35_45.dat', 'data/Ps_d5_12_34_15_25.dat', 'data/Ps_d5_12_13_34_15_25.dat', 'data/Ps_d5_12_13_14_34_15_25.dat', 'data/Ps_d5_12_13_23_34_15_25.dat', 'data/Ps_d5_12_13_23_34_15_25_35.dat']
                ]


####################
# Helper functions #
####################

rest = lambda l: l[1:]

def positive_poly(f):
    num, den = f.numerator_denominator()
    if den != 1:
        return "denominator problem:\n%s" % den
    else:
        coeffs = num.polynomial(ZZ).coefficients()
        if any(c < 0 for c in coeffs):
            return "numerator problem:\n%s" % num
    return ""


cpdef int cfact(int n):
    cdef int i, ret
    ret = 1
    for i in xrange(1, n+1):
        ret *= i
    return ret

cpdef int cfallfact(int n, int k):
    cdef int i, ret
    ret = 1
    for i in xrange(n, n-k, -1):
        ret *= i
    return ret

cpdef int cprod(list l):
    cdef int ret
    ret = 1
    for n in l:
        ret *= n
    return ret

cpdef int cbinomial(int n, int k):
    return cfallfact(n, k)/cfact(k)

cpdef int cmultinomial(list l):
    cdef int total, den, i
    den = 1
    total = 0
    for i in l:
        den *= cfact(i)
        total += i
    return cfact(total)/den

# give the number of partitions of a set size sum(L) into parts with sizes
# given by the elements of l
cpdef int num_set_partitions(list l):
    cdef int ret
    counts = Counter(l)
    ret = cmultinomial(l)
    for i in counts.itervalues():
        ret /= cfact(i)
    return ret   
