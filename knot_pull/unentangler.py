from __future__ import print_function
import os,tempfile
import subprocess

from .vector_ops import *
from .kpclasses import chainDeepCopy
from .dowker_code_new import dowker_code
from .config import VERBOSE


def run_through_hashes_wanda(atoms_lists):
    """For each atom list """
    out = []
    for x, atoms in enumerate(atoms_lists):
        if len(atoms) < 5:
            out.append("01")
            continue
        _kfp1 = check_with_wanda(atoms)
        if int(_kfp1.strip("-").strip("+")[0]) > 6:
            _inner_kfp = list(map(str, in_depth_search(atoms)))
            if len(_inner_kfp) > 1:
                 _kfp1 += "-{}(x{})".format(_inner_kfp[0], " # ".join(_inner_kfp[1:]))
        out.append(_kfp1)
    return " # ".join(out)

def run_through_hashes_dowker(atoms_lists):
    """For each atom list """
    out = []
    out_full = []
    for x, atoms in enumerate(atoms_lists):
        if len(atoms) < 5:
            out.append("01")
            out.append("[]")
            continue
        dowker,code = dowker_code(atoms)
        out.append(dowker)
        out_full.append(str(code))
    return " # ".join(out)," # ".join(out_full)

def neigh(i, j, n):
    ns = []
    for x, y in [(i - 1, j), (i - 1, j + 1), (i, j + 1)]:
        if x >= 0 and y >= 0 and x < n and y < n and not (x == i and y == j):
            ns.append((x, y))
    return ns


def new_cell(out, i, j):
    ns = neigh(i, j, len(out))
    d = {}
    for x, y in ns:
        d[out[x][y]] = d.get(out[x][y], 0) + 1
    d[0] = 0
    new = sorted(d.items(), key=lambda _: _[1], reverse=True)[0]
    new = new[0] if new[1] == 3 else out[i][j]  # new[1]>=4
    return new

def knot_check_fnc_wanda(atoms):
    _atoms = chainDeepCopy(atoms)
    _kfp = check_with_wanda(_atoms)
    return _kfp

def knot_check_fnc_dowker(atoms):#deprecated
    _atoms = chainDeepCopy(atoms)
    _kfp = dowker_code(_atoms)
    return "{}0".format(len(_kfp))

def pretty_print(matrix):
    plen = max(max(len(str(x)) for x in y) for y in matrix)
    for row in matrix:
        print (" ".join(map(lambda x: str(x).rjust(plen),row)))


def in_depth_search(atoms,knot_check_fnc=knot_check_fnc_wanda):
    """Makes a matrix of results to find overlapping topologies"""
    out = [[0 for x in range(len(atoms))] for y in range(len(atoms))]
    out_set = set([])
    for i in range(len(atoms) - 2):
        for j in range(i, len(atoms)):
            if j - i >= 5:
                _kfp = knot_check_fnc(atoms[i:j])
                out_set.add(_kfp)
                out[i][j] = 999 if _kfp[:3] == "999" else int(_kfp) if _kfp[0] != "0" else 0
    for i in range(len(atoms) - 2):
        for j in range(len(atoms) - 1, i, -1):
            if not out[i][j]:
                out[i][j] = new_cell(out, i, j)
    largest = {}
    for i in range(len(atoms) - 2):
        for j in range(len(atoms) - 1, i, -1):
            _tmp = set([_ for x in out[:i + 1] for _ in x[j:]])
            if 0 not in _tmp and (len(_tmp) == 1 or (len(_tmp) == 2 and "#" in _tmp)):
                if not (i + 1, len(out) - j) in largest:
                    largest[(i + 1, len(out) - j)] = out[i][j]
                    for _x in range(i + 1):
                        for _y in range(j, len(out)):
                            out[_x][_y] = "#"
    found = {}
    for (i, j), k in sorted(largest.items(), key=lambda _: _[0][0] * _[0][1], reverse=True):
        for (x, y) in found.keys():
            if i <= x and j <= y and k == found[(x, y)]:
                break
        else:
            found[(i, j)] = k
    response = [x[1] for x in sorted(found.items(), key=lambda _: _[0][0] * _[0][1])]
    return response


def write_xyz(atoms,fname):
    for i,a in enumerate(atoms):
        fname.write("{} {} {} {}\n".format(i,a.vec[0],a.vec[1],a.vec[2]))

def check_with_wanda(atom_list):
    """run programKnot without KMT to determine topology of a substructure"""
    with tempfile.NamedTemporaryFile() as atom_file, tempfile.NamedTemporaryFile() as output:
        write_xyz(atom_list, atom_file)
        atom_file.flush()
        p = subprocess.Popen(
            "./programKnot_no_simpler {fname} 0 1 | grep -ve 'HEAD' -ve '^$' > {out_fname}".format(fname=atom_file.name,
                                                                                                   out_fname=output.name),
            shell=True)
        p.wait()  # _ = raw_input("wcisnij cos")
        kfp = output.read()
        try:
            int(kfp.split()[0]) #knot output should be an int
        except:
            print ("Cannot parse",kfp)
            print (atom_file.name, output.name)
            input("czekam na reakcje - przeczytaj outputy zanim je usune")
        kfp = kfp.split()[0]
    return kfp

def make_chain_neighbourhoods(chains):
    neigh = {c.chain:[_.chain for _ in c.neighbours] for c in chains}
    trans = {c.chain:c for c in chains}
    all_chains = list(trans.keys())
    neighbourhoods = []
    while all_chains:
        queue = [all_chains.pop(0)]
        n = []
        while queue:
            cur = queue.pop(0)
            n.append(cur)
            for _ in neigh[cur]:
                if _ not in n and _ not in queue:
                    queue.append(_)
        neighbourhoods.append([trans[x] for x in n])
        for x in n[1:]:
            all_chains.pop(all_chains.index(x))
    return neighbourhoods

def unentangle(chains,outfile=0):
    """Finds all non-interlocking segments, then runs knot detection"""
    for chain in chains:
        #chain.wanda = run_through_hashes_wanda(chain.atom_lists)
        chain.dowker,chain.dowker_code = run_through_hashes_dowker(chain.atom_lists)
        if VERBOSE: print ("Dowker for chain {}: {}".format(chain.chain,chain.dowker))

    neighs = make_chain_neighbourhoods(chains)
    if VERBOSE:
        print ("Chain neighbourhoods:", neighs)
    #W = " U ".join("  #  ".join(x.wanda for x in neigh) for neigh in neighs)
    D = "  U  ".join(" # ".join("[{}]{}".format(x.dowker,x.chain_name) for x in neigh) for neigh in neighs)
    C = "  U  ".join(" # ".join("[{}]{}".format(x.dowker_code,x.chain_name) for x in neigh) for neigh in neighs)
    return D,C
