from __future__ import print_function
from .vector_ops import *
from .config import NUMBER_PRECISION_FUNCTION,VERBOSE
from .kpclasses import Bead
from .downloader import get_from_afar
from numpy import array as Vector
import re


#def download_from_pdb(pdbid,savename=''):
#    url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={}".format(pdbid)
#    urllib.urlretrieve(url,savename if savename else "{}.pdb".format(pdbid))

def read_from_web(pdbid,selected_chains='',begin=None,end=None):
    coords,err = get_from_afar(pdbId=pdbid,chain=selected_chains)
    if err:
        return [],[],err
    a=[]
    curc = None
    cnames = []
    bid = bnum = eid = enum = None
    if begin is not None:
        if 'i' in begin:
            bid = int(begin.replace("i", ""))
        else:
            bnum = int(begin)
    if end is not None:
        if 'i' in end:
            eid = int(end.replace("i", ""))
        else:
            enum = int(end)
    count = 0

    first_rid = None
    last_rid = None

    for atom in coords:
        ch,rid,pos = atom
        if first_rid is None:
            first_rid = rid
        #if selected_chains and ch not in selected_chains: # Done already in get_from_afar
        #    continue
        count += 1
        if eid is not None:
            if rid > eid:
                break
        if enum is not None:
            if count > enum:
                break
        if bid is not None:
            if rid < bid:
                continue
        if bnum is not None:
            if count < bnum:
                continue
        if curc is not None and ch != curc:
            if begin is not None or end is not None:
                #last_rid = rid
                break
            a[-1].end = True
        if ch not in cnames:
            cnames.append(ch)
        pos = list(map(NUMBER_PRECISION_FUNCTION,pos))
        new_vec = Vector(pos)
        new_atom = Bead(new_vec,"CA")
#        nv = get_middlepoint(a[-1].vec,new_vec)
        if a and (a[-1].end != True) and point_distance(a[-1].vec,new_vec)>4.:
            pd = point_distance(a[-1].vec,new_vec)
            #cd = pd/4
            l = a[-1].vec
            middle = get_middlepoint(l,new_vec)
            if pd>12.:
                lm = get_middlepoint(l,middle)
                a.append(Bead(Vector(lm),"CA"))
                a.append(Bead(Vector(middle),"CA"))
                rm = get_middlepoint(new_vec,middle)
                a.append(Bead(Vector(rm),"CA"))
            else:
                a.append(Bead(Vector(middle),"CA"))
        a.append(new_atom)
        curc = ch
    last_rid = rid
    for i,at in enumerate(a):
        at.setId(i)
        if i>0:
            at.setNhand(a[i-1])
        if i<len(a)-1:
            at.setChand(a[i+1])
    atoms = a
    err = ""
    if not a:
        err = "Begin/end values ({},{}) are incorrect for chain {} with length {} and residue numbers ({},{})".format(
            begin,end,cnames[0],count,first_rid,last_rid)
    return atoms,cnames,err


def read_from_pdb(filename,selected_chains='',begin=None,end=None,rna=False):
    a = []
    last = None
    last_chain = None
    Calpha = []
    cnames = []
    if rna:
        atom_to_be = 'P '
    else:
        atom_to_be = 'CA'

    if selected_chains:#len(selected_chain)>1:
        for c in selected_chains:
            Calpha.append(re.compile("ATOM  .{7}"+atom_to_be+".{6}"+c+"([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8})"))#).{23}C"))
    else:
        Calpha.append(re.compile("ATOM  .{7}"+atom_to_be+".{7}([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8})"))#.{23}C"))

    bid = bnum = eid = enum = None
    if begin is not None:
        if 'i' in begin:
            bid = int(begin.replace("i", ""))
        else:
            bnum = int(begin)
    if end is not None:
        if 'i' in end:
            eid = int(end.replace("i", ""))
        else:
            enum = int(end)
    count = 0

    first_rid = None
    last_rid = None

    found_chain = False

    _avg_dist = 0
    _avg_cnt = 0

    with open(filename) as input:
        for line in input:
            if "ENDMDL" in line:
                break
            for Ca in Calpha:
                if Ca.match(line):
                    found_chain = True
                    if line[21] not in cnames: cnames.append(line[21])
                    if last_chain is not None and last_chain != line[21]:
                        if begin is not None or end is not None:
                            break
                        if selected_chains:
                            selected_chains.pop(selected_chains.index(last_chain))
                            if not selected_chains:
                                break
                        else:
                            a[-1].end = True
                    g = Ca.findall(line)[0]
                    rid = int(g[0])
                    if last != rid:
                        last_rid = rid
                        if first_rid is None:
                            first_rid = rid
                        count += 1
                        if eid is not None:
                            if rid > eid:
                                break
                        if enum is not None:
                            if count > enum:
                                break
                        if bid is not None:
                            if rid < bid:
                                continue
                        if bnum is not None:
                            if count < bnum:
                                continue
                        new_vec = Vector(list(map(NUMBER_PRECISION_FUNCTION,g[1:])))
                        last = int(g[0])
                        new_atom = Bead(new_vec,"CA")
                        if a and (a[-1].end != True):
                            _avg_dist += point_distance(a[-1].vec, new_vec)
                            _avg_cnt += 1
                        if a and (a[-1].end != True) and point_distance(a[-1].vec, new_vec) > 4.:
                            pd = point_distance(a[-1].vec, new_vec)
                            # cd = pd/4
                            l = a[-1].vec
                            middle = get_middlepoint(l, new_vec)
                            if pd > 12.:
                                lm = get_middlepoint(l, middle)
                                a.append(Bead(Vector(lm), "CA"))
                                a.append(Bead(Vector(middle), "CA"))
                                rm = get_middlepoint(new_vec, middle)
                                a.append(Bead(Vector(rm), "CA"))
                            else:
                                a.append(Bead(Vector(middle), "CA"))
                        a.append(new_atom)
                    last_chain = line[21]

    if VERBOSE:
        print ("Average distance in file is: {}".format(_avg_dist / _avg_cnt))

    err = ""
    if not a:
        if found_chain:
            err = "Begin/end values ({},{}) are incorrect for chain {} with length {} and residue numbers ({},{})".format(
            begin,end,cnames[0],count,first_rid,last_rid)
        else:
            err = "No chain(s) {} found!".format(selected_chains)
    for i,at in enumerate(a):
        at.setId(i)
        if i>0:
            at.setNhand(a[i-1])
        if i<len(a)-1:
            at.setChand(a[i+1])
    atoms = a
    return atoms,cnames,err

def guess_format(fname):
    Ca = re.compile("ATOM  .{7}CA.{7}([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8}).{23}C")
    xyz = re.compile("^\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+$")
    with open(fname) as input:
        for line in input:
            if Ca.match(line):
                return "pdb"
            if xyz.match(line):
                return "xyz"
        return "pdb"


def read_from_xyz(filename, begin=None, end = None):#,start,stop):
    a=[]
    ccount = 1
    bid = bnum = eid = enum = None
    if begin is not None:
        if 'i' in begin:
            bid = int(begin.replace("i", ""))
        else:
            bnum = int(begin)
    if end is not None:
        if 'i' in end:
            eid = int(end.replace("i", ""))
        else:
            enum = int(end)
    count = 0

    first_rid = None
    last_rid = None

    proper_line = re.compile("^\s*([0-9-]+)\s+(-*[0-9.]+)\s+(-*[0-9.]+)\s+(-*[0-9.]+)")

    found_smthing = False

    with open(filename) as input:
        for line in input:
            if not line: continue
            if line.strip() == "END":
                if begin is not None or end is not None:
                    break
                a[-1].end = True
                ccount += 1
                continue
            if proper_line.match(line):
                found_smthing = True
                line = proper_line.findall(line)[0]
                rid = int(line[0])
                pos = list(map(NUMBER_PRECISION_FUNCTION, line[1:]))
                if first_rid is None:
                    first_rid = rid
                # if selected_chains and ch not in selected_chains: # Done already in get_from_afar
                #    continue
                count += 1
                if eid is not None:
                    if rid > eid:
                        break
                if enum is not None:
                    if count > enum:
                        break
                if bid is not None:
                    if rid < bid:
                        continue
                if bnum is not None:
                    if count < bnum:
                        continue
                #if line.strip() == "SOFTEND":
                #    continue
                last_rid = rid
                new_vec = Vector(pos)
                new_atom = Bead(new_vec,"CA")
        #        nv = get_middlepoint(a[-1].vec,new_vec)
                if a and (a[-1].end!=True) and point_distance(a[-1].vec,new_vec)>4.:
                    pd = point_distance(a[-1].vec,new_vec)
                    #cd = pd/4
                    l = a[-1].vec
                    middle = get_middlepoint(l,new_vec)
                    if pd>12.:
                        lm = get_middlepoint(l,middle)
                        a.append(Bead(Vector(lm),"CA"))
                        a.append(Bead(Vector(middle),"CA"))
                        rm = get_middlepoint(new_vec,middle)
                        a.append(Bead(Vector(rm),"CA"))
                    else:
                        a.append(Bead(Vector(middle),"CA"))
                a.append(new_atom)
    err = ''
    if not a:
        if found_smthing:
            err = "Begin/end values ({},{}) are incorrect for chain with length {} and residue numbers ({},{})".format(
            begin,end,count,first_rid,last_rid)
        else:
            err = "Improper XYZ file format!"

    for i,at in enumerate(a):
        at.setId(i)
        if i>0:
            at.setNhand(a[i-1])
        if i<len(a)-1:
            at.setChand(a[i+1])
    atoms = a
    return atoms, range(ccount),err
