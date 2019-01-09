import urllib
from .vector_ops import *
from .config import NUMBER_PRECISION_FUNCTION
from .kpclasses import Bead
from numpy import array as Vector
import re


def download_from_pdb(pdbid,savename=''):
    url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={}".format(pdbid)
    urllib.urlretrieve(url,savename if savename else "{}.pdb".format(pdbid))


def read_from_pdb(filename,selected_chain=''):
    a = []
    last = None
    last_chain = None
    Calpha = []
    if selected_chain: selected_chain = list(selected_chain)
    if len(selected_chain)>1:
        for c in selected_chain:
            Calpha.append(re.compile("ATOM  .{7}CA.{6}"+c+"([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8}).{23}C"))
    elif selected_chain:
        Calpha.append(re.compile("ATOM  .{7}CA.{6}"+selected_chain[0]+"([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8}).{23}C"))
    else:
        Calpha.append(re.compile("ATOM  .{7}CA.{7}([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8}).{23}C"))

    with open(filename) as input:
        for line in input:
#            if "END":
#                print "end"
            if "ENDMDL" in line:
#                print "breakin on endlmld"
                break
            for Ca in Calpha:
                if Ca.match(line):
                    if last_chain is not None and last_chain != line[21]:
                        selected_chain.pop(selected_chain.index(last_chain))
                        if not selected_chain:
                            break
                        else:
                            a[-1].end = True
                    g = Ca.findall(line)[0]
                    if last != int(g[0]):
                        new_vec = Vector(map(NUMBER_PRECISION_FUNCTION,g[1:]))
                        last = int(g[0])
                        new_atom = Bead(new_vec,"CA")
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
    for i,at in enumerate(a):
        at.setId(i)
        if i>0:
            at.setNhand(a[i-1])
        if i<len(a)-1:
            at.setChand(a[i+1])
    atoms = a
    return atoms

def guess_format(fname):
    Ca = re.compile("ATOM  .{7}CA.{7}([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8}).{23}C")
    xyz = re.compile("^\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+$")
    with open(fname) as input:
        for line in file:
            if Ca.match(line):
                return "pdb"
            if xyz.match(line):
                return "xyz"
        return "pdb"

def read_from_xyz(filename,save=False):#,start,stop):
    a=[]
    with open(filename) as input:
        for line in input:
            if not line: continue
            if line.strip() == "END":
                a[-1].end = True
                continue
            if line.strip() == "SOFTEND":
                continue
            line = map(NUMBER_PRECISION_FUNCTION,line.split()[1:])
            new_vec = Vector(line)
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
    #    a = a[start:stop]
    for i,at in enumerate(a):
        at.setId(i)
        if i>0:
            at.setNhand(a[i-1])
        if i<len(a)-1:
            at.setChand(a[i+1])
    atoms = a
    return atoms
