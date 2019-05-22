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

def read_from_cif(name, config, online=True):
    def getCifIndices(cif_read):
        mylines = [line.strip() for line in cif_read.split("\n") if re.search("^_atom_site\.", line)]
        fields = ["_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
                  "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_entity_id",
                  "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
                  "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.auth_asym_id",
                  "_atom_site.auth_seq_id"]
        return [mylines.index(i) for i in fields]

    if online:
        data = get_from_afar(pdbId=name)
    else:
        with open(name) as input:
            data = input.read()

    atom_name = config['atom_name']
    selected_chains = config['chains']


    cif_indices = getCifIndices(data)
    # html_coords = html.split("_atom_site.pdbx_PDB_model_num")[1].split("#")[0]
    _coords = re.compile("_atom_site.pdbx_PDB_model_num([^#]+?)#")
    _ligands = re.compile("([A-Z0-9']{3}) non-polymer +. [0-9A-Z' ]+")
    html_coords = _coords.findall(data)[0]
    atom_coords = [line for line in html_coords.split("\n") if len(line.split()) > 2 and line.split()[-1] == "1"
                   and line.split()[cif_indices[7]].strip() != "."]
    # ligands = get_ligands(pdbId,(chain if chain and chain != "." else ""))
    ligands = _ligands.findall(data)

    atom_coords = [x for x in atom_coords if x.split()[cif_indices[5]] not in ligands]

    prev = False
    acnt = 1
    rcnt = 0
    new_coords = []
    for line in atom_coords:
        cs = line.split()
        head, snum, elem, aname, alt, rname, ent, seqid, ins, X, Y, Z, occ, Bfac, chain, resid = [
            cs[i].strip().strip('"') for i in cif_indices]
        if seqid != prev:
            rcnt += 1
            prev = seqid
        acnt += 1
        if aname == atom_name and not (atom_name == "CA" and elem != "C") and alt in "A.":
            new_coords.append((chain, int(resid), (X, Y, Z)))
    err = ''
    _chains = set(_[0] for _ in new_coords)
    if not selected_chains:
        coords = new_coords
    elif not any(c in _chains for c in selected_chains):
        return [],[], "No such chain(s): {} in the structure (only {} present)".format(selected_chains, _chains)
    else:# selected_chains != '.':
        coords = filter(lambda x: x[0] in selected_chains, new_coords)


    begin,end = config['begin'],config['end']
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


def read_from_pdb(filename,config):#selected_chains='',begin=None,end=None,rna=False,atom_to_be=None):
    selected_chains, begin,end,rna,atom_to_be = config['chains'],config['begin'],config['end'],config['rna'],config['atom_name']
    a = []
    last = None
    last_chain = None
    cnames = []

    if atom_to_be is None:
        if rna:
            atom_to_be = 'P '
        else:
            atom_to_be = 'CA'
    Calpha = re.compile("ATOM  .{5} ([A-Z0-9' ]{4}).{5}(.)([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8})")

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
            if Calpha.match(line):
                atom, chain, resid, x,y,z = Calpha.findall(line)[0]
                if atom.strip() == atom_to_be.strip() and (not selected_chains or chain in selected_chains):
                    found_chain = True
                    if chain not in cnames: cnames.append(chain)
                    if last_chain is not None and last_chain != chain:
                        if begin is not None or end is not None:
                            break
                        if selected_chains:
                            selected_chains.pop(selected_chains.index(last_chain))
                            if not selected_chains:
                                break
                        else:
                            a[-1].end = True
                    rid = int(resid)
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
                        new_vec = Vector(list(map(NUMBER_PRECISION_FUNCTION,[x,y,z])))
                        last = rid
                        new_atom = Bead(new_vec,"CA",original_id=rid)
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
                    last_chain = chain

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
    Ca = re.compile("ATOM  .{7}CA.{7}([0-9 ]{4}).{4}([0-9\. -]{8})([0-9\. -]{8})([0-9\. -]{8})")
    xyz = re.compile("^\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+\s+-*\d+\.*\d+$")
    cif = re.compile("^data_[A-Z0-9]{4}")
    with open(fname) as input:
        for line in input:
            if cif.match(line):
                return "cif"
            if Ca.match(line):
                return "pdb"
            if xyz.match(line):
                return "xyz"
        return "pdb"


def read_from_xyz(filename, config):#begin=None, end = None):#,start,stop):
    a=[]
    ccount = 1
    begin, end = config['begin'], config['end']
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
                new_atom = Bead(new_vec,"CA",original_id=rid)
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
