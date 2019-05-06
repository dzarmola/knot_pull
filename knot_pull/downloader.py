from __future__ import print_function
#from future.standard_library import install_aliases
#install_aliases()
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError, URLError
except ImportError:
    from urllib2 import urlopen, URLError, HTTPError

import re

def _urlopen(address):
    x = None
    for i in range(3):
        try:
            response = urlopen(address)
            return response
        except HTTPError as e:
            x = e
            continue
        except URLError as e:
            x = e
            continue
    else:
        raise URLError("Error trying to connect to {}:\n{}".format(address,x))


def check_if_in_PDB(pdbId,chain):
    address="http://www.rcsb.org/pdb/rest/describeMol?structureId={PDBID}.{CHAIN}".format(CHAIN=chain,PDBID=pdbId)
    response = _urlopen(address)
    html = response.read()
    return "polymer" in html

def is_a_large_structure(pdbId):
    address="http://www.rcsb.org/pdb/rest/describeMol?structureId={PDBID}".format(PDBID=pdbId)
    response = _urlopen(address)
    html = response.read()
    return 'largeStructure="true"' in html

def get_ligands(pdbId,chain):
    address = "http://www.rcsb.org/pdb/rest/ligandInfo?structureId={PDBID}.{CHAIN}".format(CHAIN=chain, PDBID=pdbId)
    response = _urlopen(address)
    html = response.read()
    ligands = re.findall('chemicalID="([A-Z0-9]*)" type="non-polymer"',html.decode('utf-8'))
    return ligands

def getCifIndices(cif_read):
    mylines = [line.strip() for line in cif_read.split("\n") if re.search("^_atom_site\.",line)]
    fields = ["_atom_site.group_PDB","_atom_site.id","_atom_site.type_symbol","_atom_site.label_atom_id",
              "_atom_site.label_alt_id","_atom_site.label_comp_id","_atom_site.label_entity_id",
              "_atom_site.label_seq_id","_atom_site.pdbx_PDB_ins_code","_atom_site.Cartn_x","_atom_site.Cartn_y",
              "_atom_site.Cartn_z","_atom_site.occupancy","_atom_site.B_iso_or_equiv","_atom_site.auth_asym_id"]
    return [mylines.index(i) for i in fields]

def multiline_cif2pdb(coords, indices):
    prev = False
    acnt = 1
    rcnt = 0
    new_coords = []
    for line in coords:
        cs = line.split()
        if not indices:
            head,snum,elem,aname,alt,rname,_,ent,seqid,ins,X,Y,Z,occ,Bfac = cs[:15]
            chain = cs[23]
        else:
            head, snum, elem, aname, alt, rname, ent, seqid, ins, X, Y, Z, occ, Bfac,chain = [cs[i].strip().strip('"') for i in indices]
        #snum = str(acnt)
        if seqid!=prev:
            rcnt+=1
            prev=seqid
        #seqid = str(rcnt)
        acnt += 1
        #new = head.ljust(6) + snum.rjust(5) + " " + aname.rjust(4) + alt.replace(".", " ") + rname.rjust(3) + chain.rjust(
        #    2) + seqid.rjust(4) \
        #      + " " * 4 + X.rjust(8) + Y.rjust(8) + Z.rjust(8) + occ.rjust(6) + Bfac.rjust(6) + " " * 10 + elem.rjust(2)
        #new = map(float,(X,Y,Z))
        if elem == "C" and aname == "CA" and alt in "A.":
            new_coords.append((chain,(X,Y,Z)))
    return new_coords

def get_all_chains(pdbId):
    return get_particular_chain(pdbId,".")

def get_from_afar(pdbId, chain=''):
    if chain:
        return get_particular_chain(pdbId,chain)
    else:
        return get_all_chains(pdbId)

def get_particular_chain(pdbId,chain):
    address = "https://files.rcsb.org/view/{PDBID}.cif".format(PDBID=pdbId)
    response = _urlopen(address)
    html = response.read().decode('utf-8')
    cif_indices = getCifIndices(html)

    html_coords = html.split("_atom_site.pdbx_PDB_model_num")[1].split("#")[0]

    atom_coords = [line for line in html_coords.split("\n") if len(line.split())>2 and line.split()[-1] == "1"  and line.split()[cif_indices[7]].strip()!="."] #wywalam and line.split()[0]=="ATOM" and (chain=="." or line.split()[cif_indices[-1]]==chain)
    #####removal of ligands
    ligands = get_ligands(pdbId,(chain if chain and chain!="." else ""))
    atom_coords = [x for x in atom_coords if x.split()[cif_indices[5]] not in ligands]
    atom_coords = multiline_cif2pdb(atom_coords,cif_indices)
    if chain!=".":
        return filter(lambda x:x[0] in chain, atom_coords)
    return atom_coords

def get_chain_list(pdbId):
    address = "https://files.rcsb.org/view/{PDBID}.cif".format(PDBID=pdbId)
    response = _urlopen(address)
    html = response.read()
    html = html.split("_pdbx_poly_seq_scheme")[-1].split("#")[0].split("\n")[1:]
    chains = [line.split()[-3] for line in html if line]
    out = []
    for i in chains:
        if not i in out:
            out.append(i)
    return out

#def get_crystal_data(pdbId,path=''):
#        address = "https://files.rcsb.org/view/{PDBID}.pdb".format(PDBID=pdbId)
#        response = _urlopen(address)
#        crystal = []
#        for i in response.readlines():
#                if re.search('^(CRYST1|ORIGX.|SCALE.)' , i):
#                        crystal.append(i.strip())
#        if path and crystal:
#                with open("{}/CRYSTAL.pdb".format(path),"w",0) as out:
#                        out.write("\n".join(crystal))
#                return True
#        return False

if __name__=="__main__":
    if check_if_in_PDB("4mcb","A"):
        #print get_ligands("4mcb","A")
        print (get_particular_chain("4mcb","A"))
        print (get_all_chains("4mcb"))
    #print check_if_in_PDB("4ola", "X")