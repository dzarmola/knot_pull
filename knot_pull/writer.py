from .config import CHAINZ


def print_out_last(handle,frames):
    with open(handle,"a",0) as out:
        frame = frames[-1]
        out.write("t {}\n".format(len(frames)))
        for j,a in enumerate(frame):
            out.write("%d\t%3.5f\t%3.5f\t%3.5f\n" % (j+1,a[0],a[1],a[2]))
        out.write("\n")

def print_out_last_pdb(handle,frames,atoms):

    x=0
    #with open(handle,"a",0) as out:
    out = handle
    frame = frames[-1]
    out.write("MODEL {}\n".format(len(frames)))
    prev_id = None
    for j,(id,a,end) in enumerate(frame):
        id=int(id)
        chain = CHAINZ[x]
        if prev_id is not None:
            for _ in range(prev_id+1,id):
                out.write("ATOM  % 5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f                       C\n" % (_, chain ,_,a[0],a[1],a[2]))
        out.write("ATOM  % 5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f                       C\n" % (id, chain ,id,a[0],a[1],a[2]))
        prev_id = id
        if end:
            x+=1
    out.write("ENDMDL\n")

def print_out_one_frame(out,atoms,len_fr,chain,rna=False):
    out.write("MODEL {}\n".format(len_fr))
    if rna:
        atom_to_be,elem = 'P ','P'
    else:
        atom_to_be,elem = 'CA','C'
    for j,atom in enumerate(atoms):
        a=atom.vec
        out.write("ATOM  % 5d  %s  ALA %s%4d    %8.3f%8.3f%8.3f                       %s\n" % (j+1, atom_to_be,chain ,j+1,a[0],a[1],a[2],elem))
#            if atoms[j].end:
#                x+=1
    out.write("ENDMDL\n")


def print_out_all(handle,frames):
    with open(handle,"w",0) as out:
        for t,frame in enumerate(frames):
            out.write("t {}\n".format(t+1))
            for j,a in enumerate(frame):
                out.write("%d\t%3.5f\t%3.5f\t%3.5f\n" % (j+1,a[0],a[1],a[2]))
            out.write("\n")

def write_xyz(out,atoms,soft_end=False):
#    with open(handle,"w",0) as out:
    for i,bead in enumerate(atoms):
        out.write("{} {} {} {}\n".format(i+1,bead.x,bead.y,bead.z))
        if bead.end:
            out.write("END\n")
    if soft_end:
        out.write("SOFTEND\n")