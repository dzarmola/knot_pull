from __future__ import print_function
from numpy import array as Vector
from collections import namedtuple

from .vector_ops import *
from .kpclasses import Bead,chainDeepCopy
from .writer import *
from .reader import *
from .config import EPSILON,VERBOSE
from knot_pull.finder import find_frame


def too_far(p1 ,p2 ,p3 ,i0 ,i1):
    """Checks if crossing is possible - 9 A from all triangle points is impossible due to single connection
    length restrictions"""
    return all([point_distance(p1 ,i0 )>=9. ,point_distance(p2 ,i0 ) >9. ,point_distance(p3 ,i0 ) >9.]) \
            or all([point_distance(p1 ,i1 )>=9. ,point_distance(p2 ,i1 ) >9. ,point_distance(p3 ,i1 ) >9.])


def inspectingOtherLines(p1, p2, p3, i0, i1, greedy=0, log=0,but_not=[]):
    """Looks for connections going through the triangle. Starts with second connection from each triangle end"""

    #but_not = [x.vec for x in but_not]
    potentials = []
    potids = []

    p123normal = get_normal(p1, p2, p3)

    inspected = i0
    if inspected and not inspected.Nhand:
        if found_crossing(p1,p2,p3, inspected.vec, l1=None, p123normal=None) is not False:
            potentials.append(inspected.vec)

    while inspected and inspected.Nhand:
        if inspected.Nhand.end or (inspected in but_not and inspected.Nhand in but_not):#Skip "Virtual" connection - between chain ends
            inspected = inspected.Nhand
            continue
        l0, l1 = inspected.vec, inspected.Nhand.vec
        #greedy can have distances larger than allowed in too_far()
        i = False if (greedy == 0 and too_far(p1, p2, p3, l0, l1)) else found_crossing(p1, p2, p3, l0, l1, p123normal)

        if i is not False:
            potentials.append(i)
        inspected = inspected.Nhand

    inspected = i1
    if inspected and not inspected.Chand:
        if found_crossing(p1,p2,p3, inspected.vec, l1=None, p123normal=None) is not False:
            potentials.append(inspected.vec)
    while inspected and inspected.Chand:
        if inspected.end or (inspected in but_not and inspected.Chand in but_not):
            inspected = inspected.Chand
            continue
        l0, l1 = inspected.vec, inspected.Chand.vec
        #greedy can have distances larger than allowed in too_far()
        i = False if (greedy == 0 and too_far(p1, p2, p3, l0, l1)) else found_crossing(p1, p2, p3, l0, l1, p123normal)
        if i is not False:
            potentials.append(i)
        inspected = inspected.Chand
    return potentials

def filter_greedy(atoms, outfile='', len_frames=0):
    """Removes all superfluous atoms(kept for PyMOL purposes, to better detect topology"""
    current = atoms[0]
    latom = 0 #num of removed atoms
    changed = False
    while current.Chand and current.Chand.Chand:
        if current.Chand.end: #Skips "virtual" connections at chain ends
            current = current.Chand.Chand
            continue
        if current.end:
            current = current.Chand
            continue

        i0 = current.Nhand
        i1 = current.Chand.Chand.Chand if current.Chand and current.Chand.Chand else False

        potentials = inspectingOtherLines(current.vec, current.Chand.vec, current.Chand.Chand.vec, i0, i1, greedy=True)

        if potentials:
            current = current.Chand
            continue
        else:
            changed = changed or (not current.Chand.recent)
            current.setChand(current.Chand.Chand)
            current.Chand.setNhand(current)

            latom += 1
            continue
            ## TODO specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?

    while current.Chand: #renumbering needs to start from the C terminus
        current = current.Chand

    atoms = []
    while current:
        current.recent = max(0, current.recent - 1)
        atoms.append(current)
        current = current.Nhand
    atoms.reverse()

    for i, at in enumerate(atoms):
        at.setId(i)

    return atoms, changed

def prefilter_with_adding(atoms):
    """Removes atoms to reduce the backbone where possible. Adding new joints to long connections
    allows for reducing bumps"""

    first_bead_in_structure = atoms[0]

    next_chain_handle = first_bead_in_structure
    changed = False
    cnt=0
    while next_chain_handle is not None:
        cnt+=1
        currentN = next_chain_handle
        currentC = currentN.Chand
        if currentC is None: #to bypass that currentN can be an end
            break
        while currentC.Chand and not currentC.end:
            currentC = currentC.Chand
        next_chain_handle = currentC.Chand
        assert currentC.end or next_chain_handle is None
        #print(cnt,(currentN.Nhand,(currentN.id,currentN.end),currentN.Chand),(currentC.Nhand,(currentC.id,currentC.end),currentC.Chand))
        verbose = 0

        Ntriangle = [currentN,currentN.Chand,currentN.Chand.Chand if currentN.Chand else None]
        Ctriangle = [currentC,currentC.Nhand,currentC.Nhand.Nhand if currentC.Nhand else None]
        #if 1 or verbose:
            #print("starting with",len(atoms))
        #    print([(a.id,a.end) for a in atoms])
         #   print("handle",next_chain_handle,next_chain_handle.id if next_chain_handle is not None else None)

        while None not in Ntriangle and None not in Ctriangle:# and not any(x.end for x in Ntriangle+Ctriangle): #maybe just one of them?
            if verbose: print(currentN.id,currentN,currentC.id,currentC)
            if verbose:  print(Ntriangle,Ctriangle)
            if Ntriangle[0] == Ctriangle[0]: #just passed the overlap (for odd) or lack of it (even)
                break
            # TODO: multichain
            """
                    if current.Chand.end: #skip "Virtual" connections at chain ends
                current = current.Chand.Chand
                continue
            if current.end:
                current = current.Chand
                continue"""
            overlap = sum(nt in Ctriangle for nt in Ntriangle)
            #if overlap: print("OVERLAP",overlap)

            #TODO dodac overalp 2 gdy a' w a,(b,c'),(c,b')

            if overlap == 3:
                if verbose: print("overlap=3")
                potentials = inspectingOtherLines(Ntriangle[0].vec, Ntriangle[1].vec, Ntriangle[2].vec,
                                                  Ntriangle[0].Nhand,#.vec,
                                                  Ntriangle[2].Chand)#.vec)
                if potentials:
                    if verbose: print("have potetnails")
                    if point_distance(Ntriangle[0].vec, Ntriangle[1].vec) > 3.: #new join on long connection
                        #breaks join Ntraingle[0],Ntriangle[1] (same as Ctirange[1], Ctriangle[2])
                        new_N = get_middlepoint(Ntriangle[0].vec, Ntriangle[1].vec)
                        if verbose: print("adding nowy at", new_N)
                        nowy_id = (Ntriangle[0].original_id + Ntriangle[1].original_id) / 2.
                        if int(nowy_id) != nowy_id:
                            if int(nowy_id) != Ntriangle[0].original_id:
                                nowy_id = int(nowy_id)
                            else:
                                nowy_id = int(nowy_id) + 1
                        nowy = Bead(Vector(new_N), "CA", nowy_id)
                        nowy.recent = 2
                        nowy.setId((Ntriangle[0].id + Ntriangle[1].id) / 2.)
                        nowy.setNhand(Ntriangle[0])
                        nowy.setChand(Ntriangle[1])
                        Ntriangle[1].setNhand(nowy)
                        Ntriangle[0].setChand(nowy)
                    if point_distance(Ntriangle[1].vec, Ntriangle[2].vec) > 3.:  # new join on long connection
                        # breaks join Ntraingle[1],Ntriangle[2] (same as Ctriangle[0], Ctriangle[1])
                        new_N = get_middlepoint(Ntriangle[1].vec, Ntriangle[2].vec)
                        if verbose: print("adding nowy at", new_N)
                        nowy_id = (Ntriangle[1].original_id + Ntriangle[2].original_id) / 2.
                        if int(nowy_id) != nowy_id:
                            if int(nowy_id) != Ntriangle[1].original_id:
                                nowy_id = int(nowy_id)
                            else:
                                nowy_id = int(nowy_id) + 1
                        nowy = Bead(Vector(new_N), "CA", nowy_id)
                        nowy.recent = 2
                        nowy.setId((Ntriangle[1].id + Ntriangle[2].id) / 2.)
                        nowy.setNhand(Ntriangle[1])
                        nowy.setChand(Ntriangle[2])
                        Ntriangle[2].setNhand(nowy)
                        Ntriangle[1].setChand(nowy)
                    currentN = currentN.Chand
                    Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                    currentC = currentC.Nhand
                    Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                    continue
                else:
                    if verbose: print("dont have potetnails")
                    if point_distance(Ntriangle[0].vec, Ntriangle[2].vec) < 4.:  # TODO may be too far?
                        # if can be remove, and the distance is not so great - remove the atom
                        changed = changed or (not Ntriangle[1].recent)
                        Ntriangle[0].setChand(Ntriangle[2])
                        Ntriangle[2].setNhand(Ntriangle[0])
                        ## specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?
                    else:
                        position = get_middlepoint(Ntriangle[0].vec, Ntriangle[2].vec)
                        if verbose: print("adding nowy at", position)
                        # if we can reduce but the distance is too great - remove the atom, put a new one on the connection
                        # reduces traingles to straight line
                        changed = changed or point_distance(Ntriangle[1].vec, position) > EPSILON
                        nowy_id = (Ntriangle[0].original_id + Ntriangle[2].original_id) / 2.

                        nowy = Bead(Vector(position), "CA", nowy_id)
                        nowy.recent = 2
                        nowy.setId(Ntriangle[1].id)  # replaces current atom
                        nowy.setNhand(Ntriangle[0])
                        nowy.setChand(Ntriangle[2])
                        Ntriangle[2].setNhand(nowy)
                        Ntriangle[0].setChand(nowy)
                        currentN = nowy
                        Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                    break #removed middle from overlappig triangle, nothing more to do
            else:
                Npotentials = inspectingOtherLines(Ntriangle[0].vec, Ntriangle[1].vec, Ntriangle[2].vec,
                                                   Ntriangle[0].Nhand,#.vec if Ntriangle[0].Nhand else None,
                                                   Ntriangle[2].Chand,#.vec if Ntriangle[2].Chand else None,
                                                   but_not=Ctriangle)
                Cpotentials = inspectingOtherLines(Ctriangle[0].vec, Ctriangle[1].vec, Ctriangle[2].vec,
                                                   Ctriangle[2].Nhand,#.vec if Ctriangle[2].Nhand else None,
                                                   Ctriangle[0].Chand,#.vec if Ctriangle[0].Chand else None,
                                                   but_not=Ntriangle)
                if verbose: print("overlap < 3")
                if Npotentials:
                    if verbose: print("have Nportenials")
                    if Cpotentials:
                        if verbose: print("have Cportenials")
                        #nothing can be done
                        currentC=currentC.Nhand
                        Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                        currentN = currentN.Chand
                        Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                        continue # it should break at the beginning of next loop pass if overlap was 2
                    else:
                        C123normal = get_normal(Ctriangle[0].vec, Ctriangle[1].vec, Ctriangle[2].vec)
                        # check if they ARE crossing each othe
                        c0v, c1v, c2v = [c.vec for c in Ctriangle]
                        n0v, n1v, n2v = [n.vec for n in Ntriangle]
                        if any(map(lambda x: x is not False, [found_crossing(c0v, c1v, c2v, n0v, n1v, C123normal),
                                                              found_crossing(c0v, c1v, c2v, n1v, n2v, C123normal)])):
                            if verbose: print("N jest w C")
                            currentC = currentC.Nhand
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            currentN = currentN.Chand
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            continue
                        if verbose: print("dont have other Cportenials")
                        if point_distance(Ctriangle[0].vec, Ctriangle[2].vec) < 4.:  # TODO may be too far?
                            # if can be remove, and the distance is not so great - remove the atom
                            changed = changed or (not Ctriangle[1].recent)
                            Ctriangle[0].setNhand(Ctriangle[2])
                            Ctriangle[2].setChand(Ctriangle[0])
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            currentN = currentN.Chand
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            continue
                             ## specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?
                        else:
                            position = get_middlepoint(Ctriangle[0].vec, Ctriangle[2].vec)
                            if verbose: print("adding nowy at", position)
                            # if we can reduce but the distance is too great - remove the atom, put a new one on the connection
                            # reduces traingles to straight line
                            changed = changed or point_distance(Ctriangle[1].vec, position) > EPSILON
                            nowy_id = (Ctriangle[0].original_id + Ctriangle[2].original_id) / 2.
                            nowy = Bead(Vector(position), "CA", nowy_id)
                            nowy.recent = 2
                            nowy.setId(Ctriangle[1].id)  # replaces current atom
                            nowy.setChand(Ctriangle[0])
                            nowy.setNhand(Ctriangle[2])
                            Ctriangle[2].setChand(nowy)
                            Ctriangle[0].setNhand(nowy)
                            currentC = nowy
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            currentN = currentN.Chand
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            continue
                else:
                    if verbose: print("dont have other Nportenials")
                    if Cpotentials:
                        if verbose: print("have Cportenials")
                        N123normal = get_normal(Ntriangle[0].vec, Ntriangle[1].vec, Ntriangle[2].vec)
                        # check if they ARE crossing each othe
                        c0v, c1v, c2v = [c.vec for c in Ctriangle]
                        n0v, n1v, n2v = [n.vec for n in Ntriangle]
                        if any(map(lambda x: x is not False, [found_crossing(n0v, n1v, n2v, c0v, c1v, N123normal),
                                                              found_crossing(n0v, n1v, n2v, c1v, c2v, N123normal)])):
                            if verbose: print("C jest w N")
                            currentN = currentN.Chand
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            currentC = currentC.Nhand
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            continue
                        if point_distance(Ntriangle[0].vec, Ntriangle[2].vec) < 4.:  # TODO may be too far?
                            # if can be remove, and the distance is not so great - remove the atom
                            changed = changed or (not Ntriangle[1].recent)
                            Ntriangle[0].setChand(Ntriangle[2])
                            Ntriangle[2].setNhand(Ntriangle[0])
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            currentC = currentC.Nhand
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            continue
                             ## specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?
                        else:
                            position = get_middlepoint(Ntriangle[0].vec, Ntriangle[2].vec)
                            if verbose: print("adding nowy at", position)
                            # if we can reduce but the distance is too great - remove the atom, put a new one on the connection
                            # reduces traingles to straight line
                            changed = changed or point_distance(Ntriangle[1].vec, position) > EPSILON
                            nowy_id = (Ntriangle[0].original_id + Ntriangle[2].original_id) / 2.

                            nowy = Bead(Vector(position), "CA", nowy_id)
                            nowy.recent = 2
                            nowy.setId(Ntriangle[1].id)  # replaces current atom
                            nowy.setNhand(Ntriangle[0])
                            nowy.setChand(Ntriangle[2])
                            Ntriangle[2].setNhand(nowy)
                            Ntriangle[0].setChand(nowy)
                            currentN = nowy
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            currentC = currentC.Nhand
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            continue
                    else:
                        if verbose: print("dont have any portenials")
                        if overlap == 2:
                            if verbose: print("overlap==2")
                            break

                            new_BC = get_middlepoint(Ntriangle[1].vec, Ctriangle[1].vec)
                            if verbose: print("adding nowy at", new_BC)
                            changed = True
                            nowy_id = (Ntriangle[1].original_id +Ctriangle[1].original_id)/2.
                            nowy = Bead(Vector(new_BC), "CA", nowy_id)
                            nowy.recent = 2
                            nowy.setId(Ntriangle[1].id)  # replaces current atom

                            Ntriangle[0].setChand(nowy)
                            Ctriangle[0].setNhand(nowy)
                            nowy.setNhand(Ntriangle[0])
                            nowy.setChand(Ctriangle[0])
                            Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                            Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                            continue #tera bedzie overlap N.A,nowy,C.A
                        else: #separate triangles
                            if verbose: print("overlap<2")
                            N123normal = get_normal(Ntriangle[0].vec, Ntriangle[1].vec, Ntriangle[2].vec)
                            C123normal = get_normal(Ctriangle[0].vec, Ctriangle[1].vec, Ctriangle[2].vec)
                            #check if they ARE crossing each othe
                            c0v,c1v,c2v = [c.vec for c in Ctriangle]
                            n0v,n1v,n2v = [n.vec for n in Ntriangle]
                            if verbose: print("checking",[found_crossing(n0v,n1v,n2v,c0v,c1v,N123normal),
                                    found_crossing(n0v,n1v,n2v,c1v,c2v,N123normal),
                                    found_crossing(c0v, c1v, c2v, n0v, n1v, C123normal),
                                    found_crossing(c0v, c1v, c2v, n1v, n2v, C123normal)])
                            if any(map(lambda x: x is not False,[found_crossing(n0v,n1v,n2v,c0v,c1v,N123normal),
                                    found_crossing(n0v,n1v,n2v,c1v,c2v,N123normal),
                                    found_crossing(c0v, c1v, c2v, n0v, n1v, C123normal),
                                    found_crossing(c0v, c1v, c2v, n1v, n2v, C123normal)])):
                                if verbose: print("traingles cross each other")
                                currentN = currentN.Chand
                                currentC = currentC.Nhand
                                Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                                Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                                #to be determined if move them
                                continue

                            # check if they would cross each other
                            xt1 = found_crossing(Ntriangle[0].vec, Ntriangle[1].vec, Ntriangle[2].vec,Ctriangle[0].vec,Ctriangle[2].vec,N123normal)
                            xt2 = found_crossing(Ctriangle[0].vec, Ctriangle[1].vec, Ctriangle[2].vec,Ntriangle[0].vec,Ntriangle[2].vec,C123normal)
                            #print("xts:",xt1,xt2)
                            #changed to or below, bo mzoe czegos nei widze
                            if xt1 is not False or xt2 is not False: #should always be bool(xt1)==bool(xt2)
                                currentN = currentN.Chand
                                currentC = currentC.Nhand
                                Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                                Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                                # to be on the safe side
                                continue

                                X = get_middlepoint(xt1,xt2) #moving lines should meet here, last safe point
                                changed = True
                                new_N_B = get_middlepoint(xt1,X)
                                if verbose: print("adding nowy at", new_N_B)
                                nowy_N = Bead(Vector(new_N_B), "CA", Ntriangle[1].original_id)
                                nowy_N.recent = 2
                                nowy_N.setId(Ntriangle[1].id)
                                Ntriangle[0].setChand(nowy_N)
                                Ntriangle[2].setNhand(nowy_N)
                                nowy_N.setNhand(Ntriangle[0])
                                nowy_N.setChand(Ntriangle[2])
                                new_C_B = get_middlepoint(xt2,X)
                                if verbose: print("adding nowy at", new_C_B)
                                nowy_C = Bead(Vector(new_C_B), "CA", Ctriangle[1].original_id)
                                nowy_C.recent = 2
                                nowy_C.setId(Ctriangle[1].id)
                                Ctriangle[0].setNhand(nowy_C)
                                Ctriangle[2].setChand(nowy_C)
                                nowy_C.setNhand(Ctriangle[2])
                                nowy_C.setChand(Ctriangle[0])
                                #move on, wont move further here
                                currentN = currentN.Chand
                                Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                                currentC = currentC.Nhand
                                Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]

                            else: #safe to do both triangles
                                if verbose: print("no xts,distances are",point_distance(Ntriangle[0].vec, Ntriangle[2].vec),point_distance(Ctriangle[0].vec, Ctriangle[2].vec))
                                if point_distance(Ntriangle[0].vec, Ntriangle[2].vec) < 4.:  # TODO may be too far?
                                    # if can be remove, and the distance is not so great - remove the atom
                                    changed = changed or (not Ntriangle[1].recent)
                                    Ntriangle[0].setChand(Ntriangle[2])
                                    Ntriangle[2].setNhand(Ntriangle[0])
                                    Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                                    ## specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?
                                else:
                                    position = get_middlepoint(Ntriangle[0].vec, Ntriangle[2].vec)
                                    if verbose: print("adding nowy at",position)
                                    # if we can reduce but the distance is too great - remove the atom, put a new one on the connection
                                    # reduces traingles to straight line
                                    changed = changed or point_distance(Ntriangle[1].vec, position) > EPSILON
                                    nowy_id = (Ntriangle[0].original_id + Ntriangle[2].original_id) / 2.

                                    nowy = Bead(Vector(position), "CA", nowy_id)
                                    nowy.recent = 2
                                    nowy.setId(Ntriangle[1].id)  # replaces current atom
                                    nowy.setNhand(Ntriangle[0])
                                    nowy.setChand(Ntriangle[2])
                                    Ntriangle[2].setNhand(nowy)
                                    Ntriangle[0].setChand(nowy)
                                    currentN = nowy
                                    Ntriangle = [currentN, currentN.Chand, currentN.Chand.Chand if currentN.Chand else None]
                                if point_distance(Ctriangle[0].vec, Ctriangle[2].vec) < 4.:  # TODO may be too far?
                                    # if can be remove, and the distance is not so great - remove the atom
                                    changed = changed or (not Ctriangle[1].recent)
                                    Ctriangle[0].setNhand(Ctriangle[2])
                                    Ctriangle[2].setChand(Ctriangle[0])
                                    Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                                    continue
                                    ## specjalnie nie zamieniam current - watpliwe, ale a noz widelec jeszcze jeden?
                                else:
                                    position = get_middlepoint(Ctriangle[0].vec, Ctriangle[2].vec)
                                    if verbose: print("adding nowy at", position)
                                    # if we can reduce but the distance is too great - remove the atom, put a new one on the connection
                                    # reduces traingles to straight line
                                    changed = changed or point_distance(Ctriangle[1].vec, position) > EPSILON
                                    nowy_id = (Ctriangle[0].original_id + Ctriangle[2].original_id) / 2.
                                    nowy = Bead(Vector(position), "CA", nowy_id)
                                    nowy.recent = 2
                                    nowy.setId(Ctriangle[1].id)  # replaces current atom
                                    nowy.setChand(Ctriangle[0])
                                    nowy.setNhand(Ctriangle[2])
                                    Ctriangle[2].setChand(nowy)
                                    Ctriangle[0].setNhand(nowy)
                                    currentC = nowy
                                    Ctriangle = [currentC, currentC.Nhand, currentC.Nhand.Nhand if currentC.Nhand else None]
                                    continue
        if verbose: print("changed is",changed)


    #print(list((x.id,x.original_id) for x in atoms))

    #while currentN.Chand:
    #    currentN = currentN.Chand

    currentN = first_bead_in_structure
    atoms = []
    while currentN:
        currentN.recent = max(0, currentN.recent - 1)
        atoms.append(currentN)
        currentN = currentN.Chand
    #atoms.reverse()

    #print("after filtering atoms count:",len(atoms),atoms[0],atoms[0].Chand)

    for i, at in enumerate(atoms):
        at.setId(i)
    #exit()
    #exit()

    return atoms, changed

def run_through_filtering(atoms, config, greedy=0):
    """Governs reduction process and when it should end (when there is no improvement)"""

    frames = []
    frames.append([(x.original_id, x.vec, x.end) for x in atoms])
    o_and_t = config.get('outfile',0) and config.get('trajectory',0)

    if o_and_t:
        #print_out_last_pdb(config, frames, atoms)
        print_out_last_frame(config, atoms, len(frames))

    if not config or not config['no_pulling']:
        atoms, changed = prefilter_with_adding(atoms)

        frames_without_shortening = 0
        atoms_num = len(set(_.tuple() for _ in atoms))

        while changed:  # TODO should check RMSD up to 2 steps back
            frames.append([(x.original_id, list(x.vec), x.end) for x in atoms])
            if len(frames) >= 3 and len(frames[-1]) == len(frames[-3]):
                if all([all(i == j for i, j in zip(x[1], y[1])) for x, y in zip(frames[-1], frames[-3])]):
                    # if two frames back we had the exact same positions
                    break
                if len(set(_.tuple() for _ in atoms)) == atoms_num:
                    frames_without_shortening += 1
                else:
                    atoms_num = len(set(_.tuple() for _ in atoms))
                    frames_without_shortening = 0
                if frames_without_shortening > 2*len(atoms):
                    if VERBOSE: print ("Breaking filtering due to to not enough change")
                    break

            if o_and_t:
                #print_out_last_pdb(config, frames, atoms)
                print_out_last_frame(config, atoms, len(frames))

            atoms, changed = prefilter_with_adding(atoms)

        #    frames.append([(x.original_id, x.vec) for x in atoms])

    if VERBOSE:
        print ("Finished filtering, have {} atoms".format(len(atoms)))

    if greedy: #reduce number of atoms (improves topology detection)
        atoms, _ = filter_greedy(atoms, config, len(frames))
        frames.append([(x.original_id, x.vec, x.end) for x in atoms])
        if VERBOSE:
            print ("Finished greedy filtering, have {} atoms".format(len(atoms)))

    if o_and_t:
        #print_out_last_pdb(config, frames, atoms)
        print_out_last_frame(config, atoms, len(frames))

    return frames,atoms

def run_through_division(atoms):
    """Finds connections which can be safely cut (connected segments can't be simplified even separately) """
    if VERBOSE: print ("Starting with {} atoms".format(len(atoms)))
    cuts = []
    for i in range(1, len(atoms) - 2):  # No sense in checking last edges
        j = i + 1
        if VERBOSE: print("Cutting edge {}".format(j))
        _middle = get_middlepoint(atoms[i].vec, atoms[j].vec)
        _tmp1 = chainDeepCopy(atoms[:j] + [Bead(_middle, "CA")])
        _, _tmp1a = run_through_filtering(_tmp1, {}, greedy=1)
        if VERBOSE: print("Left from {} to {}".format(len(_tmp1),len(_tmp1a)))
        _tmp2 = chainDeepCopy([Bead(_middle, "CA")] + atoms[j:])  #
        _, _tmp2a = run_through_filtering(_tmp2, {}, greedy=1)
        if VERBOSE: print("Right from {} to {}".format(len(_tmp2),len(_tmp2a)))
        if _tmp1 == _tmp1a and _tmp2 == _tmp2a:  # no reduction
            # print i,
            cuts.append(i)

    cuts = [0] + cuts + [len(atoms) - 1]  # fix to get proper segment indices
    atom_lists = []
    if VERBOSE: print("Cuts: {}".format(cuts))

    for i in range(0, len(cuts) - 1):
        cutl = cuts[i]
        cutr = cuts[i + 1]
        lista = []
        if cutl:  # add middlepoints of our breakable edges
            _mid_l = get_middlepoint(atoms[cutl].vec, atoms[cutl + 1].vec)
            lista.append(Bead(_mid_l, "CA"))
        lista += atoms[cutl + bool(cutl):cutr + 1]
        if cutr < len(atoms) - 2:
            _mid_r = get_middlepoint(atoms[cutr].vec, atoms[cutr + 1].vec)
            lista.append(Bead(_mid_r, "CA"))
        lista = chainDeepCopy(lista)
        atom_lists.append(lista)

    return atom_lists  # returns all separateable segments


def pull(atoms, config, chain_names, get_representative=0,greedy=0):#, greedy_file="",trajectory=True,quiet=False, chain_names=(),rna=False):
    outfile = config['outfile']

    frames, atoms = run_through_filtering(atoms, config, greedy=greedy)
    quiet = config['quiet']
    repr = []
    if get_representative:
        Atom = namedtuple("Atom","original_id x y z end")
        repr = find_frame(frames)
        repr = [Atom(x[0],x[1][0],x[1][1],x[1][2],x[2]) for x in repr]

    chains = []
    model_num = len(frames)+1
    separate_chains = divide_into_bead_chains(atoms)
    for c, chain in enumerate(separate_chains):
        ch = Chain(c,chain)
        if chain_names: ch.setChainName(chain_names[c])
        if outfile:
            if not quiet: model_num = ch.print2file(config,model_num)
        chains.append(ch)

    for _,c1 in enumerate(chains):
        for c2 in chains[_+1:]:
            if are_chains_linked(c1,c2):
                c1.neighbours.append(c2)
                c2.neighbours.append(c1)

    return chains,repr

def same_vecs(l1,l2):
    return all(np.array_equal(x.vec,y.vec) for x,y in zip(l1,l2))

def are_chains_linked(ch1,ch2):
    _tmp_both = chainDeepCopy(ch1.atoms+ch2.atoms)
    _tmp_both[len(ch1.atoms)-1].end = True
    _tmp_ch1 = chainDeepCopy(ch1.atoms)
    _tmp_ch2 = chainDeepCopy(ch2.atoms)

    _,_both = run_through_filtering(_tmp_both, {}, greedy=1)
    _,_ch1 = run_through_filtering(_tmp_ch1, {}, greedy=1)
    _,_ch2 = run_through_filtering(_tmp_ch2, {}, greedy=1)
    return not same_vecs(_both,_ch1+_ch2)


class Chain:
    def __init__(self,chain,atoms):
        self.chain = chain
        self.atoms = atoms
        self.atom_lists = run_through_division(atoms)
        self.knots = []
        self.neighbours = []
        self.closing_path = closeTheCurve(atoms)
        self.closing_paths = [closeTheCurve(a) for a in self.atom_lists]
        self.wanda = None
        self.dowker = None
        self.dowker_code = None
        self.chain_name = ""


    def setChainName(self,c):
        self.chain_name = "({})".format(c)

    def print2file(self,config,model_num):
        #print_out_one_frame(outfile, self.atoms, model_num,self.chain)
        model_num+=1
        for alist in self.atom_lists:
            #print_out_one_frame(outfile, alist, model_num,self.chain,self.rna)
            print_out_last_frame(config, alist, model_num, cur_chain=self.chain_name[1:-1], final=True)
            model_num+=1
        return model_num



if __name__ == "__main__":
    import sys
    with open(sys.argv[1]) as input:
        atoms = input.readlines()
    if len(sys.argv) > 2:
        atoms = atoms[int(sys.argv[2]):int(sys.argv[3]) + 1]
    print ("Result is", run_through_filtering(atoms, sys.argv[1] + (
        "{}_{}".format(sys.argv[2], sys.argv[3]) if len(sys.argv) > 2 else "") + ".outall.pdb"))


