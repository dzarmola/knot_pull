from __future__ import print_function
import math

from .vector_ops import *
from .config import VERBOSE
from .kpclasses import Code, Crossing, CodeHistory, DowkerError



#### Helper functions #########

def find_z_old(l0, l1, p):
    """Finds z axis position"""
    print(l0,l1,p)
    diffx = l1[0] - l0[0]
    diffpx = p[0] - l0[0]
    diffz = l1[2] - l0[2]
    print(diffx,diffpx,diffz)
    diffpz = diffpx * diffz
    if diffpz:
        diffpz /= diffx
    print(diffpz)
    return l0[2] + diffpz

def find_z(l0, l1, p):
    """Finds z axis position"""
    x0,y0,z0 = l0
    x1,y1,z1 = l1
    full_pd = point_distance_2d(l0,l1)
    cross_pd = point_distance_2d(l0,p)
    dz = z1-z0
    return z0 + dz*(cross_pd/full_pd)



def fix_numbering(code,start):
    """Start is the lowest number from which we should remove 2"""
    for i,l in enumerate(code):
        if abs(l[0])>=start:
            l[0] = l[0]/abs(l[0]) * (abs(l[0])-2)
        if abs(l[1])>=start:
            l[1] = l[1]/abs(l[1]) * (abs(l[1])-2)
        code[i] = l
    return code


def fix_minus(loop):
    """Removes "-" from uneven number"""
    if uneven(loop)<0:
        loop = map(abs,loop)
    return loop


def shortest_loop(code):
    """Sorts crossing by shortest loop described"""
    return sorted(code,key=lambda x:abs(abs(x[0])-abs(x[1])))


def same_side_cross(v1,v2):
    """Checks if one line crosses both crossings from the same side"""
    if abs(v1[0]-v2[0]) == 1:
        return v1[0].top == v2[0].top
    elif abs(v1[1]-v2[0]) == 1:
        return v1[1].top == v2[0].top
    elif abs(v1[1] - v2[1]) == 1:
        return v1[1].top == v2[1].top
    elif abs(v1[0] - v2[1]):
        return v1[0].top == v2[1].top
    else:
        return False


def flatten(x):
    return list(x[0]) + list(x[1])


def sign(x):
    return -1 if x < 0 else 1


def add(x,y):
    return x - y if sign(x) == 1 else -(abs(x) - y)


def index(lista, co):
    """.index with on absolute values"""
    for i in range(len(lista)):
        if abs(lista[i]) == co:
            return i
    else:
        raise IndexError


def uneven(x):
    return x[0] if x[0]%2 else x[1]


def even(x):
    return x[1] if x[0]%2 else x[0]


def shorten_code(lista):
    lista = sorted(lista,key=lambda x:uneven(x))
    return tuple([even(x) for x in lista])


### HELPERS END HERE ####


def get_dt_code(atoms):
    """Calculates DT code based on protein chain"""
    cnt = 1
    crosses = {}
    code = Code()
    for l in range(len(atoms) - 1):
        this_cross = []
        for k in range(len(atoms) - 1):
            if l == k or abs(l - k) == 1: continue
            line1 = [atoms[l], atoms[l + 1]]
            line2 = [atoms[k], atoms[k + 1]]
            if intersect_2d(*(line1 + line2)) and not len(set(tuple(x) for x in line1 + line2)) < 4:
            #if intersect_2d(*(line1 + line2)) and not len(set(tuple(x[:2]) for x in line1 + line2)) < 4:
                cross = get_crossing_2d(atoms[l], atoms[l + 1], atoms[k], atoms[k + 1])
                this_cross.append((point_distance_2d(atoms[l], cross), k, cross))
        for _, k, cross in sorted(this_cross):
            val = 1
            skein = ''
            if find_z(atoms[l], atoms[l + 1], cross) > find_z(atoms[k], atoms[k + 1], cross):
                if cnt % 2 == 0:
                    val = -1
                skein = get_skein(atoms[l], atoms[l+1], atoms[k], atoms[k+1])
            crosses[tuple(sorted([l, k]))] = crosses.get(tuple(sorted([l, k])), []) + [cnt * val,skein]

            cnt += 1
    if not crosses:
        return None
    crossings = list(zip(*sorted(crosses.items(), key=lambda x:x[1])))[1]
    if VERBOSE: print("Reading crosses:",crossings)
    code.read_in(crossings)
    #code.read_in(crosses.values())
    return code


def renumber_dt_list(lista):
    nums = sorted(lista)
    code = [(nums.index(x)+1)*2 for x in lista]
    return code


def renumber_dt_code(codeobj):
    mod_code = codeobj.mod_dowker_code()
    code = codeobj.dowker_code()

    code = renumber_dt_list(code)
    mod_code = [x*(mod_code[i]/abs(mod_code[i])) for i,x in enumerate(code)]
    codeobj.set_mod_dowker_code(mod_code)


def issublistof(overlist,sublist):
    l = len(sublist)
    for i in range(len(overlist)-l):
        if sublist == overlist[i:i+l]:
            return True
    return False


def is_permuted(lista,code_len):
    sl = sorted(lista)
    if issublistof(lista+lista,sl):
        return all(sl[i+1]-sl[i] in [2,(code_len-(len(sl)-1))*2] for i in range(len(sl)-1))
    return False


def splitz(seq, smallest):
    group = []
    for num in seq:
        if num != smallest:
            group.append(num)
        elif group:
            yield group
            group = []
    if group:
        yield group


def find_permutations_in(code):
    # TODO instead of "None not in" => filter(x is not None)
    ###just like blast - find a seed of len 3 and try to expand
    dbcode = code + code
    subcodes = []
    lc = len(code)
    i = 0
    while i<lc:
        j=i+3
        if None not in dbcode[i:j]:
            if is_permuted(dbcode[i:j],max(filter(lambda x:x is not None,dbcode))/2):
                while j-(i-1)<=len(code) and None not in dbcode[i-1:j] and is_permuted(dbcode[i - 1:j],max(dbcode)/2):
                    i -= 1
                while j+1-i<=len(code) and None not in dbcode[i:j+1] and is_permuted(dbcode[i:j + 1],max(dbcode)/2):
                    j += 1
                if len(list(filter(lambda x: x is not None, dbcode))) - (j-i)*2 in [1,2]:
                    i=j-1
                    continue
                subcodes.append(dbcode[i:j])
                for _ in range(i,j):
                    dbcode[_] = None
                    dbcode[(_+lc)%(lc*2)] = None
                i=0
                continue
        i+=1
    dbcode = dbcode[:lc]
    if None in dbcode and dbcode[-1] != None:
        while dbcode[0]!=None:
            dbcode.append(dbcode.pop(0))
    rest = list(splitz(dbcode, None))
    subcodes += rest
    subcodes = list(map(renumber_dt_list,subcodes))
    return subcodes


def translate_dt_list(code):
    DT_CODES = {(4, 6, 2): "3_1",
                (4, 6, 8, 2): "4_1",
                (6, 8, 10, 2, 4): "5_1",
                (4, 8, 10, 2, 6): "5_2",
                (4, 8, 12, 10, 2, 6): "6_1",
                (4, 8, 10, 12, 2, 6): "6_2",
                (4, 8, 10, 2, 12, 6): "6_3",
                (8,10,12,14,2,4,6) : "7_1",
                (4,10,14,12,2,8,6) : "7_2",
                (6, 10, 12, 14, 2, 4, 8) : "7_3",
                (6, 10, 12, 14, 4, 2, 8) : "7_4",
                (4, 10, 12, 14, 2, 8, 6) : "7_5",
                (4, 8, 12, 2, 14, 6, 10) : "7_6",
                (4, 8, 10, 12, 2, 14, 6): "7_7",
                }
    out = []
    for c in code:
        for dc in DT_CODES:
            if issublistof(list(dc)+list(dc),c) and issublistof(c+c,list(dc)):
                out.append(DT_CODES[dc])
                break
        else:
            if len(c) == 5:
                out.append("5_2")
            elif len(c) == 3:
                out.append("3_1")
            elif len(c) == 4:
                out.append("4_1")
            else:
                out.append("{}_0".format(len(c)))
    return "#".join(out)


def translate_dt_code(code):
    DT_CODES = {(4, 6, 2): "3_1",
                (4, 6, 8, 2): "4_1",
                (6, 8, 10, 2, 4): "5_1",
                (4, 8, 10, 2, 6): "5_2",
                (4, 8, 12, 10, 2, 6): "6_1",
                (4, 8, 10, 12, 2, 6): "6_2",
                (4, 8, 10, 2, 12, 6): "6_3",
                (8, 10, 12, 14, 2, 4, 6): "7_1",
                (4, 10, 14, 12, 2, 8, 6): "7_2",
                (6, 10, 12, 14, 2, 4, 8): "7_3",
                (6, 10, 12, 14, 4, 2, 8): "7_4",
                (4, 10, 12, 14, 2, 8, 6): "7_5",
                (4, 8, 12, 2, 14, 6, 10): "7_6",
                (4, 8, 10, 12, 2, 14, 6): "7_7",
                }
    code_str = code.dowker_str()
    code_str_r = code.dowker_rstr()
    for k, v in DT_CODES.items():
        k_new = "-".join(map(str, k + k))
        if code_str in k_new or code_str_r in k_new:
            return v
    return code


def unfurl(code):
    """Takes dowker code sorted by length
    unfurls crossed loops which are passed through another line"""
    changed = False
    for _,loop in enumerate(code.gen_short_loop()):
        i,a = loop.min().val, loop.max().val
        if loop.len_loop() == 3: #@TODO doesnt work on wraparound - doesnt have to, added renumbering
            to_del = [code.index(loop)]
            rg_i_a = range(i+1,a)
            fr,fz = None,None #first in range, first_zip_line
            sr,sz = None,None #second in range, second_zip_line
            for _,l in enumerate(code):
                if l[0].val == rg_i_a[0]:
                    fr,fz = l
                    l1 = l
                    to_del.append(_)
                elif l[0].val == rg_i_a[1]:
                    sr,sz = l
                    l2 = l
                    to_del.append(_)
                if l[1].val == rg_i_a[0]:
                    fz,fr = l
                    l1 = l
                    to_del.append(_)
                elif l[1].val == rg_i_a[1]:
                    sz,sr = l
                    l2 = l
                    to_del.append(_)
            if fr is not None and sr is not None and abs(fr.val-sr.val) == abs(fz.val-sz.val) == 1:
                if same_side_cross(Crossing(fr.val, fr.top, fz.val, fz.top), loop) or same_side_cross(
                        Crossing(sr.val, sr.top, sz.val, sz.top), loop):
                    if fz - fr > 0: #najpierw petla
                        if fz.val-sz.val == -1: # kierunki zgodne
                            _1 = Crossing(fr.val-1,fr.top,fz.val-1,fz.top)
                            _2 = Crossing(sr.val - 1, sr.top, sz.val - 3, sz.top)
                        elif fz.val-sz.val == 1: #kierunki przeciwne
                            _1 = Crossing(fr.val - 1, fr.top, fz.val - 3, fz.top)
                            _2 = Crossing(sr.val - 1, sr.top, sz.val - 1, sz.top)
                    else: #petla pozniej
                        if fz.val - sz.val == -1:  # kierunki zgodne
                            _1 = Crossing(fr.val - 1, fr.top, sz.val, fz.top) #zipy zamiana!
                            _2 = Crossing(sr.val - 1, sr.top, fz.val, sz.top)
                        elif fz.val - sz.val == 1:  # kierunki przeciwne
                            _1 = Crossing(fr.val - 1, fr.top, sz.val, fz.top)  # zipy zamiana!
                            _2 = Crossing(sr.val - 1, sr.top, fz.val, sz.top)
                    _1.skein = l1.skein
                    _2.skein = l2.skein
                    code.remove(loop)
                    code.pop(l1)
                    code.pop(l2)
                    code.add(_1)
                    code.add(_2)
                    changed = True
                    break
    return changed


def neighbours_val(a,b,len_c):
    if b-a==1: return True
    if b == len_c*2 and a==1:
        return True
    return False


def both_neighbours(v1,v2,len_c):
    return (abs(v1.even()-v2.uneven())==1 or v1.even()+v2.uneven()==len_c*2) and (abs(v2.even()-v1.uneven())==1 or v2.even()+v1.uneven()==len_c*2)


def dowker_dbl_loops(code):
    """Finds a loop crossing either over or under a line, brings it back
    (both crossing on the same side of the line)"""
    trash = set([])
    changed = 0

    found = 1
    while found:
        found = 0
        for k1, v1 in enumerate(code):
            if found: break
            for k2, v2 in enumerate(code[k1 + 1:]):
                if same_side_cross(v1, v2):
                    if both_neighbours(v1, v2, len(code)):
                        if VERBOSE: print("Removing {} {} from {}: ".format(v1, v2, code))
                        code.remove(v1,v2)
                        if VERBOSE: print("{}".format(code))
                        changed = True
                        found = True
                        break

    found = 1
    while found:
      found = 0
      for k1, v1 in enumerate(code):
        if found:
            break
        for k2, v2 in enumerate(code[k1+1:]):
            #if v1 in trash or v2 in trash:
            #    continue
            if same_side_cross(v1, v2) and min(map(abs,[v1[0]-v2[0],v1[0]-v2[1],v1[1]-v2[0],v1[1]-v2[1]]))==1: #are neighbours
                v1_neigh = code.different_neighbour(v1,v2)
                v2_neigh = code.different_neighbour(v2, v1)
                v1_loop_val = v1.min() if abs(v1.min()-v2.min())==1 or abs(v1.min()-v2.max())==1 else v1.max()
                v2_loop_val = v2.min() if abs(v1.min()-v2.min())==1 or abs(v1.max()-v2.min())==1 else v2.max()
                if v2_loop_val < v1_loop_val:
                    v1,v2 = v2,v1
                    v1_loop_val, v2_loop_val = v2_loop_val, v1_loop_val
                if not same_side_cross(v2,v2_neigh): #tylko nastepnik? # not (same_side_cross(v1, v1_neigh) or same_side_cross(v2, v2_neigh)):
                    sstop = v1_loop_val.top
                    ssn = code.different_neighbour(v2,v1)
                    if ssn.min() < v1_loop_val < v2_loop_val < ssn.max():
                        sswhich = 0 if abs(ssn[0]-v2_loop_val) == 1 else 1
                        try:
                            progeny = code.pullulate()
                            pssn = progeny.different_neighbour(v2, v1)
                            progeny.find_and_remove(v1,v2)
                            #progeny.find_and_remove(v2)
                            pnewss = pssn[sswhich].val
                            pnewss2 = progeny.find_error(pnewss)
                            progeny.make_way(pnewss)
                            progeny.make_way(pnewss2)
                            if pnewss2 <= pnewss: pnewss += 1
                            new = Crossing(pnewss, sstop, pnewss2, not sstop)
                            new.skein = v2.skein
                            progeny.add(new)

                            code.remove(v1)
                            code.remove(v2)
                            newss = ssn[sswhich].val
                            newss2 = code.find_error(newss)
                            code.make_way(newss)
                            code.make_way(newss2)
                            if newss2 <= newss: newss += 1
                            new = Crossing(newss,sstop,newss2,not sstop)
                            new.skein = v2.skein
                            code.add(new)
                            found = 1
                            changed = 1
                            break
                        except DowkerError as e:
                            if VERBOSE: print("407",e)
                            continue
    if VERBOSE: print("after progenys", code)
    '''found = 1
    while found:
        found =0
        for k1, v1 in enumerate(code):
            if found: break
            for k2, v2 in enumerate(code[k1+1:]):
                #if v1 in trash or v2 in trash:
                #    continue
                if same_side_cross(v1, v2):
                    if not both_neighbours(v1,v2,len(code)):
                        v1_neigh = code.different_neighbour(v1,v2)
                        v2_neigh = code.different_neighbour(v2, v1)
                        if v2.min() - v1.min() == v2.max()-v1.max() and (v2.min()- v1.max() ==1 or v1.min()-v2.max() == 1) \
                            and not (v1.min() < v2.min() < v1.max() or v2.min() < v1.min() < v2.max() ): #TODO add wraparound
                            if same_side_cross(v1,v1_neigh) != same_side_cross(v2,v2_neigh): #if just one - can do this
                                pass
                                """trash.add(v1)
                                trash.add(v2)
                                b=0
                                print "adding", v1_neigh, v1, v2, v2_neigh
                                print v1.min(),v1.max(),v2.min(),v2.max()
                                break"""
                                #wywalam bo psulo (-4,11),(3,-12),(-2,7),(1,-8),(6,.9),(-10,13),(.5,14)
                            elif same_side_cross(v1,v1_neigh) and same_side_cross(v2,v2_neigh):
                                trash.add(v1)
                                trash.add(v2)
                                trash.add(v1_neigh)
                                trash.add(v2_neigh)
                                px=v1
                                py=v2
                                x=v1_neigh
                                y=v2_neigh
                                nx = code.different_neighbour(x,px)
                                ny = code.different_neighbour(y,py)
                                while same_side_cross(x,nx) and same_side_cross(y,ny): ### @TODO check how long
                                    trash.add(nx)
                                    trash.add(ny)
                                    px=x
                                    x=nx
                                    py=y
                                    y=ny
                                    nx = code.different_neighbour(x, px)
                                    ny = code.different_neighbour(y, py)

                                code.remove(*trash)
                                found = 1
                                changed = 1
                                break
                            else: #both neighbours differ
                                a,b = v1.even(),v2.uneven()
                                c,d = v1.uneven(),v2.even()
                                if abs(a-b)!=1:
                                    a,b = v1.uneven(),v2.even()
                                    c,d = v1.even(),v2.uneven()
                                #now we find the previous corssing
                                for cr in code:
                                    if c-1 in cr and d-1 in cr:
                                        break
                                else:
                                    break #there is no previous crossing
                                _1 = Crossing(a.val,a.top,cr.even().val if a%2 else cr.uneven().val,not a.top)
                                _2 = Crossing(b.val,b.top,cr.uneven().val if a%2 else cr.even().val,not b.top)
                                e,f = cr
                                if e>f:
                                    if c>d:
                                        _3 = Crossing(c.val,e.top,d.val,f.top)
                                    else:
                                        _3 = Crossing(c.val,f.top,d.val,e.top)
                                else:
                                    if c>d:
                                        _3 = Crossing(c.val,f.top,d.val,e.top)
                                    else:
                                        _3 = Crossing(c.val,e.top,d.val,f.top)
                                code.pop(v1)
                                code.pop(v2)
                                code.pop(cr)
                                code.add(_1)
                                code.add(_2)
                                code.add(_3)
                                changed = 1
                                found = 1
                                break
            #else:
            #    b=1
    if VERBOSE: print("after multiunfold", code)'''
    #trash = sorted(trash)

    return changed

def third_redei(code):
    changed = False
    for k1, v1 in enumerate(code):
        for k2, v2 in enumerate(code[k1+1:]):
            if same_side_cross(v1, v2):
                if not both_neighbours(v1,v2,len(code)) and (abs(v1.uneven()-v2.even())==1 or abs(v1.even()-v2.uneven())==1):
                    #v1_neigh = code.different_neighbour(v1,v2)
                    #v2_neigh = code.different_neighbour(v2, v1)
                    a, b = v1.even(), v2.uneven()
                    c, d = v1.uneven(), v2.even()
                    if abs(a - b) != 1:
                        a, b = v1.uneven(), v2.even()
                        c, d = v1.even(), v2.uneven()
                    cross_to_cross = (min(c,d)+1,max(c,d)-1)
                    #print "cross to cross",cross_to_cross,"for",v1,v2,a,b
                    our_cross = [cr.has_values(cross_to_cross) for cr in code]
                    #print our_cross,code
                    if sum(our_cross)==1:
                        our_cross = [cr for cr in code if cr.has_values(cross_to_cross)][0]
                        e,f = our_cross.even(),our_cross.uneven()
                        if a == v1.uneven():
                            e, f = our_cross.uneven(), our_cross.even()

                        _1 = Crossing(a.val,a.top,f.val,not a.top)
                        _2 = Crossing(b.val,b.top,e.val,not b.top)
                        if e > f:
                            if c > d:
                                _3 = Crossing(c.val, e.top, d.val, f.top)
                            else:
                                _3 = Crossing(c.val, f.top, d.val, e.top)
                        else:
                            if c > d:
                                _3 = Crossing(c.val, f.top, d.val, e.top)
                            else:
                                _3 = Crossing(c.val, e.top, d.val, f.top)
                        if VERBOSE: print( "was",v1,v2,our_cross, "will be",_1,_2,_3)
                        changed = True
                        code.pop(v1)
                        code.pop(v2)
                        code.pop(our_cross)
                        _1.skein = v1.skein
                        _2.skein = v2.skein
                        _3.skein = our_cross.skein
                        code.add(_1)
                        code.add(_2)
                        code.add(_3)
                        return changed # added 11.02
    return changed

def fix_dbl_loop(code, dbl_loops):
    """Corrects numbering after removing dbl loop"""
    m = (len(code) + len(dbl_loops) * 2) * 2
    ad = [0 for _ in range(m)]
    for loop in dbl_loops:
        for p in flatten(loop):
            for _ in range(p, m + 1):
                ad[_ - 1] += 1
    for k, v in enumerate(code):
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        code[k] = v


def just_one_in(loop,l):
    return (l[0] in range(loop.min().val,loop.max().val)) != (l[1] in range(loop.min().val, loop.max().val))

def both_in(loop,l):
    return (l[0] in range(loop.min().val,loop.max().val)) and (l[1] in range(loop.min().val, loop.max().val))


def can_be_untwisted(code):
    trash = []
    for i,loop in enumerate(code.gen_short_loop()):
        if not any(just_one_in(loop,x) for x in code if x!=loop):
            trash.append(loop)
    ch = bool(trash)
    while trash:
        loop = trash.pop()
        code.remove(loop)
        for cross in code:
            if both_in(loop,cross):
                cross.reverse_topo()
    return ch

def minimize_full_code(code):
    maks = len(code)*2
    num_possibilities = maks
    possible_codes = []

    tmp_code = sorted([((l[0].top,l[0].val),(l[1].top,l[1].val)) for l in code],key=lambda x:x[0][1] if x[0][1]%2 else x[1][1])
    possible_codes.append(tmp_code)
    while num_possibilities:
        num_possibilities -= 1
        next_code = []
        for i,c in enumerate(possible_codes[-1]):
            (s1,v1),(s2,v2) = c
            #s1,v1 = l1.top,l1.val
            #s2,v2 = l2.top,l2.val
            if v1 > 1:
                v1 -= 1
            else:
                v1 = maks
            if v2 > 1:
                v2 -= 1
            else:
                v2 = maks
            next_code.append(((s1,v1),(s2,v2)))
        next_code = sorted(next_code,key=lambda x:x[0][1] if x[0][1]%2 else x[1][1])
        possible_codes.append(next_code)
    #for p in possible_codes:
    #    print(p)
    #for y in possible_codes:
    #    print([x[0][1] if x[1][1]%2 else x[1][1] for x in y])
    possible_codes = sorted(possible_codes,key=lambda y: [x[0][1] if x[1][1]%2 else x[1][1] for x in y])
    #print("min",possible_codes[0])
    return possible_codes[0]

def minimize_sub_code(code):
    maks = len(code)*2
    num_possibilities = maks
    possible_codes = []
    tmp_code = [(2*i+1,code[i]) for i in range(len(code))]
    #tmp_code = sorted([((l[0].top,l[0].val),(l[1].top,l[1].val)) for l in code],key=lambda x:x[0][1] if x[0][1]%2 else x[1][1])
    possible_codes.append(tmp_code)
    while num_possibilities:
        num_possibilities -= 1
        next_code = []
        for i,c in enumerate(possible_codes[-1]):
            v1,v2 = c
            #s1,v1 = l1.top,l1.val
            #s2,v2 = l2.top,l2.val
            if v1 > 1:
                v1 -= 1
            else:
                v1 = maks
            if v2 > 1:
                v2 -= 1
            else:
                v2 = maks
            next_code.append((v1,v2))
        next_code = sorted(next_code,key=lambda x:x[0] if x[0]%2 else x[1])
        possible_codes.append(next_code)
    #for p in possible_codes:
    #    print(p)
    #for y in possible_codes:
    #    print([x[0] if x[1]%2 else x[1] for x in y])
    possible_codes = sorted(possible_codes,key=lambda y: [x[0] if x[1]%2 else x[1] for x in y])
    #print("min",possible_codes[0])
    return [x[0] if x[1]%2 else x[1] for x in possible_codes[0]]


def dowker_loop(code):
    """Finds single loops (twisted chain with no crossings)"""
#    loops = []
    trashes = []
    for i,v in enumerate(code):
        if neighbours_val(v.min().val,v.max().val,len(code)):
            trashes.append(i)
    ch = bool(trashes)
    while trashes:
        code.remove(trashes.pop())
    return ch


def fix_loop(code, loops):
    """Corrects numbering after removing single loop"""
    m = (len(code) + len(loops)) * 2
    ad = [0 for _ in range(m)]
    for loop in loops:
        for p in loop:
            for _ in range(p, m + 1):
                ad[_ - 1] += 1
    for i, v in enumerate(code):
        if ad[abs(v[0]) - 1]%2:
            if even(v) > 0:
                v = [even(v),-uneven(v)]
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        code[i] = v


def find_axis(atoms):
    atoms = [_.vec for _ in atoms]
    code = get_dt_code(atoms)
    if code is None:
        return None
    """r=1
    while any(c[0]%2==c[1]%2 for c in code):
        if VERBOSE: print("Changing axis")
        if r==3:
            #print "Need more axes changes"
            break
        for i, e in enumerate(atoms):
            e = [e[r % 3], e[(r + 1) % 3], e[(r + 2) % 3]]
            atoms[i] = e
        code = get_dt_code(atoms)
        if VERBOSE: print("got",code)
        if code is None:
            return None
        r+=1"""

    #rotate y axis

    sinT = math.sin(15)
    cosT = math.cos(15)
    r = 1
    while any(c[0] % 2 == c[1] % 2 for c in code):
        if VERBOSE: print("Changing axis by 15deg")
        if r == 23:
            break
        for i, e in enumerate(atoms):
            x,y,z = e
            nx  = x*cosT - z *sinT
            nz = z*cosT + x*sinT
            e = [nx,y,nz] #[e[r % 3], e[(r + 1) % 3], e[(r + 2) % 3]]
            atoms[i] = e
        code = get_dt_code(atoms)
        if code is None:
            return None
        r+=1
    return code


def from_tuples(atoms):
    code=Code()
    code.read_in(atoms)
    return code

def dowker_code(atoms, from_atoms=True):

    if from_atoms:
        code = find_axis(atoms)
    else:
        code = from_tuples(atoms)
    if code is None:
        return "01",[]
    changed = 1
    dc = code
    #print("will check",dc)
    dc.check_yo()
    #print(dc)
    dc.mark_first()
    history = CodeHistory()
    history.add(dc)
    if VERBOSE: print("code",code)
    while changed and dc.dowker_code():
        changed = dowker_loop(dc)
        dc.check_yo()
        if VERBOSE: print( "single",dc,changed)
        changed = dowker_dbl_loops(dc) or changed  # if not loops else False #bo po co liczyc na zapas
        dc.check_yo()
        if VERBOSE: print ("double",dc,changed)
        changed = can_be_untwisted(dc) or changed  # if not loops else False #bo po co liczyc na zapas
        dc.check_yo()
        #        fix_dbl_loop(dc, dbl_loops)
        if VERBOSE: print ("untwist", dc,changed)
        changed2 = unfurl(dc)
        dc.check_yo()
        if VERBOSE: print ("unfurl",dc,changed2)
        changed = changed or changed2
        dc.check_yo()
        while changed2:
            changed2 = unfurl(dc)
            dc.check_yo()
            if VERBOSE: print ("unfurl",dc,changed2)
        changed = changed or third_redei(dc)
        dc.check_yo()
        if VERBOSE: print("3rd redei", dc, changed)
        if not history.add(dc):
            if VERBOSE: print("Nothing changed, changing perspective")
            break

    if VERBOSE: print ("Finally",dc.mod_dowker_code())

    if not dc.dowker_code():
        return "01",[]

    #if VERBOSE: print("before changing perspective",dc)
    dc.start_later_by(6)

    changed = 1
    dc = code
    while changed and dc.dowker_code():
        changed = dowker_loop(dc)
        dc.check_yo()
        if VERBOSE: print ("single",dc)
        changed = dowker_dbl_loops(dc) or changed  # if not loops else False #bo po co liczyc na zapas
        dc.check_yo()
        if VERBOSE: print ("double",dc)
        changed = can_be_untwisted(dc) or changed  # if not loops else False #bo po co liczyc na zapas
        dc.check_yo()
        #        fix_dbl_loop(dc, dbl_loops)
        if VERBOSE: print ("untwist", dc)
        changed2 = unfurl(dc)
        dc.check_yo()
        if VERBOSE: print ("unfurl",dc)
        changed = changed or changed2
        while changed2:
            changed2 = unfurl(dc)
            dc.check_yo()
            if VERBOSE: print ("unfurl",dc)
        changed = changed or third_redei(dc)
        dc.check_yo()
        if not history.add(dc):
            break
    if VERBOSE: print ("Finally",dc)

    if not dc.dowker_code():
        return "01",[]

    sub_codes = find_permutations_in(dc.dowker_code())
    if VERBOSE: print (sub_codes)
    min_sub_codes = []
    for _s in sub_codes:
        min_sub_codes.append(minimize_sub_code(_s))

    translated = translate_dt_list(min_sub_codes)
    if VERBOSE: print (translated)

    #print("Koncowy",dc)
    #print("Koncowy", dc.dowker_code())
    #print(minimize_code(dc))


    return translated,dc.mod_dowker_code()#sub_codes

if __name__ == "__main__":
    #import sys
    #lista = eval("[" + sys.argv[1].replace('.',"") + "]")
    #print lista
    lista = [(-4,11),(3,-12),(-2,7),(1,-8),(6,-9),(-10,13),(-5,14)]
    #code = Code()
    #code.read_in(lista)
    dowker_code(lista,from_atoms=False)
#    print find_permutations_in([8,6,10,2,4,14,16,12])
