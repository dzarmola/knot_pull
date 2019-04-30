from .vector_ops import *


def find_z(l0, l1, p):
    diffx = l1[0] - l0[0]
    diffpx = p[0] - l0[0]
    diffz = l1[2] - l0[2]
    diffpz = diffpx * diffz / diffx
    return l0[2] + diffpz


def get_dt_code(atoms):
    cnt = 1
    crosses = {}
    for l in xrange(len(atoms) - 1):
        this_cross = []
        for k in xrange(len(atoms) - 1):
            if l == k or abs(l - k) == 1: continue
            line1 = [atoms[l], atoms[l + 1]]
            line2 = [atoms[k], atoms[k + 1]]
            if intersect_2d(*(line1 + line2)) and not len(set(tuple(x[:2]) for x in line1 + line2)) < 4:
                #                cross = get_crossing((atoms[l].vec, atoms[l+1].vec),([atoms[k].vec,atoms[k+1].vec]))
                #                print get_crossing((atoms[l].vec, atoms[l+1].vec),([atoms[k].vec,atoms[k+1].vec])),get_crossing((atoms[l].vec, atoms[l+1].vec),([atoms[k+1].vec,atoms[k].vec]))
                cross = get_crossing_2d(atoms[l], atoms[l + 1], atoms[k], atoms[k + 1])
                #                print atoms[l].vec, atoms[l+1].vec,atoms[k].vec,atoms[k+1].vec
                #                print "##########"
                this_cross.append((point_distance_2d(atoms[l], cross), k, cross))
        #        print sorted(this_cross)
        for _, k, cross in sorted(this_cross):
            #                print "cross",l,k
            val = 1
            if cnt % 2 == 0 and find_z(atoms[l], atoms[l + 1], cross) > find_z(atoms[k], atoms[k + 1], cross):
                val = -1
            crosses[tuple(sorted([l, k]))] = crosses.get(tuple(sorted([l, k])), []) + [cnt * val]
            cnt += 1
    return crosses


def translate_dt_code(code):
    DT_CODES = {(4, 6, 2): "31",
                (4, 6, 8, 2): "41",
                (6, 8, 10, 2, 4): "51",
                (4, 8, 10, 2, 6): "52",
                (4, 8, 12, 10, 2, 6): "61",
                (4, 8, 10, 12, 2, 6): "62",
                (4, 8, 10, 2, 12, 6): "63"}
    #    DT_CODES_STR = {}
    code_str = "-".join(map(lambda x: str(abs(x)), code))
    for k, v in DT_CODES.items():
        k_new = "-".join(map(str, k + k))
        #        DT_CODES_STR[k_new] = v
        if code_str in k_new:
            return v
    return "{}0:{}".format(len(code),code)
    #return "999:%s" % code_str

def shortest_loop(code):
    return sorted(code,key=lambda x:abs(x[0]-x[1]))

def same_side_cross(v1,v2):
    p1 = v1[0] if v1[0] % 2 == 0 else v1[1] #even
    p2 = v2[0] if v2[0] % 2 == 0 else v2[1] #even
    return p1 / abs(p1) != p2 / abs(p2) # rozne znaki = oba crossingi z tej same strony

def unfurl(code): 
    """Takes dowker code sorted by length
    unfurls crossed loops which are passed through another line"""
    changed = False
    for _,loop in enumerate(code):
        if abs(loop[1]-loop[0]) == 3:
            to_del = [_]
            i,a = min(loop),max(loop)
            rg_i_a = range(i+1,a)
            fr,fz = None,None #first in range, first_zip_line
            sr,sz = None,None #second in range, second_zip_line
            for _,l in enumerate(code):
                if abs(l[0]) in rg_i_a:
                    to_del.append(_)
                    if abs(l[0]) == rg_i_a[0]:
                        fr,fz = l
                    else:
                        sr,sz = l
                if abs(l[1]) in rg_i_a:
                    to_del.append(_)
                    if abs(l[1]) == rg_i_a[0]:
                        fz,fr = l
                    else:
                        sz,sr = l
            if abs(fr)-abs(sr) == abs(fz)-abs(sz): #both are increasing
                if same_side_cross((fz,fz),loop) or same_side_cross((sr,sz),loop):
                    # there are two consecutive with one direction
                    print "will change",[code[_] for _ in to_del]
                    print fr,fz,sr,sz
                    for _ in sorted(to_del,reverse=True):
                        code.pop(_)
                    code.append(fix_minus([-sz,i]))
                    code.append(fix_minus([-fz,i+1]))
    #                print code
                    code = fix_numbering(code,i+2)
    #                print code
    #                exit()
                    changed = True
                    break
    return code,changed

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
    if uneven(loop)<0:
        loop = map(abs,loop)
    return loop

def dowker_dbl_loops_old(atoms, crosses):
    cnt = len(crosses) * 2
    dbl_loops = []
    for k1, v1 in crosses.items():
        for k2, v2 in crosses.items():
            if k1 != k2:
                if abs(v1[0]) - abs(v2[0]) == -1 * (abs(v1[1]) - abs(v2[1])):
                    all = sorted(map(abs, v1 + v2))
                    crossA = get_crossing_2d(atoms[k1[0]], atoms[k1[0] + 1], atoms[k1[1]], atoms[k1[1] + 1])
                    crossB = get_crossing_2d(atoms[k2[0]], atoms[k2[0] + 1], atoms[k2[1]], atoms[k2[1] + 1])
                    if all[2] - all[1] == 1:
                        allK = sorted(map(abs, k1 + k2))
                        if not allK[1] in k1:
                            k1, k2 = k2, k1
                            ###### problem jest jak sa crossy 1,4 i 1,6 <-- 1,1,4,6, trzeba jakos oznaczy, ze chodzi raczej o cos z czym crossujemy
                        print k1, k2, allK, crosses
                        _ = k1[(k1.index(allK[1]) + 1) % 2]
                        d = find_z(atoms[allK[1]], atoms[allK[1] + 1], crossA) - find_z(atoms[_], atoms[_ + 1], crossA)
                        dd = d / abs(d)
                        _ = k2[(k2.index(allK[2]) + 1) % 2]
                        e = find_z(atoms[allK[2]], atoms[allK[2] + 1], crossA) - find_z(atoms[_], atoms[_ + 1], crossB)
                        ee = e / abs(e)
                        if ee == dd:
                            dbl_loops.append([allK[1], allK[2], [crossA, d], [crossB, e]])
                            break
                    elif all[0] - all[3] % cnt == 1:
                        raise IndexError("dbl loop na N/C koncu, halp.")

    return dbl_loops


def dowker_dbl_loops_old(crosses):
    dbl_loops = set([])
    trash = set([])
    for k1, v1 in crosses.items():
        for k2, v2 in crosses.items():  ## assumuje ze v2 sa dalej
            if k1 == k2: continue
            a1 = map(abs, v1)
            a2 = map(abs, v2)
            if (abs(a1[0] - a2[0]) == 1 and abs(a1[1] - a2[1]) == 1) or (
                    abs(a1[0] - a2[1]) == 1 and abs(a1[1] - a2[0]) == 1):
                # mam dwa sasiadujace punkty
                #            if min(map(abs,v2))-max(map(abs,v1))==1 or max(map(abs,v2))
                p1 = v1[0] if v1[0] % 2 == 0 else v1[1]
                p2 = v2[0] if v2[0] % 2 == 0 else v2[1]
                if p1 / abs(p1) != p2 / abs(p2):  # rozne znaki = oba crossingi z tej same strony
                    dbl_loops.add(tuple(sorted((tuple(a1), tuple(a2)))))
                    trash.add(k1)
                    trash.add(k2)
                    print "found dbl"
    for t in trash:
        del crosses[t]
    return dbl_loops

def dowker_dbl_loops(code):
    dbl_loops = set([])
    trash = set([])
    for k1, v1 in enumerate(code):
        for k2, v2 in enumerate(code[k1+1:]):  ## assumuje ze v2 sa dalej
            a1 = map(abs, v1)
            a2 = map(abs, v2)
            if (abs(a1[0] - a2[0]) == 1 and abs(a1[1] - a2[1]) == 1) or (
                    abs(a1[0] - a2[1]) == 1 and abs(a1[1] - a2[0]) == 1):
                # mam dwa sasiadujace punkty
                #            if min(map(abs,v2))-max(map(abs,v1))==1 or max(map(abs,v2))
                p1 = v1[0] if v1[0] % 2 == 0 else v1[1]
                p2 = v2[0] if v2[0] % 2 == 0 else v2[1]
                if p1 / abs(p1) != p2 / abs(p2):  # rozne znaki = oba crossingi z tej same strony
                    dbl_loops.add(tuple(sorted((tuple(a1), tuple(a2)))))
                    trash.add(k1)
                    trash.add(k2+k1+1)
                    #print "found dbl"
    #print "in dbl loop",code, dbl_loops,trash
    trash = sorted(trash)
    while trash:
        code.pop(trash.pop())
    return dbl_loops


flatten = lambda x: list(x[0]) + list(x[1])
sign = lambda x: -1 if x < 0 else 1
add = lambda x, y: x - y if sign(x) == 1 else -(abs(x) - y)


def fix_dbl_loop(code, dbl_loops):
    m = (len(code) + len(dbl_loops) * 2) * 2
    ad = [0 for _ in xrange(m)]
    for loop in dbl_loops:
        for p in flatten(loop):
            for _ in xrange(p, m + 1):
                ad[_ - 1] += 1
    #    print ad
    #    print crosses
    for k, v in enumerate(code):
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        #        v = [v[0]-ad[abs(v[0])-1],v[1]-ad[abs(v[1])-1]]
        code[k] = v

def fix_dbl_loop_old(crosses, dbl_loops):
    m = (len(crosses) + len(dbl_loops) * 2) * 2
    ad = [0 for _ in xrange(m)]
    for loop in dbl_loops:
        for p in flatten(loop):
            for _ in xrange(p, m + 1):
                ad[_ - 1] += 1
    #    print ad
    #    print crosses
    for k, v in crosses.items():
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        #        v = [v[0]-ad[abs(v[0])-1],v[1]-ad[abs(v[1])-1]]
        crosses[k] = v


#    print crosses

def dowker_loop_old(crosses):
    loops = []
    trashes = []
    for k, v in crosses.items():
        if abs(abs(v[0]) - abs(v[1])) == 1:
            print "found loop"
            loops.append(v)
            trashes.append(k)
    for k in trashes:
        del crosses[k]
    return loops

def dowker_loop(code):
    loops = []
    trashes = []
    for i,v in enumerate(code):
        if abs(abs(v[0]) - abs(v[1])) == 1:
            loops.append(v)
            trashes.append(i)
    #print "in dowker loop",code,loops,trashes
    while trashes:
        code.pop(trashes.pop())
    #print "in dowker loop", code, loops, trashes
    return loops

def fix_loop(code, loops):
    m = (len(code) + len(loops)) * 2
    ad = [0 for _ in xrange(m)]
    for loop in loops:
        for p in loop:
            for _ in xrange(p, m + 1):
                ad[_ - 1] += 1
    for i, v in enumerate(code):
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        code[i] = v

def fix_loop_old(crosses, loops):
    m = (len(crosses) + len(loops)) * 2
    ad = [0 for _ in xrange(m)]
    for loop in loops:
        for p in loop:
            for _ in xrange(p, m + 1):
                ad[_ - 1] += 1
    #    print ad
    #    print crosses
    for k, v in crosses.items():
        v = [add(v[0], ad[abs(v[0]) - 1]), add(v[1], ad[abs(v[1]) - 1])]
        #        v = [v[0]-ad[abs(v[0])-1],v[1]-ad[abs(v[1])-1]]
        crosses[k] = v


def fix_dbl_loop_old(atoms, loop):
    i, j, cA, cB = loop
    atoms[i] = cA[0][:2] + [cA[1]]
    atoms[j] = cB[0][:2] + [cB[1]]
    for e in xrange(j - 1, i, -1):
        atoms.pop(e)


def plot_atoms(atoms):
    xs = [_[0] for _ in atoms]
    ys = [_[1] for _ in atoms]
    zs = [_[2] for _ in atoms]
    zs = map(lambda x: zs.index(x), zs)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(xs, ys, '-')  # ,linewidth=2.0,color="white")
    plt.scatter(xs, ys, c=zs)
    plt.show()


def same_side(c1, c2, i):
    #	g1 = (c1[i]%2==0 and c1[i]<0) or (c1[i]%2 and c1[(i+1)%2]>0)
    #	g2 = (c2[i]%2==0 and c2[i]<0) or (c2[i]%2 and c2[(i+1)%2]>0)
    #	return g1==g2
    p2 = c2[1] if c2[0] % 2 else c2[0]
    # dodac dol
    #	print "ss",c1,i,c1[(i+1)%2],p2
    return (c1[i] % 2 and c1[(i + 1) % 2] > 0 and p2 < 0) or (c1[i] % 2 == 0 and c1[i] < 0 and p2 > 0)


#    return (c1[0]*c1[1]>0 and c2[0]*c2[1]<0) or (c1[0]*c1[1]<0 and c2[0]*c2[1]>0)

def od(val, d=1):
    return (abs(val) + d) * (-1 if val < 0 else 1)


def index(lista, co):
    for i in xrange(len(lista)):
        if abs(lista[i]) == co:
            return i
    else:
        raise IndexError


def remove_overlay(clist, idx=0):
    print clist
    from copy import copy
    drabinka = sorted(clist, key=lambda x: x[0] if x[0] % 2 else x[1])
    drabinka = [sorted(x, key=lambda x: ((x % 2) + 1) % 2) for x in drabinka]
    kolejka = [False for x in drabinka + drabinka]
    for elem in drabinka:
        kolejka[abs(elem[0]) - 1] = elem
        kolejka[abs(elem[1]) - 1] = elem
    #    print drabinka
    #    print kolejka

    pary = []
    print "kolejka", kolejka
    for i in xrange(1, len(kolejka) - 1):
        #	print kolejka[i],kolejka[i+1]
        if same_side(kolejka[i], kolejka[i + 1], index(kolejka[i], i + 1)):
            pary.append(i)
    #	    print kolejka[i],kolejka[i+1]
    print "pary", pary
    if idx >= len(pary): return clist, None
    i = pary[-1 - idx]
    next = kolejka[(i + 2) % len(kolejka)]
    s1 = kolejka[i]
    s2 = kolejka[i + 1]
    print "wybrane", s1, s2
    i = 0 if s2[0] - s1[0] == 0 else 1
    dc = [copy(x) for x in drabinka]
    dc.pop(dc.index(s1))
    dc.pop(dc.index(s2))
    print "dc", dc
    next = dc.index(next)
    print "next", next
    for p in xrange(len(dc)):
        for _ in [0, 1]:
            for s in sorted(map(abs, s1 + s2), reverse=True):
                if abs(dc[p][_]) > s:
                    dc[p][_] = od(dc[p][_], -1)
    #            if abs(dc[p][_]) > abs(s1[0]):
    # print "1",dc[p][_],s1[0]
    #		dc[p][_] = od(dc[p][_], -1)
    # if abs(dc[p][_]) > abs(s2[0]):
    # print "2",dc[p][_],s2[0]
    # dc[p][_] = od(dc[p][_], -1)
    # if abs(dc[p][_]) > abs(s1[1]):
    # print "3",dc[p][_],s1[1]
    # dc[p][_] = od(dc[p][_], -1)
    # if abs(dc[p][_]) > abs(s2[1]):
    # print "4",dc[p][_],s2[1]
    # dc[p][_] = od(dc[p][_], -1)
    new = [min(dc[next])]  # fix minusy!!
    print "new dc", dc  # ,new,dc[next]
    print "new", new
    for p in xrange(len(dc)):
        for _ in [0, 1]:
            if abs(dc[p][_]) >= abs(new[0]):
                # print 5,dc[p],new[0]
                dc[p][_] = od(dc[p][_], 1)
    print "newer dc", dc
    bads = [p for p in dc if p[0] % 2 == p[1] % 2]
    _ = min(max(map(abs, b)) for b in bads)
    print "bads", bads
    if not _ % 2 and ((s2[i] % 2 and s1[i] < 0) or (s1[i] % 2 and s2[i] < 0)):
        _ *= -1
    elif _ % 2 and ((s2[i] % 2 and s1[i] < 0) or (s1[i] % 2 and s2[i] < 0)):
        new[0] *= -1
    new.append(_)
    print new
    for p in xrange(len(dc)):
        for _ in [0, 1]:
            if abs(dc[p][_]) >= abs(new[1]):
                dc[p][_] = od(dc[p][_], 1)
    print "even newer", dc
    bads = [p for p in dc if p[0] % 2 == p[1] % 2]
    if bads:
        return dc, idx + 1
    dc.append(new)

    # fix -
    for e in xrange(len(dc)):
        np = 0 if dc[e][0] % 2 else 1
        p = (np + 1) % 2  # 0 if dc[e][0]%2 else 1
        if dc[e][np] < 0:
            dc[e][np] *= -1
        elif dc[e][1] > 0 and p == 0 and dc[e][p] > 0:
            dc[e][p] *= -1

    #        if (dc[e][0]%2 and dc[e][0]<0) or (dc[e][1]%2 and dc[e][0]<0):
    #            dc[e][0] = abs(dc[e][0])
    #            dc[e][1] = abs(dc[e][1])
    #        elif dc[e][0]%2 == 0 and dc[e][1]>0:
    #            dc[e][0] *= -1
    print "final", dc
    return dc, 0

uneven = lambda x: x[0] if x[0]%2 else x[1]
even = lambda x: x[1] if x[0]%2 else x[0]

def shorten_code(lista):
    lista = sorted(lista,key=lambda x:uneven(x))
    return tuple([even(x) for x in lista])


def dowker_code(atoms, rotated=False):

    if not rotated:
        atoms = [_.vec for _ in atoms]  # + [atoms[0].vec]
    else:
        for i, e in enumerate(atoms):
            e = [e[rotated % 3], e[(rotated + 1) % 3], e[(rotated + 2) % 3]]
            atoms[i] = e

    ##### DT code #######
    crosses = get_dt_code(atoms)
    dt_code_long = sorted(crosses.values(), key=lambda x: x[0] if x[0] % 2 else x[1])
    dt_code_short = [x[0] if x[0] % 2 == 0 else x[1] for x in dt_code_long]
    ### Fixing based on DT code ####
    cnt = len(dt_code_long) * 2
    to_pop = []
    changed = 1
    dc = crosses.values()
    print dc
    while changed:
        changed = 0
        loops = dowker_loop(dc)
        changed = bool(loops)
        fix_loop(dc, loops)
        print "single",dc,loops
        dbl_loops = dowker_dbl_loops(dc)  # if not loops else False #bo po co liczyc na zapas
        changed = changed or bool(dbl_loops)
        fix_dbl_loop(dc, dbl_loops)
        print "double",dc,loops
        sort_dc = shortest_loop(dc)
        dc,changed2 = unfurl(sort_dc)
        print "unfrl",dc
        while changed2:
            sort_dc = shortest_loop(dc)
            dc,changed2 = unfurl(sort_dc)
        changed = changed or changed2
        dc = sorted(dc, key=lambda x:min(map(abs,x)))

    short_dc = shorten_code(dc)

    translated = translate_dt_code(short_dc) if len(short_dc) >= 3 else "0"
    #print "Final code:", short_dc, translated  # translate_dt_code(dt_code)
    return translated

    exit()
    ret = 0
    while ret is not None:
        print "DC IS", dc
        dc, ret = remove_overlay(dc, ret + 1)
    print dc,ret

    short_dc = shorten_code(dc)

    translated = translate_dt_code(short_dc) if len(short_dc) >= 3 else "0"
    print "Final code:", short_dc, translated  # translate_dt_code(dt_code)
    #    import matplotlib.lines as mlines
    ###### PLOT KNOT ######
    ##### END PLOT KNOT #####
    return translated

    exit()

    plot_atoms(atoms)
    if rotated:
        return
    else:
        print "rot1"
        dowker_code(atoms, 1)
        print "rot2"
        dowker_code(atoms, 1)
    exit()
    """xs = [_[0] for _ in atoms]
    ys = [_[1] for _ in atoms]
    zs = [_[2] for _ in atoms]
    zs = map(lambda x: zs.index(x), zs)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(xs,ys,'-')#,linewidth=2.0,color="white")
    plt.scatter(xs,ys,c=zs )
    plt.show()"""
    dt_code = dt_code_short
    while any(loops) or dbl_loops:  # abs(x[0]%cnt)-abs(x[1])%cnt == 1 for x in dt_code_long):
        #        print "in while",crosses,loops,[cr for cr in crosses]
        while loops:
            if crosses:
                plot_atoms(atoms)
            to_pop = []
            #        print "Have a loop!",loops
            for cr in crosses:
                if crosses[cr] in loops:
                    l, k = cr
                    cross = get_crossing_2d(atoms[l], atoms[l + 1], atoms[k], atoms[k + 1])
                    #                print l,k,cross
                    new_l = [cross[0], cross[1], find_z(atoms[l], atoms[l + 1], cross)]
                    new_k = [cross[0], cross[1], find_z(atoms[k], atoms[k + 1], cross)]
                    #                print new_l,new_k
                    atoms[l] = new_l
                    atoms[k + 1] = new_k
                    if k - l != 1:
                        to_pop += range(l + 1, k)
                    #                print "merguje",l,k
                    break

            for _ in sorted(list(set(to_pop)), reverse=True):
                atoms.pop(_)

            crosses = get_dt_code(atoms)
            dt_code_long = sorted(crosses.values(), key=lambda x: x[0] if x[0] % 2 else x[1])
            dt_code_short = [x[0] if x[0] % 2 == 0 else x[1] for x in dt_code_long]
            cnt = len(dt_code_long) * 2
            loops = [x for x in dt_code_long if abs((abs(x[0]) % cnt) - (abs(x[1]) % cnt)) == 1]

        crosses = get_dt_code(atoms)
        dbl_loops = dowker_dbl_loops(crosses)

        while dbl_loops:
            if dbl_loops:
                fix_dbl_loop(atoms, dbl_loops[0])
            crosses = get_dt_code(atoms)
            dbl_loops = dowker_dbl_loops(crosses)
        dt_code_long = sorted(crosses.values(), key=lambda x: x[0] if x[0] % 2 else x[1])
        dt_code_short = [x[0] if x[0] % 2 == 0 else x[1] for x in dt_code_long]
        cnt = len(dt_code_long) * 2
        loops = [x for x in dt_code_long if abs((abs(x[0]) % cnt) - (abs(x[1]) % cnt)) == 1]

        dt_code = dt_code_short  # [x[0] if x[0]%2==0 else x[1] for x in sorted(crosses.values(),key=lambda x: x[0] if x[0]%2 else x[1])]
    #        dbl_loops = dowker_dbl_loops(atoms,crosses)

    #    else:
    #        dt_code = dt_code_short
    #       print "No loops"
    #    print "koniec",crosses,dt_code,len(atoms)
    translated = translate_dt_code(dt_code) if len(dt_code) >= 3 else "0"
    print "Final code:", dt_code, translated  # translate_dt_code(dt_code)
    #    import matplotlib.lines as mlines
    ###### PLOT KNOT ######
    ##### END PLOT KNOT #####
    return translated
