from __future__ import print_function
import numpy as np
from math import sqrt, acos
tmp = 0


def add_v3v3(v0, v1):
    return (
        v0[0] + v1[0],
        v0[1] + v1[1],
        v0[2] + v1[2],
    )


def sub_v3v3(v0, v1):
    return (
        v0[0] - v1[0],
        v0[1] - v1[1],
        v0[2] - v1[2],
    )


def dot_v3v3(v0, v1):
    return (
            (v0[0] * v1[0]) +
            (v0[1] * v1[1]) +
            (v0[2] * v1[2])
    )


def len_squared_v3(v0):
    return dot_v3v3(v0, v0)


def mul_v3v3(v0, v1):
    x = v0[1] * v1[2] - v0[2] * v1[1]
    y = v0[2] * v1[0] - v0[0] * v1[2]
    z = v0[0] * v1[1] - v0[1] * v1[0]
    return (x, y, z)


def get_normal(p1, p2, p3):
    a, b, c = np.cross((p2 - p1), (p3 - p1))
    d = sum(np.array([a, b, c]) * p1) * -1
    return (a, b, c, d)


def mul_v3_fl(v0, f):
    #    return (
    #        v0[0] * f,
    #        v0[1] * f,
    #        v0[2] * f,
    #        )
    return np.array(v0) * f


def get_outside_points(Nend, Cend, dist, sec_vec=(1., 1., 1.)):
    os_bialka = Nend - Cend  # sub_v3v3(Nend,Cend)
    perpendicular = np.cross(os_bialka, sec_vec)  # mul_v3v3(os_bialka,sec_vec)
    if sum(perpendicular) == 0:
        perpendicular = np.cross(os_bialka, np.array([1., 2., 3.]))  # mul_v3v3(os_bialka,(1.,2.,3.))
    per_len = point_distance(np.array([0., 0., 0.]), perpendicular)
    if per_len < dist:
        perpendicular = perpendicular * dist / (per_len)  # mul_v3_fl(perpendicular, dist/per_len)
    return [Cend + perpendicular, Nend + perpendicular]


def check_if_all_acutes(Nend, Cend, point):
    a = point_distance(Nend, Cend)
    b = point_distance(Cend, point)
    c = point_distance(point, Nend)
    #conditions = [(a ** 2 + b ** 2 > c ** 2), (a ** 2 + c ** 2 > b ** 2)]  # ,(c**2+b**2 > a**2)] #TODO changed 13.12.18
    #based on https://stackoverflow.com/questions/25877875/triangle-type-given-three-points-acute-obtuse-right
    conditions = [(a ** 2 + b ** 2 > c ** 2), (a ** 2 + c ** 2 > b ** 2), (c ** 2 + b ** 2 > a ** 2)]
    return all(conditions), max(a, b)

def check_if_end_angles_acute(Nend,Cend,point):
    """from https://math.stackexchange.com/questions/361412/finding-the-angle-between-three-points"""
    # angle at Nend
    abN = sub_v3v3(Nend, Cend)
    bcN = sub_v3v3(point, Nend)
    acN = dot_v3v3(abN, bcN) > 0
    # angle at Cend
    ab = sub_v3v3(Cend,Nend)
    bc = sub_v3v3(point,Cend)
    acC = dot_v3v3(ab,bc)>0
    return acN and acC, max(bcN,bc)


def calc_area(p1, p2, p3):
    return get_triangle_area(point_distance(p1, p2), point_distance(p2, p3), point_distance(p3, p1))


def checkIfLED(atoms):
    try:
        max_dist = 0
        Nend = atoms[0].vec  # A
        Cend = atoms[-1].vec  # B
        for p in atoms[1:-1]:
            a, b, c, d = get_normal(Nend, Cend, p.vec)
            abc = np.array([a, b, c])
            sign = sum(abc * atoms) + d >= 0  # (a*atoms[0].x + b*atoms[0].y + c*atoms[0].z +d) >= 0
            for o in atoms[1:-1]:
                ov = o.vec
                if o != p:
                    s = (sum(abc * ov) + d > 0) or (sum(abc * ov) + d == 0 and not is_on_line(ov, Nend, Cend))
                    if sign != s:
                        break
            else:
                return True
        return False
    except:
        print( "wwyalilo sie w checkIfLED")
        return False

def checkIfProperLED(atoms):
    Nend,Cend = atoms[0],atoms[-1]

    points = atoms[1:-1]
    for i,p1 in enumerate(points):
        for j,p2 in enumerate(points[i+1:]):
            for p3 in points[i+j:]:
                if found_crossing(p1, p2, p3, Nend, Cend, get_normal(p1, p2, p3)) is not False:
                    return False,None
    return True,[]

def divide_into_chains(atoms):
    atom_lists = []
    alist = []
    for a in atoms:
        alist.append(a.vec)
        if a.end:
            atom_lists.append(alist)
            alist  = []
    atom_lists.append(alist)
    return atom_lists

def divide_into_bead_chains(atoms):
    atom_lists = []
    alist = []
    for a in atoms:
        alist.append(a)
        if a.end:
            atom_lists.append(alist)
            alist  = []
    atom_lists.append(alist)
    return atom_lists

def closeTheCurve_old(atoms):
    single_chain = (not any([a.end for a in atoms]))
    if single_chain:
        works,joining = checkIfHalfGlobe([a.vec for a in atoms])
        result = [joining if works else None]
    else:
        alists = divide_into_chains(atoms)
        result = []
        for al in alists:
            w,j = checkIfHalfGlobe(al) #### TODO !! Should direct joining path away from other chains
            result.append(j if w else None)
    return result

def closeTheCurve(atoms):
    works,joining = checkIfHalfGlobe([a.vec for a in atoms])
    result = joining if works else None
    return result


def checkIfHalfGlobe(atoms):
    ## TODO add possibility a LED like design ( all atoms on one hemisphere in relation to line between ends)
    """Checks if in relation to both ends, other atoms are on one hemisphere """
    max_dist = 0
    Nend = atoms[0]  # A
    Cend = atoms[-1]  # B
    for i, point in enumerate(atoms[1:-1]):
        #res, md = check_if_all_acutes(Nend, Cend, point)
        res,md = check_if_end_angles_acute(Nend, Cend, point)
        if res:
            max_dist = max(max_dist, md)
        else:
            return checkIfProperLED(atoms)

    outsiders = get_outside_points(Nend, Cend, max_dist) ### TODO here should direct joining away from other chains
    return True, outsiders

def isect_line_plane_v3_4d(p0, p1, plane, epsilon=1e-6):
    plane = np.array(plane)
    u = p1 - p0  # sub_v3v3(p1, p0)
    dot = plane[:3].dot(u)  # dot_v3v3(plane, u)

    if abs(dot) > epsilon:
        # calculate a point on the plane
        # (divide can be omitted for unit hessian-normal form).
        p_co = plane[:3] * (-plane[3] / len_squared_v3(plane))
        w = p0 - p_co  # sub_v3v3(p0, p_co)
        fac = -1 * (plane[:3].dot(w)) / dot  # -dot_v3v3(plane, w) / dot
        u = u * fac  # mul_v3_fl(u, fac)
        return p0 + u  # add_v3v3(p0, u)
    else:
        return None


def point_distance(p1, p2):
    return np.linalg.norm(p1 - p2)


def point_distance_2d(p1, p2):
    return sqrt(((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2))  # + (p1[2]-p2[2])**2))


def get_closest(p1, p2, points):  # first give the first moving point - larger movement on this side
    dists = []
    for p in points:
        a = point_distance(p1, p)
        b = point_distance(p2, p)
        dists.append((a + b, a, b, p))
    return sorted(dists)[0]  # closest point to first moved atom ##CHANGE


def intersect(p123normal, l0, l1):
    #    print get_normal(p1,p2,p3)
    return isect_line_plane_v3_4d(l0, l1, p123normal)


def get_triangle_area(p1p2dist, p2odist,
                      op1dist):  # traingle formed by obstacle, with defualt pulling vector at the bottom
    hperim = sum([p1p2dist, p2odist, op1dist]) / 2.
    area = sqrt(hperim * (hperim - p1p2dist) * (hperim - p2odist) * (hperim - op1dist))
    return area


def get_triangle_height(p1p2dist, p2odist,
                        op1dist):  # traingle formed by obstacle, with defualt pulling vector at the bottom
    area = get_triangle_area(p1p2dist, p2odist, op1dist)
    return area * 2 / p1p2dist


def get_point_on_line_closest_to_point(l0, l1, p):
    os_bialka = l0 - l1  # sub_v3v3(l0,l1)
    perpendicular = np.cross(os_bialka, p)  # mul_v3v3(os_bialka,p)
    # hopefully get_crossing computes crossing of two 3d lines
    pp = p + perpendicular  # add_v3v3(p,perpendicular)

    cross = get_crossing([l0, l1], [p, pp])
    return cross


def in_triangle(p1, p2, p3, x):
    global tmp
    c1 = p2 - p1  # sub_v3v3(p2,p1)
    c2 = p3 - p1  # sub_v3v3(p3,p1)
    c3 = x - p1  # sub_v3v3(x,p1)
    a = np.array([
        [c1[0], c2[0]],
        [c1[1], c2[1]]])
    b = np.array([[c3[0]], [c3[1]]])
    try:
        s = np.linalg.solve(a, b)
        tmp = s
    except np.linalg.linalg.LinAlgError:
        s = np.linalg.lstsq(a, b)[0]
        print ("LinAlgError:",tmp, s)
    return s[0] >= 0 and s[1] >= 0 and s[0] + s[1] <= 1


def found_crossing(p1, p2, p3, l0, l1, p123normal):
    inter = intersect(p123normal, l0, l1)
    if inter is not None and in_triangle(p1, p2, p3, inter) and is_on_line(inter, l0, l1):
        return inter
    else:
        return False


def get_list_of_points_on_line(p1, p2):
    line_vec = sub_v3v3(p2, p1)
    llen = point_distance(p1, p2)  # sqrt(len_squared_v3())
    unit_line_vec = mul_v3_fl(line_vec, llen)
    out = []
    for i in range(1, 10):
        out.append(add_v3v3(p1, mul_v3_fl(unit_line_vec, i / 10.)))
    return out


def get_triangle_height_unit_vector(p1, p2, p3):
    l0 = p2
    l1 = get_point_on_line_closest_to_point(p1, p3, p2)
    line = sub_v3v3(l1, l0)
    llen = 1 / point_distance(l0, l1)
    line = [x * llen for x in line]
    return line


def point_outward_moved_by_one(p1, p2, p3, p):
    return sum(p, get_triangle_height_unit_vector(p1, p2, p3))


def point_to_line_distance(p1, l0, l1):
    return point_distance(p1, get_point_on_line_closest_to_point(l0, l1, p1))


def is_on_line(point, l0, l1):
    lx = min(l0[0], l1[0])
    rx = max(l0[0], l1[0])
    ly = min(l0[1], l1[1])
    ry = max(l0[1], l1[1])
    lz = min(l0[2], l1[2])
    rz = max(l0[2], l1[2])

    return lx <= point[0] <= rx and ly <= point[1] <= ry and lz <= point[2] <= rz


def get_point_on_line(p1, p3, dist):
    p = p1 + ((p3 - p1) * dist / point_distance(p1, p3))
    return p


def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])

    # Return true if line segments AB and CD intersect


def intersect_2d(A, B, C, D):
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def get_middlepoint(p1, p2):
    mid = (p1 + p2) * .5
    return mid


def get_somepoint(p1, p2, mul):
    mid = mul_v3_fl(add_v3v3(p1, p2), mul)
    return mid


def get_point_over_obstacle(obstacle, p2, k=False):  # k = dist from second onne # TODO brakuje absoluta?!
    dist = k
    if dist:
        p = add_v3v3(p2, mul_v3_fl(sub_v3v3(p2, obstacle), dist / point_distance(p2, obstacle)))
    else:
        p = mul_v3_fl(add_v3v3(obstacle, p2), .5)
    return p


def trilateration(s1, s2, s3):
    P1P2 = plane_from_spheres(s1, s2)
    P2P3 = plane_from_spheres(s2, s3)
    P3P1 = plane_from_spheres(s3, s1)

    # teraz przeciecie tych plaszczyzna da mi linie

    p1 = s1[:3]
    p2 = s2[:3]
    p3 = s3[:3]

    Q = get_normal(p1, p2, p3)

    H = four_plane_intersect(P1P2, P2P3, P3P1, Q)

    U = Q[:3]  # = normal z Q?

    t = sqrt(s1[3] ** 2 - point_distance(H, p1) ** 2)

    Point1 = sub_v3v3(H, mul_v3_fl(U, t))
    Point2 = add_v3v3(H, mul_v3_fl(U, t))

    return (Point1, Point2)


def four_plane_intersect(p1, p2, p3, p4):
    a = np.array([
        [p1[0], p2[0], p3[0], p4[0]],
        [p1[1], p2[1], p3[1], p4[1]],
        [p1[2], p2[2], p3[2], p4[2]]])
    b = np.array([[p1[3]], [p2[3]], [p3[3]], [p4[3]]])  # ,[c3[2]]])
    return np.linalg.solve(a, b)


def plane_from_spheres(s1, s2):
    x0, y0, z0, r0 = s1
    x1, y1, z1, r1 = s2
    A = 2 * x1 - 2 * x0
    B = 2 * y1 - 2 * y0
    C = 2 * z1 - 2 * z0
    D = x0 ** 2 + y0 ** 2 + z0 ** 2 - x1 ** 2 - y1 ** 2 - z1 ** 2 - r0 ** 2 - r1 ** 2
    return (A, B, C, -D)  # chyba -


def line_sphere_insect(line, sphere):
    p0, p1 = line
    c, r = sphere

    A = (p0[0] - c[0]) ** 2 + (p0[1] - c[0]) ** 2 + (p0[2] - c[2]) ** 2 - r ** 2
    C = (p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2
    B = (p1[0] - c[0]) ** 2 + (p1[1] - c[0]) ** 2 + (p1[2] - c[2]) ** 2 - A - C - r ** 2
    D = B ** 2 - 4 * A * C
    t = False
    if D == 0:
        t = -B / (2 * A)
    elif D > 0:
        x0 = (-B - sqrt(D)) / (2 * A)
        x1 = (-B + sqrt(D)) / (2 * A)
        t = min(x0, x1)
    if t:
        x = p0[0] * (1 - t) + p1[0] * t
        y = p0[1] * (1 - t) + p1[1] * t
        z = p0[2] * (1 - t) + p1[2] * t
        return (x, y, z)


def get_crossing_2d(p, pr, q, qs):
    def sub(p, q):
        return (p[0] - q[0], p[1] - q[1])

    def cross(p, q):
        return p[0] * q[1] - p[1] * q[0]

    def mul(p, s):
        return list(map(lambda x: x * s, p))

    def sum(p, q):
        return (p[0] + q[0], p[1] + q[1])

    r = (pr[0] - p[0], pr[1] - p[1])
    s = (qs[0] - q[0], qs[1] - q[1])

    t = cross(sub(q, p), s) / cross(r, s)
    assert 0 <= t and t <= 1
    return sum(p, mul(r, t))


def get_crossing(l0, l1):
    M = [(l0[1] - l0[0]), (l1[1] - l1[0]) * -1]
    #    R = sub_v3v3(l1[0],l0[0])
    R = l1[0] - l0[0]
    a = np.array(M)
    b = np.array(R)
    s = np.linalg.lstsq(a.T, b)[0]
    c = l0[0] + (l0[1] - l0[0]) * s[0]
    return c


def get_kite_triangle(p1, p2, p3, p):
    c = get_crossing((p1, p), (p2, p3))
    return (p, c, p3)


if __name__ == "__main__":
    p1 = (0., 0., 0.)
    p2 = (0., 3., 0.)
    p3 = (3., 0., 0.)
    l0 = (1., 0.5, 1.)
    l1 = (1., 0.5, -1.)
    inter = intersect(p1, p2, p3, l0, l1)
    a1 = [0, 0, 0]
    a2 = [5, 5, 5]
    b1 = [0, 5, 5]
    b2 = [5, 0, 0]
    g = get_crossing((a1, a2), (b1, b2))
