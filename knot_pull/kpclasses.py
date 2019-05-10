from __future__ import print_function
from .vector_ops import point_distance
from numpy import ndarray, array as Vector


class Bead(object):
    def __init__(self, vec, type="CA", original_id=None):
        self.x = vec[0]
        self.y = vec[1]
        self.z = vec[2]
        if not isinstance(vec,ndarray):
            vec = Vector(vec)
        self.vec = vec
        self.prev_vec = None
        self.type = type
        self.movement = None # acceleration
        self.history = []
        self.Nhand = None
        self.Chand = None
        self.Ndist = None
        self.Cdist = None
        self.prevId = None
        self.id = None
        self.end = False
        self.recent = False
        self.original_id = original_id

    def __repr__(self):
        return "%r,%r,%r" % (self.vec[0],self.vec[1],self.vec[2])

    def __str__(self):
        return "%r,%r,%r" % (self.vec[0],self.vec[1],self.vec[2])

    def tuple(self):
        return self.x, self.y, self.z

    def isCa(self):
        return self.type=="CA"

    def setId(self,n):
        if self.original_id == None:
            self.original_id = n
        if self.id: self.prevId = self.id
        self.id = n

    def setNhand(self,N):
        self.Nhand = N
        self.Ndist = point_distance(self.vec,N.vec) if N else None

    def setChand(self,C):
        self.Chand = C
        self.Cdist = point_distance(self.vec,C.vec) if C else None


def chainDeepCopy(atom_list):
    """Deep copy of a Bead list"""
    new = []
    for a in atom_list:
        new.append(Bead(a.vec ,"CA" ,a.original_id))
    for i ,e in enumerate(new):
        e.setId(i)
        if i> 0:
            e.setNhand(new[i - 1])
        if i < len(new) - 1:
            e.setChand(new[i + 1])
    return new

class Line(object):
    def __init__(self, num, top):
        self.val = num
        self.top = top # boolean

    def topsign(self):
        return "-" if self.val%2==0 and self.top else ""

    def __repr__(self):
        return "{}{}".format(("-" if self.val%2==0 and self.top else "." if self.top else ""),self.val)

    def __gt__(self,l):
        if type(l) in [type(1),type(1.)]:
            return self.val
        else:
            return self.val > l.val

    def __lt__(self,l):
        if type(l) in [type(1),type(1.)]:
            return self.val < l
        else:
            return self.val < l.val

    def __eq__(self,l):
        if type(l) in [type(1),type(1.)]:
            return self.val == l
        else:
            return self.val == l.val

    def __add__(self,l):
        if type(l) in [type(1),type(1.)]:
            return self.val+l
        else:
            return self.val+l.val

    def __sub__(self,l):
        if type(l) in [type(1),type(1.)]:
            return self.val-l
        else:
            return self.val-l.val

    def __int__(self):
        return self.val

    def __mod__(self,l):
        return self.val%l

    def __abs__(self):
        return self.val

    def __trunc__(self):
        return self.val


class Crossing(object):
    def __init__(self,num1,top1,num2,top2):
        self.l1 = Line(num1,top1)
        self.l2 = Line(num2,top2)
        self._pair = (self.l1,self.l2)

    def reverse_topo(self):
        self.l1.top = not self.l1.top
        self.l2.top = not self.l2.top

    def __getitem__(self,key):
        return self._pair[key]

    def pair(self):
        return self._pair

    def min(self):
        return self._pair[0] if self._pair[0].val <self._pair[1].val else self._pair[1]

    def max(self):
        return self._pair[1] if self.l1.val <self.l2.val else self._pair[0]

    def even(self):
        return self._pair[1] if self.l1.val%2 else self._pair[0]

    def uneven(self):
        return self._pair[0] if self.l1.val%2 else self._pair[1]

    def reverse(self):
        self.l1.top = not self.l1.top
        self.l2.top = not self.l2.top

    def __repr__(self):
        return "({},{})".format(str(self.min()),str(self.max()))

    def remove_previous(self,*cnts):
        for c in sorted(cnts,reverse=True):
            if self.l1.val>c:
                self.l1.val -= 1
            if self.l2.val>c:
                self.l2.val -= 1

    def remove_previous_one(self,idx,*cnts):
        for c in sorted(cnts,reverse=True):
            if self._pair[idx].val>c:
                self._pair[idx].val -= 1

    def values(self):
        return [x.val for x in self._pair]

    def len_loop(self):
        return abs(self.l1.val-self.l2.val)

    def sum_loop(self):
        return abs(self.l1.val+self.l2.val)

    def to_code(self):
        return (self.uneven().val,self.even().val)

    def to_mod_code(self):
        return (self.uneven().val,self.even().val * (-1 if self.even().top else 1))

    def has_values(self,tup):
        return self.l1.val!=self.l2.val and self.l1.val in tup and self.l2.val in tup

    def __eq__(self, cr):
        return self.even().val == cr.even().val and self.uneven().val == cr.uneven().val and \
            self.even().top == cr.even().top and self.even().top == cr.even().top

class DowkerError(ValueError):
    """ Raised when Dowker code cannot be corrected """
    pass

class Code(object):
    def __init__(self):
        self.crossings = []
        self._mod_dowker_code = None
        self.history = set([])

    def __repr__(self):
        return ",".join(str(x) for x in self.crossings)

    def __iter__(self):
        return self.crossings.__iter__()

    def __len__(self):
        return len(self.crossings)

    def __getitem__(self,key):
        return self.crossings[key]

    def check_yo(self):
        if any(cr.l1.val%2==cr.l2.val%2 or cr.l1.top==cr.l2.top for cr in self.crossings):
            raise ValueError("Dowker broke: {}".format(self))

    def set_mod_dowker_code(self,lista):
        self._mod_dowker_code = lista

    def start_later_by(self,val):
        for c in self.crossings:
            c.l1.val += val
            c.l1.val %= (len(self.crossings)*2)
            c.l2.val += val
            c.l2.val %= (len(self.crossings)*2)
            if c.l1.val == 0: c.l1.val = (len(self.crossings)*2)
            if c.l2.val == 0: c.l2.val = (len(self.crossings)*2)

    def index(self,k):
        return self.crossings.index(k)

    def add(self,cr):
        self.crossings.append(cr)

    def make_way(self,i):
        for c in self.crossings:
            if c.l1.val >= i: c.l1.val += 1
            if c.l2.val >= i: c.l2.val += 1

    def pop(self,cr):
        try:
            id = cr if type(cr) == type(1) else self.crossings.index(cr)
            self.crossings.pop(id)
        except ValueError as e:
            if self.crossings:
                raise e
            else:
                pass


    def find_error(self, omc_inserted):
        smax = len(self.crossings)*2
        values = list(range(1,smax+2+1))
        for c1,c2 in self.crossings:
            _c1,_c2 = c1.val,c2.val
            if _c1>=omc_inserted: _c1+=1
            if _c2>=omc_inserted: _c2+=1
            if _c1%2!=_c2%2: #correct
                values = list(filter(lambda x: x <= _c1 or  x > _c2, values))
            else:
                values = list(filter(lambda x: _c1 < x <= _c2, values))
        values = list(values)
        if len(values) == 1:
                return values[0]
        else:
                raise DowkerError("cannot correct the code", self.crossings, values)

    def remove(self,*cr):
        vals = []
        indices = []
        for c in cr:
            id = c if type(c) == type(1) else self.crossings.index(c)
            vals += self.crossings[id].values()
            indices.append(id)
        for id in sorted(indices,reverse=True):
            self.crossings.pop(id)
        for cr in self.crossings:
            cr.remove_previous(*vals)

    def find_and_remove(self,*cr):
        vals = []
        for c in cr:
            for id,_ in enumerate(self.crossings):
                if set(c.values()) == set(_.values()):
                    break
            else:
                raise ValueError("{} is not in {}".format(c,self.crossings))
            vals += self.crossings[id].values()
            self.crossings.pop(id)
        for cr in self.crossings:
            cr.remove_previous(*vals)

    def gen_short_loop(self):
        def gen_short():
            for i in sorted(self.crossings, key=lambda x:x.len_loop()):
                yield i
        return gen_short()

    def read_in(self,code):
        for c in code:
            e = 1 if c[0]%2 else 0
            if c[e]>0:
                cr = Crossing(abs(c[0]),c[0]!=c[e],abs(c[1]),c[1]!=c[e])
            else:
                cr = Crossing(abs(c[0]),c[0]==c[e],abs(c[1]),c[1]==c[e])
            self.add(cr)

    def dowker_code(self): #absolute values are implied
        if self._mod_dowker_code is not None:
            return list(map(abs,self._mod_dowker_code))
        cd = sorted([x.to_code() for x in self.crossings])
        return [x[1] for x in cd]

    def mod_dowker_code(self):
        if self._mod_dowker_code is not None:
            return self._mod_dowker_code
        cd = sorted([x.to_mod_code() for x in self.crossings])
        return [x[1] for x in cd]

    def dowker_str(self):
        d = self.dowker_code()
        return "".join(map(str,d))

    def dowker_rstr(self):
        d = self.dowker_code()
        d.reverse()
        return "".join(map(str,d))

    def different_neighbour(self,for_what,diff_than):
        fw = for_what
        dt = diff_than
        if abs(fw[0]-dt[0]) == 1:
            val = fw[0] + (fw[0]-dt[0])
        elif abs(fw[1]-dt[0]) == 1:
            val = fw[1] + (fw[1]-dt[0])
        elif abs(fw[1] - dt[1]) == 1:
            val = fw[1] + (fw[1] - dt[1])
        elif abs(fw[0] - dt[1]) == 1:
            val = fw[0] + (fw[0] - dt[1])
        else:
            return None
        val = val%(len(self.crossings)*2)
        if not val: val = len(self.crossings)*2
        for c in self.crossings:
            if val in c:
                return c
        return None

    def pullulate(self):
        _copy = Code()
        for l1,l2 in self.crossings:
            _copy.add(Crossing(l1.val,l1.top,l2.val,l2.top))
        return _copy

class CodeHistory(object):
    def __init__(self):
        self.history = set([])
    def add(self,code):
        cr = code.dowker_str()# tuple(code.crossings)
        if cr not in self.history:
            self.history.add(cr)
            return True
        else:
            return False
