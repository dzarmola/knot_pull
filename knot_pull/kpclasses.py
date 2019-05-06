
from .vector_ops import point_distance

class VectorN(object):
    def __init__(self,*args):
        if len(args)==1:
            v = args[0]
            self.x = v.x
            self.y = v.y
            self.z = v.z
        else:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
    def __repr__(self):
        return "\t".join(self.x,self.y,self.z)

    def __str__(self):
        return "\t".join(self.x,self.y,self.z)



class Bead(object):
    def __init__(self, vec, type, original_id=None):
        self.x = vec[0]
        self.y = vec[1]
        self.z = vec[2]
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
        if any(cr.l1.val%2==cr.l2.val for cr in self.crossings):
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

    def pop(self,cr):
        id = cr if type(cr) == type(1) else self.crossings.index(cr)
        self.crossings.pop(id)

    def remove_old(self,cr):
        id = cr if type(cr) == type(1) else self.crossings.index(cr)
        vals = self.crossings[id].values()
        self.crossings.pop(id)
        for cr in self.crossings:
            cr.remove_previous(*vals)

    def remove(self,*cr):
        vals = []
        for c in cr:
            id = c if type(c) == type(1) else self.crossings.index(c)
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
        elif abs(fw[0] - dt[1]):
            val = fw[0] + (fw[0] - dt[1])
        else:
            return None
        val = val%(len(self.crossings)*2)
        for c in self.crossings:
            if val in c:
                return c
        return None

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