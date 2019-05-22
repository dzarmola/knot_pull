from __future__ import print_function
from numpy import array as Vector
from .dowker_code_new import get_dt_code

def find_frame(frames):
    length = None
    cur = None
    for _,frame in enumerate(frames[::-1]):
        if len(frame) <= 3:
            return frame
        ldc = 0
        atoms = []
        for x in frame:
            atoms.append(Vector(x[1]))
            if x[2]:
                ldc += len(get_dt_code(atoms))
                atoms = []
        ldc += len(get_dt_code(atoms))
        if ldc < 3:
            return frame
        if length is not None and length < ldc:
            return cur
        length = ldc
        cur = frame
    return cur

