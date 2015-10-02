#cython: embedsignature=True
# distutils: language = c++
cimport numpy as np
import numpy as np
from cython.operator cimport dereference as deref

class RocheError(Exception):
    pass

# helper functions
cdef rconv(obj, Vec3* p):
    """Create a Subs::Vec3 from a subs.Vec3"""
    p.set(obj.x, obj.y, obj.z)



def bspot(q, rad, accel=1.0e-7):
    """returns position and stream velocity on stream at radius rad
    
    (x,y,vx,vy) = bspot(q, rad, acc=1.e-7)
    
    Args:
        q (float): mass ratio = M2/M1
        rad (float): radius to aim for
        acc (float[optional]): computationa accuracy parameter
    """
    cdef Vec3 r, v, rs, vs
    roche_strinit(q,r,v)
    try:
        roche_stradv(q,r,v,rad,accel,1.0e-2)
        return r.x(), r.y(), v.x(), v.y()
    except:
        raise RocheError("roche.stradv: never achieved desired radius")
        
 
def face(q,spin,dirn,rref,pref,star=2,acc=1.e-5):
    """
    p,d,r,g = face(q, spin, dirn, rref, pref, star=2, acc=1.e-5), returns position and direction of element of specific Roche potential.
    
     q      -- mass ratio = M2/M1
     spin   -- ratio spin/orbital frequency
     dirn   -- direction (a Vec3) to take from centre of mass of star in question.
     rref   -- reference radius greater than any radius of potential in question.
     pref   -- the potential to aim for.
     star   -- 1 or 2 for primary or secondary star. acc    -- accuracy in terms of separation of location.

    Returns p = position, d = direction perpendicular to face, r = radius from centre of mass, g = gravity.
    """    
    
    assert q>0, "roche.face: q <= 0"
    assert rref > 0, "roche.face: rref must be > 0"
    assert (acc > 0) or (acc < 0.1), "roche.face: acc <= 0 or acc > 0.1"
    assert star in [1,2], "roche.face: star must be either 1 or 2"
    
    # compute roche lobe
    cdef Vec3 pvec, dvec
    cdef Vec3 cdirn
    cdef double r, g
    rconv(dirn, &cdirn)
    
    cdef STAR rstar
    rstar = PRIMARY if star==1 else SECONDARY
    try:
        roche_face(q,rstar,spin,cdirn,rref,pref,acc,pvec,dvec,r,g)
    except Exception as err:
        raise RocheError(err)
    return ((pvec.x(),pvec.y(),pvec.z()), (dvec.x(),dvec.y(),dvec.z()), r, g)
    
def xl1(q):
    """Calculate the inner Lagrangian point distance.
    
    Args:
        q (float): mass ratio = M2/M1
    """
    return roche_xl1(q)