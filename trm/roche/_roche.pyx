#cython: embedsignature=True
# distutils: language = c++
cimport numpy as np
import numpy as np
from cython.operator cimport dereference as deref

from croche cimport *
from csubs cimport Vec3

class RocheError(Exception):
    pass

# helper functions
cdef Vec3 rconv(obj):
    """Create a Subs::Vec3 from a subs.Vec3"""
    cdef Vec3 p
    p.set(obj.x, obj.y, obj.z)
    return p

cdef STAR star_enum(star):
    cdef STAR rstar
    rstar = PRIMARY if star==1 else SECONDARY
    return rstar
    
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
    cdef Vec3 cdirn = rconv(dirn)
    cdef double r, g
    #rconv(dirn, &cdirn)
    
    rstar = star_enum(star)
    try:
        roche_face(q,rstar,spin,cdirn,rref,pref,acc,pvec,dvec,r,g)
    except Exception as err:
        raise RocheError(err)
    return ((pvec.x(),pvec.y(),pvec.z()), (dvec.x(),dvec.y(),dvec.z()), r, g)

def ref_sphere(q,spin,ffac,star=2):
    """
    (rref,pref) = ref_sphere(q, spin, ffac, star=2), returns reference radius and potential needed for face.
    
    q      -- mass ratio = M2/M1
    spin   -- ratio spin/orbital frequency
    ffac   -- linear filling factor of star in question, defined as the radius of the star along the line of
              centres towards its companion divided by the Roche lobe radius in that direction. For spin = 1
              the latter is simply the distance to the L1 point, but otherwise you need to use modified L1
              radii as returned by xl11 or xl12.
    star   -- 1 or 2 for primary or secondary star.
    """
    cdef double rref, pref
    rstar = star_enum(star)
    roche_ref_sphere(q,rstar,spin,ffac,rref,pref)
    return rref, pref

def findi(q, deltaphi, acc=1.e-4, di=1.0e-5):
    """
    findi(q, deltaphi, acc=1.e-4, di=1.e-5), computes inclination for a given mass ratio and phase width
    """
    assert q>0, "roche.findi: q must be > 0"
    assert deltaphi > 0 and deltaphi <= 0.25, "roche.findi: pwidth out of range 0 to 0.25"
    assert acc > 0 and acc <= 0.1, "roche.findi: acc <= 0 or acc > 0.1"
    assert di > 0 and di <= 10, "roche.findi: di <= 0 or di > 10."
    
    cdef double ilo=65., ihi=90.
    cdef double phi = deltaphi/2.
    cdef Vec3 earth1 = roche_set_earth(ilo,phi)
    cdef Vec3 earth2 = roche_set_earth(ihi,phi)
    cdef Vec3 r
    cdef bool elo = roche_fblink(q, SECONDARY, 1.0, 1.0, acc, earth1, r)
    cdef bool ehi = roche_fblink(q, SECONDARY, 1.0, 1.0, acc, earth2, r)
    
    cdef double iangle = 0
    if elo and ehi:
        iangle = -2
    elif (not elo) and (not ehi):
        iangle = -1
    else:
        while ihi - ilo > di:
            iangle = (ilo+ihi)/2.
            earth1 = roche_set_earth(iangle,phi)
            if (roche_fblink(q,SECONDARY,1.0,1.0,acc,earth1,r)):
                ihi = iangle
            else:
                ilo = iangle
        iangle = (ilo+ihi)/2
    return iangle

def findq(i, deltaphi, acc=1.e-4, dq=1.0e-5, qlo=0.001, qhi=2.):
    """
    findq(i, deltaphi, acc=1.e-4, dq=1.e-5, qlo=0.001, qhi=2.), computes mass ratio q for a given phase width and inclination
    """
    assert i>0 and i<=90, "roche.findq: iangle out of range 0 to 90"
    assert deltaphi > 0 and deltaphi <= 0.25, "roche.findq: pwidth out of range 0 to 0.25"
    assert acc > 0 and acc <= 0.1, "roche.findq: acc <= 0 or acc > 0.1"
    assert dq > 0 and dq <= 0.1, "roche.findq: dq <= 0 or di > 0.1"
    

    cdef double phi = deltaphi/2.
    cdef Vec3 r
    cdef Vec3 earth = roche_set_earth(i,phi)

    cdef bool elo = roche_fblink(qlo, SECONDARY, 1.0, 1.0, acc, earth, r)
    cdef bool ehi = roche_fblink(qhi, SECONDARY, 1.0, 1.0, acc, earth, r)
    
    cdef double q=0
    if elo and ehi:
        q = -2
    elif (not elo) and (not ehi):
        q = -1
    else:
        while qhi - qlo > dq:
            q = (qlo+qhi)/2.
            if (roche_fblink(q,SECONDARY,1.0,1.0,acc,earth,r)):
                qhi = q
            else:
                qlo = q
        q = (qlo+qhi)/2
    return q
 
def findphi(q,iangle,delta=1.0e-6):
    """
    findphi(q, i, delta=1.e-6), computes deltaphi for a given mass ratio and inclination
    """
    assert q>0, "roche.findphi: q <= 0"
    assert iangle>0 and iangle <= 90, "roche.findphi: iangle out of range 0 to 90"
    assert 0 < delta <= 0.001, "roche.findphi: delta <= 0 or delta > 0.001" 
    
    cdef Vec3 r = Vec3(0,0,0)
    cdef double ingress, egress
    cdef bool status = roche_ingress_egress(q,SECONDARY,1.0,1.0,iangle,delta,r,ingress,egress)
    if not status: 
        raise RocheError("roche.findphi: the centre of mass of the white dwarf is not eclipsed")
    return egress-ingress
          
def fblink(q,iangle,phi,r,ffac=1.,acc=1.0e-4,star=2,spin=1):
    """fblink(q, i, phi, r, ffac=1., acc=1.e-4, star=2, spin=1), computes whether a point is eclipsed or not
    
    Args:
        r (subs.Vec3): position to check
    """
    assert 0 < iangle <= 90, "roche.fblink: iangle out of range 0 to 90"
    assert 0 < ffac <= 1, "roche.fblink: ffac out of range 0 to 1"
    assert 0 < acc < 0.1, "roche.fblink: acc <= 0 or acc > 0.1"
    assert star==1 or star==2, "roche.fblink: star must be either 1 or 2"
    
    cdef STAR rstar = star_enum(star)
    cdef Vec3 earth = roche_set_earth(iangle,phi)
    cdef Vec3 rvec = rconv(r)
    if roche_fblink(q,rstar,spin,ffac,acc,earth,rvec):
        return True
    else:
        return False

def ineq(q,iangle,x,y,z=0,ffac=1.,delta=1.0e-7,star=2,spin=1):
    """ in,out = ineg(q, i, x, y, z=0, ffac=1., delta=1.e-7, star=2, spin=1), computes ingress and egress phase of a point"""
    assert 0 < iangle <= 90, "roche.ineg: iangle out of range 0 to 90"
    assert q>0, "roche.ineg: q must be > 0"
    assert 0 < ffac <= 1, "roche.ineg: ffac out of range 0 to 1"
    assert delta > 0, "roche.ineg: delta <= 0"
    assert star==1 or star==2, "roche.ineg: star must be either 1 or 2"
    
    cdef Vec3 r = Vec3(x,y,z)
    cdef double ingress, egress
    cdef STAR rstar = star_enum(star)
    cdef bool status = roche_ingress_egress(q,rstar,spin,ffac,iangle,delta,r,ingress,egress)
    if not status: 
        raise RocheError("roche.ineg: point is not eclipsed")
    return ingress, egress
    
def xl1(q):
    """Calculate the inner Lagrangian point distance.
    
    Args:
        q (float): mass ratio = M2/M1
    """
    return roche_xl1(q)