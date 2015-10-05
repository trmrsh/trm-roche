from libcpp.string cimport string
from libcpp cimport bool

from csubs cimport Vec3

cdef extern from "trm/roche.h" namespace "Roche":
    enum STAR:
        PRIMARY, SECONDARY
        
    void roche_face "Roche::face" \
        (double q, STAR star, double spin, \
        Vec3& dirn, double rref, double pref, double acc,\
        Vec3& pvec, Vec3& dvec, double& r, double& g) except+
        
    double roche_stradv "Roche::stradv" \
        (double q, Vec3&, Vec3&, double, double, double) except+
        
    void roche_strinit "Roche::strinit" (double q, Vec3& r, Vec3& v)
    
    void roche_ref_sphere "Roche::ref_sphere" \
        (double q, STAR star, double spin, double ffac,\
        double& rref, double& pref)
        
    Vec3 roche_set_earth "Roche::set_earth" (double iangle, double phase)
    
    Vec3 roche_set_earth "Roche::set_earth" \
        (double cosi, double sini, double phase)
        
    bool roche_fblink "Roche::fblink" \
        (double q, STAR star, double spin, double ffac, double acc,\
        const Vec3& earth, const Vec3& p)
        
    bool roche_ingress_egress "Roche::ingress_egress" \
        (double q, STAR star, double spin, double ffac, \
        double iangle, double delta, const Vec3& r,double& \
        ingress, double& egress)

    double roche_xl1 "Roche::xl1" (double q)
    double roche_xl2 "Roche::xl2" (double q)
    double roche_xl3 "Roche::xl3" (double q)
    double roche_xl11 "Roche::xl11" (double q)
    double roche_xl12 "Roche::xl12" (double q)
    
    double roche_rpot "Roche::rpot" (double q, Vec3& p)
    
    cdef cppclass Roche_Error:
        Roche_Error()
        Roche_Error(string)
        
