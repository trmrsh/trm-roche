from libcpp.string cimport string
from libcpp cimport bool

from csubs cimport Vec3

    
cdef extern from "trm/roche.h" namespace "Roche":
    enum STAR:
        PRIMARY, SECONDARY
    void roche_face "Roche::face" (double, STAR, double, Vec3&, double, double, double, Vec3&, Vec3&, double&, double&) except+
    double roche_stradv "Roche::stradv" (double q, Vec3&, Vec3&, double, double, double) except+
    void roche_strinit "Roche::strinit" (double, Vec3, Vec3)
    void roche_ref_sphere "Roche::ref_sphere" (double, STAR, double, double, double&, double&);
    Vec3 roche_set_earth "Roche::set_earth" (double, double)
    Vec3 roche_set_earth "Roche::set_earth" (double, double, double)
    bool roche_fblink "Roche::fblink" (double, STAR, double, double, double, const Vec3&, const Vec3&)
    bool roche_ingress_egress "Roche::ingress_egress" (double, STAR, double , double,\
                        double, double, const Vec3&,\
                        double&, double&)
    double roche_xl1 "Roche::xl1" (double)
    cdef cppclass Roche_Error:
        Roche_Error()
        Roche_Error(string)
        
