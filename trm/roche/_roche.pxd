from libcpp.string cimport string

cdef extern from "trm/vec3.h" namespace "Subs":
    cdef cppclass Vec3:
        Vec3()
        Vec3(double, double, double)
        Vec3(double*)
        double& x()
        double& y()
        double& z()
        void set(double, double, double)
        void set(double*)
        void get(double*)

    
cdef extern from "trm/roche.h" namespace "Roche":
    enum STAR:
        PRIMARY, SECONDARY
    void roche_face "Roche::face" (double, STAR, double, Vec3&, double, double, double, Vec3&, Vec3&, double&, double&) except+
    double roche_stradv "Roche::stradv" (double q, Vec3&, Vec3&, double, double, double) except+
    void roche_strinit "Roche::strinit" (double, Vec3, Vec3)
    double roche_xl1 "Roche::xl1" (double)
    cdef cppclass Roche_Error:
        Roche_Error()
        Roche_Error(string)
        
