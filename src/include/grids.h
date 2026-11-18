#ifndef GRIDS_H
#define GRIDS_H
#include <complex>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>


class Grids
{
public:
    Grids();
    Grids(int Nr,double rmin,double rmax,int NY,double Ymin,double Ymax,int NT,double Tmin,double Tmax,double x0=0.01);
    ~Grids();

    int getNr(){return Nr_;}
    // r-points written to mv_dipoles.dat: the first Nr_out of the Nr evolved points (default all).
    // The remaining points are a buffer against the rmax edge (daughters beyond rmax are not integrated).
    void setNrOut(int n){ Nr_out_ = (n > 0 && n < Nr_) ? n : Nr_; }
    int getNrOut(){return Nr_out_;}
    double getrmaxOut(){return rmin_*std::pow(rmax_/rmin_, (Nr_out_-1.)/(Nr_-1.));}  // grid node, log spacing
    int getNY(){return NY_;}
    int getNT(){return NT_;}
    
    double getrmin(){return rmin_;}
    double getrmax(){return rmax_;}

    double getYmin(){return Ymin_;}
    double getYmax(){return Ymax_;}

    double getTmin(){return Tmin_;}
    double getTmax(){return Tmax_;}

    double getdT(){return dT_;}
    double getdY(){return dY_;}
    double getx0(){return x0_;}

    void setYskip(int yskip){YSkip_=yskip;}

    double T_at(int i){return Tmin_+dT_*i; }
    double Y_at(int i){return Ymin_+dY_*i; }
    

    std::string write_out_config();
private:
    // r 
    int Nr_;
    int Nr_out_;
    double rmin_;
    double rmax_; 
    // Y
    int NY_;
    double Ymin_;
    double Ymax_; 
    double dY_;
    double x0_;     // x at Y = 0 (Y = ln(x0/x)); BK does not use it, it labels the Y grid for McDIPPER
    int YSkip_;
    // T
    int NT_;
    double Tmin_;
    double Tmax_; 
    double dT_;
};

inline Grids::Grids(/* args */)
{
}

inline Grids::Grids(int Nr,double rmin,double rmax,int NY,double Ymin,double Ymax,int NT,double Tmin,double Tmax,double x0){
    Nr_=Nr;
    Nr_out_=Nr;
    rmin_=rmin;
    rmax_=rmax;
    // Y
    NY_=NY;
    Ymin_=Ymin;
    Ymax_=Ymax; 
    dY_ = (NY > 1) ? (Ymax-Ymin)/(NY-1.) : 0.0;   // NY=1: initial condition only
    x0_ = x0;
    // T
    NT_=NT;
    Tmin_=Tmin;
    Tmax_=Tmax; 
    dT_ = (NT > 1) ? (Tmax-Tmin)/(NT-1.) : 0.0;   // NT=1: single T = Tmin
    std::cerr<< "[ Grids ]: Grid created! " <<std::endl;
}

inline std::string Grids::write_out_config(){
    std::ostringstream config_str;
    config_str << "Grid:\n";
    config_str << "    r-grid:\n";
    config_str << "        Nr:   " << Nr_out_<<"\n";
    config_str << "        rmin: " << rmin_<<"# in fm \n";
    config_str << "        rmax: " << getrmaxOut()<<"\n";
    config_str << "        dr: Evolved in log-scale\n";
    config_str << "        Nr_evolved:   " << Nr_<<"   # buffer: evolved up to rmax_evolved, written up to rmax\n";
    config_str << "        rmax_evolved: " << rmax_<<"\n";
    config_str << "    Y-grid:\n";
    config_str << "        x0:   " << x0_<<"\n";
    config_str << "        NY:   " << NY_<<"\n";
    config_str << "        Ymin: " << Ymin_<<"\n";
    config_str << "        Ymax: " << Ymax_<<"\n";
    config_str << "        dY:   " << dY_<<"\n";
    config_str << "    T-grid:\n";
    config_str << "        NT:   " << NT_<<"\n";
    config_str << "        Tmin: " << Tmin_<<" # in fm^{-2}\n";
    config_str << "        Tmax: " << Tmax_<<" # in fm^{-2}\n";
    config_str << "        dT:   " << dT_<<"\n";
    return config_str.str();    
}


inline Grids::~Grids()
{
}




#endif


