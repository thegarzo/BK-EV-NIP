#ifndef GRIDS_H
#define GRIDS_H
#include <complex>
#include <cmath>


class Grids
{
public:
    Grids();
    Grids(int Nr,double rmin,double rmax,int NY,double Ymin,double Ymax,int NT,double Tmin,double Tmax);
    ~Grids();

    int getNr(){return Nr_;}
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

    void setYskip(int yskip){YSkip_=yskip;}

    double T_at(int i){return Tmin_+dT_*i; }
    double Y_at(int i){return Ymin_+dY_*i; }
    

    std::string write_out_config();
private:
    // r 
    int Nr_;
    double rmin_;
    double rmax_; 
    // Y
    int NY_;
    double Ymin_;
    double Ymax_; 
    double dY_;
    int YSkip_;
    // T
    int NT_;
    double Tmin_;
    double Tmax_; 
    double dT_;
};

Grids::Grids(/* args */)
{
}

Grids::Grids(int Nr,double rmin,double rmax,int NY,double Ymin,double Ymax,int NT,double Tmin,double Tmax){
    Nr_=Nr;
    rmin_=rmin;
    rmax_=rmax;
    // Y
    NY_=NY;
    Ymin_=Ymin;
    Ymax_=Ymax; 
    dY_ = (Ymax-Ymin)/(NY-1.);
    // T
    NT_=NT;
    Tmin_=Tmin;
    Tmax_=Tmax; 
    dT_ = (Tmax-Tmin)/(NT-1.);
    std::cerr<< "[ Grids ]: Grid created! " <<std::endl;
}

std::string Grids::write_out_config(){
    std::ostringstream config_str;
    config_str << "Grid:\n";
    config_str << "    r-grid:\n";
    config_str << "        Nr:   " << Nr_<<"\n";
    config_str << "        rmin: " << rmin_<<"# in fm \n";
    config_str << "        rmax: " << rmax_<<"\n";
    config_str << "        dr: Evolved in log-scale\n";
    config_str << "    Y-grid:\n";
    config_str << "        x0:   " << 0.01<<"\n";
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


Grids::~Grids()
{
}




#endif


