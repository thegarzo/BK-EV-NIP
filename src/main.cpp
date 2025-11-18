// file: main.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>



#include "include/grids.h"
#include "include/params.h"
#include "include/dipole.h"
#include "include/BK.h"



int main() {
    // simple example: r in [0.1, 1.0] (avoid r=0 if your data is singular) with 5 points

    // const double Qs = 1.0;         // saturation scale
    // const double r_min = 1e-6;     // smallest r
    // const double r_max = 10.0;     // largest r
    // const int n_points = 300;      // number of grid points

    
    MV_parameters PARS;
    BKParameters BKPars;
    

    // Lattice
    double rmin = 0.001; // in fm 
    double rmax = 10.0 ; // In fm !!
    int NR = 141;
  
    
    //Dyn values
    // double T0 = 4.0;// in fm^{-2}
    BKAux::fill_in_all();

    // Basic MV according to PRD 88, 114020 (2013)
    PARS = BKAux::MVPars;
    BKPars = BKAux::BK_MV;
    // MVgamma according to PRD 88, 114020 (2013)
    // PARS = BKAux::MVgPars;
    // BKPars = BKAux::BK_MVg;
    // MVe according to PRD 88, 114020 (2013)
    // PARS = BKAux::MVePars;
    // BKPars = BKAux::BK_MVe;
    std::cerr<< "[ main ]: Parameters Set! " <<std::endl;

    double Ymin=0;double Ymax=12;
    int NY=121;

    double Tmin=7.1;double Tmax=10.0;
    int NT=30;

    Grids GRID( NR, rmin,rmax, NY,Ymin, Ymax,NT, Tmin, Tmax); 
    std::cerr<< "[main]: Parameters Set! " <<std::endl;
    BK::RunModel("MV_new", GRID, PARS, BKPars );
   
    return 0;
}



    // std::vector<double> r(n_points);
    // std::vector<double> S(n_points);

    // // log-spaced grid
    // const double log_rmin = std::log10(r_min);
    // const double log_rmax = std::log10(r_max);
    // for (int i = 0; i < n_points; ++i) {
    //     double log_r = log_rmin + i * (log_rmax - log_rmin) / (n_points - 1);
    //     r[i] = std::pow(10.0, log_r);
    // }

    // // compute S(r)
    // for (int i = 0; i < n_points; ++i) {
    //     double rq = r[i] * Qs;
    //     S[i] = sin(3*rq)*exp(-0.5*rq)+0.5;
    // }

    // // write to file
    // std::ofstream out("dipole_test_hybrid.txt");
    // if (!out.is_open()) {
    //     std::cerr << "Error: cannot open output file.\n";
    //     return 1;
    // }

    // out << "# Test dipole S(r) = exp(-(r*Qs)^1.2 * log(1 + 1/(r*Qs)))\n";
    // out << "# Qs = " << Qs << "\n";
    // out << "# Columns: r    S\n";
    // out << std::scientific << std::setprecision(12);

    // for (int i = 0; i < n_points; ++i) {
    //     out << r[i] << "    " << S[i] << "\n";
    // }

    // out.close();
    // std::cout << "Wrote dipole_test_hybrid.txt with " << n_points << " points.\n";

    // Dipole dip;
    // dip.init_from_file("dipole_test_hybrid.txt", /*interp_in_log=*/true);

    // double r0= 1e-6;
    // double rM =1e1;
    // int Nr=1000;
    // double dR = std::log(rM/r0)/(Nr-1.);

    // for (size_t i = 0; i < Nr; i++)
    // {

    //     double r = r0*exp(i*dR);
    //     std::cout<<  r << "\t" << dip.eval(r)<<std::endl;
    // }
    


    // Build it 
    
    // std::vector<double> r_tV; 
    // std::vector<double> q_tV;
    // int NR2 =160; 
    // double rmin2 = 0.001; double rmax2 = 20;
    

    // // dip.create_r_grid(NR2, double r_max, double r_min,std::vector<double>& r_grid,std::vector<double>& q_grid)

    // // BK::Test_Kernels( NR2,rmax2,rmin2, 0.1, 0.0,BKPars);
    // // BK::Test_Kernels( NR2,rmax2,rmin2, 0.1, M_PI/4.,BKPars);
    // // BK::Test_Kernels( NR2,rmax2,rmin2, 0.1, M_PI/2.,BKPars);
    // // BK::Test_Kernels( NR2,rmax2,rmin2, 0.1, M_PI,BKPars);
    
    // double res, err;
    // for (size_t i = 0; i < 100; i++)
    // {
    //     double u = i*log(rmax/rmin)/(100-1.);
    //     BKPars.r = rmin*exp(u);
    //     BK::RHS_integral (BKPars,res, err);
    //     std::cout<< BKPars.r << "\t"<< res<< "\t"<< err << std::endl; 
    // }
    