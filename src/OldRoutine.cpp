// file: main.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>



// #include "include/grids.h"
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
    int NR = 120;
    double rmax=5.0 ; // In fm !!
    
    //Dyn values
    // double T0 = 4.0;// in fm^{-2}
    BKAux::fill_in_all();

    // Basic MV according to PRD 88, 114020 (2013)
    // PARS = BKAux::MVPars;
    // BKPars = BKAux::BK_MV;
    // MVgamma according to PRD 88, 114020 (2013)
    PARS = BKAux::MVgPars;
    BKPars = BKAux::BK_MVg;
    // MVe according to PRD 88, 114020 (2013)
    // PARS = BKAux::MVePars;
    // BKPars = BKAux::BK_MVe;

    BKPars.mp=0.2*GeV_to_fmm1;

    
    Dipole * SF = new Dipole();
    SF->init_MV(NR, rmax, rmin, PARS);

    std::ofstream dipole_B;
    std::ostringstream dipole_B_name;
    dipole_B_name << "fig3/S_MVg_0_Y_0_nBK_"<<BK::nGC<<"_BKm_"<<BKPars.mp*fmm1_to_GeV<<"GeV.txt";
    dipole_B.open (dipole_B_name.str());
    for (size_t ir = 0; ir < SF->getNr(); ir++){ dipole_B << SF->r_at(ir) << "\t" << SF->S_at(ir) << std::endl;}
    dipole_B.close();

    std::ofstream QS_file;
    std::ostringstream QS_file_name;
    QS_file_name << "fig3/S_MVg_Qs_nBK_"<<BK::nGC<<"_BKm_"<<BKPars.mp*fmm1_to_GeV<<"GeV.txt";
    QS_file.open (QS_file_name.str());
    
    Dipole * SF1 = new Dipole();
    SF1->init_MV(NR, rmax, rmin, PARS);
    double rS, Qs2;
    SF1->find_effective_saturation_scale(rS, Qs2);
    QS_file<< 0<< "\t" << 0<< "\t"<< rS<< "\t"<<  Qs2*fmm2_to_GeV2 <<std::endl;
    Dipole * SF2 = new Dipole();
    SF2->init_MV(NR, rmax, rmin, PARS);
    Dipole * SF3 = new Dipole();
    SF3->init_MV(NR, rmax, rmin, PARS);
    
    double Y=0;
    int NY= 140;
    std::cout<< "Ymax = "<< BK::dY*(NY-1.) << "_xmin_"<< 0.01 * exp(-BK::dY*(NY-1.)) << std::endl;
    for (size_t iy = 1; iy < NY; iy++)
    {

        Y=BK::getY(iy);
        BK::make_one_step( BKPars,SF,SF1,SF2,SF3);
        if (iy%10==0 || iy ==138 ){
            std::ofstream dipole;
            std::ostringstream dipole_name;
            dipole_name << "fig3/S_MVg_"<<iy<<"_Y_"<< Y <<"_nBK_"<<BK::nGC<<"_BKm_"<<BKPars.mp*fmm1_to_GeV<<"GeV.txt";
            dipole.open (dipole_name.str());
            for (size_t ir = 0; ir < SF->getNr(); ir++){ dipole << SF->r_at(ir) << "\t" << SF->S_at(ir) << std::endl;}
            dipole.close();
        }
        SF1->find_effective_saturation_scale(rS, Qs2);
        QS_file<< iy<< "\t" << Y<< "\t"<< rS<< "\t"<<  Qs2*fmm2_to_GeV2 << std::endl;
    }
    QS_file.close();
    
    SF->clear();
    SF1->clear();
    SF2->clear();
    SF3->clear();
    
   
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
    