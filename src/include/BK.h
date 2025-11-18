
#ifndef BK_H
#define BK_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

#include "params.h"
#include "output.h"

enum class KernelType : int { LO = 0,  rcParent =1,  rcDaughter =2,  rcBalitsky =3,  rcSmallest =4};


struct BKParameters
{
    // Kernel values
    double Lambda;
    double mu2;
    double C2;
    double zeta;
    double mp;
    KernelType KType;
    //Dipole values
    double r;
    Dipole * Dip;
};

namespace BKAux{
    MV_parameters MVPars;
    MV_parameters MVgPars;
    MV_parameters MVePars;
    BKParameters BK_MV;
    BKParameters BK_MVg;
    BKParameters BK_MVe;

    // These functions fill in the parameters obtained in PHYSICAL REVIEW D 88, 114020 (2013)
    void fill_in_MV(){
        // BASIC MV:  MV_parameters= { Qs02;gamma; Lambda;T;sigma0;ec;}
        MVPars.Qs02 = 0.104 *GeV2_to_fmm2 ;
        MVPars.gamma = 1.0;
        MVPars.Lambda = LQCD * GeV_to_fmm1 ;
        MVPars.sigma0 = 18.81 *mb_to_fm2 ;
        MVPars.T = 1./MVPars.sigma0 ; 
        MVPars.ec = 1;
        /// ALL MUST BE IN FM! 
        // BK Parameters= {}
        // Kernel values
        BK_MV.Lambda = LQCD * GeV_to_fmm1;
        BK_MV.mu2 = pow(0.28*GeV_to_fmm1,2.0); 
        BK_MV.C2 = 14.5;
        BK_MV.zeta= 0.2;
        BK_MV.mp = 0.4*GeV_to_fmm1; 
        BK_MV.KType = KernelType::rcDaughter;
        //Dipole values
        BK_MV.r = 0.2;
    };
    void fill_in_MVe(){
        // MVe:  MV_parameters= { Qs02;gamma; Lambda;T;sigma0;ec;}
        MVePars.Qs02 = 0.060 *GeV2_to_fmm2 ;
        MVePars.gamma = 1.0;
        MVePars.Lambda = LQCD * GeV_to_fmm1 ;
        MVePars.sigma0 = 16.36 *mb_to_fm2 ;
        MVePars.T = 1./MVePars.sigma0;
        MVePars.ec = 18.9;  /// ALL MUST BE IN FM! 
        // BK Parameters= {}
        // Kernel values
        BK_MVe.Lambda = LQCD * GeV_to_fmm1;
        BK_MVe.mu2 = pow(0.28*GeV_to_fmm1,2.0); 
        BK_MVe.C2 = 7.2;
        BK_MVe.zeta= 0.2;
        BK_MVe.mp = 0.4*GeV_to_fmm1; 
        BK_MVe.KType = KernelType::rcDaughter;
        //Dipole values
        BK_MVe.r = 0.2;
    };
    void fill_in_MVg(){
        // MVe:  MV_parameters= { Qs02;gamma; Lambda;T;sigma0;ec;}
        MVgPars.Qs02 = 0.165 *GeV2_to_fmm2 ;
        MVgPars.gamma = 1.135;
        MVgPars.Lambda = LQCD * GeV_to_fmm1 ;
        MVgPars.sigma0 = 16.45 *mb_to_fm2 ;
        MVgPars.T =1./MVgPars.sigma0 ;
        MVgPars.ec = 1;
        /// ALL MUST BE IN FM! 
        //g BK Parameters= {}
        // Kernel values
        BK_MVg.Lambda = LQCD * GeV_to_fmm1;
        BK_MVg.mu2 = pow(0.28*GeV_to_fmm1,2.0); 
        BK_MVg.C2 = 6.35;
        BK_MVg.zeta= 0.2;
        BK_MVg.mp = 0.4*GeV_to_fmm1; 
        BK_MVg.KType = KernelType::rcDaughter;
        //Dipole values
        BK_MVg.r = 0.2;
    };

    void fill_in_all(){
        fill_in_MV();
        fill_in_MVe();
        fill_in_MVg();
    };
}

namespace BK{
    double beta0 = (11*NC-2*NF)/3.0;
    double pref = NC/(2*M_PI*M_PI);
    int nGC= 40;
    double wGC= 2*M_PI/nGC;
    double dY=0.1;

    double BKEpsAbs=0.0;
    double BKEpsRel=1e-4;
    double BKLimit=1e6;

    double NaNLimit= 1e-25;
    double BKClipTol=1e-25;

    bool no_handlebars=true;

    double getY(int iY){return iY*dY;}

    double  alphaS(double r, BKParameters &BKP){
        double IR = BKP.mu2*BKP.C2 *pow(BKP.Lambda,-2.);
        double UV = 4*BKP.C2*pow(r*BKP.Lambda, -2.0);
        return (4*M_PI) / (beta0*BKP.zeta* log (  pow(IR, 1.0/BKP.zeta)  + pow(UV, 1.0/BKP.zeta) ) )  ;
    }

    // Here come the kernels 
    double KLO(double r, double rp, double phi, BKParameters &BKP){
        //Normal Leading Order Kernel, AlphaS is fixed 
        double xdotz=r*rp*cos(phi) ;
        double x2 =sqrt(fabs(r*r + rp*rp -2*xdotz));
        double kernel_sum = pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp),2.0);
        kernel_sum +=  pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2),2.0);
        kernel_sum +=  2*( BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp)  ) * (BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2)) * (xdotz-rp*rp)/(rp*x2); 
        return pref*alpha_S*kernel_sum;
    }


    double rcKDaughter(double r, double rp, double phi, BKParameters &BKP){
        //rcBK Leading Order Kernel, daughter dipoles set the scale 
        double xdotz=r*rp*cos(phi) ;
        double x2 =sqrt(fabs(r*r + rp*rp -2*xdotz));
        double aSz=alphaS(rp,BKP);
        double aSx2=alphaS(x2,BKP);
        double kernel_sum =  aSz*pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp),2.0);
        kernel_sum += aSx2 * pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2),2.0);
        kernel_sum += 2*sqrt( aSx2 * aSz) * ( BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp)  ) * (BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2)) * (xdotz-rp*rp)/(rp*x2); 
        return pref*kernel_sum; 
    }

    double rcKParent(double r, double rp, double phi, BKParameters &BKP){
        double xdotz=r*rp*cos(phi) ;
        double x2 =sqrt(fabs(r*r + rp*rp -2*xdotz));
        double kernel_sum = pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp),2.0);
        kernel_sum +=  pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2),2.0);
        kernel_sum +=  2*( BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp)  ) * (BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2)) * (xdotz-rp*rp)/(rp*x2); 
        return pref*alphaS(r,BKP)*kernel_sum; 
    }

    double rcKSmallest(double r, double rp, double phi, BKParameters &BKP){
        double xdotz=r*rp*cos(phi) ;
        double x2 =sqrt(fabs(r*r + rp*rp -2*xdotz));
        double alphaS_smallest=alphaS(std::min({r, rp, x2}),BKP);
        double kernel_sum = pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp),2.0);
        kernel_sum +=  pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2),2.0);
        kernel_sum +=  2*( BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp)  ) * (BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2)) * (xdotz-rp*rp)/(rp*x2); 
        return pref*alphaS_smallest * kernel_sum; 
    }

    double rcKBalitsky(double r, double rp, double phi, BKParameters &BKP){
        double xdotz=r*rp*cos(phi) ;
        double x2 =sqrt(fabs(r*r + rp*rp -2*xdotz));
        double aSx = alphaS(r,BKP);
        double aSz=alphaS(rp,BKP);
        double aSx2=alphaS(x2,BKP);
        double kernel_sum = pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *rp),2.0) *(aSz/aSx2 -1.0);
        kernel_sum +=  pow(BKP.mp * gsl_sf_bessel_K1(BKP.mp *x2),2.0) * (aSx2/aSz -1.0);
   
        return rcKParent(r, rp, phi, BKP) + pref * aSx * kernel_sum; 
    }
    
    double EvaluateKernel(double r, double rp, double phi, BKParameters &BKP){
        double Kernel;
        switch (BKP.KType)
        {
        case KernelType :: LO : 
            Kernel= KLO(r, rp, phi,BKP);
            break;
        case KernelType ::rcParent : 
            Kernel= rcKParent(r, rp, phi,BKP);
            break;
        case KernelType ::rcDaughter: 
            Kernel= rcKDaughter(r, rp, phi,BKP);
            break;
        case KernelType ::rcSmallest: 
            Kernel= rcKSmallest(r, rp, phi,BKP);
            break;
        case KernelType ::rcBalitsky: 
            Kernel= rcKBalitsky(r, rp, phi,BKP);
            break;
        default:
            break;
        }
        return Kernel;
    }

    void Test_Kernels(int Nr_t, double r_max_t, double r_min_t, double rp, double phi, BKParameters &BKP)
    {
        std::vector<double> r_tmp(Nr_t);
        double dq= std::log(r_max_t/r_min_t)/(Nr_t-1.);
        for (int i = 0; i < Nr_t; ++i) {
            r_tmp[i] = r_min_t*std::exp(dq*i);
        }
        std::ofstream dipole_B;
        std::ostringstream dipole_B_name;
        dipole_B_name << "test_kernels_rp_"<< rp << "_phi_"<< phi <<"_Lambda_"<<BKP.Lambda <<"_mu2_"<< BKP.mu2<<"_C2_"<<BKP.C2<<"_zeta_"<<BKP.zeta<<"_mp_"<<BKP.mp <<".txt";
        dipole_B.open (dipole_B_name.str());
        for (size_t i = 0; i < Nr_t; i++){ 
            BKP.KType = KernelType::LO;
            dipole_B << r_tmp[i] << "\t" << EvaluateKernel(r_tmp[i], rp, phi,BKP);
            BKP.KType = KernelType::rcParent;
            dipole_B << "\t" << EvaluateKernel(r_tmp[i], rp, phi,BKP);
            BKP.KType = KernelType::rcDaughter;
            dipole_B << "\t" << EvaluateKernel(r_tmp[i], rp, phi,BKP);
            BKP.KType = KernelType::rcSmallest;
            dipole_B << "\t" << EvaluateKernel(r_tmp[i], rp, phi,BKP);
            BKP.KType = KernelType::rcBalitsky;
            dipole_B << "\t" << EvaluateKernel(r_tmp[i], rp, phi,BKP)<< std::endl;  
        }
        dipole_B.close();
    }
    
    double RHS_integrand (double u, void * params) {
        BKParameters pars = *(BKParameters *) params;
        // double x = exp(-pars.Y);
        double Irrp=0.0;
        double rp =pars.Dip->r_min() *exp(u);

        double phi_i, Wi, Ki,x2_i;
        // Here comes the Gauss-Chebyshev quadrature
        for (size_t i = 1; i <= nGC; i++)
        {
            // get the angle 
            phi_i = (2.0*i-1.0)*M_PI/(2.0*nGC);
            x2_i = sqrt( fabs( pow(pars.r ,2. ) + pow(rp,2.) - 2*rp* pars.r *cos(phi_i) ) );
            // get the kernel 
            Ki = EvaluateKernel(pars.r, rp, phi_i, pars);   
            // get the Dipoles 
            Wi=( pars.Dip->eval(rp))*(pars.Dip->eval(x2_i) ) - pars.Dip->eval(pars.r);
            // Mix together
            Irrp+= Ki*Wi;
        }
        return  wGC*rp*rp*Irrp;
        /////// RECHECK THE DEFINITION WITH Y/YBAR!
    }

    void RHS_integral (BKParameters &bkpars ,double &result,double &error ) {
            gsl_function F;
            F.function = &RHS_integrand;
            F.params = &bkpars;

            double uMax= log(bkpars.Dip->r_max()/bkpars.Dip->r_min() );

            gsl_integration_workspace * w= gsl_integration_workspace_alloc (BKLimit);
            int status = gsl_integration_qags (&F, 0.0, uMax, BKEpsAbs,BKEpsRel, BKLimit, w, &result, &error); 
            if (status) {
                status = gsl_integration_qags (&F, 0.0,  uMax,BKEpsAbs,1e-3, BKLimit ,w, &result, &error);
                if (status) { 
                    status = gsl_integration_qags (&F, 0.0,  uMax, BKEpsAbs,1e-2, BKLimit ,w, &result, &error);
                    if (status) {
                        std::cerr << "[BK-ERROR]: integration failure! status = " << gsl_strerror (status )<< ". " << result << "+_" << error << std::endl; exit(EXIT_FAILURE);
                    }else {
                        std::cerr << "[BK-ERROR]: reduced relative and absolute error (1e-2) integration. Result = " << result << "+-" << error <<std::endl;
                    }
                } else {
                    std::cerr << "\n[BK-ERROR]: reduced relative error (1e-3) integration. Result = " << result << "+-" << error <<std::endl; }
            }
            gsl_integration_workspace_free (w);

   // Check status of integration and reduce accuracy if failing
            if (result!=result){
                std::cerr  << "[BK-ERROR]: Found NaN at " << bkpars.r  <<"\t" << std::endl; exit(0);  
            }

    }

    void clip_evolution(std::vector<double> &arr , int Narr, double clipTol) {
        for (std::size_t i = 0; i < Narr; ++i) {
            if (arr[i] < clipTol)  arr[i] = 0;
        }
    }

    void make_one_step(BKParameters BKProxy, Dipole * S, Dipole * S1, Dipole * S2, Dipole * S3){
        std::vector<double> k1(S->getNr());
        std::vector<double> k2(S->getNr());
        std::vector<double> k3(S->getNr());
        std::vector<double> k4(S->getNr());

        BKParameters BKPars={BKProxy.Lambda,BKProxy.mu2,BKProxy.C2,BKProxy.zeta,BKProxy.mp,BKProxy.KType,BKProxy.r,S};
        BKParameters BKPars1={BKProxy.Lambda,BKProxy.mu2,BKProxy.C2,BKProxy.zeta,BKProxy.mp,BKProxy.KType,BKProxy.r,S1};
        BKParameters BKPars2={BKProxy.Lambda,BKProxy.mu2,BKProxy.C2,BKProxy.zeta,BKProxy.mp,BKProxy.KType,BKProxy.r,S2};
        BKParameters BKPars3={BKProxy.Lambda,BKProxy.mu2,BKProxy.C2,BKProxy.zeta,BKProxy.mp,BKProxy.KType,BKProxy.r,S3};

        std::vector<double> S_proxy (S->getNr());
        double res,err;
    

        for (size_t ir= 0; ir < S->getNr(); ir++){ 
            // std::cerr<< "In k1 "<< i <<std::endl;
            BKPars.r=BKPars.Dip->r_at(ir);
            RHS_integral (BKPars,res,err );
            k1[ir]=res;
            S_proxy[ir] = BKPars.Dip->S_at(ir) + 0.5* dY *k1[ir];
        }
        BKPars1.Dip->evolve_dipole(S_proxy);

        for (size_t ir= 0; ir < S->getNr(); ir++){ 
            // std::cerr<< "In k1 "<< i <<std::endl;
            BKPars1.r=BKPars1.Dip->r_at(ir);
            RHS_integral (BKPars1,res,err );
            k2[ir]=res;
            S_proxy[ir] = BKPars.Dip->S_at(ir) + 0.5* dY *k2[ir];
        }
        BKPars2.Dip->evolve_dipole(S_proxy);

        for (size_t ir= 0; ir < S->getNr(); ir++){ 
            // std::cerr<< "In k1 "<< i <<std::endl;
            BKPars2.r=BKPars2.Dip->r_at(ir);
            RHS_integral (BKPars2,res,err );
            k3[ir]=res;
            S_proxy[ir] = BKPars.Dip->S_at(ir) + dY *k3[ir];
        }
        BKPars3.Dip->evolve_dipole(S_proxy);

        for (size_t ir= 0; ir < S->getNr(); ir++){ 
            // std::cerr<< "In k1 "<< i <<std::endl;
            BKPars3.r=BKPars3.Dip->r_at(ir);
            RHS_integral (BKPars3,res,err );
            k4[ir]=res;
        }

        for (size_t ir= 0; ir < S->getNr(); ir++){ 
            // BKPars.r=BKPars.Dip->r_at(ir);
            S_proxy[ir]= BKPars.Dip->S_at(ir) + dY *(k1[ir] + 2*k2[ir] + 2*k3[ir] + k4[ir] )/6.;
            // std::cout << Y_t<< "\t" <<  glue->xgx(Y_t) << "\t" << S_proxy[i] << std::endl;
        }
        
        clip_evolution(S_proxy,S->getNr(), BKClipTol);
        // std::cout<< "Now it is mu2 = " << glue->get_mu2()<<std::endl;
        BKPars.Dip->evolve_dipole(S_proxy);
    }


bool findFolder(std::string& PATH) {

    std::filesystem::path dir(PATH);

    while (true) {
        if (std::filesystem::exists(dir)) {
            std::cout << "Folder \"" << dir.string() << "\" already exists. Overwrite? (y/n): ";
            char choice;
            std::cin >> choice;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // clear input buffer

            if (choice == 'y' || choice == 'Y') {
                std::filesystem::remove_all(dir);
                if (std::filesystem::create_directory(dir)) {
                    std::cout << "Folder overwritten successfully.\n";
                    return true;
                } else {
                    std::cerr << "Error overwriting folder.\n";
                    return false;
                }
            } 
            else if (choice == 'n' || choice == 'N') {
                std::cout << "Please enter an alternative folder path: ";
                std::getline(std::cin, PATH);
                dir = std::filesystem::path(PATH);  // update path reference and recheck
                continue;
            } 
            else {
                std::cout << "Invalid input. Please enter 'y' or 'n'.\n";
                continue;
            }
        } 
        else {
            if (std::filesystem::create_directories(dir)) {
                std::cout << "Folder created successfully at \"" << dir.string() << "\".\n";
                PATH = dir.string();  // ensure PATH reflects the final valid folder
                return true;
            } else {
                std::cerr << "Error creating folder.\n";
                return false;
            }
        }
    }
}



void RunModel(std::string PATH, Grids grid, MV_parameters MVp, BKParameters BKp ){

    if(no_handlebars){gsl_set_error_handler_off();}
    dY=grid.getdY();
    //Here we write out the config file
    findFolder(PATH);
    std::ofstream config_f;
    std::ostringstream config_f_name;
    config_f_name << PATH<< "/config.yaml";
    config_f.open (config_f_name.str());
    config_f << "General:\n";
    config_f << "    Nc:        " << NC<< "\n";
    config_f << "    Nf:        " << NF<< "\n";
    config_f << "    LambdaQCD: " << LQCD << " # GeV \n"; // GeV
    config_f << std::endl;
    config_f << "MV-Parameters:\n";
    config_f << "    Qs02:   " << MVp.Qs02 * fmm2_to_GeV2 << " # GeV^2\n"; // GeV
    config_f << "    gamma:  " << MVp.gamma<< "\n";
    config_f << "    Lambda: " << MVp.Lambda * fmm1_to_GeV<< " # GeV \n";// GeV
    config_f << "    sigma0: " << MVp.sigma0 << "\n";
    config_f << "    ec:     " << MVp.ec<< "\n";
    config_f << std::endl;
    config_f << "BK-Parameters:\n";
    config_f << "    Lambda: " << BKp.Lambda * fmm1_to_GeV<<" # GeV \n";
    config_f << "    mu2:    " << BKp.mu2* fmm2_to_GeV2<<" # GeV^2\n";
    config_f << "    C2:     " << BKp.C2<<"\n";
    config_f << "    zeta:   " << BKp.zeta<<"\n";
    config_f << "    mp:     " << BKp.mp* fmm1_to_GeV <<" # GeV \n";
    config_f << "    dY:     " << dY << "\n";
    config_f << "    nGC:    " << nGC << "\n";
    switch (BKp.KType)
    {
    case KernelType::LO :
        config_f << "    Ktype:  LO \n";break;
    case KernelType::rcDaughter :
        config_f << "    Ktype:  Daughter \n";break;
    case KernelType::rcParent :
        config_f << "    Ktype:  Parent \n";break;
    case KernelType::rcSmallest :
        config_f << "    Ktype:  Smallest \n";break;
    case KernelType::rcBalitsky :
        config_f << "    Ktype:  Balitsky \n";break;
    default:
        std::cerr<< "[ERROR] KType not found!";
        exit(0);
        break;
    }
    config_f << std::endl;
    config_f << grid.write_out_config();
    config_f << std::endl;
    config_f.close();
    std::cerr<< "[BK]: Config File written out! " <<std::endl;
    
    // Ok.... let's go!
    std::ofstream qs2_f;
    std::ostringstream qs2_f_name;
    qs2_f_name << PATH<< "/mv_Q2.dat";
    qs2_f.open (qs2_f_name.str());
    
    std::ofstream dipole_f;
    std::ostringstream dipole_f_name;
    dipole_f_name << PATH<< "/mv_dipoles.dat";
    dipole_f.open (dipole_f_name.str());

    Dipole * SF = new Dipole();
    Dipole * SF1 = new Dipole();
    Dipole * SF2 = new Dipole();
    Dipole * SF3 = new Dipole();
    Dipole * SA = new Dipole();

    for (size_t iT = 0; iT < grid.getNT(); iT++)
    {
        // Create the dipoles in fundamental and adjoint reps + the proxies. 
        double Ti = grid. T_at(iT);
        if(Ti>0){

            MVp.T=Ti;

            SF->init_MV(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            SF1->init_MV(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            SF2->init_MV(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            SF3->init_MV(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            
            // Adjoint 
            SA->init_MVA(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            // Dipole * SA1 = new Dipole(); SA1->init_MVA(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            // Dipole * SA2 = new Dipole(); SA2->init_MVA(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            // Dipole * SA3 = new Dipole(); SA3->init_MVA(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);

            for (size_t iY = 0; iY < grid.getNY(); iY++)
            {
                /* Here goes the evolution with respect to Y */
                double Y=BK::getY(iY);
                if(iY>0) BK::make_one_step(BKp,SF,SF1,SF2,SF3);
                // if(iY>0) BK::make_one_step(BKp,SA,SA1,SA2,SA3);
                /* After step goes the extraction of Qs */
                double rSA, Qs2A;
                double rSF, Qs2F;
                std::vector<double> SA_proxy (SA->getNr());
                for (size_t ir = 0; ir < SA->getNr(); ir++){SA_proxy[ir] = pow(SF->S_at(ir),CA/CF);}
                SA->evolve_dipole(SA_proxy);
            
                
                SA->find_effective_saturation_scale(rSA, Qs2A , std::exp(-0.5),1e-5, 10000);
                SF->find_effective_saturation_scale(rSF, Qs2F , std::exp(-0.5),1e-5, 10000);
                                                   
                qs2_f << Ti << "\t"<< Y << "\t" << Qs2F*fmm2_to_GeV2 << "\t" << Qs2A*fmm2_to_GeV2  << std::endl;

                /* And output to general array */
                for (size_t ir = 0; ir < SF->getNr(); ir++){ 
                    dipole_f << Ti << "\t"<< Y << "\t" <<  SF->r_at(ir) << "\t" << SF->S_at(ir)<< "\t" << SA->S_at(ir) <<std::endl ;
                }
                // Now some output
                if (iY%2 ==0){   
                    double percentage_done1 = double(iT)/double(grid.getNT());
                    double percentage_done2 = double(iY)/double(grid.getNY());
                    Output::printProgress2(percentage_done1,percentage_done2);
                }
                
            }
        }
        else{
            Dipole * SF = new Dipole(); SF->init_MV(grid.getNr(), grid.getrmax(), grid.getrmin(), MVp);
            for (size_t iY = 0; iY < grid.getNY(); iY++)
            {
                double Y=BK::getY(iY);
                for (size_t ir = 0; ir < SF->getNr(); ir++){ 
                    dipole_f << 0.0 << "\t"<< Y << "\t" <<  SF->r_at(ir) << "\t" << 0.0<< "\t" << 0.0 <<std::endl ;
                }
                qs2_f << 0.0 << "\t"<< Y << "\t" << 0.0 << "\t" << 0.0  << std::endl;
            }
        }
        
    }
    qs2_f.close();
    dipole_f.close();

    SF->clear();
    SF1->clear();
    SF2->clear();
    SF3->clear();

    SA->clear();
    // SA1->clear();
    // SA2->clear();
    // SA3->clear();
    
    }
}

#endif // MV_HPP
