
#ifndef DIPOLE_H
#define DIPOLE_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "params.h"
struct MV_parameters
{ 
    double Qs02;
    double gamma;
    double Lambda;
    double T;
    double sigma0;
    double ec;
};

class Dipole {
public:
    Dipole() noexcept;
    // disallow copy to avoid double-free of gsl resources; allow move
    Dipole(const Dipole&) = delete;
    Dipole& operator=(const Dipole&) = delete;
    Dipole(Dipole&& other) noexcept;
    Dipole& operator=(Dipole&& other) noexcept;
    ~Dipole() noexcept;

    // initialize from std::vector (r and S must have same size >=2)
    // interp_in_log: if true, will build interpolator in log(r) space
    void init(const std::vector<double>& r_grid, const std::vector<double>& S_grid,
              bool interp_in_log = false);

    // initialize from raw arrays with size n
    void init_from_arrays(const double* r_grid, const double* S_grid, std::size_t n,
                          bool interp_in_log = false);

    // initialize from file with two columns: r  S (optional log-space)
    void init_from_file(const std::string& filename, bool interp_in_log = false);
    
    void create_r_grid(int Nr_t, double r_max, double r_min,std::vector<double>& r_grid,std::vector<double>& q_grid);
    void init_MV(int Nr_t, double r_max, double r_min, MV_parameters &params);
    void init_MVA(int Nr_t, double r_max, double r_min, MV_parameters &params);
    void evolve_dipole(const std::vector<double>& S_grid);
    // access S by index (no bounds check in release, add check if you want)
    double operator[](std::size_t i) const noexcept;

    // return r at index
    double r_at(std::size_t i) const noexcept;
    double S_at(std::size_t i) const noexcept;

    // number of grid points
    std::size_t size() const noexcept;

    // evaluate interpolator at arbitrary r_query.
    // If r_query is outside the stored r-range, we clamp to the endpoint S value.
    double eval(double r_query)const;
    bool find_effective_saturation_scale(double& r_sat_out,
                                     double& Qs2_out,
                                     double target = std::exp(-0.5),
                                     double rel_tol = 1e-6,
                                     int    max_iter = 1000) const;


    // Rebuild interpolator from current r_vec and S_vec (useful after modifying arrays)
    void rebuildInterpolator();

    // Clear internal data and free GSL objects
    void clear() noexcept;
    

    // minimal getters
    double r_min() const noexcept;
    double r_max() const noexcept;
    int getNr();
    

    // other tools
    void dump_dipole();
    void dump_dipole_interpolated(int Nr_t, double r_max_t, double r_min_t);


private:
    std::vector<double> r_vec_;
    std::vector<double> q_vec_;
    std::vector<double> S_vec_;

    double rmin_;
    double rmax_;

    bool interp_in_log_ = false;   // if true, build spline in log(r) space

    // GSL interpolation objects
    gsl_interp_accel* acc_;
    gsl_spline* spline_;

    // points on the grid;
    int Nr;

    // helper to free gsl objects if allocated
    void free_gsl() noexcept;

    // Tools 
    double MV_Dipole(double r, MV_parameters &pars){
        double exponent = 0.25*pow( r*r * pars.Qs02 , pars.gamma)  * log( pow(r*pars.Lambda, -1) +  E1*pars.ec ) ; 
        return exp( - (pars.sigma0 * pars.T)* exponent);
    }
    double MVA_Dipole(double r, MV_parameters &pars){
        double exponent = 0.25*pow( r*r * pars.Qs02 , pars.gamma)  * log( pow(r*pars.Lambda, -1) +  E1*pars.ec ) ; 
        return exp( - (CA/CF)*(pars.sigma0 * pars.T)* exponent);
    }
};

#endif // DIPOLE_S_HPP
