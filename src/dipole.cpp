// file: dipole.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <stdexcept>

#include "include/dipole.h"

Dipole::Dipole() noexcept
    : r_vec_(),q_vec_(), S_vec_(), acc_(nullptr), spline_(nullptr) {}

Dipole::Dipole(Dipole&& other) noexcept
    : r_vec_(std::move(other.r_vec_)),
      q_vec_(std::move(other.q_vec_)),
      S_vec_(std::move(other.S_vec_)),
      acc_(other.acc_),
      spline_(other.spline_)
{
    other.acc_ = nullptr;
    other.spline_ = nullptr;
}


double Dipole::operator[](std::size_t i) const noexcept {
    return S_vec_[i];
}

double Dipole::r_at(std::size_t i) const noexcept {
    return r_vec_[i];
}
double Dipole::S_at(std::size_t i) const noexcept {
    return S_vec_[i];
}

std::size_t Dipole::size() const noexcept {
    return S_vec_.size();
}

double Dipole::r_min() const noexcept {
    return r_vec_.empty() ? 0.0 : r_vec_.front();
}
double Dipole::r_max() const noexcept {
    return r_vec_.empty() ? 0.0 : r_vec_.back();
}
int Dipole::getNr(){
    return r_vec_.size();
}


Dipole& Dipole::operator=(Dipole&& other) noexcept {
    if (this != &other) {
        free_gsl();
        r_vec_ = std::move(other.r_vec_);
        S_vec_ = std::move(other.S_vec_);
        acc_ = other.acc_;
        spline_ = other.spline_;
        other.acc_ = nullptr;
        other.spline_ = nullptr;
    }
    return *this;
}

Dipole::~Dipole() noexcept {
    free_gsl();
}

void Dipole::free_gsl() noexcept {
    if (spline_) {
        gsl_spline_free(spline_);
        spline_ = nullptr;
    }
    if (acc_) {
        gsl_interp_accel_free(acc_);
        acc_ = nullptr;
    }
}


//// ROUTINES


void Dipole::create_r_grid(int Nr_t, double r_max, double r_min,std::vector<double>& r_grid,std::vector<double>& q_grid){
    // log-spaced grid
    double dq= std::log(r_max/r_min)/(Nr_t-1.);
    for (int i = 0; i < Nr_t; ++i) {
        q_grid[i] = dq*i;
        r_grid[i] = r_min*std::exp(q_grid[i]);
    }
}

void Dipole::init_MV(int Nr_t, double r_max, double r_min, MV_parameters &params){
    
    Nr=Nr_t;
    std::vector<double> r_tmp(Nr);
    std::vector<double> q_tmp(Nr);
    std::vector<double> S_tmp(Nr);
    create_r_grid(Nr,r_max,r_min, r_tmp, q_tmp);
    for (int i = 0; i < Nr_t; ++i) {
        S_tmp[i] = MV_Dipole(r_tmp[i], params);
    }
    r_vec_=r_tmp;
    q_vec_=q_tmp;
    S_vec_=S_tmp;

    // build interpolator
    rebuildInterpolator();
}   

void Dipole::init_MVA(int Nr_t, double r_max, double r_min, MV_parameters &params){
    
    Nr=Nr_t;
    std::vector<double> r_tmp(Nr);
    std::vector<double> q_tmp(Nr);
    std::vector<double> S_tmp(Nr);
    create_r_grid(Nr,r_max,r_min, r_tmp, q_tmp);
    for (int i = 0; i < Nr_t; ++i) {
        S_tmp[i] = MVA_Dipole(r_tmp[i], params);
    }
    r_vec_=r_tmp;
    q_vec_=q_tmp;
    S_vec_=S_tmp;

    // build interpolator
    rebuildInterpolator();
}   



void Dipole::rebuildInterpolator() {
    // free old gsl objects if present
    free_gsl();

    const std::size_t n = q_vec_.size();
    if (n < 2)
        throw std::runtime_error("Not enough points to build interpolator.");

    // allocate new GSL objects
    acc_ = gsl_interp_accel_alloc();
    if (!acc_)
        throw std::runtime_error("Failed to allocate gsl_interp_accel.");

    spline_ = gsl_spline_alloc(gsl_interp_cspline, static_cast<size_t>(n));
    if (!spline_) {
        gsl_interp_accel_free(acc_);
        acc_ = nullptr;
        throw std::runtime_error("Failed to allocate gsl_spline.");
    }

    // prepare x-array: either r_vec_ or log(r_vec_)
    std::vector<double> xvec(n);
    for (std::size_t i = 0; i < n; ++i) xvec[i] = q_vec_[i];

    // initialize spline: gsl copies the pointers (so xvec must persist until this call completes)
    int status = gsl_spline_init(spline_, xvec.data(), S_vec_.data(), static_cast<size_t>(n));
    if (status != GSL_SUCCESS) {
        free_gsl();
        throw std::runtime_error(std::string("gsl_spline_init failed: ") + gsl_strerror(status));
    }
}

double Dipole::eval(double r_query) const {
    if (S_vec_.empty() || r_vec_.empty())
        throw std::runtime_error("Interpolator not initialized (empty data).");

    if (!(r_query > 0.0))
        throw std::invalid_argument("r_query must be > 0 for log interpolation.");

    double xq = std::log(r_query/r_vec_.front());
    if (xq <= q_vec_.front()) return S_vec_.front();
    if (xq >= q_vec_.back()) return S_vec_.back();
    if (!spline_ || !acc_)
        throw std::runtime_error("Interpolator not built. Call init/rebuildInterpolator first.");
    return gsl_spline_eval(spline_, xq, acc_);
    // } else {
    //     if (r_query <= r_vec_.front()) return S_vec_.front();
    //     if (r_query >= r_vec_.back())  return S_vec_.back();
    //     if (!spline_ || !acc_)
    //         throw std::runtime_error("Interpolator not built. Call init/rebuildInterpolator first.");
    //     return gsl_spline_eval(spline_, r_query, acc_);
    // }
}



void Dipole::clear() noexcept {
    r_vec_.clear();
    S_vec_.clear();
    free_gsl();
}



void Dipole::dump_dipole(){
    std::ofstream dipole_B;
    std::ostringstream dipole_B_name;
    dipole_B_name << "dipole_dump.txt";
    dipole_B.open (dipole_B_name.str());
   std::size_t ND = q_vec_.size();
   for (size_t i = 0; i < ND; i++){ dipole_B << r_vec_[i] << "\t" << q_vec_[i] << "\t" << S_vec_[i] << std::endl;}
    dipole_B.close();
}

// void Dipole::dump_dipole(){
//     std::ofstream dipole_B;
//     std::ostringstream dipole_B_name;
//     dipole_B_name << "dipole_dump.txt";
//     dipole_B.open (dipole_B_name.str());
//    std::size_t ND = q_vec_.size();
//    for (size_t i = 0; i < ND; i++){ dipole_B << r_vec_[i] << "\t" << q_vec_[i] << "\t" << S_vec_[i] << std::endl;}
//     dipole_B.close();
// }

void Dipole::dump_dipole_interpolated(int Nr_t, double r_max_t, double r_min_t){

    std::vector<double> r_tmp(Nr_t);
    std::vector<double> q_tmp(Nr_t);
    // std::vector<double> S_tmp(Nr_t);
    create_r_grid(Nr_t,r_max_t,r_min_t, r_tmp, q_tmp);

    std::ofstream dipole_B;
    std::ostringstream dipole_B_name;
    dipole_B_name << "dipole_dump_interpolated.txt";
    dipole_B.open (dipole_B_name.str());
   for (size_t i = 0; i < Nr_t; i++){ dipole_B << r_tmp[i] << "\t" << q_tmp[i] << "\t" << eval( r_tmp[i])  << std::endl;}
    dipole_B.close();
}



// void Dipole::evolve_dipole(const double* r_grid, const double* S_grid, std::size_t n,
//                                bool interp_in_log) {
//     if (!r_grid || !S_grid)
//         throw std::invalid_argument("Null pointer passed to init_from_arrays.");
//     if (n < 2)
//         throw std::invalid_argument("Need at least two grid points.");
//     std::vector<double> r_tmp(r_grid, r_grid + n);
//     std::vector<double> S_tmp(S_grid, S_grid + n);
//     init(r_tmp, S_tmp, interp_in_log);
// }

///// THIS IS FOR LATER 


void Dipole::evolve_dipole(const std::vector<double>& S_grid) {
    if (r_vec_.size() != S_grid.size())
        throw std::invalid_argument("r_grid and S_grid must have same size.");
    if (r_vec_.size() < 2)
        throw std::invalid_argument("Need at least two grid points to interpolate.");
    // check monotonic increase of r_grid and validate positivity if using log
    for (std::size_t i = 1; i < r_vec_.size(); ++i) {
        if (!(r_vec_[i] > r_vec_[i-1]))
            throw std::invalid_argument("r_grid must be strictly increasing.");
    }
    // copy data
    S_vec_ = S_grid;
    // build interpolator
    rebuildInterpolator();
}

bool Dipole::find_effective_saturation_scale(double& r_sat_out,
                                             double& Qs2_out,
                                             double target,
                                             double rel_tol,
                                             int    max_iter) const
{
    const std::size_t n = size();
    if (n < 2)
        throw std::runtime_error("Need at least two points to find saturation scale.");

    // Helper: sign function
    auto sgn = [](double x) { return (x > 0) - (x < 0); };

    // Evaluate S on the native r-grid to find a bracket [i, i+1]
    // where (S - target) changes sign.
    double f_prev = eval(r_at(0)) - target;
    if (std::abs(f_prev) == 0.0) {
        r_sat_out = r_at(0);
        Qs2_out = 2. *pow(r_sat_out,-2.);
        return true;
    }

    std::size_t i_lo = 0, i_hi = 0;
    bool bracketed = false;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        double f_curr = eval(r_at(i + 1)) - target;
        if (f_prev == 0.0) {
            r_sat_out = r_at(i);
            Qs2_out = 2. *pow(r_sat_out,-2.);
            return true;
        }
        if (sgn(f_prev) != sgn(f_curr)) {
            i_lo = i;
            i_hi = i + 1;
            bracketed = true;
            break;
        }
        f_prev = f_curr;
    }

    if (!bracketed) {
        // No crossing on the provided grid: either always above or always below target
        // -> the "solution" lies outside [r_min, r_max].
        return false;
    }

    // Bisection in log-space between r_lo and r_hi (more stable across decades)
    double r_lo = r_at(i_lo);
    double r_hi = r_at(i_hi);
    double x_lo = std::log(r_lo/r_vec_.front());
    double x_hi = std::log(r_hi/r_vec_.front());

    double f_lo = eval(r_lo) - target;
    double f_hi = eval(r_hi) - target;

    // Standard bisection
    for (int it = 0; it < max_iter; ++it) {
        double x_mid = 0.5 * (x_lo + x_hi);
        double r_mid = r_vec_.front()*std::exp(x_mid);
        double f_mid = eval(r_mid) - target;

        // Convergence: relative bracket size
        if (std::abs(x_hi - x_lo) <= std::max(1e-15, rel_tol * std::abs(x_mid))) {
            r_sat_out = r_mid;
            Qs2_out = 2. *pow(r_sat_out,-2.);
            return true;
        }

        // Keep the side with the sign change
        if ((f_lo > 0 && f_mid > 0) || (f_lo < 0 && f_mid < 0)) {
            x_lo = x_mid;
            r_lo = r_mid;
            f_lo = f_mid;
        } else {
            x_hi = x_mid;
            r_hi = r_mid;
            f_hi = f_mid;
        }
    }

    // If we exit the loop, return the midpoint as best effort
    double x_mid = 0.5 * (x_lo + x_hi);
    r_sat_out = r_vec_.front()*std::exp(x_mid);
    Qs2_out = 2. *pow(r_sat_out,-2.);
    return true;
}


// void Dipole::init(const std::vector<double>& r_grid, const std::vector<double>& S_grid,
//                    bool interp_in_log) {
//     if (r_grid.size() != S_grid.size())
//         throw std::invalid_argument("r_grid and S_grid must have same size.");
//     if (r_grid.size() < 2)
//         throw std::invalid_argument("Need at least two grid points to interpolate.");
//     // check monotonic increase of r_grid and validate positivity if using log
//     for (std::size_t i = 1; i < r_grid.size(); ++i) {
//         if (!(r_grid[i] > r_grid[i-1]))
//             throw std::invalid_argument("r_grid must be strictly increasing.");
//     }
//     if (interp_in_log) {
//         for (double rv : r_grid) {
//             if (!(rv > 0.0))
//                 throw std::invalid_argument("All r values must be > 0 for log interpolation.");
//         }
//     }

//     // copy data
//     r_vec_ = r_grid;
//     S_vec_ = S_grid;
//     interp_in_log_ = interp_in_log;

//     // build interpolator
//     rebuildInterpolator();
// }


// void Dipole::init_from_arrays(const double* r_grid, const double* S_grid, std::size_t n,
//                                bool interp_in_log) {
//     if (!r_grid || !S_grid)
//         throw std::invalid_argument("Null pointer passed to init_from_arrays.");
//     if (n < 2)
//         throw std::invalid_argument("Need at least two grid points.");
//     std::vector<double> r_tmp(r_grid, r_grid + n);
//     std::vector<double> S_tmp(S_grid, S_grid + n);
//     init(r_tmp, S_tmp, interp_in_log);
// }


// void Dipole::init_from_file(const std::string& filename, bool interp_in_log) {
//     std::ifstream infile(filename);
//     if (!infile.is_open())
//         throw std::runtime_error("Cannot open file: " + filename);

//     std::vector<double> r_vals;
//     std::vector<double> S_vals;

//     std::string line;
//     while (std::getline(infile, line)) {
//         // trim leading spaces for robust comment detection
//         std::string::size_type pos = line.find_first_not_of(" \t\r\n");
//         if (pos == std::string::npos) continue;
//         if (line[pos] == '#') continue;

//         std::istringstream iss(line);
//         double r, S;
//         if (!(iss >> r >> S)) continue; // skip malformed line

//         r_vals.push_back(r);
//         S_vals.push_back(S);
//     }

//     infile.close();

//     if (r_vals.size() < 2)
//         throw std::runtime_error("File must contain at least two valid data lines.");

//     // check monotonic r
//     for (std::size_t i = 1; i < r_vals.size(); ++i) {
//         if (!(r_vals[i] > r_vals[i - 1]))
//             throw std::runtime_error("r column in file must be strictly increasing.");
//     }

//     if (interp_in_log) {
//         for (double rv : r_vals) {
//             if (!(rv > 0.0))
//                 throw std::runtime_error("All r values must be > 0 for log interpolation.");
//         }
//     }

//     init(r_vals, S_vals, interp_in_log);
// }


// double Dipole::eval(double q_query) const {
//     if (S_vec_.empty() || r_vec_.empty())
//         throw std::runtime_error("Interpolator not initialized (empty data).");

//     if (!(q_query >= 0.0))
//         throw std::invalid_argument("q_query must be > 0 for log interpolation.");

//     // double xq = std::log(r_query/r_vec_.front());
//     if (q_query <= q_vec_.front()) return S_vec_.front();
//     if (q_query >= q_vec_.back()) return S_vec_.back();
//     if (!spline_ || !acc_)
//         throw std::runtime_error("Interpolator not built. Call init/rebuildInterpolator first.");
//     return gsl_spline_eval(spline_, q_query, acc_);
// }