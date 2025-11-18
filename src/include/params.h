#ifndef PARAMS_H
#define PARAMS_H
#include <complex>
#include <cmath>

const std::complex<double> M_I = sqrt(std::complex<double>(-1.));
const double E1 = exp(1.0);

const int NC = 3;
const int dA = NC*NC-1;
const double CA = double(NC);
const double CF = double(dA)/(2.*double(NC));
const double CBAR = (NC*NC-4)/(2.*double(NC));
const int NF =3;

const double qB = 1./3.;
const double alpha_S = 0.3;

const double LQCD= 0.241; // In GEV!

const double hbarc= 0.1973;
const double fm_to_GeVm1 =1/hbarc;
const double fm2_to_GeVm2 =fm_to_GeVm1*fm_to_GeVm1;

const double GeV_to_fmm1 =1/hbarc;
const double GeV2_to_fmm2 =GeV_to_fmm1*GeV_to_fmm1;

const double GeVm1_to_fm =hbarc;
const double GeVm2_to_fm2 =GeVm1_to_fm*GeVm1_to_fm;

const double fmm1_to_GeV =hbarc;
const double fmm2_to_GeV2 =fmm1_to_GeV*fmm1_to_GeV;

const double mb_to_fm2 = 0.1;
const double fm2_to_mb = 1/mb_to_fm2;
#endif