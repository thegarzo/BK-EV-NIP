# BK-EV-NIP
Balitsky-Kovchegov evolution of impact-parameter-independent dipoles at LO, with running coupling
(several prescriptions), as a function of rapidity Y and local thickness T.

Please cite [XXX] when using this code.

## Requisites

### GSL
This evolution code uses the GSL library for integration routines, as well as for the evaluation of
special functions. The GSL library can be found [here](https://www.gnu.org/software/gsl/).
The Makefile finds it through `gsl-config`. A C++17 compiler is required.

## How to compile
```
make
```
and to run
```
./compute
```
If the output folder already exists, the program asks whether to overwrite it.

## Setting up a run (`src/main.cpp`)
Select an initial condition and an evolution setup by (un)commenting a pair of lines:

| Initial condition | Evolution     | Description                                                        |
|-------------------|---------------|--------------------------------------------------------------------|
| `MVPars`          | `BK_MV`       | MV-Logarithm IC [PRD 88, 114020], massive kernel, daughter rc       |
| `MVgPars`         | `BK_MVg`      | MV-Logarithm with anomalous dimension γ                             |
| `MVePars`         | `BK_MVe`      | MV-Logarithm with infrared constant e_c                             |
| `MVPars` (etc.)   | `BK_MV_Fit` (`_MVg_Fit`, `_MVe_Fit`) | evolution as in the fits of PRD 88, 114020: massless kernel, Balitsky rc, α_s frozen at 0.7, S_A = S_F² |
| `MVRegPars`       | `BK_MVReg`    | MV-Regulated IC and evolution of PRD 111, 054033 (Eq. 8)            |

Grids: `Grids(NR, rmin, rmax, NY, Ymin, Ymax, NT, Tmin, Tmax)`, with r in fm and T in fm⁻² (T > 0).
`NY = 1` or `NT = 1` is allowed.

### Options (fields of `BKParameters`)
- `KType`: `LO`, `rcParent`, `rcDaughter`, `rcSmallest`, `rcBalitsky`
- `mp`: infrared regulator of the kernel; `mp = 0` gives the massless kernel
- `aSType`: `Smooth` (μ₀, ζ regularisation) or `Frozen` (one-loop, frozen at `alpha_freeze`)
- `Adjoint`: `Casimir` (S_A = S_F^{C_A/C_F}) or `Squared` (S_A = S_F²)

The initial condition type is set by `MV_parameters::IC` (`MVLog` or `MVReg`).

## Output
Written to the folder given in `RunModel`:
- `config.yaml`: all parameters and grids
- `mv_dipoles.dat`: columns `T  Y  r[fm]  S_F  S_A`
- `mv_Q2.dat`: columns `T  Y  Qs²_F[GeV²]  Qs²_A[GeV²]`, with S(r_s) = e^{-1/2} and Qs² = 2/r_s²

## Validation
- The MV-Regulated IC reproduces Fig. 1 of PRD 111, 054033.
- The `_Fit` setups reproduce Table I (Qs²), Fig. 2 (σ_r) and Fig. 3 (N(r, x)) of PRD 88, 114020.
