# Dam3D-DRBEM: 3D Arch Dam–Reservoir Interaction (PhD research code, 1997–1998)

Fortran 90 research code written by **Yasin Fahjan** for the PhD thesis
*Arch Dam–Reservoir Interaction Modeling and Application to Karakaya Dam*
(Boğaziçi University, İstanbul, 1998),
supervised by Assoc. Prof. Osman S. Börekçi and Prof. Mustafa Erdik.

The code is published **as an archive**: it documents the method used in the thesis
and in the papers listed below. It is not maintained and is not meant for design work.

## What it does

**Dam3D** couples a finite element model of the dam with a boundary element model of the reservoir:

- **Dam (FEM):** 3D thick-shell and 8-node solid elements; free vibration (eigenvalue) analysis.
- **Reservoir (DRBEM):** dual reciprocity boundary elements (linear 4-node and quadratic 8-node) for the hydrodynamic pressure on the upstream face (Helmholtz equation, water compressibility, reservoir-bottom absorption through the reflection coefficient α).
- **Radiation condition** at the far end of the reservoir, with three options:
  - `_D`: 2D DRM radiation matrix (2-node linear / 3-node quadratic)
  - `_F`: 2D FEM radiation matrix (4-node linear / 8-node quadratic)
  - `_A`: analytical expression (2-node linear) combined with the continuum solution
- **Analysis types:** free vibration; dam only; full dam–reservoir interaction; reservoir with deformable or rigid dam motion.
- **Response:** frequency response functions, response spectrum and time-history analysis using FFT.

## Repository layout

| Folder | Content |
|---|---|
| `dam3d/src/general` | Main program (`DArch3D.f90`), global data module, shape functions, modal and frequency-domain solution drivers |
| `dam3d/src/structure` | 3D FEM of the dam (shell and solid elements, eigenvalue problem) |
| `dam3d/src/hydro` | 3D DRBEM of the reservoir (`ABem3D.f90` module) |
| `dam3d/src/radiation` | Radiation-boundary options A, D and F |
| `dam3d/src/graphics` | Output of results for plotting |
| `dam3d/tools/karakaya-mesh` | Mesh generator used for the Karakaya Dam model (dam, reservoir and radiation boundary) |
| `drbem3d-helmholtz` | Stand-alone 3D DRBEM solver for the Helmholtz equation (linear and quadratic elements), Jan 1998 |
| `bem3d-helmholtz` | Stand-alone 3D BEM solver for the Helmholtz equation, with benchmark inputs in `examples/` |
| `bem3d-laplace` | 3D BEM solver for the Laplace equation (linear elements), with examples |
| `mesh3d` | Box-type surface mesh generator for the 3D BEM programs (Aug 1997), with examples |

Files named `*1.f90` (`HModal1`, `RModal_D1`, `RRad_D1`) are alternative trial versions that were not part of the final build.

## Building

The code was written for **Microsoft Fortran PowerStation / Digital Visual Fortran** on Windows and uses the
**IMSL** numerical library (`use msimsl`) and the Microsoft `portlib` module. It will not compile with
gfortran or Intel Fortran without changes.

The original Dam3D build (DRM radiation option) compiled these files, with the modules first:

1. Modules: `ADam3D`, `AShape`, `AModal`, `AFem3D`, `ABem3D`, `ARad_D`, `AGraphics`
2. Everything else in `general`, `structure`, `hydro`, `graphics`, plus `RInput`, `RModal_D`, `RRad_D`, `DModalRD`

For the FEM or analytical radiation option, use the `_F` or `_A` files (`ARad_F`, `RModal_F`, `RRad_F`, `DModalRF`, ...) instead of the `_D` files.

### External routines a port needs to replace

| IMSL routine | Purpose | LAPACK equivalent |
|---|---|---|
| `LINRG`, `DLINRG` | Real matrix inverse | `?GETRF` + `?GETRI` |
| `LINCG`, `DLINCG` | Complex matrix inverse | `ZGETRF` + `ZGETRI` |
| `LSLRG` | Real linear system | `?GESV` |
| `LSLCG`, `DLSLCG` | Complex linear system | `ZGESV` |
| `LINDS` | Symmetric positive-definite inverse | `?POTRF` + `?POTRI` |
| `EVESF` | Some eigenpairs of a real symmetric matrix | `?SYEVX` |
| `GVCSP` | Generalized symmetric eigenproblem | `?SYGV` |
| `EVCCG`, `DEVCCG` | Complex general eigenproblem | `ZGEEV` |
| `GQRUL` | Gauss quadrature points and weights | any Gauss–Legendre routine |

Other non-standard features: DEC variable format expressions (`<n>` inside `FORMAT`), and the `date()` / `clock()` functions from `portlib`.

### Running

Each program asks for an input file name and an output file name at the console. The input format can be read from
the `*Input*.f90` routines of each program; the examples in `examples/` show the layout for the stand-alone BEM programs.
The original Dam3D input files are not included.

## References

- Fahjan, Y. M., Börekçi, O. S., Erdik, M. (2003). Earthquake-induced hydrodynamic pressures on a 3D rigid dam–reservoir system using DRBEM and a radiation matrix. *International Journal for Numerical Methods in Engineering*, 56, 1511–1532. https://doi.org/10.1002/nme.623
- Fahjan, Y. M., Börekçi, O. S. (2000). 3-D DRM dam–reservoir interaction modeling. *4th International Conference on Advances in Civil Engineering (ACE 2000)*, Eastern Mediterranean University, Gazimağusa.
- Fahjan, Y. M. (1998). *Arch Dam–Reservoir Interaction Modeling and Application to Karakaya Dam*. PhD thesis, Boğaziçi University, İstanbul.

The DRBEM formulation follows Partridge, P. W., Brebbia, C. A., Wrobel, L. C. (1992). *The Dual Reciprocity Boundary Element Method*. Computational Mechanics Publications.

## Citation

If you use this code, please cite the 2003 IJNME paper above (see `CITATION.cff`).

## Author

Prof. Dr. Yasin Fahjan, Istanbul Technical University, Department of Civil Engineering · ORCID [0000-0003-1254-4526](https://orcid.org/0000-0003-1254-4526)

## License

MIT, see `LICENSE`.
