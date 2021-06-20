[![Build Status](https://github.com/Normaliz/Normaliz/workflows/Run%20tests/badge.svg)](https://github.com/Normaliz/Normaliz/actions)
[![Code Coverage](https://codecov.io/github/Normaliz/Normaliz/coverage.svg)](https://codecov.io/gh/Normaliz/Normaliz)

# Normaliz - a tool for discrete convex geometry

Normaliz is a open source tool for computations in affine monoids, vector configurations, rational polyhedra and rational cones. Normaliz now computes rational and algebraic polyhedra, i.e., polyhedra defined over real algebraic extensions of QQ.

## Computation goals

- convex hulls and dual cones
- conversion from generators to constraints and vice versa
- projections of cones and polyhedra
- triangulations, disjoint decompositions and Stanley decompositions
- Hilbert bases of rational, not necessarily pointed cones
- normalizations of affine monoids
- lattice points of polytopes and (unbounded) polyhedra
- automorphism groups (euclidean, integral, rational/algebraic, combinatorial)
- face lattices and f-vectors
- Euclidean and lattice normalized volumes of polytopes
- Hilbert (or Ehrhart) series and (quasi) polynomials under Z-gradings (for example, for rational polytopes)
- generalized (or weighted) Ehrhart series and Lebesgue integrals of - polynomials over rational polytopes

Normaliz offers the API libnormaliz that allows the user to access Normaliz computations from any C++ program.

The frontend Normaliz reads input files and writes output files. There is a wide variety of input types to specify polyhedra and lattices by generators (vertices, extreme rays, bases) or by constraints (inequalities, equations and congruences). The user sets computation goals and chooses algorithmic variants through command line options or the input file.

Online exploration of Normaliz: <https://mybinder.org/v2/gh/Normaliz/NormalizJupyter/master>

## Sample input and output

The file 2cone.in from the directory example contains

    amb_space 2
    cone 2
    1 3
    2 1

It defines a cone in two-dimensional real space by its extreme rays. 
![2-dimensional cone](https://github.com/Normaliz/Normaliz/blob/master/doc/2cone.jpg)

The command

    normaliz example/2cone

runs Normaliz with its default computation goals. It produces the output file 2cone.out (here typeset in two columns):

    4 Hilbert basis elements          embedding dimension = 2
    2 extreme rays                    rank = 2 (maximal)
    2 support hyperplanes             external index = 1
                                      internal index = 5
                                      original monoid is not integrally closed
    
    size of triangulation   = 1       rank of class group = 0
    resulting sum of |det|s = 5       finite cyclic summands:
                                      5: 1  
    No implicit grading found
    
    ***********************************************************************
    
    4 Hilbert basis elements:         2 extreme rays:
     1 1                               1 3
     1 2                               2 1
     1 3
     2 1                              2 support hyperplanes:
                                       -1  2
                                        3 -1

The main point was the computation of the Hilbert basis (encircled in red in the figure).
                                        
## Platforms
Each [release](https://github.com/Normaliz/Normaliz/releases) contains executables for Linux 64, MacOS X and MS Windows 64.

## Interfaces
Normaliz can be called from several other systems:
- [CoCoA](http://cocoa.dima.unige.it)
- [GAP](https://github.com/gap-packages/NormalizInterface)
- [Macaulay2](http://www2.macaulay2.com/Macaulay2/)
- [polymake](https://polymake.org/doku.php)
- [Singular](https://www.singular.uni-kl.de/)
- [SageMath](https://www.sagemath.org/)

The Python package [`PyNormaliz`](https://github.com/Normaliz/PyNormaliz) by Sebastian Gutsche provides an environment for interactive access. It is contained in the source package of Normaliz.

`jNormaliz` by Vicinius Almendra and Bogdan Ichim provides a GUI to Normaliz

![Normaliz.jpg](https://github.com/Normaliz/Normaliz/blob/master/doc/jNormaliz.jpg)

## Optional packages
For its basic functionality Normaliz needs only GMP. Parallelization is based on OpenMP. For the computation of integrals [CoCoALib](http://cocoa.dima.unige.it) is used.

For algebraic polyhedra Normaliz needs [Flint](https://www.flintlib.org/), [arb](https://arblib.org/), [antic](https://github.com/wbhart/antic/) and [e-antic](https://github.com/flatsurf/e-antic)

The computation of automorphism groups uses [nauty](https://users.cecs.anu.edu.au/~bdm/nauty).


## Installation

All files can be found at https://github.com/Normaliz/Normaliz/releases.

For the binary package (ready made executable program) download and extract
- the executable for your system (`normaliz-x.y.zLinux64.zip`, `normaliz-x.y.zMacOS.zip` or `normaliz-x.y.zLinux64.zip`).

For the source package download  and extract
-  `normaliz-x.y.z.zip` (or tar.gz)

From the source one can compile Normaliz oneself on Linux or MacOS by one of the installation scripts
- `install_normaliz_with_opt.sh` (only rational polyhedra)
- `install_normaliz_with_eantic.sh` (with algebraic polyhedra)

## Docker image

available from https://hub.docker.com/r/normaliz/normaliz/

## Distributions

Normaliz is available as a Debian, Gentoo and Ubuntu package, as well as from Conda (Linux, MacOS, MS Windows).
