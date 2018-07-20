
# Normaliz -- a tool for discrete convex geometry

Normaliz is a open source tool for computations in affine monoids, vector configurations, rational polyhedar and rational cones. The variant QNormaliz now computes algebraic polyhedra, i.e., polyhedra defined over real algebraic extensions of QQ.

## Computation goals

- convex hulls and dual cones
- conversion from generators to constraints and vice versa
- projections of cones and polyhedra
- triangulations, disjoint decompositions and Stanley decompositions
- Hilbert basis of rational, not necessarily pointed cones
- normalization of affine monoids
- lattice points of rational polytopes and (unbounded) polyhedra
- Euclidean and lattice normalized volumes of polytopes
- Hilbert (or Ehrhart) series and (quasi) polynomials under Z-gradings (for example, for rational polytopes)
- generalized (or weighted) Ehrhart series and Lebesgue integrals of - polynomials over rational polytopes

Normaliz offers the API libnormaliz that allows the user to access tNormaliz computations from any C++ program.

The frontend Normaliz reads input files and writes output files. There is a wide variety of input types to specify polyhedra and lattices by generators (vertices, extreme rays) or by constraints (inequalities, equations and congruences). The user sets computation goals and chooses algorithmic variants through comman line options.

Online exploration of Normaliz: https://mybinder.org/v2/gh/Normaliz/NormalizJupyter/master

## Sample input and output

The file 2cone.in from the directory example contains

    amb_space 2
    cone 2
    1 3
     2 1

It defines a cone in two-dimensional real space by its extreme rays. 
![2-dimensional cone](https://github.com/Normaliz/Normaliz/blob/master/doc/2cone.pdf)

The command

    normaliz example/2cone

runs Normaliz with its default computation goals. Ot produces the output file 2cone.out:

## Platforms
Each releae contains executables for Linux 64, MacOS X and MS Windows 64.

## Interfaces
Normaliz can be called from several other systems:
- [CoCoA](http://cocoa.dima.unige.it)
- [GAP](https://github.com/gap-packages/NormalizInterface)
- [Macaulay2](http://www2.macaulay2.com/Macaulay2/)
- [polymake](https://polymake.org/doku.php)
- [Singular](https://www.singular.uni-kl.de/)
- [SageMath](http://www.sagemath.org/)

The Python packages `PyNormaliz` and `PyQNormaliz` provide an envirinment for interactive access. They are part of the Normaliz repository.

`jNormaliz` provides a GUI to Normaliz

![Normaliz.jpg](https://github.com/Normaliz/Normaliz/blob/master/doc/jNormaliz.jpg)

## Optional packages
For its basic functionality Normaliz needs only GMP and Boost. Pararllelization is based on OpenMP. For the computation of integrals [CoCoALib](http://cocoa.dima.unige.it) is used.

QNormaliz needs [Flint](http://www.flintlib.org/), [antic](https://github.com/wbhart/antic), [arb](http://arblib.org/) and [e-antic](https://github.com/videlec/e-antic).

## Installation

Download  and decompress
- the source basic package `normaliz-x.y.z.zip` (or tar.gz) or the extended source package `normaliz-.y.z-full.zip` (`x.y.z` denotes the version) from the release page of this repository. 

Download and decompress
- the executable for your system (`normaliz-x.y.zLinux64.zip`, `normaliz-x.y.zMacOS.zip` or `normaliz-x.y.zLinux64.zip`).

Or compile Normaliz yourself on Linux or MacOS by one of the installation scripts
- `install_normaliz_with_opt.sh` (only Normaliz)
- `install_normaliz_with_qnormaliz_eantic.sh` (Normaliz and QNormaliz)

## Docker image

available from https://hub.docker.com/r/normaliz/normaliz/








