## Building Normaliz yourself

Normaliz offers the luxury of three build systems:

1. autotools,
2. cmake
3. “classical” make using the “handwritten” Makefile.classic.

The basic steps are the same for all three systems, namely

* configuration
* compilation 
* installation.

The main difference is the way how the build system is configured whereas compilation and installation are essentially identical. In all cases the compilation of NmzIntegrate is included.

In the following we describe the basic steps of autotools and cmake for Linux 64 bit under certain standard assumptions. The file INSTALL in the directory source contains more information.

The autotools scripts have been written by [Matthias Köppe](https://github.com/mkoeppe). The Normaliz team thanks him cordially for his generous help.

## Compiler prerequisites

We require some C++11 features (e.g. `std::exception_ptr`), supported by:

* GNU g++ 4.4 
* clang++ 2.9 
* Intel icpc 12.0

See [here](https://github.com/Normaliz/Normaliz/issues/26) for a more detailed discussion. The mentioned compilers are also able to handle OpenMP 3.0, with the exception of clang++, there the first OpenMP support was introduced in 3.7.

## Required libraries

For compiling Normaliz the following libraries are needed:

* GMP including the C++ wrapper (libgmpxx and libgmp)
* Boost (headers only)

We will only discuss the basic use of cmake for compilation, see the file `source/INSTALL` for additional information, especially on how to use customized paths. Also the use of the autotools system is explained in this file. The “classical” Makefile is meant for development. Therefore you should use autotools or cmake.

The installation will store the files in standard locations, and we assume in the following that they do not need individual include paths.

## Optional packages

As discussed in the manual, Normaliz can profit from the use of SCIP. If you want to use it, SCIP must be installed before the compilation of Normaliz, independently of the method used for building Normaliz.

To build SCIP download the scipoptsuite at http://scip.zib.de/. Notice that SCIP is not distributed under GPL, but the ZIB Academic License (http://scip.zib.de/academic.txt). Unpack it and then compile it with

```make ZLIB=false GMP=false READLINE=false scipoptlib```

Another optional package is CoCoALib, but only to the extent to that NmzIntegrate is optional. If you want to compile NmzIntegrate, CoCoALib must be installed first. The following sequence of commands will install it in the subdirectory CoCoA of your home directory.

```mkdir ~/CoCoA/
cd ~/CoCoA/
wget http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-0.99543.tgz
tar xvf CoCoALib-0.99543.tgz
cd CoCoALib-0.99543
./configure --threadsafe-hack
make library -j2
```

If CoCoALib-0.99543 is no longer available, replace it by a newer version.

## Getting Started

Using code from the GitHub repository, the first thing the user should do is run

```./bootstrap.sh```

## autotools

To build Normaliz with the autotools system, navigate to the Normaliz directory and issue the following sequence of commands:

```
./configure
make
```

This will compile Normaliz, but most likely without SCIP and NmzIntegrate since they need the optional libraries mentioned above and these must be found. If they are not located at standard places, you must specify their paths. Examples (on the machine of a Normaliz team member):

```./configure --with-scipoptsuite-src=$HOME/SCIP/cipoptsuite-3.2.0/```

or

```./configure --with-cocoalib=$HOME/CoCoA/CoCoALib-0.99543```

or with both paths. If the libraries are found, Normaliz will be compiled with SCIP and NmzIntegrate will be made, respectively, by the `make` command. Check the terminal output of `./configure` for success.

The next step is 

```make```

After this step you will find `normaliz` (and `nmzIntegrate`) in the directory source (and `maxsimplex` in its directory). 

The last, optional step is

```sudo make install```

It copies the header files, the library `libnormaliz` and the executables (except `maxsimplex`) into subdirectories of `/usr/local`. It is of course possible to specify another installation path in the call of `./configure`.

Note: Unfortunately, the paths for SCIP are version dependent. We have tested versions 3.2.0 and 3.2.1.

## cmake

You may need to install cmake first:

```sudo apt-get cmake cmake-curses-gui```

To build Normaliz with cmake, start by creating a build directory within the Normaliz directory, say BUILD. Then change the working directory to BUILD.

The basic configuration (equivalent to configure of autotools) is

```cmake ../source```

Then `make` and `make install` will complete the basic installation. For the inclusion of SCIP, use (for example)

```SCIP_DIR=$HOME/SCIP/scipoptsuite-3.2.0/ cmake ../source``` 
     
replacing `$HOME/SCIP/scipoptsuite-3.2.0/` with your own path to SCIP if necessary. Similarly,

```COCOA_DIR=$HOME/CoCoA/CoCoALib-0.99543 cmake ../source/```

Then `make` and `make install` will complete the work. After make the executables can be found in BUILD and its subdirectories `genEhrhart` and `maxsimplex`.

The main advantage of cmake is the existence of a GUI in which you can change most settings originally chosen by cmake. Call `ccmake ../source` (2 times c) or, for a more sophisticated version, `cmake-gui ../source`.

Note: Unfortunately, the paths for SCIP are version dependent. The configuration files for SCIP presently can find the versions 3.2.0 and 3.2.1. For another version you must edit the file `FindSCIP.cmake` in `source/cmake/Modules`.

## Mac OS X

Currently Apple does not supply a compiler which supports OpenMP. We recommend the use of a current gcc.

To install this, and other prerequisites, we recommend using a package manager such as [Fink](http://www.finkproject.org/), [Homebrew](http://brew.sh/), or [MacPorts](https:// www.macports.org/).

You can then follow the instructions for Linux.

Note: Do not try to compile Normaliz with static libraries for Mac OS X.

## Windows

One can compile Windows executables with the Cygwin port of GCC. Unfortunately it is not compatible to OpenMP.

Using Visual Studio is a bit tricky. Microsoft’s C++ compiler does not support OpenMP 3.0. Creating a Normaliz Visual Studio project via cmake is currently not fully supported. The executables that are offered in the Normaliz distribution have been compiled with icpc and a manually created project. Please contact us if you want to build Normaliz on Windows.

## Copyright and how to cite

Normaliz 3.1 is free software licensed under the GNU General Public License, version 3. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with the program. If not, see `http://www.gnu.org/licenses/`.
Please refer to Normaliz in any publication for which it has been used:

    W. Bruns, B. Ichim, T. Römer, R. Sieg and C. Söger: Normaliz. Algorithms for rational cones and affine monoids. Available at `http://normaliz.uos.de` 

The corresponding `\bibitem`:

```\bibitem{Normaliz} W. Bruns, B. Ichim, T. R\"omer, R. Sieg and C. S\"oger:
Normaliz. Algorithms for rational cones and affine monoids.
Available at \url{http://normaliz.uos.de}.
```

A BibTeX entry:

```
@Misc{Normaliz,
author = {W. Bruns and B. Ichim and T. R\"omer and R. Sieg and C. S\"oger},
title = Normaliz. Algorithms for rational cones and affine monoids,
howpublished ={Available at \url{http://normaliz.uos.de}}
```

It is now customary to evaluate mathematicians by such data as numbers of publications, citations and impact factors. The data bases on which such dubious evaluations are based do not list mathematical software. Therefore we ask you to cite the article

```W. Bruns, B. Ichim and C. Söger. The power of pyramid decomposition in Normaliz. J. Symb. Comp. 74 (2016), 513–536.```

in addition. This is very helpful for the younger members of the team.