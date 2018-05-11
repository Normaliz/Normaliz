# PyQNormaliz - An interface to QNormaliz


## What is PyQQNormaliz

PyQQNormaliz provides an interface to QNormaliz (https://www.normaliz.uni-osnabrueck.de) via libQNormaliz.
It offers the complete functionality of QNormaliz, and can be used interactively from python. For a first example,
see [this introduction](examples/PyQNormaliz_Tutorial.ipynb) by Richard Sieg.

## Requirements

* python 2.7 or higher or python 3.4 or higher
* Recent version of QNormaliz (https://github.com/Normaliz/Normaliz/tree/enfnormaliz2018)

## Installation

You need to have QNormaliz installed on your system via the following commands
```
git clone https://github.com/Normaliz/Normaliz.git
cd Normaliz
git checkout enfnormaliz2018
NORMALIZ_PATH=${PWD}
./install_normaliz_with_qnormaliz_eantic.sh
make install
```
is enough. If you prefer or are not able to install it, you need to set CPATH and
LD_LIBRARY_PATH accordingly.

After that, you can install PyQNormaliz via
```
git clone https://github.com/sebasguts/PyQNormaliz.git
cd PyQNormaliz
python setup.py build_ext --include-dirs=${NORMALIZ_PATH}/nmz_opt_lib/include --library-dirs=${NORMALIZ_PATH}/nmz_opt_lib/lib
python setup.py install
```

## Usage

See this [notebook](https://nbviewer.jupyter.org/github/sebasguts/PyQNormaliz/blob/master/examples/Dodecahedron.ipynb) as an example.