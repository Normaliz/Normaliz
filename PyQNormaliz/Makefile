all: module2 module3

module3: NormalizModule.cpp setup.py
	python3 setup.py build_ext --inplace

module2: NormalizModule.cpp setup.py
	python2 setup.py build_ext --inplace

install: install3

install3: module3
	python3 setup.py install

install2: module2
	python2 setup.py install

clean:
	rm PyNormaliz*.so
