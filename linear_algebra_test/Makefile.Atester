LIBNMZ = ../source/libnormaliz/libnormaliz.a

all: LAtester LAtester-mpz

LAtester: LAtester.o $(LIBNMZ)
	g++ -O3 LAtester.o $(LIBNMZ) -lgmpxx -lgmp -fopenmp -o LAtester

LAtester.o: LAtester.cpp
	g++ -O3 LAtester.cpp -c -o LAtester.o 

LAtester-mpz: LAtester-mpz.o $(LIBNMZ)
	g++ -O3 LAtester-mpz.o $(LIBNMZ) -lgmpxx -lgmp -fopenmp -o LAtester-mpz

LAtester-mpz.o: LAtester.cpp
	g++ -O3 LAtester.cpp -DUSE_MPZ -c -o LAtester-mpz.o 

clean:
	rm -f LAtester LAtester.o LAtester-mpz LAtester-mpz.o
