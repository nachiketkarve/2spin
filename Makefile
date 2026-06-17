all: 2spin cyclicPerturbation spectralFunction classical2spin classicalSpectralFunction

2spin: src/2spin.cpp
	g++ -c ./src/2spin.cpp -fopenmp -o ./lib/2spin.o -O3 -I./include -I./../Libraries/eigen

cyclicPerturbation: src/cyclicPerturbation.cpp
	g++ -c ./src/cyclicPerturbation.cpp -fopenmp -o ./lib/cyclicPerturbation.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbation.out -fopenmp -O3 ./lib/cyclicPerturbation.o ./lib/2spin.o

spectralFunction: src/spectralFunction.cpp
	g++ -c ./src/spectralFunction.cpp -fopenmp -o ./lib/spectralFunction.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _spectralFunction.out -fopenmp -O3 ./lib/spectralFunction.o ./lib/2spin.o

classical2spin: src/classical2spin.cpp
	g++ -c ./src/classical2spin.cpp -fopenmp -o ./lib/classical2spin.o -O3 -I./include -I./../Libraries/eigen

classicalSpectralFunction: src/classicalSpectralFunction.cpp
	g++ -c ./src/classicalSpectralFunction.cpp -fopenmp -o ./lib/classicalSpectralFunction.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _classicalSpectralFunction.out -fopenmp -O3 ./lib/classicalSpectralFunction.o ./lib/classical2spin.o