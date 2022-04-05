# makefile generate poisson mat

CXX=g++


INCEIGEN=-I/usr/include/eigen3
ARMADILLO=-larmadillo

generate_poisson_mat: generate_poisson_mat.cpp
	$(CXX) $^ $(INCEIGEN) $(ARMADILLO) -o $@



clean:
	rm -rf generate_poisson_mat
