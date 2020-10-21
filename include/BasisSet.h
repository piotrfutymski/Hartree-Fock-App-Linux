#pragma once
#include <vector>
#include <iostream>
#include <armadillo>
#include <chrono>
#include <future>
#include <fstream>
#include "Nucleon.h"
#include "ContractedGTO.h"

class BasisSetLoadException : public std::exception{
    const char * what () const throw () {
      return "Unable to load AO file";
   }
};


class BasisSet
{
public:

	BasisSet() {delete[] _D_Matrix;};

	BasisSet(const BasisSet & b) = delete;

	void createBasisSet(const std::vector<Nucleon> & nucleons);
	void createTestBasisSet(const std::vector<Nucleon> & nucleons, double l, double s);

	void calculateIntegrals();

	std::vector<std::vector<std::tuple<double, double, Position>>> loadAOs(int n);
	//void saveAOs(const std::vector<std::vector<std::tuple<double, double, Position>>>&, int n);

	int getSize();
    int getPrimitiveCount();

	double getH(int i, int j);
	double getS(int i, int j);
	double getDI(int i, int j, int k, int l);

	double getFunctionValue(int n, const Position& p)const;

private:

	std::vector<ContractedGTO> _basisFucntion;

    std::vector<std::vector<double>> _SP_Matrix;

	std::vector<std::vector<double>> _S_Matrix;

	std::vector<std::vector<double>> _H_Matrix;

	double * _D_Matrix;

	std::vector<Nucleon> _nucleons;

    int _size;
	int _allPrimitives;

private:

    void calculatePrimitiveOverlapIntegrals();
	void calculateOneElectronHamiltonians();
	void calculateOverlapIntegrals();
	void calculateTwoElectronIntegrals();

	void calculateDContracted();
	double calculateDPrimitive();

};

