#pragma once
#include <vector>
#include <string>
#include <chrono>
#include <exception>
#include <armadillo>
#include <cmath>
#include "BasisSet.h"
#include "EigenSolver.h"
#include "Nucleon.h"
#include "Logger.h"

class BasisSetException : public std::exception{
    const char * what () const throw () {
      return "Unable to load basis Set, check basis_set folder";
   }
};


class Mol
{
public:

	Mol(){};
	Mol(const std::vector<Nucleon>& n);
	Mol(int charge);

    // init

	void setMOcount(int c);
    void setApproxIntegrals(bool approx);

    // the most important function - it makes all computation for molecule

    void HFComputation(Logger & logger);

    //

    bool computationDone()const;

    //results

	double getElectronicEnergy()const;
	double getTotalEnergy()const;
    double getRepulsionEnergy()const;
    double getOribtalEnergy(int m)const;
    arma::vec getMolecularCoeficents(int m);
	double countMolecularFunction(int m, const Position & p)const;
    double countMolecularFunction(int m, double x)const;
    double countDencity(const Position & p)const;
    double getTotalTime()const;

    bool doesNotConv()const;

    int getMOcount()const;

private:

    //molecule settings

	std::vector<Nucleon> _nucleons;
	int _MOcount{0};
    bool _approxIntegrals{ false };

    //molecule state

    bool _computationDone{ false };

    // set

    BasisSet _basisSet;

    // matrix for computation

	arma::mat _FMatrix;
    arma::mat _SMatrix;
    arma::mat _PMatrix;
	arma::mat _oldPMatrix;
    arma::mat _orbitalCoeficients;

    // for computations

    int _baseSize{ 0 };

    // results

    double _nuclearRepulsion;
	double _electronicEnergy;
    arma::vec _orbitalEnergies;

    bool _doesNotConv{false};

    //timing

    double _totalTime;

private:

    //loading set from files

	void initBasisSet(Logger & logger);

    //calculating Integrals

	void calculateIntegrals(Logger & logger);

    //Self Consistent Field Method

	void HFSelfConsistentFieldMethod(Logger & logger);


    //Help functions

    void HFProcedure();
    void recalculateP();
    void recalculateF();
    void recalculateEnergy();

    void calculateNuclearRepulsion();

};