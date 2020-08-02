#include "Mol.h"

Mol::Mol(const std::vector<Nucleon>& n)
{
	_nucleons = n;
	for (int i = 0; i < n.size(); i++)
	{
		_MOcount += n[i].charge;
		_nucleons[i].p.x/=BOHR_TO_ANGSTROM;
		_nucleons[i].p.y/=BOHR_TO_ANGSTROM;
		_nucleons[i].p.z/=BOHR_TO_ANGSTROM;
		_nucleons[i].p.reloadSpherical();
	}

	_MOcount = (_MOcount + 1) / 2;

	_nuclearRepulsion = 0.0;
	for (int i = 0; i < _nucleons.size(); i++)
	{
		for (int j = 0; j < i; j++)
		{
			_nuclearRepulsion += _nucleons[i].charge * _nucleons[j].charge / (_nucleons[i].p - _nucleons[j].p).r;
		}
	}
	
}

Mol::Mol(int charge):Mol({ {Position{0,0,0}, charge } })
{
}

void Mol::setMOcount(int c)
{
	_MOcount = c;
	_orbitalCoeficients.resize(_baseSize, c);
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < _baseSize; j++)
		{
			_orbitalCoeficients(j, i) = 0.0;
		}
	}
}

void Mol::setApproxIntegrals(bool approx)
{
    _approxIntegrals = approx;
}

void Mol::HFComputation(Logger & logger)
{
    logger.say("Starting HF Computations for the molecule\n");
	auto start = std::chrono::steady_clock::now();
	try{
		 this->initBasisSet(logger);
	}catch (const std::exception & e){
            throw BasisSetException();
			return;
    }

	this->calculateIntegrals(logger);
	this->HFSelfConsistentFieldMethod(logger);

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	_totalTime = elapsed_seconds.count();
	_computationDone = true;

	logger.say("Computations for the molecule finished\n");
	logger.say("Computation time: ");
	logger.say(std::to_string(_totalTime));
	logger.say("s\n");
}

bool Mol::computationDone()const
{
	return _computationDone;
}

void Mol::initBasisSet(Logger & logger)
{
	logger.say("Loading basis set...\n");

	_basisSet.createBasisSet(_nucleons);
	_baseSize = _basisSet.getSize();
	this->setMOcount(_MOcount);
	_FMatrix.resize(_baseSize, _baseSize);
	_SMatrix.resize(_baseSize, _baseSize);
	_PMatrix.resize(_baseSize, _baseSize);

	logger.say("Basis Set loaded\n");
	logger.say("It contains: ");
	logger.say(std::to_string(_baseSize));
	logger.say(" functions and in total ");
	logger.say(std::to_string(_basisSet.getPrimitiveCount()));
	logger.say(" primitives\n");

}

void Mol::calculateIntegrals(Logger & logger)
{
	logger.say("Calculating integrals...\n");
	_basisSet.calculateIntegrals();
	logger.say("Integrals calculated\n");
	for (int r = 0; r < _baseSize; r++)
	{
		for (int s = 0; s <= r; s++)
		{
			_SMatrix(r, s) = _basisSet.getS(r, s);
			_SMatrix(s, r) = _basisSet.getS(r, s);
		}
	}
			
}

void Mol::HFProcedure()
{
	_oldPMatrix = _PMatrix;
	EigenSolver::solve(_FMatrix, _SMatrix, _orbitalCoeficients, _orbitalEnergies);
	this->recalculateP();
	this->recalculateF();
	this->recalculateEnergy();
}


void Mol::HFSelfConsistentFieldMethod(Logger & logger)
{
	logger.say("Starting SCF Method...\n");
	double oldE;
	double newE = 1000000;
	int i = 0;
	while(true){
		HFProcedure();
		oldE = newE;
		newE = this->getTotalEnergy();

		i++;

		if(fabs((double)(oldE-newE)) < 0.0001)
			break;

		if (i > 150)
		{
			logger.say("Energy doesn't convarges!\n");
			logger.say("Braking SCF\n");
			_doesNotConv = true;
			break;
		}			
	} 

	if(abs(oldE-newE) < 0.00001)
		logger.say("Lack of significant changes in energy\nEnd of SCF\n");
}

void Mol::recalculateP()
{
	for (int i = 0; i < _baseSize; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double s = 0.0;
			for (int k = 0; k < _MOcount; k++)
			{
				s += _orbitalCoeficients(i, k) * _orbitalCoeficients(j, k);
			}
			_PMatrix(i, j) =  2 * s;
			_PMatrix(i, j) = _PMatrix(i, j) * 0.5 + _oldPMatrix(i, j) * 0.5;
			_PMatrix(j, i) = _PMatrix(i, j);

		}
	}

}

void Mol::recalculateF()
{
	for (int r = 0; r < _baseSize; r++)
	{
		for (int s = 0; s <= r; s++)
		{
			double sum = _basisSet.getH(r, s);
			for (int p = 0; p < _baseSize; p++)
			{
				for (int q = 0; q < _baseSize; q++)
				{
					sum += _PMatrix(q, p) * (_basisSet.getDI(r, p, s, q) - _basisSet.getDI(r, p, q, s) / 2.0);
				}
			}
			_FMatrix(r, s) = sum;
			_FMatrix(s, r) = sum;
		}
	}
}

void Mol::recalculateEnergy()
{
	_electronicEnergy = 0.0;
	for (int r = 0; r < _baseSize; r++)
	{
		for (int s = 0; s < _baseSize; s++)
		{
			_electronicEnergy += _PMatrix(s, r) * (_basisSet.getH(s, r) + _FMatrix(r, s));
		}
	}
	_electronicEnergy /= 2.0;
}

double Mol::getTotalEnergy()const
{
	return _electronicEnergy + _nuclearRepulsion;
}

double Mol::getRepulsionEnergy()const
{
	return _nuclearRepulsion;
}

double Mol::getOribtalEnergy(int m)const
{
	return _orbitalEnergies(m);
}

double Mol::getElectronicEnergy()const
{
	return _electronicEnergy;
}

arma::vec Mol::getMolecularCoeficents(int m)
{
	arma::vec res;
	res.resize(_baseSize);

	for (int i = 0; i < _baseSize; i++)
	{
		res(i) = _orbitalCoeficients(i,m);
	}
	return res;
}
double Mol::countMolecularFunction(int m, const Position & p)const
{
	double res = 0.0;
	auto pos = p;
	pos.x/=BOHR_TO_ANGSTROM;
	pos.y/=BOHR_TO_ANGSTROM;
	pos.z/=BOHR_TO_ANGSTROM;
	pos.reloadSpherical();

	for (int i = 0; i < _baseSize; i++)
	{
		res += _basisSet.getFunctionValue(i, pos) * _orbitalCoeficients(i, m);
	}
	return res;
}

double Mol::countMolecularFunction(int m, double x)const
{
	return this->countMolecularFunction(m, {x,0,0});
}
double Mol::countDencity(const Position & p)const
{
	double res = 0.0;
	for (int i = 0; i < _MOcount; i++)
	{
		double v = this->countMolecularFunction(i,p);
		res += v*v;
	}

	return 2.0*res;
	
}
double Mol::getTotalTime()const
{
	return _totalTime;
}

bool Mol::doesNotConv()const
{
	return _doesNotConv;
}

int Mol::getMOcount()const
{
	return _MOcount;
}