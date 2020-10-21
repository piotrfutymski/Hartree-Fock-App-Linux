#include "BasisSet.h"
#include "Mol.h"


void BasisSet::createBasisSet(const std::vector<Nucleon>& nucleons)
{
	_nucleons = nucleons;
	_size = 0;
	_allPrimitives = 0;
	for (int i = 0; i < nucleons.size(); i++)
	{
		auto date = this->loadAOs(_nucleons[i].charge);
		for (int j = 0; j < date.size(); j++)
		{
			_basisFucntion.push_back({});
			for (int k = 0; k < date[j].size(); k++)
			{
				(_basisFucntion.end() - 1)->addPrimitive(
					std::get<1>(date[j][k]), {
						 nucleons[i].p + std::get<2>(date[j][k]),
						 std::get<0>(date[j][k])
					}
				);
				_allPrimitives++;
			}
			_size++;

		}
	}
}

void BasisSet::createTestBasisSet(const std::vector<Nucleon>& nucleons, double b, double s)
{
	_nucleons = nucleons;
	for (int i = 0; i < nucleons.size(); i++)
	{
		auto date = this->loadAOs(_nucleons[i].charge);
		//if (nucleons[i].charge == 8)
		{
			double tab[] = { 4.72797,  1.19034 ,  0.359412, 0.126751 };

			double c = sqrt(1 / (2 - 2 * exp(-s*b*b)));
			date.push_back({
				{s * tab[0], c , {b / sqrt(2 * tab[0]),0,0}},
				{s * tab[0], -c, {-b / sqrt(2 * tab[0]),0,0}}
				});
			date.push_back({
				{s * tab[1], c , {b / sqrt(2 * tab[1]),0,0}},
				{s * tab[1], -c, {-b / sqrt(2 * tab[1]),0,0}}
				});
			date.push_back({
				{s * tab[2], c , {b / sqrt(2 * tab[2]),0,0}},
				{s * tab[2], -c, {-b / sqrt(2 * tab[2]),0,0}}
				});
			date.push_back({
				{s * tab[3], c , {b / sqrt(2 * tab[3]),0,0}},
				{s * tab[3], -c, {-b / sqrt(2 * tab[3]),0,0}}
				});
			date.push_back({
				{s * tab[0], c , {0,b / sqrt(2 * tab[0]),0}},
				{s * tab[0], -c, {0,-b / sqrt(2 * tab[0]),0}}
				});
			date.push_back({
				{s * tab[1], c , {0,b / sqrt(2 * tab[1]),0}},
				{s * tab[1], -c, {0,-b / sqrt(2 * tab[1]),0}}
				});
			date.push_back({
				{s * tab[2], c , {0,b / sqrt(2 * tab[2]),0}},
				{s * tab[2], -c, {0,-b / sqrt(2 * tab[2]),0}}
				});
			date.push_back({
				{s * tab[3], c , {0,b / sqrt(2 * tab[3]),0}},
				{s * tab[3], -c, {0,-b / sqrt(2 * tab[3]),0}}
				});
			date.push_back({
				{s * tab[0], c , {0,0,b / sqrt(2 * tab[0])}},
				{s * tab[0], -c, {0,0,-b / sqrt(2 * tab[0])}}
				});
			date.push_back({
				{s * tab[1], c , {0,0,b / sqrt(2 * tab[1])}},
				{s * tab[1], -c, {0,0,-b / sqrt(2 * tab[1])}}
				});
			date.push_back({
				{s * tab[2], c , {0,0,b / sqrt(2 * tab[2])}},
				{s * tab[2], -c, {0,0,-b / sqrt(2 * tab[2])}}
				});
			date.push_back({
				{s * tab[3], c , {0,0,b / sqrt(2 * tab[3])}},
				{s * tab[3], -c, {0,0,-b / sqrt(2 * tab[3])}}
				});
		}	
		for (int j = 0; j < date.size(); j++)
		{
			_basisFucntion.push_back({});
			for (int k = 0; k < date[j].size(); k++)
			{
				(_basisFucntion.end() - 1)->addPrimitive(
					std::get<1>(date[j][k]), {
						 nucleons[i].p + std::get<2>(date[j][k]),
						 std::get<0>(date[j][k])
					}
				);
			}		
		}

		//this->saveAOs(date, _nucleons[i].charge);
	}

}

void BasisSet::calculateOneElectronHamiltonians()
{
	_H_Matrix.resize(_basisFucntion.size());
	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_H_Matrix[i].resize(i+1);
		for (int j = 0; j <= i; j++)
		{
			_H_Matrix[i][j] = _basisFucntion[i].calculateHIntegral(_basisFucntion[j], _nucleons);
		}
	}

}

void BasisSet::calculateOverlapIntegrals()
{
	_S_Matrix.resize(_basisFucntion.size());

	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_S_Matrix[i].resize(i + 1);
		for (int j = 0; j <= i; j++)
		{
			_S_Matrix[i][j] = _basisFucntion[i].calculateOverlapIntegral(_basisFucntion[j]);
		}
	}
}

void BasisSet::calculateTwoElectronIntegrals()
{

	auto siz = _basisFucntion.size();

	_D_Matrix = (double *)malloc(siz*siz*siz*siz*sizeof(double));
	

	double * gausianData = (double *)malloc(5 * _allPrimitives * sizeof(double));
	int * info = (int*)malloc(siz * sizeof(int));
	
	int v = 0;
	int sum = 0;
	for(int i = 0; i <siz; ++i)
	{
		sum+=_basisFucntion[i].getPrimitives().size();
		info[i] = sum;
		for (auto & p: _basisFucntion[i].getPrimitives())
		{
			gausianData[5*v] = p.second.getalfa();
			gausianData[5*v+1] = p.first;
			auto pos = p.second.getR();
			gausianData[5*v+2] = pos.x;
			gausianData[5*v+3] = pos.y;
			gausianData[5*v+4] = pos.z;
			v+=5;
		}
		
	}

	calculateDIntegrals(gausianData , info, _D_Matrix, siz, _allPrimitives);

	free(gausianData);
	free(info);

}

void BasisSet::calculateIntegrals()
{
	this->calculateOverlapIntegrals();
	this->calculateOneElectronHamiltonians();
	this->calculateTwoElectronIntegrals();
}

std::vector<std::vector<std::tuple<double, double, Position>>> BasisSet::loadAOs(int n)
{
	std::vector<std::vector<std::tuple<double, double, Position>>> res;
	std::fstream file("basis_set/AO" + std::to_string(n) +".txt", std::ios::in);
	if (!file.is_open())
		throw BasisSetLoadException();
	int na;
	file >> na;
	res.resize(na);
	for (int i = 0; i < res.size(); i++)
	{
		file >> na;
		res[i].resize(na);
		for (int j = 0; j < res[i].size(); j++)
		{
			double a;
			double b;
			Position p;

			file >> a >> b >> p.x >> p.y >> p.z;
			p.reloadSpherical();
			res[i][j] = std::make_tuple(a, b, p);
		}
	}
	file.close();
	return res;
}

int BasisSet::getSize()
{
    return _size;
}

int BasisSet::getPrimitiveCount()
{
    return _allPrimitives;
}

double BasisSet::getH(int i, int j)
{
	if (i < j)
		return _H_Matrix[j][i];
	return _H_Matrix[i][j];
}

double BasisSet::getS(int i, int j)
{
	if (i < j)
		return _S_Matrix[j][i];
	return _S_Matrix[i][j];
}

double BasisSet::getDI(int i, int j, int k, int l)
{
	auto siz = _basisFucntion.size();
	return _D_Matrix[i*siz*siz*siz + j*siz*siz + k*siz + l];
}

double BasisSet::getFunctionValue(int n, const Position& p)const
{
	return _basisFucntion[n].F_X(p).real();
}
