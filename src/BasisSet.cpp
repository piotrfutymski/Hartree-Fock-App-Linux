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
		std::vector<std::future<void>> futures;
		for (int j = 0; j <= i; j++)
		{
			auto wsk = &_H_Matrix[i][j];
			auto bs = &_basisFucntion;
			auto n = &_nucleons;
			futures.push_back(std::async(std::launch::async, [wsk, i, j, bs, n]() {
				*wsk = ((*bs)[i]).calculateHIntegral((*bs)[j], *n); }));
		}
		for (auto& e : futures)
			e.wait();
	}

}

void BasisSet::calculateOverlapIntegrals()
{
	_S_Matrix.resize(_basisFucntion.size());

	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_S_Matrix[i].resize(i + 1);
		std::vector<std::future<void>> futures;
		for (int j = 0; j <= i; j++)
		{
			auto wsk = &_S_Matrix[i][j];
			auto bs = &_basisFucntion;
			futures.push_back(std::async(std::launch::async, [wsk, i, j, bs]() {
				*wsk = (*bs)[i].calculateOverlapIntegral((*bs)[j]); }));
		}
		for (auto& e : futures)
			e.wait();
	}
}

void BasisSet::calculateTwoElectronIntegrals()
{
	auto siz = _basisFucntion.size();

	_D_Matrix = new double[siz*siz*siz*siz];
	double * S_Matrix = new double[_allPrimitives*_allPrimitives];

	double * gausianData = new double[5 * _allPrimitives];
	int * info = new int[siz];

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
			v++;
		}		
	}

	for (int i = 0; i < _allPrimitives; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double res = calculateSPrimitive(i, j, gausianData, _allPrimitives);
			S_Matrix[i *_allPrimitives + j] = res;
			S_Matrix[j * _allPrimitives + i] = res;
		}	
	}
	


	for (int p = 0; p < siz; p++)
	{
		for (int q = 0; q <= p; q++)
		{
			for (int r = 0; r < siz; r++)
			{
				for (int s = 0; s <= r; s++)
				{
					double val = this->calculateDContracted(p, q, r, s, info, gausianData, S_Matrix, siz, _allPrimitives);
					int siz2 = siz*siz;
					int siz3 = siz*siz2;
					_D_Matrix[p*siz3 + r*siz2 + q*siz + s] = val;
					_D_Matrix[q*siz3 + r*siz2 + p*siz + s] = val;
					_D_Matrix[p*siz3 + s*siz2 + q*siz + r] = val;
					_D_Matrix[q*siz3 + s*siz2 + p*siz + r] = val;
				}			
			}		
		}		
	}

	delete[] gausianData;
	delete[] info;
	delete[] S_Matrix;

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

double BasisSet::calculateDContracted(int p, int q, int r, int s, int *info, double * primitives, double * SMat, int num, int pnum)
{
	int startP = 0;
    if(p!=0)
        startP = info[p - 1];
    int endP = info[p];

    int startQ = 0;
    if(q!=0)
        startQ = info[q - 1];
	int endQ = info[q];
	
	int startR = 0;
    if(r!=0)
        startR = info[r- 1];
	int endR = info[r];
	
	int startS = 0;
    if(s!=0)
        startS = info[s - 1];
    int endS = info[s];

	double res = 0.0;
    for(int i = startP; i < endP; i++)
    {
		for(int j = startQ; j < endQ; j++)
		{
			for(int k = startR; k < endR; k++)
			{
				for(int l = startS; l < endS; l++)
				{					
					double r = calculateDPrimitive(i, j, k, l, pnum, SMat, primitives);
					res +=primitives[i*5 + 1] * primitives[j*5 + 1] * primitives[k*5 + 1] * primitives[l*5 + 1] * r;
					
				}
			}
		}          
	}
	return res;
}

double BasisSet::calculateDPrimitive(int p, int q, int r, int s, int pnum, double * SMat, double * primitives)
{
	double Spq = SMat[p*pnum + q];
	double Srs = SMat[r*pnum + s];
	double Spqrs = Spq * Srs;
	double ap = primitives[p*5];
	double aq = primitives[q*5];
	double ar = primitives[r*5];
	double as = primitives[s*5];
	double A = ap + aq;
	double B = ar + as;
	double AB = A + B;
	double _A = 1.0/A;
	double _B = 1.0/B;
	double kx = _A*(ap*primitives[p*5 + 2]+ aq*primitives[q*5 + 2]);
	double ky = _A*(ap*primitives[p*5 + 3]+ aq*primitives[q*5 + 3]);
	double kz = _A*(ap*primitives[p*5 + 4]+ aq*primitives[q*5 + 4]);
	double lx = _B*(ar*primitives[r*5 + 2]+ as*primitives[s*5 + 2]);
	double ly = _B*(ar*primitives[r*5 + 3]+ as*primitives[s*5 + 3]);
	double lz = _B*(ar*primitives[r*5 + 4]+ as*primitives[s*5 + 4]);
	double xx = kx - lx;
	double yy = ky - ly;
	double zz = kz - lz;
	double r2 = xx*xx + yy*yy + zz*zz;

	return 2.0/M_SQRTPI * (sqrt(A*B))/sqrt(AB) * BoysCalculator::boys(A*B/AB * r2)*Spqrs;
}

double BasisSet::calculateSPrimitive(int row, int collumn, double * primitives, int primitivesNum)
{
	 double b = primitives[5*row] + primitives[5*collumn];
    double B = primitives[5*row] * primitives[5*collumn];
    double A = 4.0 * B / (b*b);
    double x = primitives[5*row+2] - primitives[5*collumn + 2];
    double y = primitives[5*row+3] - primitives[5*collumn + 3];
    double z = primitives[5*row+4] - primitives[5*collumn + 4];
    double R2 = x*x + y*y + z*z;
    return pow(A, 0.75)* exp(-B*R2/b);
}