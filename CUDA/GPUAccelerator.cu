#include "GPUAccelerator.h"
#define BLOCK_SIZE 16

__device__ double boysSmall(double x)
{
    static double BV[] = 
{
1.000000000000,
0.967643312160,
0.937150027892,
0.908392875737,
0.881254033291,
0.855624389925,
0.831402868958,
0.808495804396,
0.786816367708,
0.766284040554,
0.746824129652,
0.728367320368,
0.710849265820,
0.694210208622,
0.678394632588,
0.663350941956,
0.649031165907,
0.635390686316,
0.622387986848,
0.609984421694,
0.598144002338,
0.586833200921,
0.576020768860,
0.565677569502,
0.555776423688,
0.546291967195,
0.537200519122,
0.528479960322,
0.520109621122,
0.512070177558,
0.504343555480,
0.496912841897,
0.489762202988,
0.482876808268,
0.476242760420,
0.469847030352,
0.463677397081,
0.457722392049,
0.451971247559,
0.446413848984,
0.441040690473,
0.435842833889,
0.430811870723,
0.425939886757,
0.421219429278,
0.416643476636,
0.412205409980,
0.407898987001,
0.403718317527,
0.399657840853,
0.395712304644,
0.391876745326,
0.388146469828,
0.384517038595,
0.380984249771,
0.377544124458,
0.374192892992,
0.370926982140,
0.367743003168,
0.364637740707,
0.361608142366,
0.358651309035,
0.355764485825,
0.352945053614,
0.350190521132,
0.347498517576,
0.344866785694,
0.342293175312,
0.339775637284,
0.337312217822,
0.334901053181,
0.332540364685,
0.330228454060,
0.327963699058,
0.325744549352,
0.323569522681,
0.321437201231,
0.319346228237,
0.317295304790,
0.315283186830,
0.313308682323,
0.311370648606,
0.309467989882,
0.307599654866,
0.305764634567,
0.303961960194,
0.302190701188,
0.300449963358,
0.298738887130,
0.297056645891,
0.295402444421,
0.293775517417,
0.292175128099,
0.290600566886,
0.289051150151,
0.287526219036,
0.286025138332,
0.284547295425,
0.283092099289,
0.281658979534,
0.280247385507,
0.278856785434,
0.277486665612,
0.276136529634,
0.274805897662,
0.273494305728,
0.272201305079,
0.270926461544,
0.269669354940,
0.268429578504,
0.267206738352,
0.266000452963,
0.264810352693,
0.263636079304,
0.262477285521,
0.261333634607,
0.260204799960,
0.259090464723,
0.257990321423,
0.256904071611,
0.255831425529,
0.254772101794,
0.253725827081,
0.252692335841,
0.251671370013,
0.250662678757,
0.249666018201,
0.248681151189,
0.247707847053,
0.246745881381
};

static double BV1[] = 
{
-0.333083484239,
-0.314029396598,
-0.296048129678,
-0.279291029637,
-0.263667620959,
-0.249093674909,
-0.235492604639,
-0.222793065663,
-0.210929655631,
-0.199841234672,
-0.189472361134,
-0.179771155739,
-0.170689554158,
-0.162183863228,
-0.154213386949,
-0.146740132601,
-0.139729490312,
-0.133149181125,
-0.126969225689,
-0.121162026412,
-0.115702200598,
-0.110565888223,
-0.105731370437,
-0.101177983214,
-0.096887184942,
-0.092841323629,
-0.089024372318,
-0.085421220763,
-0.082017661884,
-0.078801080661,
-0.075759388651,
-0.072881190472,
-0.070156201792,
-0.067574825957,
-0.065127973736,
-0.062807057222,
-0.060604619989,
-0.058513348888,
-0.056526338665,
-0.054637396552,
-0.052840626041,
-0.051130555734,
-0.049501813537,
-0.047950205750,
-0.046470626569,
-0.045059331544,
-0.043712351849,
-0.042425847890,
-0.041196780765,
-0.040021474523,
-0.038897433219,
-0.037821587695,
-0.036791361624,
-0.035804259113,
-0.034858060327,
-0.033950585676,
-0.033079933681,
-0.032243904768,
-0.031440989201,
-0.030669361550,
-0.029927438749,
-0.029213706843,
-0.028527046225,
-0.027865725005,
-0.027228778232,
-0.026615053272,
-0.026023175448,
-0.025452424279,
-0.024901543272,
-0.024369771034,
-0.023856424392,
-0.023360330337,
-0.022880700612,
-0.022417102286,
-0.021968485531,
-0.021534407646,
-0.021114203053,
-0.020707328408,
-0.020313304483,
-0.019931175432,
-0.019560677023,
-0.019201678192,
-0.018853195283,
-0.018515059763,
-0.018186979533,
-0.017868090341,
-0.017558622278,
-0.017257696772,
-0.016965208183,
-0.016680907100,
-0.016404384851,
-0.016135353664,
-0.015873695303,
-0.015618694393,
-0.015370675122,
-0.015129040354,
-0.014893623779,
-0.014664282233,
-0.014440722951,
-0.014222639566,
-0.014010089644,
-0.013802692655,
-0.013600457371,
-0.013403029742,
-0.013210458893,
-0.013022146195,
-0.012838569485,
-0.012659123416,
-0.012483746092,
-0.012312425524,
-0.012144877104,
-0.011981333758,
-0.011821411072,
-0.011664726629,
-0.011511637421,
-0.011361829921,
-0.011215356754,
-0.011071901285,
-0.010931455852,
-0.010794075300,
-0.010659485239,
-0.010527584082,
-0.010398323924,
-0.010271908768,
-0.010148013097,
-0.010026364905,
-0.009907257403,
-0.009790487282,
-0.009676031837,
-0.009563681208
};

static double BV2[] = 
{
0.2,
0.186267510056,
0.173525020480,
0.161702737212,
0.150888219476,
0.140711545944,
0.131408646703,
0.122674599290,
0.114654257894,
0.107168734074,
0.100303336978,
0.093838199973,
0.087862208486,
0.082297787070,
0.077191233635,
0.072350539267,
0.067897558212,
0.063754886389,
0.059819526970,
0.056213244796,
0.052783697844,
0.049825832248,
0.046847946942,
0.044178150594,
0.041610062122,
0.039291560650,
0.037115328014,
0.035007491708,
0.033096984029,
0.031301990151,
0.029569044709,
0.027992703021,
0.026483975351,
0.025121778250,
0.023822709918,
0.022554196417,
0.021430604160,
0.020408257842,
0.019436269999,
0.018461458385,
0.017525754869,
0.016680255532,
0.015884079039,
0.015182934701,
0.014487169683,
0.013783015311,
0.013151951134,
0.012614540756,
0.012009918690,
0.011494055390,
0.010999843478,
0.010523445904,
0.010058835149,
0.009696774185,
0.009282074869,
0.008880652487,
0.008529700339,
0.008201830089,
0.007851548493,
0.007562644780,
0.007313996553,
0.006974913180,
0.006760381162,
0.006470605731,
0.006211631000,
0.006053820252,
0.005838274956,
0.005636543036,
0.005381099880,
0.005235254765,
0.005030728877,
0.004868466407,
0.004792410880,
0.004556208849,
0.004373945296,
0.004247587174,
0.004115894437,
0.003968000412,
0.003903977573,
0.003759298474,
0.003627311438,
0.003508456051,
0.003459598869,
0.003338880837,
0.003252103925,
0.003159917891,
0.003065265715,
0.002996873111,
0.002863284200,
0.002791732550,
0.002723764628,
0.002649065107,
0.002632074058,
0.002488613129,
0.002448838204,
0.002396479249,
0.002363391221,
0.002239506692,
0.002223722637,
0.002126548439,
0.002061165869,
0.002042524517,
0.002003166825,
0.001977387816,
0.001901797950,
0.001861613244,
0.001794632524,
0.001783378422,
0.001734465361,
0.001681677997,
0.001621190459,
0.001608580351,
0.001593299210,
0.001528590918,
0.001497741789,
0.001457422972,
0.001463249326,
0.001391887665,
0.001410230994,
0.001362800598,
0.001330845058,
0.001318227500,
0.001276485622,
0.001248061657,
0.001223806292,
0.001206837595,
0.001193787903,
0.001179732382,
0.001132018864,
0.001088831574,
};
    int a = int(x/0.1);

    double d = x - a * 0.1;

    double D3 = (BV2[a+1]-BV2[a])*10.0;

    if(d > 0.05)
    {
        d = 0.1-d;
        return BV[a+1] - BV1[a+1]*d + BV2[a+1]/2.0*d*d - D3/6.0*d*d*d;
    }
    return BV[a] + BV1[a]*d + BV2[a]/2.0*d*d + D3/6.0*d*d*d;
}

__device__ double boysBig(double x)
{
    if(x > 50000)
        return M_SQRTPI/(2*sqrt(x));
    return M_SQRTPI/(2*sqrt(x)) - 1.0e-08*log(x);   // dla x > 13
}


__device__ double boys(double x) 
{   
    if(x == 0)
        return 1.0;

    if(x >= 12.89)
        return boysBig(x);
    else if(x > 0.05)
        return boysSmall(x);
    else
        return boysSmall(x) - x*2.5e-04;
}

__global__ void countSElement(double * primitives, double * SMat, int primitivesNum)
{
    int row = blockIdx.x*blockDim.x + threadIdx.x;
    int collumn = blockIdx.y*blockDim.y + threadIdx.y;
    if(row > collumn)
        return;
    if(row >= primitivesNum || collumn >= primitivesNum)
        return;
    int matPos1 = row * primitivesNum + collumn;
    int matPos2 = collumn * primitivesNum + row;
    double b = primitives[5*row] + primitives[5*collumn];
    double B = primitives[5*row] * primitives[5*collumn];
    double A = 4.0 * B / (b*b);
    double x = primitives[5*row+2] - primitives[5*collumn + 2];
    double y = primitives[5*row+3] - primitives[5*collumn + 3];
    double z = primitives[5*row+4] - primitives[5*collumn + 4];
    double R2 = x*x + y*y + z*z;
    double res = pow(A, 0.75)* exp(-B*R2/b);
    SMat[matPos1] = res;
    SMat[matPos2] = res;
}

__device__ double countDPrimElement(double * primitives, double *SMat, int pnum, int p, int q, int r, int s)
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

	return 2.0/M_SQRTPI * (sqrt(A*B))/sqrt(AB) * boys(A*B/AB * r2)*Spqrs;

}

__global__ void countDElement(int * info, double * primitives, double *SMat, double * Mat, int num, int pnum)
{
    int index = blockIdx.x*blockDim.x + threadIdx.x;
    int num2 = num*num;
    int num3 = num2*num;
    int num4 = num3*num;
    
    if( index > num4)
        return;

    int value = index;
    int s = value % num; 
    value /= num;
    int q = value % num;
    value /= num;
    int r = value % num;
    value /= num;
    int p = value % num;
    if(p < q || s < r)
        return;

    int index2 = q*num3 + r*num2+ p*num + s;
    int index3 = p*num3 + s*num2+ q*num + r;
    int index4 = q*num3 + s*num2+ p*num + r;

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
					double r = countDPrimElement(primitives, SMat, pnum, i, j, k, l);
					res +=primitives[i*5 + 1] * primitives[j*5 + 1] * primitives[k*5 + 1] * primitives[l*5 + 1] * r;
					
				}
			}
		}          
	}
	Mat[index] = res;
    Mat[index2] = res;
    Mat[index3] = res;
    Mat[index4] = res;
}


__host__
void calculateDIntegrals(double * orbitalData, int * contractedData, double * integralData, unsigned int baseSize, unsigned int primitiveSize)
{
    double * dOribtalData = nullptr;
    int * dContractedData = nullptr;
    double * d_S_Mat = nullptr;
    double * dIntegralData = nullptr;

    cudaMalloc((void**)&dOribtalData, primitiveSize * sizeof(double));
    cudaMalloc((void**)&dContractedData, baseSize * sizeof(int));
    cudaMalloc((void**)&d_S_Mat, primitiveSize*primitiveSize*sizeof(double));
    cudaMalloc((void**)&dIntegralData, baseSize * baseSize * baseSize * baseSize *sizeof(double));

    cudaMemcpy(dOribtalData, orbitalData, primitiveSize * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dContractedData, contractedData, baseSize * sizeof(int), cudaMemcpyHostToDevice);

    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid((primitiveSize + BLOCK_SIZE - 1)/BLOCK_SIZE, (primitiveSize + BLOCK_SIZE - 1)/BLOCK_SIZE);
    countSElement<<<dimGrid, dimBlock>>>(dOribtalData, d_S_Mat, primitiveSize);
    cudaDeviceSynchronize();

    int numb = BLOCK_SIZE * BLOCK_SIZE;
    int numg = (baseSize * baseSize * baseSize * baseSize + BLOCK_SIZE * BLOCK_SIZE - 1)/(BLOCK_SIZE * BLOCK_SIZE);
    countDElement<<<numg,numb>>>(dContractedData, orbitalData, d_S_Mat, dIntegralData, baseSize, primitiveSize);
    cudaDeviceSynchronize();

    cudaMemcpy(integralData, dOribtalData, baseSize * baseSize * baseSize * baseSize  * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dOribtalData);
    cudaFree(dContractedData);
    cudaFree(d_S_Mat);
    cudaFree(dIntegralData);
}