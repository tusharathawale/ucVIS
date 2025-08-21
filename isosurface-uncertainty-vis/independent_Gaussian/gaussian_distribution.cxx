#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <eigen3/Eigen/Dense>
//#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "gaussian_distribution.h"

# define infinity 10000
# define minus_infinity -10000
#define _USE_MATH_DEFINES

/*
  We need a functor that can pretend it's const,
  but to be a good random number generator 
  it needs mutable state.
*/
namespace Eigen {
namespace internal {
template<typename Scalar> 
struct scalar_normal_dist_op 
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

using namespace std;

double gaussian_distribution::normalCumulativeDistribution(double t, double mu, double sig)
{
	return 0.5*(1 + erf((double)((t-mu)/(sqrt(2)*sig))));
}

double gaussian_distribution::hinkleyGaussianRatioIndependent(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double alphaVal)
{
	// if sigNumerator or sigDenominator are zero, a divide by zero error can occur, so add small epsilon
	if (sigNumerator < 0.0001)
		sigNumerator = 0.0001;
	if (sigDenominator < 0.0001)
		sigDenominator = 0.0001;

	double az = sqrt((double)(pow(alphaVal,2)/pow(sigNumerator,2)) + (double)(1.0/pow(sigDenominator,2)));
	//cout<<"az:"<<az<<"\n";
	double bz = (double)(muNumerator*alphaVal/pow(sigNumerator,2)) + (double)(muDenominator/pow(sigDenominator,2));
	//cout<<"bz:"<<bz<<"\n";
	double c = (double)(pow(muNumerator,2)/pow(sigNumerator,2)) + (double)(pow(muDenominator,2)/pow(sigDenominator,2));
	//cout<<"c:"<<c<<"\n";
	double dz = exp((pow(bz,2) - c*pow(az,2))/(2*pow(az,2)));
	//cout<<"dz:"<<dz<<"\n";

	double t1 = normalCumulativeDistribution((double)(bz/az),0,1);
	double t2 = normalCumulativeDistribution((double)(-bz/az),0,1);

	double alphaDensity = (double)((bz*dz)/pow(az,3))*(double)(1.0/(sqrt(2*M_PI)*sigNumerator*sigDenominator))*(t1-t2) + (double)(1.0/(pow(az,2)*M_PI*sigNumerator*sigDenominator)*exp((double)(-c/2))); 

	return alphaDensity;

}

double gaussian_distribution::hinkleyGaussianRatioCorrelated(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double rhoNumDenom, double alphaVal)
{
	//Compute parameters of transformed random variable that make numertor and denominator uncorrelated
	double muNumeratorDash = muNumerator - (double)(rhoNumDenom*muDenominator*sigNumerator/sigDenominator);
	//cout<<"muNumeratorDash:"<<muNumeratorDash<<"\n";
	double sigNumeratorDash = sigNumerator*sqrt(1.0 - rhoNumDenom*rhoNumDenom);
	//cout<<"sigNumeratorDash:"<<sigNumeratorDash<<"\n";
	double offset = (double)(rhoNumDenom*sigNumerator/sigDenominator);
	//cout<<"offset:"<<offset<<"\n";
	double alphaDensity = hinkleyGaussianRatioIndependent(muNumeratorDash, muDenominator, sigNumeratorDash, sigDenominator, alphaVal -offset);
	//cout<<"alphaDensity:"<<alphaDensity<<"\n";
	return alphaDensity;
}

void gaussian_distribution::computeAlphaUncertainty(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double rhoNumDenom, double* expt, double* cross_prob, double* var)
{

	// Compute alpha probability at 100 points
	double* alphaDensity = new double[100];
	double alphaSum = 0;	
	double alphaExpected = 0;
	double alphaVar = 0;	
	
	for (int i=0; i<100; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		double alphaVal = (double)i*(1.0/100);
		double density = hinkleyGaussianRatioCorrelated(muNumerator, muDenominator, sigNumerator, sigDenominator, rhoNumDenom, alphaVal);
		//cout<<"analytical density:"<<density<<"\n";		
		alphaDensity[i] = density;
		// total density over [0,1] for normalization
		alphaSum = alphaSum + density;
	}

	// Compute expected value
	for (int i=0; i<100; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		double alphaVal = (double)(i)*(1.0/100);
		// Use normalized density
		alphaExpected = alphaExpected + (double)(alphaDensity[i]/alphaSum)*alphaVal;
	}
	*expt = alphaExpected;
	//cout<<"alpha expected:"<<alphaExpected<<"\n";
	
	// Compute variance
	for (int i=0; i<100; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		double alphaVal = (double)(i)*(1.0/100);
		alphaVar = alphaVar + (double)(alphaDensity[i]/alphaSum)*pow((alphaVal - *expt),2);
	}
	//cout<<"alphaVar:"<<alphaVar<<"\n";
	*expt = alphaExpected;
	// double check cross prob. Currently, a placeholder.
	*cross_prob = alphaSum;
	*var = alphaVar;
}


// Piecewise function returned assuming data is sampled from a kernel density estimation
void gaussian_distribution::gaussian_alpha_pdf(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var)
{
	//cout<<"Isovalue:"<<c<<"\n";

	//cout<<"muX:"<<muX<<"muY:"<<muY<<"\n";

        //cout<<"varX:"<<varX<<"covarXY:"<<covarXY<<"\n";
	// Standard devition c-X
	double sigC_X=sqrt(varX);
	//cout<<"sig C-X:"<<sigC_X<<"\n";

	// Standard deviation Y-X
	double sig_Y_X = sqrt(varX + varY -2*covarXY);
	//cout<<"sig Y-X:"<<sig_Y_X<<"\n";

	// Compute correlation between c-X and Y-X
	double pearsonCorrelation = (double)((varX - covarXY)/(sigC_X*sig_Y_X));
        //cout<<"rhoXY:"<<pearsonCorrelation<<"\n";
	//cout<<"rhoXY:"<<(double)(covarXY/(sqrt(varX)*sqrt(varY)))<<"\n";

	// The Pearson's correlation can go above 1 because of numerical instability if there is perfect correlation between X and Y
	if (pearsonCorrelation > 1)
		pearsonCorrelation = 1;

	if (pearsonCorrelation < -1)
		pearsonCorrelation = -1;
	
	// compute expected value, crossing probability, and variance of alpha distribution using the Hinkley's approach
	computeAlphaUncertainty(c - muX, muY - muX, sigC_X, sig_Y_X, pearsonCorrelation, expt, cross_prob, var);

        //cout<<"(Closed) expected val:"<<*expt<<"var:"<<*var<<"\n";
}



// Piecewise function returned assuming data is sampled from a kernel density estimation
void gaussian_distribution::gaussian_alpha_pdf_MonteCarlo(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var, int nn)
{

  int size = 2; // Dimensionality (rows)
  //int nn=10;     // How many samples (columns) to draw
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng

  // Define mean and covariance of the distribution
  Eigen::VectorXd mean(size);       
  Eigen::MatrixXd covar(size,size);

  mean  <<  muX,  muY;
  covar <<  varX, covarXY,
           covarXY,  varY;

  Eigen::MatrixXd normTransform(size,size);

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if 
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors() 
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform 
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise() 
                           + mean;

  //std::cout << "Mean\n" << mean << std::endl;
  //std::cout << "Covar\n" << covar << std::endl;
  //std::cout << "Samples\n" << samples << std::endl;
  //std::cout << "Samples[0][0]\n" << samples.coeff(0,0) << std::endl;
  //std::cout << "Samples[1][0]\n" << samples.coeff(1,0) << std::endl;

  // Calculate inverse linear interpolation 
  std::vector<double> alphaArr;
  for (int i=0; i<nn; i++) 
  {
	double alphaVal;

	// Avoid division by 0
	if (abs(samples.coeff(1,i)-samples.coeff(0,i)) >= 0.001)
	{
		alphaVal = (double)(c-samples.coeff(0,i))/(samples.coeff(1,i)-samples.coeff(0,i));

		// Consider [0,1] domain
		if ((alphaVal >= 0) && (alphaVal<=1))
		{
			alphaArr.push_back(alphaVal);
		}
  	}
  }

  //for (int i=0; i<alphaArr.size(); i++)
  //	std::cout << alphaArr[i] <<" ";
  //std::cout<<"\n";

 // Create histogram  of inverse linear interpolation (following code picked from web)

  std::sort(alphaArr.begin(), alphaArr.end());

  map<double, int> histogram;

  double bin = 0; //Choose your starting bin
  const double bin_width = 0.001; //Choose your bin interval
  for (int e=0; e<alphaArr.size(); e++)
  {
  	while (alphaArr[e] >= (bin + bin_width)) 
		bin += bin_width;
	++histogram[bin];
  }

  map<double, int>::iterator it;
  // Histogram bin centers
  std::vector<double> alphaBin;
  // Histogram bin count
  std::vector<double> alphaBinProb;

  double totalCount = 0;

  for (it=histogram.begin(); it!=histogram.end(); it++)
  {
    //printf("[%.2f,%.2f[ : %d\n", it->first, (it->first) + bin_width, it->second);

    alphaBin.push_back(it->first + (double)(bin_width/2.0));
    alphaBinProb.push_back(it->second);	
    totalCount+=it->second;
  }

  // Normalize the histogram
  for (int i=0; i<alphaBinProb.size(); i++)
  {
  	alphaBinProb[i] = (double)(alphaBinProb[i]/totalCount);
  }

  //for (int i=0; i<alphaBinProb.size(); i++)
  	//std::cout << alphaBinProb[i] <<" ";
  //std::cout<<"\n";
	
  double expt_val = 0;
  //std::cout << "For expectation: "; 
  // Compute the expected alpha value
  for (int i=0; i<alphaBinProb.size(); i++) 
  {	
        //std::cout << alphaBin[i] <<" ";
	//std::cout << alphaBinProb[i] <<"\n";
  	expt_val+=alphaBinProb[i]*alphaBin[i];
  }

  double var_val = 0; 
  // Compute the expected alpha value
  for (int i=0; i<alphaBinProb.size(); i++)
  	var_val+=alphaBinProb[i]*pow(alphaBin[i] - expt_val, 2);

  *expt = expt_val;
  *var = var_val;

  //cout<<"(MC) expected val:"<<expt_val<<"var:"<<var_val<<"\n";
  // Crossing probability yet to be written
}


/*
int main()
{
	//mu1: 2.85714 delta1: 57.1429 mu2: 11.4286 delta2: 114.286
	
	// ratio distribution assuming data is uniformly distributed
	alpha_distribution* d = new alpha_distribution();

	// ratio distribution assuming data is sampled from kde
	alpha_distribution* kde_d = new alpha_distribution();*/
	
/*	// Allocate memory for mu and delta arrays
	double* mu1 = new double[2];
	double* mu2 = new double[2];
	double* delta1 = new double[2];
	double* delta2 = new double[2];

	// kernel 1
	mu1[0] = 5;
	delta1[0] = 3;
	mu1[1] = 9;
	delta1[1] = 1;

	// kernel 2
	mu2[0] = 23;
	delta2[0] = 3;
	mu2[1] = 27;
	delta2[1] = 1;*/

	/*piecewise p1, p2;

	// Assuming data is uniformly distributed
	p1 = d->alpha_pdf(8,7, 9, 4, 12.3);
	p2 = d->getPdfOver0To1(p1);

	double expected, crossProb, var;
	d->Compute0To1(p2, &expected, &crossProb, &var);
	cout<<"Expected Value is:"<<expected<<"\n";
	cout<<"Crossing Probability is:"<<crossProb<<"\n";
	cout<<"Variance is:"<<var<<"\n";

	//p1 = d->alpha_pdf(mu1[0], delta1[0], mu2[0], delta2[0], 10);*/
	/*cout<<p1.numPieces<<"\n";
	for(int i=0;i<=p1.numPieces;i++)
		cout<<p1.limits[i]<<" ";
	cout<<"\n";

	cout<<p2.numPieces<<"\n";
	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.limits[i]<<" ";
	cout<<"\n";*/

	/*// Assuming data is kde sampled
	p2 = kde_d->kde_alpha_pdf(mu1, delta1, mu2, delta2, 10);
	cout<<p2.numPieces<<"\n";
	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.limits[i]<<" ";
	cout<<"\n";

	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.pc[i].type<<" ";
	cout<<"\n ";*/

//	return 0;
//}
