//---------------------------------------------------------------------------------------
/*
	Title : Uncertainty Quantification in Linear Interpolation for Isosurface Extraction	
       Authors: Tushar Athawale, Alireza Entezari
        Date  : Jun 27, 2013.
*/
//---------------------------------------------------------------------------------------

#include <viskores/Math.h>
#include <cmath>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <eigen3/Eigen/Dense>


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



class gaussian_distribution{

	private :
		

	public :
		
		// alpha distribution parameters returned assuming data is sampled from a gaussian distribution
		VISKORES_EXEC void gaussian_alpha_pdf(viskores::Float64 muX, viskores::Float64 varX, viskores::Float64 muY, viskores::Float64 varY, viskores::Float64 covarXY, viskores::Float64 c, viskores::Float64* expt, viskores::Float64* cross_prob, viskores::Float64* var);
    
        // Piecewise function returned assuming data is sampled from a kernel density estimation
        VISKORES_EXEC void gaussian_alpha_pdf_MonteCarlo(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var, int nn);

		// helper function for gaussian_alpha_pdf
		VISKORES_EXEC void computeAlphaUncertainty(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 rhoNumDenom, viskores::Float64* expt, viskores::Float64* cross_prob, viskores::Float64* var);

		// helper function for gaussian_alpha_pdf: utilization of the Hinkley's derivation
		VISKORES_EXEC viskores::Float64 hinkleyGaussianRatioCorrelated(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 rhoNumDenom, viskores::Float64 alphaVal);
		
		// helper function for gaussian_alpha_pdf: utilization of the Hinkley's derivation	
		VISKORES_EXEC viskores::Float64 hinkleyGaussianRatioIndependent(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 alphaVal);

		// helper function for gaussian_alpha_pdf: normal cumulative distribution	
		VISKORES_EXEC viskores::Float64 normalCumulativeDistribution(viskores::Float64 t, viskores::Float64 mu, viskores::Float64 sig);

};


// Piecewise function returned assuming data is sampled from a kernel density estimation
VISKORES_EXEC inline void gaussian_distribution::gaussian_alpha_pdf(viskores::Float64 muX, viskores::Float64 varX, viskores::Float64 muY, viskores::Float64 varY, viskores::Float64 covarXY, viskores::Float64 c, viskores::Float64* expt, viskores::Float64* cross_prob, viskores::Float64* var)
{
	//cout<<"Isovalue:"<<c<<"\n";

	//cout<<"muX:"<<muX<<"muY:"<<muY<<"\n";

        //cout<<"varX:"<<varX<<"covarXY:"<<covarXY<<"\n";
	// Standard devition c-X
	viskores::Float64 sigC_X=viskores::Sqrt(varX);
	//cout<<"sig C-X:"<<sigC_X<<"\n";

	// Standard deviation Y-X: This can lead to nan if value to square root is negative
  // Avoid square root of negative value (which should be highly unlikely)
    
  viskores::Float64 sig_Y_X;
    
  if (varX + varY -2*covarXY > 0)
    sig_Y_X = viskores::Sqrt(varX + varY -2*covarXY);
  else
    sig_Y_X = viskores::Sqrt(0.00000001);
	//cout<<"sig Y-X:"<<sig_Y_X<<"\n";

	// Compute correlation between c-X and Y-X
	viskores::Float64 pearsonCorrelation = (varX - covarXY)/(sigC_X*sig_Y_X);
    //cout<<"rhoXY:"<<pearsonCorrelation<<"\n";
	//cout<<"rhoXY:"<<(double)(covarXY/(viskores::Sqrt(varX)*viskores::Sqrt(varY)))<<"\n";

	// The Pearson's correlation can go above 1 because of numerical instability if there is perfect correlation between X and Y
	if (pearsonCorrelation > 1)
		pearsonCorrelation = 0.99;

	if (pearsonCorrelation < -1)
		pearsonCorrelation = -0.99;
	
	// compute expected value, crossing probability, and variance of alpha distribution using the Hinkley's approach
	computeAlphaUncertainty(c - muX, muY - muX, sigC_X, sig_Y_X, pearsonCorrelation, expt, cross_prob, var);

        //cout<<"(Closed) expected val:"<<*expt<<"var:"<<*var<<"\n";
}


VISKORES_EXEC inline void gaussian_distribution::computeAlphaUncertainty(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 rhoNumDenom, viskores::Float64* expt, viskores::Float64* cross_prob, viskores::Float64* var)
{

	// Compute alpha probability at 100 points
	constexpr viskores::IdComponent NUM_POINTS = 100;
	viskores::Float64* alphaDensity = new viskores::Float64[NUM_POINTS];
	viskores::Float64 alphaSum = 0;	
	viskores::Float64 alphaExpected = 0;
	viskores::Float64 alphaVar = 0;	
	
	for (viskores::IdComponent i=0; i<NUM_POINTS; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		viskores::Float64 alphaVal = (viskores::Float64)i*(1.0/100);
		viskores::Float64 density = hinkleyGaussianRatioCorrelated(muNumerator, muDenominator, sigNumerator, sigDenominator, rhoNumDenom, alphaVal);
		//cout<<"analytical density:"<<density<<"\n";		
		alphaDensity[i] = density;
		// total density over [0,1] for normalization
		alphaSum = alphaSum + density;
	}

	// Compute expected value
	for (viskores::IdComponent i=0; i<100; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		viskores::Float64 alphaVal = (viskores::Float64)(i)*(1.0/100);
		// Use normalized density
		alphaExpected = alphaExpected + (viskores::Float64)(alphaDensity[i]/alphaSum)*alphaVal;
	}
	*expt = alphaExpected;
	//cout<<"alpha expected:"<<alphaExpected<<"\n";
	
	// Compute variance
	for (viskores::IdComponent i=0; i<100; i++) 
	{	
		// we sample only [0,1] range for marching cubes			
		viskores::Float64 alphaVal = (viskores::Float64)(i)*(1.0/100);
		alphaVar = alphaVar + (viskores::Float64)(alphaDensity[i]/alphaSum)*viskores::Pow((alphaVal - *expt),2);
	}
	//cout<<"alphaVar:"<<alphaVar<<"\n";
	*expt = alphaExpected;
	// double check cross prob. Currently, a placeholder.
	*cross_prob = alphaSum;
	*var = alphaVar;
}

// Piecewise function returned assuming data is sampled from a kernel density estimation
VISKORES_EXEC inline void gaussian_distribution::gaussian_alpha_pdf_MonteCarlo(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var, int nn)
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
  //    std::cout << alphaArr[i] <<" ";
  //std::cout<<"\n";

 // Create histogram  of inverse linear interpolation (following code picked from web)

  std::sort(alphaArr.begin(), alphaArr.end());

  std::map<double, int> histogram;

  double bin = 0; //Choose your starting bin
  const double bin_width = 0.001; //Choose your bin interval
  for (int e=0; e<alphaArr.size(); e++)
  {
      while (alphaArr[e] >= (bin + bin_width))
        bin += bin_width;
    ++histogram[bin];
  }

  std::map<double, int>::iterator it;
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





VISKORES_EXEC inline viskores::Float64 gaussian_distribution::hinkleyGaussianRatioCorrelated(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 rhoNumDenom, viskores::Float64 alphaVal)
{
	//Compute parameters of transformed random variable that make numertor and denominator uncorrelated
	viskores::Float64 muNumeratorDash = muNumerator - (viskores::Float64)(rhoNumDenom*muDenominator*sigNumerator/sigDenominator);
	//cout<<"muNumeratorDash:"<<muNumeratorDash<<"\n";
	viskores::Float64 sigNumeratorDash = sigNumerator*viskores::Sqrt(1.0 - rhoNumDenom*rhoNumDenom);
	//cout<<"sigNumeratorDash:"<<sigNumeratorDash<<"\n";
	viskores::Float64 offset = (viskores::Float64)(rhoNumDenom*sigNumerator/sigDenominator);
	//cout<<"offset:"<<offset<<"\n";
	viskores::Float64 alphaDensity = hinkleyGaussianRatioIndependent(muNumeratorDash, muDenominator, sigNumeratorDash, sigDenominator, alphaVal -offset);
	//cout<<"alphaDensity:"<<alphaDensity<<"\n";
	return alphaDensity;
}


VISKORES_EXEC inline viskores::Float64 gaussian_distribution::hinkleyGaussianRatioIndependent(viskores::Float64 muNumerator, viskores::Float64 muDenominator, viskores::Float64 sigNumerator, viskores::Float64 sigDenominator, viskores::Float64 alphaVal)
{
	// if sigNumerator or sigDenominator are zero, a divide by zero error can occur, so add small epsilon
	if (sigNumerator < 0.0001)
		sigNumerator = 0.0001;
	if (sigDenominator < 0.0001)
		sigDenominator = 0.0001;

	viskores::Float64 az = viskores::Sqrt((viskores::Float64)(viskores::Pow(alphaVal,2)/viskores::Pow(sigNumerator,2)) + (viskores::Float64)(1.0/viskores::Pow(sigDenominator,2)));
	//cout<<"az:"<<az<<"\n";
	viskores::Float64 bz = (viskores::Float64)(muNumerator*alphaVal/viskores::Pow(sigNumerator,2)) + (viskores::Float64)(muDenominator/viskores::Pow(sigDenominator,2));
	//cout<<"bz:"<<bz<<"\n";
	viskores::Float64 c = (viskores::Float64)(viskores::Pow(muNumerator,2)/viskores::Pow(sigNumerator,2)) + (viskores::Float64)(viskores::Pow(muDenominator,2)/viskores::Pow(sigDenominator,2));
	//cout<<"c:"<<c<<"\n";
	viskores::Float64 dz = viskores::Exp((viskores::Pow(bz,2) - c*viskores::Pow(az,2))/(2*viskores::Pow(az,2)));
	//cout<<"dz:"<<dz<<"\n";

	viskores::Float64 t1 = normalCumulativeDistribution((viskores::Float64)(bz/az),0,1);
	viskores::Float64 t2 = normalCumulativeDistribution((viskores::Float64)(-bz/az),0,1);

	viskores::Float64 alphaDensity = (viskores::Float64)((bz*dz)/viskores::Pow(az,3))*(viskores::Float64)(1.0/(viskores::Sqrt(2*viskores::Pi())*sigNumerator*sigDenominator))*(t1-t2) + (viskores::Float64)(1.0/(viskores::Pow(az,2)*viskores::Pi()*sigNumerator*sigDenominator)*viskores::Exp((viskores::Float64)(-c/2))); 

	return alphaDensity;
}


VISKORES_EXEC inline viskores::Float64 gaussian_distribution::normalCumulativeDistribution(viskores::Float64 t, viskores::Float64 mu, viskores::Float64 sig)
{
	return 0.5*(1 + std::erf((viskores::Float64)((t-mu)/(viskores::Sqrt(2)*sig))));
}
