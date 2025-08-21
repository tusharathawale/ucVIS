//---------------------------------------------------------------------------------------
/*
	Title : Uncertainty Quantification in Linear Interpolation for Isosurface Extraction	
       Authors: Tushar Athawale, Alireza Entezari
        Date  : Jun 27, 2013.
*/
//---------------------------------------------------------------------------------------

#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"

# define infinity 10000
# define minus_infinity -10000

using namespace std;

class gaussian_distribution{

	private :
		

	public :
		
		// alpha distribution parameters returned assuming data is sampled from a gaussian distribution
		void gaussian_alpha_pdf(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var);

		// helper function for gaussian_alpha_pdf
		void computeAlphaUncertainty(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double rhoNumDenom, double* expt, double* cross_prob, double* var);

		// helper function for gaussian_alpha_pdf: utilization of the Hinkley's derivation
		double hinkleyGaussianRatioCorrelated(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double rhoNumDenom, double alphaVal);
		
		// helper function for gaussian_alpha_pdf: utilization of the Hinkley's derivation	
		double hinkleyGaussianRatioIndependent(double muNumerator, double muDenominator, double sigNumerator, double sigDenominator, double alphaVal);

		// helper function for gaussian_alpha_pdf: normal cumulative distribution	
		double normalCumulativeDistribution(double t, double mu, double sig);

		// Ilerp distribution using Monte Carlo Gaussian sampling
		void gaussian_alpha_pdf_MonteCarlo(double muX, double varX, double muY, double varY, double covarXY, double c, double* expt, double* cross_prob, double* var, int numSamples);

};

