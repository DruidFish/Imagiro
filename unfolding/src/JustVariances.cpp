/**
  @class JustVariances

  A class to calculate the errors on the unfolded data distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "JustVariances.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Default constructor
JustVariances::JustVariances()
{
}

//Constructor with correct arguments
JustVariances::JustVariances( UnfoldingMatrix * Unfolding, SmearingMatrix * Smearing, Distribution * Data, double TrueEventNumber, Indices * IndexCalculator ) : inputSmearing(Smearing),
	inputUnfolding(Unfolding), inputDistribution(Data), inputIntegral(TrueEventNumber)
{
	binNumber = IndexCalculator->GetBinNumber();
}

//Run the calculation
void JustVariances::Calculate()
{
	//Precalculate some stuff
	fullTerms = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );
	for ( int covarianceIndexK = 0; covarianceIndexK < binNumber; covarianceIndexK++ )
	{
		for ( int effectIndexI = 0; effectIndexI < binNumber; effectIndexI++ )
		{
			fullTerms[covarianceIndexK][effectIndexI] = UnfoldingCovarianceTerm( covarianceIndexK, effectIndexI );
		}
	}

	//Calculate each element of the matrix
	for ( int covarianceIndexK = 0; covarianceIndexK < binNumber; covarianceIndexK++ )
	{
		double element = 0.0;

		//Calculate the first part of the data contribution
		for ( int effectIndexJ = 0; effectIndexJ < binNumber; effectIndexJ++ )
		{
			double dataNumber = inputDistribution->GetBinNumber(effectIndexJ);
			double firstPartData = inputUnfolding->GetElement( covarianceIndexK, effectIndexJ );
			firstPartData *= firstPartData;
			firstPartData *= dataNumber;
			firstPartData *= ( 1 - ( dataNumber / inputIntegral ) );

			element += firstPartData;
		}


		for ( int effectIndexI = 0; effectIndexI < binNumber; effectIndexI++ )
		{
			//Retrieve values for index I
			double unfoldingKI = inputUnfolding->GetElement( covarianceIndexK, effectIndexI );
			double dataI = inputDistribution->GetBinNumber(effectIndexI);

			for ( int effectIndexJ = 0; effectIndexJ < binNumber; effectIndexJ++ )
			{
				//Retrieve values from the data distribution
				double dataMultiply = dataI * inputDistribution->GetBinNumber(effectIndexJ);

				if ( effectIndexI != effectIndexJ )
				{
					//Calculate the second part of the data contribution
					double secondPartData = unfoldingKI;
					secondPartData *= inputUnfolding->GetElement( covarianceIndexK, effectIndexJ );
					secondPartData *= dataMultiply;
					secondPartData /= inputIntegral;

					element -= secondPartData;
				}

				//Calculate the smearing matrix contribution
				double smearingPart = dataMultiply;
				smearingPart *= UnfoldingCovariance( covarianceIndexK, effectIndexI, effectIndexJ );

				element += smearingPart;
			}
		}

		variances.push_back(element);
	}
}

//Destructor
JustVariances::~JustVariances()
{
}

//Retrieve a calculated variance
double JustVariances::GetVariance( int Index )
{
	if ( variances.size() == 0 )
	{
		Calculate();
	}

	return variances[Index];
}
double JustVariances::GetStandardDeviation( int Index )
{
	if ( variances.size() == 0 )
	{
		Calculate();
	}

	return sqrt( variances[Index] );
}

//Weird summation
double JustVariances::UnfoldingCovariance( int CauseIndexK, int EffectIndexI, int EffectIndexJ )
{
	double result = 0.0;
	double unfoldingKI = inputUnfolding->GetElement( CauseIndexK, EffectIndexI );
	double unfoldingKJ = inputUnfolding->GetElement( CauseIndexK, EffectIndexJ );

	for ( int indexU = 0; indexU < binNumber; indexU++ )
	{
		//Retrieve the efficiency for this U value
		double uEfficiency = inputSmearing->GetEfficiency(indexU);

		if ( indexU == CauseIndexK )
		{
			//If U==K==L, true for all R and S
			for ( int indexR = 0; indexR < binNumber; indexR++ )
			{
				//Calculate the unfolding covariance for this R
				double unfoldingRTerm;
				if ( indexR == EffectIndexI )
				{
					unfoldingRTerm = fullTerms[CauseIndexK][EffectIndexI];
				}
				else
				{
					unfoldingRTerm = -unfoldingKI / uEfficiency;

					//Check for divide by zero errors
					if ( isinf(unfoldingRTerm) )
					{
						cerr << "VARIANCE: Divide by zero error (efficiency)" << endl;
						exit(1);
					}
					if ( isnan(unfoldingRTerm) )
					{
						unfoldingRTerm = 0.0;
					}
				}

				//Loop over all S
				for ( int indexS = 0; indexS < binNumber; indexS++ )
				{
					double term = inputSmearing->GetCovarianceElement( indexR, indexS, indexU ) * unfoldingRTerm;

					//Calculate the unfolding covariance for this S
					if ( indexS == EffectIndexJ )
					{
						term *= fullTerms[CauseIndexK][EffectIndexJ];
					}
					else
					{
						term *= -unfoldingKJ / uEfficiency;

						//Check for divide by zero errors
						if ( isinf(term) ) 
						{
							cerr << "VARIANCE: Divide by zero error (efficiency)" << endl;
							exit(1);
						}       
						if ( isnan(term) )
						{
							term = 0.0;
						}
					}

					result += term;
				}
			}
		}
		else
		{
			//If U!=K, true only for R==I and S==J
			double term = inputSmearing->GetCovarianceElement( EffectIndexI, EffectIndexJ, indexU );
			term *= unfoldingKI * inputUnfolding->GetCovarianceTerm( EffectIndexI, indexU );
			term *= unfoldingKJ * inputUnfolding->GetCovarianceTerm( EffectIndexJ, indexU );

			result += term;
		}
	}

	return result;
}


//The full expression, assuming all Kroneker deltas are true
double JustVariances::UnfoldingCovarianceTerm( int IndexK, int IndexI )
{
	double unfoldingKI = inputUnfolding->GetElement( IndexK, IndexI );
	if ( unfoldingKI != 0.0 )
	{
		//Calculate first term to sum
		double sumTerm1 = 0.0;
		double smearingKI = inputSmearing->GetElement( IndexK, IndexI );
		if ( smearingKI != 0.0 )
		{
			sumTerm1 = 1.0 / smearingKI;
		}

		//Calculate second term to sum
		double sumTerm2 = 0.0;
		double efficiencyK = inputSmearing->GetEfficiency(IndexK);
		if ( efficiencyK != 0.0 )
		{
			sumTerm2 = 1.0 / efficiencyK;
		}

		double term = sumTerm1 - sumTerm2 - inputUnfolding->GetCovarianceTerm( IndexI, IndexK );
		return term * unfoldingKI;
	}
	else
	{
		return 0.0;
	}
}

//Return vectors with all values
vector<double> JustVariances::GetVariances()
{
	if ( variances.size() == 0 )
	{
		Calculate();
	}

	return variances;
}
vector<double> JustVariances::GetStandardDeviations()
{
	if ( variances.size() == 0 )
	{
		Calculate();
	}

	vector<double> deviations;

	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		deviations.push_back( sqrt( variances[binIndex] ) );
	}

	return deviations;
}
