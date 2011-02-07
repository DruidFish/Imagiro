/**
  @class CovarianceMatrix

  A class for calculating the covariance matrix for the unfolded data distribution
  Depreciated now since it's really slow, but it might work

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "CovarianceMatrix.h"
#include <iostream>
#include <cstdlib>

//Default constructor
CovarianceMatrix::CovarianceMatrix()
{
}

//Constructor with correct arguments
CovarianceMatrix::CovarianceMatrix( UnfoldingMatrix * Unfolding, SmearingMatrix * Smearing, Distribution * Data, double TrueEventNumber, Indices * IndexCalculator ) : inputSmearing(Smearing), inputUnfolding(Unfolding)
{
	//Precache some elements of the calculation
	binNumber = IndexCalculator->GetBinNumber();
	MakeUnfoldingCovarianceTerm();

	//Initialise the matrix
	matrix = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );

	//Calculate each element of the matrix
	for ( int covarianceIndexK = 0; covarianceIndexK < binNumber; covarianceIndexK++ )
	{
		//for ( int covarianceIndexL = 0; covarianceIndexL < binNumber; covarianceIndexL++ )
		//{
		//Why not only do the variances?
		int covarianceIndexL = covarianceIndexK;	

		double element = 0.0;

		//Calculate the first part of the data contribution
		for ( int effectIndexJ = 0; effectIndexJ < binNumber; effectIndexJ++ )
		{
			double firstPartData = Unfolding->GetElement( covarianceIndexK, effectIndexJ );
			firstPartData *= Unfolding->GetElement( covarianceIndexL, effectIndexJ );
			firstPartData *= Data->GetBinNumber(effectIndexJ);
			firstPartData *= ( 1 - ( Data->GetBinNumber(effectIndexJ) / TrueEventNumber ) );

			element += firstPartData;
		}

		//Calculate the second part of the data contribution
		for ( int effectIndexI = 0; effectIndexI < binNumber; effectIndexI++ )
		{
			for ( int effectIndexJ = 0; effectIndexJ < binNumber; effectIndexJ++ )
			{
				if ( effectIndexI != effectIndexJ )
				{
					double secondPartData = Unfolding->GetElement( covarianceIndexK, effectIndexI );
					secondPartData *= Unfolding->GetElement( covarianceIndexL, effectIndexJ );
					secondPartData *= ( Data->GetBinNumber(effectIndexI) * Data->GetBinNumber(effectIndexJ) / TrueEventNumber );

					element -= secondPartData;
				}
			}
		}

		//Calculate the smearing matrix contribution
		for ( int effectIndexI = 0; effectIndexI < binNumber; effectIndexI++ )
		{
			for ( int effectIndexJ = 0; effectIndexJ < binNumber; effectIndexJ++ )
			{
				double smearingPart = Data->GetBinNumber(effectIndexI);
				smearingPart *= Data->GetBinNumber(effectIndexJ);
				smearingPart *= UnfoldingCovariance( covarianceIndexK, covarianceIndexL, effectIndexI, effectIndexJ );

				element += smearingPart;
			}
		}

		matrix[covarianceIndexK][covarianceIndexL] = element;
		cout << element << " ";
		//}
	}

	//Free up some memory
	unfoldingCovarianceTerm.clear();
	cout << endl;
}

//Destructor
CovarianceMatrix::~CovarianceMatrix()
{
}

//Retrieve an element of the matrix
double CovarianceMatrix::GetElement( int CauseIndex, int EffectIndex )
{
	return matrix[CauseIndex][EffectIndex];
}

//Weird summation
double CovarianceMatrix::UnfoldingCovariance( int CauseIndexK, int CauseIndexL, int EffectIndexI, int EffectIndexJ )
{
	double result = 0.0;

	if ( CauseIndexK == CauseIndexL )
	{
		for ( int indexU = 0; indexU < binNumber; indexU++ )
		{
			if ( indexU == CauseIndexK )
			{
				//If U==K==L, true for all R and S
				for ( int indexR = 0; indexR < binNumber; indexR++ )
				{
					for ( int indexS = 0; indexS < binNumber; indexS++ )
					{
						double term = inputSmearing->GetCovarianceElement( indexR, indexS, indexU );
						term *= unfoldingCovarianceTerm[CauseIndexK][EffectIndexI][indexR][indexU]; 
						term *= unfoldingCovarianceTerm[CauseIndexL][EffectIndexJ][indexS][indexU];
						//term *= UnfoldingCovarianceTerm( CauseIndexK, EffectIndexI, indexR, indexU ); 
						//term *= UnfoldingCovarianceTerm( CauseIndexL, EffectIndexJ, indexS, indexU );

						result += term;
					}
				}
			}
			else
			{
				//If U!=K==L, true only for R==I and S==J
				double term = inputSmearing->GetCovarianceElement( EffectIndexI, EffectIndexJ, indexU );
				term *= unfoldingCovarianceTerm[CauseIndexK][EffectIndexI][EffectIndexI][indexU]; 
				term *= unfoldingCovarianceTerm[CauseIndexL][EffectIndexJ][EffectIndexJ][indexU];
				//term *= UnfoldingCovarianceTerm( CauseIndexK, EffectIndexI, EffectIndexI, indexU ); 
				//term *= UnfoldingCovarianceTerm( CauseIndexL, EffectIndexJ, EffectIndexJ, indexU );

				result += term;
			}
		}
	}
	else
	{
		for ( int indexU = 0; indexU < binNumber; indexU++ )
		{
			if ( indexU == CauseIndexK )
			{
				//If U==K!=L, true for all R, S==J
				for ( int indexR = 0; indexR < binNumber; indexR++ )
				{
					double term = inputSmearing->GetCovarianceElement( indexR, EffectIndexJ, indexU );
					term *= unfoldingCovarianceTerm[CauseIndexK][EffectIndexI][indexR][indexU];
					term *= unfoldingCovarianceTerm[CauseIndexL][EffectIndexJ][EffectIndexJ][indexU];
					//term *= UnfoldingCovarianceTerm( CauseIndexK, EffectIndexI, indexR, indexU );
					//term *= UnfoldingCovarianceTerm( CauseIndexL, EffectIndexJ, EffectIndexJ, indexU );

					result += term;
				}
			}
			else if ( indexU == CauseIndexL )
			{
				//If U==L!=K, true for all S, R==I
				for ( int indexS = 0; indexS < binNumber; indexS++ )
				{
					double term = inputSmearing->GetCovarianceElement( EffectIndexI, indexS, indexU );
					term *= unfoldingCovarianceTerm[CauseIndexK][EffectIndexI][EffectIndexI][indexU];
					term *= unfoldingCovarianceTerm[CauseIndexL][EffectIndexJ][indexS][indexU];
					//term *= UnfoldingCovarianceTerm( CauseIndexK, EffectIndexI, EffectIndexI, indexU );
					//term *= UnfoldingCovarianceTerm( CauseIndexL, EffectIndexJ, indexS, indexU );

					result += term;
				}
			}
			else
			{
				//If U!=K!=L, true only for R==I and S==J
				double term = inputSmearing->GetCovarianceElement( EffectIndexI, EffectIndexJ, indexU );
				term *= unfoldingCovarianceTerm[CauseIndexK][EffectIndexI][EffectIndexI][indexU];
				term *= unfoldingCovarianceTerm[CauseIndexL][EffectIndexJ][EffectIndexJ][indexU];
				//term *= UnfoldingCovarianceTerm( CauseIndexK, EffectIndexI, EffectIndexI, indexU );
				//term *= UnfoldingCovarianceTerm( CauseIndexL, EffectIndexJ, EffectIndexJ, indexU );

				result += term;
			}
		}
	}

	return result;
}

//Unfolding matrix and a lot of Kroneker Deltas
void CovarianceMatrix::MakeUnfoldingCovarianceTerm()
{
	//Initialise the tensor, apologies for the horrible syntax
	unfoldingCovarianceTerm = vector< vector< vector< vector<double> > > >( binNumber, vector< vector< vector<double> > >( binNumber, vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) ) ) );

	//Loop over all k, i, r and u to populate the tensor
	for ( int indexK = 0; indexK < binNumber; indexK++ )
	{
		for ( int indexI = 0; indexI < binNumber; indexI++ )
		{
			for ( int indexR = 0; indexR < binNumber; indexR++ )
			{
				for ( int indexU = 0; indexU < binNumber; indexU++ )
				{
					//Calculate first term to sum
					double sumTerm1;
					//Evaluate the delta functions
					if ( indexK == indexU && indexR == indexI )
					{
						//Retrieve the element of the smearing matrix
						double smearingElement = inputSmearing->GetElement( indexU, indexR );

						//Check for divide-by-zero errors
						if ( smearingElement == 0.0 )
						{
							cerr << "COVARIANCE: DIVIDE by zero error (smearing) with numerator = 1.0" << endl;
							exit(1);
						}
						else
						{
							sumTerm1 = 1.0 / smearingElement;
						}
					}
					else
					{
						sumTerm1 = 0.0;
					}

					//Calculate second term to sum
					double sumTerm2;
					//Evaluate the delta
					if ( indexK == indexU )
					{
						//Retrieve the efficiency
						double efficiency = inputSmearing->GetEfficiency(indexU);

						//Check for divide-by-zero errors
						if ( efficiency == 0.0 )
						{
							cerr << "COVARIANCE: DIVIDE by zero error (efficiency) with numerator = 1.0" << endl;
							exit(1);
						}
						else
						{
							sumTerm2 = 1.0 / efficiency;
						}
					}
					else
					{
						sumTerm2 = 0.0;
					}

					//Calculate third term to sum
					double sumTerm3;
					//Evaluate the delta
					if ( indexR == indexI )
					{
						//Evaluate the numerator
						sumTerm3 = inputUnfolding->GetElement( indexU, indexI ) * inputSmearing->GetEfficiency(indexU);

						if ( sumTerm3 != 0.0 )
						{
							//Retrieve an element of the smearing matrix
							double smearingElement = inputSmearing->GetElement( indexU, indexI );

							//Check for divide-by-zero errors
							if ( smearingElement == 0.0 )
							{
								cerr << "COVARIANCE: DIVIDE by zero error (smearing) with numerator = " << sumTerm3 << endl;
								exit(1);
							}
							else
							{
								sumTerm3 /= inputSmearing->GetElement( indexU, indexI );
							}
						}
					}
					else
					{
						sumTerm3 = 0.0;
					}

					double term = sumTerm1 - sumTerm2 - sumTerm3;
					term *= inputUnfolding->GetElement( indexK, indexI );
					unfoldingCovarianceTerm[indexK][indexI][indexR][indexU] = term;
				}
			}
		}
	}
}

//Unfolding matrix and a lot of Kroneker Deltas
double CovarianceMatrix::UnfoldingCovarianceTerm( int indexK, int indexI, int indexR, int indexU )
{
	//Calculate first term to sum
	double sumTerm1;
	//Evaluate the delta functions
	if ( indexK == indexU && indexR == indexI )
	{
		//Retrieve the element of the smearing matrix
		double smearingElement = inputSmearing->GetElement( indexU, indexR );

		//Check for divide-by-zero errors
		if ( smearingElement == 0.0 )
		{
			cerr << "COVARIANCE: DIVIDE by zero error (smearing) with numerator = 1.0" << endl;
			exit(1);
		}
		else
		{
			sumTerm1 = 1.0 / smearingElement;
		}
	}
	else
	{
		sumTerm1 = 0.0;
	}

	//Calculate second term to sum
	double sumTerm2;
	//Evaluate the delta
	if ( indexK == indexU )
	{
		//Retrieve the efficiency
		double efficiency = inputSmearing->GetEfficiency(indexU);
		//Check for divide-by-zero errors
		if ( efficiency == 0.0 )
		{
			cerr << "COVARIANCE: DIVIDE by zero error (efficiency) with numerator = 1.0" << endl;
			exit(1);
		}
		else
		{
			sumTerm2 = 1.0 / efficiency;
		}
	}
	else
	{
		sumTerm2 = 0.0;
	}

	//Calculate third term to sum
	double sumTerm3;
	//Evaluate the delta
	if ( indexR == indexI )
	{
		//Evaluate the numerator
		sumTerm3 = inputUnfolding->GetElement( indexU, indexI ) * inputSmearing->GetEfficiency(indexU);

		if ( sumTerm3 != 0.0 )
		{
			//Retrieve an element of the smearing matrix
			double smearingElement = inputSmearing->GetElement( indexU, indexI );

			//Check for divide-by-zero errors
			if ( smearingElement == 0.0 )
			{
				cerr << "COVARIANCE: DIVIDE by zero error (smearing) with numerator = " << sumTerm3 << endl;
				exit(1);
			}
			else
			{
				sumTerm3 /= inputSmearing->GetElement( indexU, indexI );
			}
		}
	}
	else
	{
		sumTerm3 = 0.0;
	}

	double term = sumTerm1 - sumTerm2 - sumTerm3;
	term *= inputUnfolding->GetElement( indexK, indexI );
	return term;
}
