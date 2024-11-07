//////////////////////////////////////////////////////////////
//Module that implements some linear algebra routines
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

#include "linear_algebra.hpp"

//////////////////////////////////////////////////////////////
/// Vector Functions
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
double* AllocateVector(int n)
//Allocates a double precision vector of length n
//and sets all entries to zero
{
	double* vector = new double[n];
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}

	return vector;
}

//////////////////////////////////////////////////////////////
void ZeroVector(double* vector, int n)
//Sets all entries in a vector of length n to zero
{
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}
}

//////////////////////////////////////////////////////////////
double* PrintVector(double* vector,int n)
//Prints a vector of length n to the screen
{
	for (int i=0;i<n;i++)
	{
		std::cout << vector[i] << std::endl;
	}

	return vector;
}

//////////////////////////////////////////////////////////////
void CopyVector(double* vector,double* copied_vector,int n)
//Makes a copy of a vector of length n- assumes that the memory for the
//copy has already been allocated
{
	for (int k=0;k<n;k++)
	{
		copied_vector[k] = vector[k];
	}
}
//////////////////////////////////////////////////////////////
void SubtractVectors(double* vec1,double* vec2,int n)
//Overwrites vec1 with vec1-vec2 (both of length n)
{
    for (int k=0;k<n;k++)
    {
        vec1[k] -= vec2[k];
    }
}
//////////////////////////////////////////////////////////////
void ScaleVector(double* vec,double scaleFactor,int n)
//Overwrites vec of length n with scaleFactor*vec
{
    for (int k=0;k<n;k++)
    {
        vec[k] *= scaleFactor;
    }
}
//////////////////////////////////////////////////////////////
void CombineVectors(double* vector1,double* vector2,double scale,int n)
//Carries out the operation x <- x+scale*y, where x = vector1, y = vector2
//and both of length n
{
	for (int k=0;k<n;k++)
	{
		vector1[k] += scale*vector2[k];
	}
}

//////////////////////////////////////////////////////////////
double NormVector(double* vector,int n)
//Compute the l2-norm of a vector of length n
{
	double norm = 0.0;
	for (int k=0;k<n;k++)
	{
		norm += pow(vector[k],2);
	}

	return sqrt(norm);
}

//////////////////////////////////////////////////////////////
void FindMaximum(double* vector,int n, double& max_val, int& index)
//Function to find maximum value of a vector of length n.
//Outputs:
// max_val - the maximum value
// index - index of the maximum value
{
	max_val = vector[0];
	index = 0;
    for (int k=1;k<n;k++)
    {
        if (vector[k] > max_val)
        {   
            max_val = vector[k];
			index = k;
        }
    }
}


//////////////////////////////////////////////////////////////
void DeallocateVector(double* vector)
//Deletes storage for a vector
{
	delete[] vector;
}

//////////////////////////////////////////////////////////////
/// Matrix Functions
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
DoubleMatrix AllocateDoubleMatrix(int noRows,int noCols)
//Allocates a rectangular matrix of type DoubleMatrix
//and sets the entries to 0.
{
	DoubleMatrix matrix;

	matrix.no_rows = noRows;
	matrix.no_cols = noCols;

	matrix.matrix_entries = new double*[noRows];
	for (int i=0;i<noRows;i++)
	{
		matrix.matrix_entries[i] = new double[noCols];
		for (int j=0;j<noCols;j++)
		{
			matrix.matrix_entries[i][j] = 0.0;
		}
	}

	return matrix;
}

//////////////////////////////////////////////////////////////
SymmetricMatrix AllocateSymmetricMatrix(int noRows)
//Allocates a double precision symmetric matrix and sets all entries to zero
//Stored as a lower triangular matrix
{
	SymmetricMatrix matrix;

	matrix.no_rows = noRows;
	matrix.matrix_entries = new double*[noRows];

	for (int k=0;k<noRows;k++)
	{
		matrix.matrix_entries[k] = new double[k+1];
		for (int j=0;j<=k;j++)
		{
			matrix.matrix_entries[k][j] = 0.0;
		}
	}

	return matrix;

}

//////////////////////////////////////////////////////////////
void DeallocateMatrix(SymmetricMatrix& matrix)
//Deletes a matrix of type SymmetrixMatrix
{
	for (int i=0;i<matrix.no_rows;i++)
	{
		delete[] matrix.matrix_entries[i];
	}

	delete[] matrix.matrix_entries;

	matrix.no_rows = 0;
}

//////////////////////////////////////////////////////////////
void DeallocateMatrix(DoubleMatrix& matrix)
//Deletes a matrix of type SmmetrixMatrix
{

	for (int i=0;i<matrix.no_rows;i++)
	{
		delete[] matrix.matrix_entries[i];
	}

	delete[] matrix.matrix_entries;

	matrix.no_rows = 0;
	matrix.no_cols = 0;
}

//////////////////////////////////////////////////////////////
void PrintMatrix(SymmetricMatrix matrix)
//Prints a double precision matrix to the screen
{
	for (int i=0;i<matrix.no_rows;i++)
	{
		for (int j=0;j<=i;j++)
		{
			std::cout << matrix.matrix_entries[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////
void PrintMatrix(DoubleMatrix matrix)
//Prints a double type matrix to the screen
{

	for (int i=0;i<matrix.no_rows;i++)
	{
		for (int j=0;j<matrix.no_cols;j++)
		{
			std::cout << matrix.matrix_entries[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricForwardSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution)
//For a lower triangular matrix, performs forward subsitition to solve Ax = b
//Assumes there are 1s on the main diagonal
{
	solution[0] = rhs[0];
	for (int i=1;i<matrix.no_rows;i++)
	{
		solution[i] = rhs[i];
		for (int j=0;j<i;j++)
		{
			solution[i] -= matrix.matrix_entries[i][j]*solution[j];
		}
	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricBackSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution)
//For a lower triangular matrix, performs back subsitition to solve A^Tx = b
//Assumes there are 1s on the main diagonal
{
	solution[matrix.no_rows-1] = rhs[matrix.no_rows-1];
	for (int i=matrix.no_rows-1;i>=0;i--)
	{
		solution[i] = rhs[i];
		for (int j=i+1;j<matrix.no_rows;j++)
		{
			solution[i] -= matrix.matrix_entries[j][i]*solution[j];
		}
	}
}

//////////////////////////////////////////////////////////////
void ComputeLDLFactorisation(SymmetricMatrix& matrix)
//Finds the LDLT factorisation of a symmetrix matrix and
//overwrites the lower half of the matrix
{
	for (int i=1;i<matrix.no_rows;i++)
	{
		
		for (int j=0;j<i;j++)
		{
			for (int k=0;k<j;k++)
			{
				matrix.matrix_entries[i][j] -= 
					matrix.matrix_entries[i][k]*matrix.matrix_entries[j][k]*
					matrix.matrix_entries[k][k];
			}

			matrix.matrix_entries[i][j] /= matrix.matrix_entries[j][j];

		}

		for (int k=0;k<i;k++) //Find diagonal terms
		{

			matrix.matrix_entries[i][i] -= 
				matrix.matrix_entries[i][k]*matrix.matrix_entries[i][k]*
				matrix.matrix_entries[k][k];

		}

		if (std::fabs(matrix.matrix_entries[i][i]) < 1.0-12)
		{
			std::cout << "Warning: LDL Factorisation not possible" << std::endl;
			break;
		}

	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricSolve(SymmetricMatrix& matrix, double* rhs, double* solution)
//Solves the problem Ax = b, where A is a general symmetrix matrix
{
	//Compute LDL factorisation - in place

	ComputeLDLFactorisation(matrix);


	//Forward substitution

	PerformSymmetricForwardSubstitution(matrix,rhs,solution);

	//Diagonal Division

	for (int i=0;i<matrix.no_rows;i++)
	{
		rhs[i] = solution[i]/matrix.matrix_entries[i][i];
	}

	//Backward Substitution
	PerformSymmetricBackSubstitution(matrix,rhs,solution);

}

//////////////////////////////////////////////////////////////
void MultiplyVectorByMatrix(DoubleMatrix& matrix, double* vec, double* product)
//Computes Ax, where A is a general double precision matrix
{
	int no_rows = matrix.no_rows;
	int no_cols = matrix.no_cols;

	for (int i=0;i<no_rows;i++)
	{
		product[i] = 0.0;
		for (int j=0;j<no_cols;j++)
		{
			product[i] += matrix.matrix_entries[i][j]*vec[j];
		}
	}
}

