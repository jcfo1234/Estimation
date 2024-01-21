#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <math.h>
#include "..\..\DataTypeDefinitions.h"

#define MATRIX_NORM2_TOLERANCE ((DOUBLE) 1e-12)

class CMatrix
{
   // Default constructor
   CMatrix();
public:
   // Matrix constructor
   CMatrix(const CHAR *name, INT rows, INT cols);
   // Matrix copy constructor
   CMatrix(const CMatrix &other);
   // Matrix destructor
   ~CMatrix();
   // Set Matrix name
   void SetName(const CHAR *name) { strcpy_s(ac_name, name); }
   // Set matrix component value
   void SetComponent(const INT i_therow, const INT i_thecol, const DOUBLE dValue_);
   // Set matrix to identity
   void SetIdentity();
   // Set matrix to zero
   void SetZero();
   // Enable/Disable warning display in matrix assignment
   void SetWarningEnable(BOOLEANO bSet_) { bWarningEnable = bSet_; }
   // Get matrix component value
   DOUBLE GetComponent(const INT i_therow, const INT i_thecol);
   // Get number of columns
   INT GetNumCols() { return i_columns; }
   // Get number of rows
   INT GetNumRows() { return i_rows; }
   // Get Matrix name
   const char* GetName() const { return ac_name; }
   // Find a square matrix determinant
   DOUBLE Determinant();
   // Matrix assignment operator
   CMatrix& operator = (const CMatrix &other);
   // Find cofactor matrix
   CMatrix CoFactor();
   // Find adjoint matrix
   CMatrix Adjoint();
   // Find matrix inverse
   CMatrix Inverse();
   // Find matrix numeric inverse
   CMatrix NumericInverse();
   // Find another matrix numeric inverse
   CMatrix NumericInverse2(DOUBLE dTolerance_= MATRIX_NORM2_TOLERANCE);
   // Find the matrix pseudo-inverse
   CMatrix NumericPseudoInverse(DOUBLE dTolerance_ = MATRIX_NORM2_TOLERANCE);
   // Find matrix power
   CMatrix MatrixPower(ULONG ulDegree_);
   // Matrix addition operator
   CMatrix operator + (const CMatrix& other);
   // Matrix subtraction operator
   CMatrix operator - (const CMatrix& other);
   // Matrix multiplication operator
   CMatrix operator * (const CMatrix& other);
   // Augment dimensions of Matrix
   CMatrix AugmentDim(const CMatrix& other, INT iExtraRows_, INT iExtraColumns_);
   // Submatrix of Matrix
   CMatrix SubMatrix(const CMatrix& other, INT iRowInit_, INT iRowEnd_, INT iColInit_, INT iColEnd_);
   // Matrix multiplication to a 1x1 value
   CMatrix operator * (const DOUBLE& dValue);
   // Matrix equal to operator
   BOOLEANO operator == (const CMatrix& other);
   // Matrix transpose
   CMatrix Transpose();
   // Matrix Norm-1
   DOUBLE MatrixNorm1();
   // Matrix Norm-Infinity
   DOUBLE MatrixNormInf();
   // Matrix Norm-2 (Defined as the largest singular value)
   DOUBLE MatrixNorm2();
   // Return factorial of Pade degree
   DOUBLE PadeDegreeFactorial(ULONG ulPadeDegree_);
   // Return Pade truncation degree
   ULONG MatrixPadeDegree(DOUBLE dTolerance_ = MATRIX_NORM2_TOLERANCE);
   // Return Pade truncation degree
   ULONG MatrixPadeDegree2(DOUBLE dScale_, DOUBLE dTolerance_ = MATRIX_NORM2_TOLERANCE);
   // Return Pade scaling
   DOUBLE MatrixPadeScalig(DOUBLE dRatio_ = 1);
   // Return Matrix Pade Numerator
   CMatrix PadeNumerator(ULONG ulPadeDegree_, DOUBLE dScale_ = 1);
   // Return Matrix Pade Denominator
   CMatrix PadeDenominator(ULONG ulPadeDegree_, DOUBLE dScale_ = 1);
   // Return Matrix Exponential
   CMatrix MatrixExponential();
   // Return Matrix Exponential
   CMatrix MatrixExponential2();

private:
   INT i_rows;
   INT i_columns;
   CHAR ac_name[128];
   DOUBLE** pd_Data;
   BOOLEANO bWarningEnable;
};

#endif