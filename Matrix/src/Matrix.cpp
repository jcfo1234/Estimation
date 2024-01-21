//----------------------------------------------------------------------
// Matrix.cpp
// Source file of Matrix.h
// Basic matrix library
//----------------------------------------------------------------------
#include "..\inc\Matrix.h"

//----------------------------------------------------------------------------
CMatrix::CMatrix(const CHAR *name, INT rows, INT cols) : i_rows(rows), i_columns(cols), bWarningEnable(TRUE)
{
   strcpy_s(ac_name, name);
   pd_Data = new DOUBLE*[i_rows];
   for (INT i = 0; i < i_rows; i++)
      pd_Data[i] = new DOUBLE[i_columns];

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         pd_Data[i][j] = 0.0;
      }
   }
}

//----------------------------------------------------------------------------
CMatrix::CMatrix(const CMatrix &other)
{
   strcpy_s(ac_name, other.ac_name);
   i_rows = other.i_rows;
   i_columns = other.i_columns;
   bWarningEnable = other.bWarningEnable;

   pd_Data = new DOUBLE*[i_rows];
   for (INT i = 0; i < i_rows; i++)
      pd_Data[i] = new DOUBLE[i_columns];

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         pd_Data[i][j] = other.pd_Data[i][j];
      }
   }
}

//----------------------------------------------------------------------------
CMatrix::~CMatrix()
{
   for (INT i = 0; i < i_rows; i++)
      delete[] pd_Data[i];
   delete[] pd_Data;
   i_rows = i_columns = 0;
}

//----------------------------------------------------------------------------
void CMatrix::SetComponent(const INT i_therow, const INT i_thecol, const DOUBLE dValue_)
{
   if ( (i_therow < i_rows) && (i_thecol < i_columns) )
      pd_Data[i_therow][i_thecol] = dValue_;
}

//----------------------------------------------------------------------------
void CMatrix::SetIdentity()
{
   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         if (j == i)
         {
            pd_Data[i][j] = 1;
         }
         else
         {
            pd_Data[i][j] = 0;
         }
      }
   }
}

//----------------------------------------------------------------------------
void CMatrix::SetZero()
{
   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         pd_Data[i][j] = 0;
      }
   }
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::GetComponent(const INT i_therow, const INT i_thecol)
{
   DOUBLE dComponent = 0;
   if ((i_therow < i_rows) && (i_thecol < i_columns))
      dComponent = pd_Data[i_therow][i_thecol];
   return dComponent;
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::Determinant()
{
   DOUBLE dDet = 0;
   DOUBLE **pd = pd_Data;
   switch (i_rows)
   {
   // 2 x 2 Matrix
   case 2:
   {
      dDet = pd[0][0] * pd[1][1] - pd[0][1] * pd[1][0];
      return dDet;
   }
   break;
   // 3 x 3 Matrix
   case 3:
   {
      /***
      a b c
      d e f
      g h i

      a b c a b c
      d e f d e f
      g h i g h i

      // det (A) = aei + bfg + cdh - afh - bdi - ceg.
      ***/

      DOUBLE a = pd[0][0];
      DOUBLE b = pd[0][1];
      DOUBLE c = pd[0][2];

      DOUBLE d = pd[1][0];
      DOUBLE e = pd[1][1];
      DOUBLE f = pd[1][2];

      DOUBLE g = pd[2][0];
      DOUBLE h = pd[2][1];
      DOUBLE i = pd[2][2];

      dDet = (a*e*i + b*f*g + c*d*h);
      dDet = dDet - a*f*h;
      dDet = dDet - b*d*i;
      dDet = dDet - c*e*g;

      return dDet;
   }
   break;
   // 4 x 4 matrix
   case 4:
   {
      CMatrix *temp[4];
      for (INT i = 0; i < 4; i++)
         temp[i] = new CMatrix("", 3, 3);
      // All row iteration
      for (INT k = 0; k < 4; k++)
      {
         // Discard first row
         for (INT i = 1; i < 4; i++)
         {
            INT j1 = 0;
            for (INT j = 0; j < 4; j++)
            {
               if (k == j)
                  continue;
               temp[k]->pd_Data[i - 1][j1++] = this->pd_Data[i][j];
            }
         }
      }
      // Row 1, Column j * Determinant of 3 x 3 matrix
      dDet = this->pd_Data[0][0] * temp[0]->Determinant() -
             this->pd_Data[0][1] * temp[1]->Determinant() +
             this->pd_Data[0][2] * temp[2]->Determinant() -
             this->pd_Data[0][3] * temp[3]->Determinant();
      // Delete allocated memory
      for (INT i = 0; i < 4; i++)
         delete temp[i];

      return dDet;
   }
   break;

   case 5:
   {
      CMatrix *temp[5];
      for (INT i = 0; i < 5; i++)
         temp[i] = new CMatrix("", 4, 4);

      for (INT k = 0; k < 5; k++)
      {

         for (INT i = 1; i < 5; i++)
         {
            INT j1 = 0;
            for (INT j = 0; j < 5; j++)
            {
               if (k == j)
                  continue;
               temp[k]->pd_Data[i - 1][j1++] = this->pd_Data[i][j];
            }
         }
      }
      double det = this->pd_Data[0][0] * temp[0]->Determinant() -
                   this->pd_Data[0][1] * temp[1]->Determinant() +
                   this->pd_Data[0][2] * temp[2]->Determinant() -
                   this->pd_Data[0][3] * temp[3]->Determinant() +
                   this->pd_Data[0][4] * temp[4]->Determinant();
      // Delete allocated memory
      for (INT i = 0; i < 5; i++)
         delete temp[i];
      return det;
   }
   case 6:
   case 7:
   case 8:
   case 9:
   case 10:
   case 11:
   case 12:
   default:
   {
      INT DIM = i_rows;
      CMatrix **temp = new CMatrix*[DIM];
      for (INT i = 0; i < DIM; i++)
         temp[i] = new CMatrix("", DIM - 1, DIM - 1);

      for (INT k = 0; k < DIM; k++)
      {

         for (INT i = 1; i < DIM; i++)
         {
            INT j1 = 0;
            for (INT j = 0; j < DIM; j++)
            {
               if (k == j)
                  continue;
               temp[k]->pd_Data[i - 1][j1++] = this->pd_Data[i][j];
            }
         }
      }

      DOUBLE det = 0;
      for (INT k = 0; k < DIM; k++)
      {
         if ((k % 2) == 0)
            det = det + (this->pd_Data[0][k] * temp[k]->Determinant());
         else
            det = det - (this->pd_Data[0][k] * temp[k]->Determinant());
      }

      for (INT i = 0; i < DIM; i++)
         delete temp[i];
      delete[] temp;

      return det;
   }
   break;
   }
}

//----------------------------------------------------------------------------
CMatrix& CMatrix::operator = (const CMatrix &other)
{
   if (bWarningEnable == TRUE)
   {
      if (this->i_rows != other.i_rows ||
         this->i_columns != other.i_columns)
      {
         std::cout << "\nWARNING: Assignment is taking place with by changing the number of rows and columns of the matrix" << std::endl;
      }
   }
   for (INT i = 0; i < i_rows; i++)
      delete[] pd_Data[i];
   delete[] pd_Data;
   i_rows = i_columns = 0;

   strcpy_s(ac_name, other.ac_name);
   i_rows = other.i_rows;
   i_columns = other.i_columns;

   pd_Data = new DOUBLE*[i_rows];
   for (INT i = 0; i < i_rows; i++)
      pd_Data[i] = new DOUBLE[i_columns];

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         pd_Data[i][j] = other.pd_Data[i][j];
      }
   }

   return *this;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::CoFactor()
{
   CMatrix cofactor("COF", i_rows, i_columns);
   if (i_rows != i_columns)
      return cofactor;

   if (i_rows < 2)
      return cofactor;
   else if (i_rows == 2)
   {
      cofactor.pd_Data[0][0] = pd_Data[1][1];
      cofactor.pd_Data[0][1] = -pd_Data[1][0];
      cofactor.pd_Data[1][0] = -pd_Data[0][1];
      cofactor.pd_Data[1][1] = pd_Data[0][0];
      return cofactor;
   }
   else if (i_rows >= 3)
   {
      INT DIM = i_rows;
      CMatrix ***temp = new CMatrix**[DIM];
      for (INT i = 0; i < DIM; i++)
         temp[i] = new CMatrix*[DIM];
      for (INT i = 0; i < DIM; i++)
         for (INT j = 0; j < DIM; j++)
            temp[i][j] = new CMatrix("", DIM - 1, DIM - 1);

      for (INT k1 = 0; k1 < DIM; k1++)
      {
         for (INT k2 = 0; k2 < DIM; k2++)
         {
            INT i1 = 0;
            for (INT i = 0; i < DIM; i++)
            {
               INT j1 = 0;
               for (INT j = 0; j < DIM; j++)
               {
                  if (k1 == i || k2 == j)
                     continue;
                  temp[k1][k2]->pd_Data[i1][j1++] = this->pd_Data[i][j];
               }
               if (k1 != i)
                  i1++;
            }
         }
      }

      BOOLEANO flagPositive = TRUE;
      for (INT k1 = 0; k1 < DIM; k1++)
      {
         flagPositive = ((k1 % 2) == 0);

         for (INT k2 = 0; k2 < DIM; k2++)
         {
            if (flagPositive == TRUE)
            {
               cofactor.pd_Data[k1][k2] = temp[k1][k2]->Determinant();
               flagPositive = FALSE;
            }
            else
            {
               cofactor.pd_Data[k1][k2] = -temp[k1][k2]->Determinant();
               flagPositive = TRUE;
            }
         }

      }

      for (INT i = 0; i < DIM; i++)
         for (INT j = 0; j < DIM; j++)
            delete temp[i][j];

      for (INT i = 0; i < DIM; i++)
         delete[] temp[i];

      delete[] temp;
   }
   return cofactor;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::Adjoint()
{
   CMatrix cofactor("COF", i_rows, i_columns);
   CMatrix adj("ADJ", i_rows, i_columns);
   if (i_rows != i_columns)
      return adj;

   cofactor = this->CoFactor();

   // adjoint is transpose of a cofactor of a matrix
   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         adj.pd_Data[j][i] = cofactor.pd_Data[i][j];
      }
   }
   return adj;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::Inverse()
{
   CMatrix cofactor("COF", i_rows, i_columns);
   CMatrix inv("INV", i_rows, i_columns);
   if (i_rows != i_columns)
      return inv;

   // to find out Determinant
   DOUBLE det = Determinant();

   cofactor = this->CoFactor();

   // inv = transpose of cofactor / Determinant
   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         inv.pd_Data[j][i] = cofactor.pd_Data[i][j] / det;
      }
   }
   return inv;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::NumericInverse()
{
   // This method uses the Hotelling-Bodewig numerical inverse algorithm
   CMatrix A = *this;
   // Identity matrix
   CMatrix I("Identity", i_rows, i_columns);
   I.SetIdentity();
   // Initialize matrix inverse
   DOUBLE dScale = 1 / (A.MatrixNorm1() * A.MatrixNormInf());
   CMatrix V(A.Transpose() * dScale);
   // Ideally Test is a Zero matrix
   CMatrix Test(I - A * V);
   // Norm-2 test
   while (Test.MatrixNormInf() > MATRIX_NORM2_TOLERANCE)
   {
      // Hotelling-Bodewig algorithm for finding matrix inverse
      V = V * I * 2 - V * A * V;
      Test = I - A * V;
   }
   return V;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::NumericInverse2(DOUBLE dTolerance_)
{
   // This method uses the Hotelling-Bodewig based numerical inverse algorithm
   CMatrix A = *this;
   // Identity matrix
   CMatrix I("Identity", i_rows, i_columns);
   CMatrix W("W", i_rows, i_columns);
   I.SetIdentity();
   // Initialize matrix inverse
   DOUBLE dScale = 1 / (A.MatrixNorm1() * A.MatrixNormInf());
   CMatrix V(A.Transpose() * dScale);
   // Ideally Test is a Zero matrix
   CMatrix Test(I - A * V);
   // Norm-2 test
   while (Test.MatrixNormInf() > dTolerance_)
   {
      // Hotelling-Bodewig algorithm for finding matrix inverse
      W = I *(-147) + A * V * I * 53 - A * V * A * V * 11 + A * V * A * V * A * V;
      W = I * 259 + A * V * W;
      W = I * (-301) + A * V * W;
      W = I * 231 + A * V * W;
      W = I * (-113) + A * V * W;
      W = I * 32 + A * V * W;
      V = V * W * 0.25;
      Test = I - A * V;
   }
   return V;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::NumericPseudoInverse(DOUBLE dTolerance_)
{
   // This method uses the Ben-Israel and Cohen algorithm for pseudo-inverses
   CMatrix A = *this;
   // Identity matrix
   CMatrix I("Identity", i_rows, i_columns);
   I.SetIdentity();
   // Initialize matrix pseudo-inverse
   DOUBLE dScale;
   (A.MatrixNorm1() * A.MatrixNormInf() > dTolerance_) ? dScale = 1 / (A.MatrixNorm1() * A.MatrixNormInf()) : dScale = dTolerance_;
   CMatrix clTemp(A.Transpose() * A + I * dScale);
   CMatrix V(clTemp.NumericInverse2(1e-6) * A.Transpose());
   // Ideally Test is a Zero matrix
   CMatrix Test(A - A * V * A);
   // Norm-2 test
   while (Test.MatrixNormInf() > dTolerance_)
   {
      // Ben-Israel and Cohen algorithm for finding matrix inverse
      V = V * 2 - V * A * V;
      Test = A - A * V * A;
   }
   return V;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::MatrixPower(ULONG ulDegree_)
{
   CMatrix clPower(*this);
   CMatrix clTemp(*this);
   CMatrix clIdentity(*this);
   clIdentity.SetIdentity();
   // Power greater than 0
   if (ulDegree_ > 0)
   {
      for (ULONG ul = 1; ul < ulDegree_; ul++)
      {
         clPower = clTemp * clPower;
      }
   }
   // Power equal to 0
   else
   {
      clPower = clIdentity;
   }
   return clPower;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::operator + (const CMatrix &other)
{
   if (this->i_rows != other.i_rows ||
       this->i_columns != other.i_columns)
   {
      std::cout << "Addition could not take place because number of rows and columns are different between the two matrices";
      return *this;
   }

   CMatrix result("", i_rows, i_columns);

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         result.pd_Data[i][j] = this->pd_Data[i][j] + other.pd_Data[i][j];
      }
   }
   return result;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::operator - (const CMatrix& other)
{
   if (this->i_rows != other.i_rows ||
       this->i_columns != other.i_columns)
   {
      std::cout << "Subtraction could not take place because number of rows and columns are different between the two matrices";
      return *this;
   }

   CMatrix result("", i_rows, i_columns);

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         result.pd_Data[i][j] = this->pd_Data[i][j] - other.pd_Data[i][j];
      }
   }
   return result;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::operator * (const CMatrix& other)
{
   if (this->i_columns != other.i_rows)
   {
      std::cout << "Multiplication could not take place because number of columns of 1st Matrix and number of rows in 2nd Matrix are different";
      return *this;
   }

   CMatrix result("", this->i_rows, other.i_columns);

   for (INT i = 0; i < this->i_rows; i++)
   {
      for (INT j = 0; j < other.i_columns; j++)
      {
         for (INT k = 0; k < this->i_columns; k++)
         {
            result.pd_Data[i][j] += this->pd_Data[i][k] * other.pd_Data[k][j];
         }
      }
   }
   return result;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::AugmentDim(const CMatrix& other, INT iExtraRows_, INT iExtraColumns_)
{
   CMatrix result("", this->i_rows + iExtraRows_, this->i_columns + iExtraColumns_);

   for (INT i = 0; i < other.i_rows; i++)
   {
      for (INT j = 0; j < other.i_columns; j++)
      {
         result.SetComponent(i, j, other.pd_Data[i][j]);
      }
   }

   return result;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::SubMatrix(const CMatrix& other, INT iRowInit_, INT iRowEnd_, INT iColInit_, INT iColEnd_)
{
   CMatrix result("", iRowEnd_ - iRowInit_, iColEnd_ - iColInit_);

   for (INT i = iRowInit_; i < iRowEnd_; i++)
   {
      for (INT j = iColInit_; j < iColEnd_; j++)
      {
         result.SetComponent(i - iRowInit_, j - iColInit_, other.pd_Data[i][j]);
      }
   }

   return result;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::operator * (const DOUBLE& dValue)
{
   CMatrix result("", this->i_rows, this->i_columns);
   for (INT i = 0; i < this->i_rows; i++)
   {
      for (INT j = 0; j < this->i_columns; j++)
      {
         result.pd_Data[i][j] = this->pd_Data[i][j] * dValue;
      }
   }
   return result;
}

//----------------------------------------------------------------------------
BOOLEANO CMatrix::operator == (const CMatrix& other)
{
   if (this->i_rows != other.i_rows ||
       this->i_columns != other.i_columns)
   {
      std::cout << "Comparision could not take place because number of rows and columns are different between the two matrices";
      return false;
   }

   CMatrix result("", i_rows, i_columns);

   BOOLEANO bEqual = TRUE;
   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         if (this->pd_Data[i][j] != other.pd_Data[i][j])
            bEqual = FALSE;
      }
   }
   return bEqual;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::Transpose()
{
   CMatrix trans("TR", i_columns, i_rows);

   for (INT i = 0; i < i_rows; i++)
   {
      for (INT j = 0; j < i_columns; j++)
      {
         trans.pd_Data[j][i] = pd_Data[i][j];
      }
   }
   return trans;
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::MatrixNorm1()
{
   DOUBLE dMax = 0;
   DOUBLE dColSum = 0;
   for (INT j = 0; j < i_columns; j++)
   {
      dColSum = 0;
      // Sumation of column absolute values
      for (INT i = 0; i < i_rows; i++)
         dColSum = dColSum + abs(pd_Data[i][j]);
      // Search for greatest sum of column elements
      if (dColSum > dMax)
         dMax = dColSum;
   }
   return dMax;
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::MatrixNormInf()
{
   DOUBLE dMax = 0;
   DOUBLE dColSum = 0;
   for (INT i = 0; i < i_rows; i++)
   {
      dColSum = 0;
      // Sumation of row absolute values
      for (INT j = 0; j < i_columns; j++)
         dColSum = dColSum + abs(pd_Data[i][j]);
      // Search for greatest sum of column elements
      if (dColSum > dMax)
         dMax = dColSum;
   }
   return dMax;
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::MatrixNorm2()
{
   // Iterator to find largest eigenvalue of matrix
   CMatrix clVector1("Unp1", i_rows, 1);
   CMatrix clVector2("Un1", i_rows, 1);
   // Constant vector linearly independent from Iterator vector
   CMatrix clWeightv("Wt", 1, i_columns);
   // Current matrix which we want to find out its largest eigenvalue
   CMatrix clTemp = this->Transpose() * (*this);
   // Dot product between weight vector and iterator vectors
   CMatrix dDotProduct1("Dot1", 1, 1);
   CMatrix dDotProduct2("Dot2", 1, 1);
   // Ratios
   DOUBLE dCurrRatio = 0;
   DOUBLE dOldRatio = 1000;

   // Initialize vectors
   for (INT i = 0; i < i_rows; i++)
   {
      clVector1.SetComponent(i, 0, 1);
      clWeightv.SetComponent(0, i, 1);
   }
   clWeightv.SetComponent(0, 0, 2);
   // Find largest eigenvalue
   while (abs(dCurrRatio - dOldRatio) > MATRIX_NORM2_TOLERANCE)
   {
      dOldRatio = dCurrRatio;
      clVector2 = clTemp * clVector1;
      dDotProduct1 = clWeightv * clVector1;
      dDotProduct2 = clWeightv * clVector2;
      dCurrRatio = dDotProduct2.GetComponent(0, 0) / dDotProduct1.GetComponent(0, 0);
      clVector1 = clVector2;
   }
   // Return singular value of matrix
   return sqrt(abs(dCurrRatio));
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::PadeDegreeFactorial(ULONG ulPadeDegree_)
{
   DOUBLE dRetVal = 1.0;
   // Compute Factorial of Pade Degree
   for (ULONG ul = ulPadeDegree_; ul > 1; ul--)
   {
      dRetVal *= (DOUBLE)ul;
   }
   return dRetVal;
}

//----------------------------------------------------------------------------
ULONG CMatrix::MatrixPadeDegree(DOUBLE dTolerance_)
{
   ULONG ulPadeDegree = 0;
   DOUBLE dRatio = (this->MatrixNormInf() / PadeDegreeFactorial(ulPadeDegree)) * (1 / (1 - this->MatrixNormInf() / ((DOUBLE)(ulPadeDegree + 2))));
   
   // Return the minimum Pade degree that satisfies tolerance
   while (abs(dRatio) > dTolerance_)
   {
      ulPadeDegree++;
      dRatio = (this->MatrixNormInf() / PadeDegreeFactorial(ulPadeDegree)) * (1 / (1 - this->MatrixNormInf() / ((DOUBLE)(ulPadeDegree + 2))));
   }

   return ulPadeDegree;
}

//----------------------------------------------------------------------------
ULONG CMatrix::MatrixPadeDegree2(DOUBLE dScale_, DOUBLE dTolerance_)
{
   ULONG ulPadeDegree = 0;
   DOUBLE dFactor = 1;
   while (8 * pow(this->MatrixNormInf() * dScale_, 2 * ulPadeDegree) * 1 / dFactor > dTolerance_)
   {
      ulPadeDegree++;
      dFactor = 2 * ulPadeDegree + 1;
      for (ULONG ul = 2 * ulPadeDegree; ul > ulPadeDegree + 1; ul--)
         dFactor = dFactor * pow(ul, 2);
   }
   return ulPadeDegree;
}

//----------------------------------------------------------------------------
DOUBLE CMatrix::MatrixPadeScalig(DOUBLE dRatio_)
{
   ULONG ulj = 1;
   while (this->MatrixNormInf() / pow(2, ulj) >= 1)
      ulj++;
   return pow(2, ulj);
}

//----------------------------------------------------------------------------
CMatrix CMatrix::PadeNumerator(ULONG ulPadeDegree_, DOUBLE dScale_)
{
   // Numerator Matrix
   CMatrix clNumPade(*this);
   clNumPade.SetZero();
   // State Matrix
   CMatrix clStateMatrix(*this * dScale_);
   // Cofactor constant
   DOUBLE dNumerator = 0;
   DOUBLE dDenominator = 0;

   // Compute Pade Numerator Matrix
   for (ULONG ulj = 0; ulj <= ulPadeDegree_; ulj++)
   {
      dNumerator = PadeDegreeFactorial(ulPadeDegree_ + ulPadeDegree_ - ulj) * PadeDegreeFactorial(ulPadeDegree_);
      dDenominator = PadeDegreeFactorial(ulPadeDegree_ + ulPadeDegree_) * PadeDegreeFactorial(ulj) * PadeDegreeFactorial(ulPadeDegree_ - ulj);
      clNumPade = clNumPade + clStateMatrix.MatrixPower(ulj) * (dNumerator / dDenominator);
   }

   return clNumPade;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::PadeDenominator(ULONG ulPadeDegree_, DOUBLE dScale_)
{
   // Denominator Matrix
   CMatrix clDenPade(*this);
   clDenPade.SetZero();
   // State Matrix
   CMatrix clStateMatrix(*this * dScale_);
   clStateMatrix = clStateMatrix * (-1);
   // Cofactor constant
   DOUBLE dNumerator = 0;
   DOUBLE dDenominator = 0;

   // Compute Pade Numerator Matrix
   for (ULONG ulj = 0; ulj <= ulPadeDegree_; ulj++)
   {
      dNumerator = PadeDegreeFactorial(ulPadeDegree_ + ulPadeDegree_ - ulj) * PadeDegreeFactorial(ulPadeDegree_);
      dDenominator = PadeDegreeFactorial(ulPadeDegree_ + ulPadeDegree_) * PadeDegreeFactorial(ulj) * PadeDegreeFactorial(ulPadeDegree_ - ulj);
      clDenPade = clDenPade + clStateMatrix.MatrixPower(ulj) * (dNumerator / dDenominator);
   }

   return clDenPade;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::MatrixExponential()
{
   ULONG ulTruncation = MatrixPadeDegree();
   CMatrix clNumerator(PadeNumerator(ulTruncation));
   CMatrix clDenominator(PadeDenominator(ulTruncation));

   return clDenominator.NumericInverse2() * clNumerator;
}

//----------------------------------------------------------------------------
CMatrix CMatrix::MatrixExponential2()
{
   ULONG ulm = (ULONG)MatrixPadeScalig();
   DOUBLE dScale = 1 / MatrixPadeScalig();
   ULONG ulTruncation = MatrixPadeDegree2(dScale);
   CMatrix clNumerator(PadeNumerator(ulTruncation, dScale));
   CMatrix clDenominator(PadeDenominator(ulTruncation, dScale));
   CMatrix clRatio(clDenominator.NumericInverse2() * clNumerator);

   return clRatio.MatrixPower(ulm);
}
