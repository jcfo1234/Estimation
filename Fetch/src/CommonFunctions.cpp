//----------------------------------------------------------------------
// CommonFunctions.cpp
// Source file of CommonFunctions.h
// Source code for reading input files
//----------------------------------------------------------------------
#include "..\inc\CommonFunctions.h"

//----------------------------------------------------------------------
CMatrix CalcN(CMatrix* pvLat_, DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_)
{
   // Latitude values must be contained in a vector
   if (min(pvLat_->GetNumRows(), pvLat_->GetNumCols()) != 1)
   {
      throw 5;
   }
   // Figure out the number of elements the output vector is suppose to have
   INT iNumElements = max(pvLat_->GetNumRows(), pvLat_->GetNumCols());
   // Vector storing prime vertical radii values
   CMatrix N("PrimeVerticalRadii", 1, iNumElements);
   // Eccentricity squared
   DOUBLE dEcc2 = (pow(dSemiMajorAxis_, 2) - pow(dSemiMinorAxis_, 2)) / pow(dSemiMajorAxis_, 2);
   // Weight
   DOUBLE dWElem = 0;
   for (INT i = 0; i < iNumElements; i++)
   {
      // Row vector
      if (pvLat_->GetNumRows() >= pvLat_->GetNumCols())
      {
         dWElem = sqrt(1 - dEcc2 * pow(sin(pvLat_->GetComponent(i, 0)), 2));
         N.SetComponent(0, i, dSemiMajorAxis_ / dWElem);
      }
      // Column vector
      else
      {
         dWElem = sqrt(1 - dEcc2 * pow(sin(pvLat_->GetComponent(0, i)), 2));
         N.SetComponent(0, i, dSemiMajorAxis_ / dWElem);
      }
   }
   return N;
}

//----------------------------------------------------------------------
DOUBLE CalcdNdLatitude(DOUBLE dLatitude_, DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_)
{
   DOUBLE dEcc = sqrt(1 - pow(dSemiMinorAxis_ / dSemiMajorAxis_, 2));
   return (dSemiMajorAxis_ * pow(dEcc, 2) * cos(dLatitude_) * sin(dLatitude_) / pow((1 - pow(dEcc, 2) * pow(sin(dLatitude_), 2)), 1.5));
}

//----------------------------------------------------------------------
DOUBLE CalcEccentricity(DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_)
{
   return sqrt(1 - pow(dSemiMinorAxis_ / dSemiMajorAxis_, 2));
}

//----------------------------------------------------------------------
CMatrix CalcM(CMatrix* pvLat_, DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_)
{
   // Latitude values must be contained in a vector
   if (min(pvLat_->GetNumRows(), pvLat_->GetNumCols()) != 1)
   {
      throw 6;
   }
   // Figure out the number of elements the output vector is suppose to have
   INT iNumElements = max(pvLat_->GetNumRows(), pvLat_->GetNumCols());
   // Vector storing meridian radii curvature values
   CMatrix M("MeridianRadiiCurvature", 1, iNumElements);
   // Eccentricity squared
   DOUBLE dEcc2 = (pow(dSemiMajorAxis_, 2) - pow(dSemiMinorAxis_, 2)) / pow(dSemiMajorAxis_, 2);
   // Weight
   DOUBLE dWElem = 0;
   for (INT i = 0; i < iNumElements; i++)
   {
      // Row vector
      if (pvLat_->GetNumRows() >= pvLat_->GetNumCols())
      {
         dWElem = sqrt(1 - dEcc2 * pow(sin(pvLat_->GetComponent(i, 0)), 2));
         M.SetComponent(0, i, (dSemiMajorAxis_ * (1 - dEcc2)) / (pow(dWElem, 3)));
      }
      // Column vector
      else
      {
         dWElem = sqrt(1 - dEcc2 * pow(sin(pvLat_->GetComponent(0, i)), 2));
         M.SetComponent(0, i, (dSemiMajorAxis_ * (1 - dEcc2)) / (pow(dWElem, 3)));
      }
   }
   return M;
}

//----------------------------------------------------------------------------
void GeoidToECEF(CMatrix* pvGeoidCoord_, CMatrix* pvECEFCoord_, DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_)
{
   // Check sizes of vectors coordinates vector (must be 4x1)
   if (pvGeoidCoord_->GetNumRows() != 4 && pvGeoidCoord_->GetNumCols() != 1 && pvECEFCoord_->GetNumRows() != 4 && pvECEFCoord_->GetNumCols() != 1)
      throw 7;
   CMatrix Lat("Latitude", 1, 1);
   Lat.SetComponent(0, 0, pvGeoidCoord_->GetComponent(0, 0));
   CMatrix N(CalcN(&Lat));
   // Eccentricity squared
   DOUBLE dEcc2 = (pow(dSemiMajorAxis_, 2) - pow(dSemiMinorAxis_, 2)) / pow(dSemiMajorAxis_, 2);
   // ECEF X - Position
   DOUBLE dXReturn = (N.GetComponent(0, 0) + pvGeoidCoord_->GetComponent(2, 0)) * cos(pvGeoidCoord_->GetComponent(0, 0)) * cos(pvGeoidCoord_->GetComponent(1, 0));
   // ECEF Y - Position
   DOUBLE dYReturn = (N.GetComponent(0, 0) + pvGeoidCoord_->GetComponent(2, 0)) * cos(pvGeoidCoord_->GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(1, 0));
   // ECEF Z - Position
   DOUBLE dZReturn = (pvGeoidCoord_->GetComponent(2, 0) + (1 - dEcc2) * N.GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(0, 0));

   // ECEF coordinates vector
   pvECEFCoord_->SetComponent(0, 0, dXReturn);
   pvECEFCoord_->SetComponent(1, 0, dYReturn);
   pvECEFCoord_->SetComponent(2, 0, dZReturn);
   pvECEFCoord_->SetComponent(3, 0, pvGeoidCoord_->GetComponent(3, 0));
}

//----------------------------------------------------------------------------
void ECEFToGeoid(CMatrix* pvECEFCoord_, CMatrix* pvGeoidCoord_, DOUBLE dSemiMajorAxis_, DOUBLE dSemiMinorAxis_, DOUBLE CONVERGENCE_)
{
   // Check sizes of vectors coordinates vector (must be 4x1)
   if (pvGeoidCoord_->GetNumRows() != 4 && pvGeoidCoord_->GetNumCols() != 1 && pvECEFCoord_->GetNumRows() != 4 && pvECEFCoord_->GetNumCols() != 1)
      throw 8;
   // Longitude computation
   DOUBLE dLongitude = atan2(pvECEFCoord_->GetComponent(1, 0), pvECEFCoord_->GetComponent(0, 0));
   // Distance from center of earth to projection of position vector to XY plane
   DOUBLE dp = sqrt(pow(pvECEFCoord_->GetComponent(0, 0), 2) + pow(pvECEFCoord_->GetComponent(1, 0), 2));
   // Geocentric latitude
   DOUBLE dLatGeoc = atan2(pvECEFCoord_->GetComponent(2, 0), dp);
   // Latitude iterator
   CMatrix dLatitude("Latitude", 1, 1);
   dLatitude.SetComponent(0, 0, dLatGeoc);
   // Height above ellipsoid
   DOUBLE dHeight = 0;
   CMatrix N("Radii", 1, 1);
   DOUBLE dL = 0;
   DOUBLE dEcc2 = (pow(dSemiMajorAxis_, 2) - pow(dSemiMinorAxis_, 2)) / pow(dSemiMajorAxis_, 2);
   // Find Geodetic coordinates
   while (1)
   {
      N = CalcN(&dLatitude);
      dL = pvECEFCoord_->GetComponent(2, 0) + dEcc2 * N.GetComponent(0, 0) * sin(dLatitude.GetComponent(0, 0));
      if (abs(dHeight - (sqrt(pow(dL, 2) + pow(dp, 2)) - N.GetComponent(0, 0))) < CONVERGENCE_)
         break;
      dHeight = sqrt(pow(dL, 2) + pow(dp, 2)) - N.GetComponent(0, 0);
      dLatitude.SetComponent(0, 0, atan2(pvECEFCoord_->GetComponent(2, 0) * (N.GetComponent(0, 0) + dHeight), dp * (N.GetComponent(0, 0) - N.GetComponent(0, 0) * dEcc2 + dHeight)));
   }

   // Geoid Coordinates
   pvGeoidCoord_->SetComponent(0, 0, dLatitude.GetComponent(0, 0));
   pvGeoidCoord_->SetComponent(1, 0, dLongitude);
   pvGeoidCoord_->SetComponent(2, 0, dHeight);
   pvGeoidCoord_->SetComponent(3, 0, pvECEFCoord_->GetComponent(3, 0));
}

//----------------------------------------------------------------------------
CMatrix ECEFToLLF(CMatrix* pvGeoidCoord_)
{
   // Rotation matrix
   CMatrix Rle("ECEFToLLFRotation", 4, 4);
   // R(0,0) = -sin(longitude)
   Rle.SetComponent(0, 0, -sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(0,1) = cos(longitude)
   Rle.SetComponent(0, 1, cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(0,2) = 0
   Rle.SetComponent(0, 2, 0);
   // R(0,3) = 0
   Rle.SetComponent(0, 3, 0);
   // R(1,0) = -sin(latitude) * cos(longitude)
   Rle.SetComponent(1, 0, -sin(pvGeoidCoord_->GetComponent(0, 0)) * cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(1,1) = -sin(latitude) * sin(longitude)
   Rle.SetComponent(1, 1, -sin(pvGeoidCoord_->GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(1,2) = cos(latitude)
   Rle.SetComponent(1, 2, cos(pvGeoidCoord_->GetComponent(0, 0)));
   // R(1,3) = 0
   Rle.SetComponent(1, 3, 0);
   // R(2,0) = cos(latitude) * cos(longitude)
   Rle.SetComponent(2, 0, cos(pvGeoidCoord_->GetComponent(0, 0)) * cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(2,1) = cos(latitude) * sin(longitude)
   Rle.SetComponent(2, 1, cos(pvGeoidCoord_->GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(2,2) = sin(latitude)
   Rle.SetComponent(2, 2, sin(pvGeoidCoord_->GetComponent(0, 0)));
   // R(2,3) = 0
   Rle.SetComponent(2, 3, 0);
   // R(3,0) = 0
   Rle.SetComponent(3, 0, 0);
   // R(3,1) = 0
   Rle.SetComponent(3, 1, 0);
   // R(3,2) = 0
   Rle.SetComponent(3, 2, 0);
   // R(3,3) = 1
   Rle.SetComponent(3, 3, 1);

   return Rle;
}

//----------------------------------------------------------------------------
CMatrix ECEFToLLF2(CMatrix* pvGeoidCoord_)
{
   // Rotation matrix
   CMatrix Rle("ECEFToLLFRotation", 3, 3);
   // R(0,0) = -sin(longitude)
   Rle.SetComponent(0, 0, -sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(0,1) = cos(longitude)
   Rle.SetComponent(0, 1, cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(0,2) = 0
   Rle.SetComponent(0, 2, 0);
   // R(1,0) = -sin(latitude) * cos(longitude)
   Rle.SetComponent(1, 0, -sin(pvGeoidCoord_->GetComponent(0, 0)) * cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(1,1) = -sin(latitude) * sin(longitude)
   Rle.SetComponent(1, 1, -sin(pvGeoidCoord_->GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(1,2) = cos(latitude)
   Rle.SetComponent(1, 2, cos(pvGeoidCoord_->GetComponent(0, 0)));
   // R(2,0) = cos(latitude) * cos(longitude)
   Rle.SetComponent(2, 0, cos(pvGeoidCoord_->GetComponent(0, 0)) * cos(pvGeoidCoord_->GetComponent(1, 0)));
   // R(2,1) = cos(latitude) * sin(longitude)
   Rle.SetComponent(2, 1, cos(pvGeoidCoord_->GetComponent(0, 0)) * sin(pvGeoidCoord_->GetComponent(1, 0)));
   // R(2,2) = sin(latitude)
   Rle.SetComponent(2, 2, sin(pvGeoidCoord_->GetComponent(0, 0)));

   return Rle;
}

//----------------------------------------------------------------------------
DOUBLE DegToRad(DOUBLE dDeg_)
{
   return dDeg_ * COORDTRANS_PI / 180;
}

//----------------------------------------------------------------------------
DOUBLE RadToDeg(DOUBLE dRad_)
{
   return dRad_ * 180 / COORDTRANS_PI;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeMisclosure(CMatrix* pclPointExpansion_, vector<SatelliteData>* pvEpochSVData_, vector<ULONG>* pvPRNOrder_)
{
   DOUBLE dRange = 0;
   DOUBLE dDistX = 0;
   DOUBLE dDistY = 0;
   DOUBLE dDistZ = 0;
   DOUBLE dTimeShift = 0;
   DOUBLE dMisclosure = 0;
   ULONG ulIndex = 0;
   // Misclosure vector
   CMatrix clMisclosure("Misclosure Vector", pvEpochSVData_->size(), 1);
   pvPRNOrder_->clear();
   // Iterate over collected measurement data vector
   for (vector<SatelliteData>::iterator pIterSVData = pvEpochSVData_->begin(); pIterSVData != pvEpochSVData_->end(); pIterSVData++)
   {
      // Geometric range components
      dDistX = pIterSVData->dPRNX - pclPointExpansion_->GetComponent(0, 0);
      dDistY = pIterSVData->dPRNY - pclPointExpansion_->GetComponent(1, 0);
      dDistZ = pIterSVData->dPRNZ - pclPointExpansion_->GetComponent(2, 0);
      dRange = sqrt(pow(dDistX, 2) + pow(dDistY, 2) + pow(dDistZ, 2));
      // GPS time shift
      dTimeShift = pclPointExpansion_->GetComponent(3, 0);
      // Compute misclosure
      dMisclosure = pIterSVData->dPSR - dRange - dTimeShift;
      ulIndex = pIterSVData - pvEpochSVData_->begin();
      // Store misclosure
      clMisclosure.SetComponent(ulIndex, 0, dMisclosure);
      pvPRNOrder_->push_back(pIterSVData->ulPRN);
   }

   return clMisclosure;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeDesignMatrix(CMatrix* pclPointExpansion_, vector<SatelliteData>* pvEpochSVData_)
{
   ULONG ulRowIndex = 0;
   DOUBLE dRange = 0;
   DOUBLE dDistX = 0;
   DOUBLE dDistY = 0;
   DOUBLE dDistZ = 0;
   CMatrix clA("Design Matrix", pvEpochSVData_->size(), pclPointExpansion_->GetNumRows());

   for (vector<SatelliteData>::iterator pIterSVData = pvEpochSVData_->begin(); pIterSVData != pvEpochSVData_->end(); pIterSVData++)
   {
      // Compute row index for design matrix
      ulRowIndex = pIterSVData - pvEpochSVData_->begin();
      // Compute row components of design matrix
      dDistX = pIterSVData->dPRNX - pclPointExpansion_->GetComponent(0, 0);
      dDistY = pIterSVData->dPRNY - pclPointExpansion_->GetComponent(1, 0);
      dDistZ = pIterSVData->dPRNZ - pclPointExpansion_->GetComponent(2, 0);
      dRange = sqrt(pow(dDistX, 2) + pow(dDistY, 2) + pow(dDistZ, 2));
      // Design Matrix
      clA.SetComponent(ulRowIndex, 0, dDistX / dRange);
      clA.SetComponent(ulRowIndex, 1, dDistY / dRange);
      clA.SetComponent(ulRowIndex, 2, dDistZ / dRange);
      clA.SetComponent(ulRowIndex, 3, -1);
   }

   return clA;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeObsCovMatrix(vector<SatelliteData>* pvEpochSVData_, DOUBLE* pdConstants_)
{
   ULONG ulRowIndex = 0;
   DOUBLE dVariance = 0;
   CMatrix clR("Observation Covaraince Matrix", pvEpochSVData_->size(), pvEpochSVData_->size());
   // Build the observation covariance matrix 
   for (vector<SatelliteData>::iterator pIterSVData = pvEpochSVData_->begin(); pIterSVData != pvEpochSVData_->end(); pIterSVData++)
   {
      ulRowIndex = pIterSVData - pvEpochSVData_->begin();
      dVariance = pdConstants_[DATATYPEDEFINITIONS_APRIORIVARIANCE_INDEX] / sin(DegToRad(pIterSVData->dElevation));
      clR.SetComponent(ulRowIndex, ulRowIndex, pow(dVariance, 2));
   }

   return clR;
}
