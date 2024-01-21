//----------------------------------------------------------------------
// LSCommon.cpp
// Source file of LSCommon.h
// Source code for Least Squares positioning algorithm
//----------------------------------------------------------------------
#include "..\inc\LSCommon.h"

//------------------------------------------------------------------------------------------------
void LeastSquaresProcess(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, DOUBLE* pdConstants_, string sKalmanModel_)
{
   const ULONG ulMEASUREMENT_INPUTFILE_INDEX = 0;
   const ULONG ulTRUTHPOSITION_INPUTFILE_INDEX = 1;
   const ULONG ulOUTPUTFILES_DIRECTORY_INDEX = 2;
   vector<string> vInputFileKeys;
   //# Copy input file keys
   // Measurement file
   vInputFileKeys.push_back(pvFilePaths_->at(ulMEASUREMENT_INPUTFILE_INDEX));
   // Truth position file
   vInputFileKeys.push_back(pvFilePaths_->at(ulTRUTHPOSITION_INPUTFILE_INDEX));
   // Epoch satellite data vector
   vector<SatelliteData> vEpochSatelliteData;
   // Truth Data structure
   TruthData stTruthData;
   // Keep track of GNSS time
   DOUBLE dCurrGNSSTime = 0;
   DOUBLE dNextGNSSTime = 0;
   DOUBLE dGNSSMeasTime = 0;
   // Flags indicate to continue reading at file
   BOOLEANO bReadMeas = TRUE;
   BOOLEANO bReadTruth = TRUE;
   // Store the current position
   CMatrix clCurrPos("Geoid Point of Expansion", 4, 1);
   clCurrPos.SetComponent(0, 0, DegToRad(pdConstants_[DATATYPEDEFINITIONS_LATITUDEPOINTEXPANSION_INDEX]));
   clCurrPos.SetComponent(1, 0, DegToRad(pdConstants_[DATATYPEDEFINITIONS_LONGITUDEPOINTEXPANSION_INDEX]));
   clCurrPos.SetComponent(2, 0, pdConstants_[DATATYPEDEFINITIONS_HEIGHTPOINTEXPANSION_INDEX]);
   clCurrPos.SetComponent(3, 0, pdConstants_[DATATYPEDEFINITIONS_TIMESHIFTPOINTEXPANSION_INDEX]);
   // Ensure files are properly opened
   if (bOpenFiles(&vInputFileKeys, pmpFiles_, ios::in | ios::binary))
   {
      cout << "GPS Time, Latitude, Longitude, Height" << endl;
      while (bReadMeas && bReadTruth)
      {
         dCurrGNSSTime = dNextGNSSTime;
         bReadMeas = bReadMeasEpoch(pmpFiles_, &vInputFileKeys, &vEpochSatelliteData, dCurrGNSSTime, dNextGNSSTime);
         dGNSSMeasTime = dCurrGNSSTime;
         bReadTruth = bReadTruthEpoch(pmpFiles_, &vInputFileKeys, &stTruthData, dCurrGNSSTime, dGNSSMeasTime);
         clCurrPos = SPEpochSolution(pmpFiles_, &clCurrPos, &vEpochSatelliteData, &stTruthData, pvFilePaths_->at(ulOUTPUTFILES_DIRECTORY_INDEX), pdConstants_);
         // Print Position Coordinates to confirm thread is alive
         cout << fixed << setprecision(3) << vEpochSatelliteData.at(0).dGNSSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clCurrPos.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clCurrPos.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(3) << clCurrPos.GetComponent(2, 0) << endl;
         // Clear epoch data list
         vEpochSatelliteData.clear();
      }
   }
}

//------------------------------------------------------------------------------------------------
CMatrix SPEpochSolution(map<string, fstream*>* pmpFiles_, CMatrix* pclCurrPos_, vector<SatelliteData>* pvEpochData_, TruthData* pstTruthData_, string sOutDirPath_, DOUBLE* pdConstants_)
{
   // Output files full string
   string sErrOut = sOutDirPath_;
   string sPosOut = sOutDirPath_;
   string sTrajOut = sOutDirPath_;
   string sTruthOut = sOutDirPath_;
   sErrOut.append("\\LSErrorENU.GIT");
   sPosOut.append("\\LSPosition.GIT");
   sTrajOut.append("\\LSTrajectory.GIT");
   sTruthOut.append("\\LSTruth.GIT");
   vector<ULONG> vEpochPRN;
   vector<string> vOutFiles;
   // For plotting data trends
   const ULONG ulNUM_ERROR_TRENDS = 10;
   const ULONG ulNUM_POSITION_TRENDS = 7;
   const ULONG ulNUM_TRENDS_SET1 = 4;
   const ULONG ulLATITUDE_VS_LONGITUDE = 2;
   DOUBLE adErrorValues[ulNUM_ERROR_TRENDS];
   ULONG aulErrorPrecision[ulNUM_ERROR_TRENDS];
   DOUBLE adPositionValues[ulNUM_POSITION_TRENDS];
   ULONG aulPositionPrecision[ulNUM_POSITION_TRENDS];

   // Store the full paths to output files
   vOutFiles.push_back(sErrOut);
   vOutFiles.push_back(sPosOut);
   vOutFiles.push_back(sTrajOut);
   vOutFiles.push_back(sTruthOut);
   // Open output files
   if (bOpenFiles(&vOutFiles, pmpFiles_, ios::out, sOutDirPath_))
   {
      PrintLSErrorHeader(pmpFiles_->at(sErrOut));
      PrintLSPositionHeader(pmpFiles_->at(sPosOut));
      PrintLSTrajectoryHeader(pmpFiles_->at(sTrajOut));
      PrintLSTrajectoryHeader(pmpFiles_->at(sTruthOut));
   }

   // ECEF coordinates point of expansion and rotation matrix from ECEF to Local Level Frame (LLF)
   CMatrix clECEFCoord("ECEF Point of Expansion", pclCurrPos_->GetNumRows(), 1);
   CMatrix clECEFPosError("ECEF Position Error", 4, 1);
   CMatrix clLLFPosError("LLF Position Error", 4, 1);
   CMatrix clLLFPosCov("LLF Error Covariance", 4, 4);
   CMatrix clRLE(ECEFToLLF(pclCurrPos_));
   CMatrix clCurrPos(*pclCurrPos_);
   GeoidToECEF(pclCurrPos_, &clECEFCoord);

   // Misclosures
   CMatrix clMisc(ComputeMisclosure(&clECEFCoord, pvEpochData_, &vEpochPRN));
   // Geometry matrix
   CMatrix clA(ComputeDesignMatrix(&clECEFCoord, pvEpochData_));
   // Observation covariance matrix
   CMatrix clCL(ComputeObsCovMatrix(pvEpochData_, pdConstants_));
   // State covariance matrix
   CMatrix clCX(ComputeStateCovariance(&clCL, &clA));
   // Error correction vector
   CMatrix clDelta(ComputeStateCorrection(&clCL, &clA, &clMisc));

   DOUBLE dErrorNorm2 = sqrt(pow(clDelta.GetComponent(0, 0), 2) + pow(clDelta.GetComponent(1, 0), 2) + pow(clDelta.GetComponent(2, 0), 2) + pow(clDelta.GetComponent(3, 0), 2));

   while (dErrorNorm2 > clCX.MatrixNorm2() * LSCOMMON_ERROR_TOLERANCE)
   {
      // Update point of expansion
      clECEFCoord = clECEFCoord + clDelta;
      // Misclosures
      clMisc = ComputeMisclosure(&clECEFCoord, pvEpochData_, &vEpochPRN);
      // Geometry matrix
      clA = ComputeDesignMatrix(&clECEFCoord, pvEpochData_);
      // Observation covariance matrix
      clCL = ComputeObsCovMatrix(pvEpochData_, pdConstants_);
      // State covariance matrix
      clCX = ComputeStateCovariance(&clCL, &clA);
      // Error correction vector
      clDelta = ComputeStateCorrection(&clCL, &clA, &clMisc);
      // Compute norm of state error vector
      dErrorNorm2 = sqrt(pow(clDelta.GetComponent(0, 0), 2) + pow(clDelta.GetComponent(1, 0), 2) + pow(clDelta.GetComponent(2, 0), 2) + pow(clDelta.GetComponent(3, 0), 2));
   }

   // Recompute geoid curvilinear coordinates
   ECEFToGeoid(&clECEFCoord, pclCurrPos_);
   clRLE = ECEFToLLF(pclCurrPos_);
   clCurrPos.SetComponent(0, 0, RadToDeg(pclCurrPos_->GetComponent(0, 0)));
   clCurrPos.SetComponent(1, 0, RadToDeg(pclCurrPos_->GetComponent(1, 0)));
   clCurrPos.SetComponent(2, 0, pclCurrPos_->GetComponent(2, 0));
   clCurrPos.SetComponent(3, 0, pclCurrPos_->GetComponent(3, 0));
   // ECEF position error
   clECEFPosError.SetComponent(0, 0, pstTruthData_->dXCoord - clECEFCoord.GetComponent(0, 0));
   clECEFPosError.SetComponent(1, 0, pstTruthData_->dYCoord - clECEFCoord.GetComponent(1, 0));
   clECEFPosError.SetComponent(2, 0, pstTruthData_->dZCoord - clECEFCoord.GetComponent(2, 0));
   clECEFPosError.SetComponent(3, 0, clECEFCoord.GetComponent(3, 0));
   // LLF position error
   clLLFPosError = clRLE * clECEFPosError;
   clLLFPosCov = clRLE * clCX * clRLE.Transpose();
   // Truth trajectory matrices
   CMatrix clTruthGeoid("Truth LLF Position", 4, 4);
   CMatrix clTruthECEF("Truth ECEF Position", 4, 4);
   // Truth Data
   clTruthECEF.SetComponent(0, 0, pstTruthData_->dXCoord);
   clTruthECEF.SetComponent(1, 0, pstTruthData_->dYCoord);
   clTruthECEF.SetComponent(2, 0, pstTruthData_->dZCoord);
   ECEFToGeoid(&clTruthECEF, &clTruthGeoid);
   // Set Geoid curvilinear position vector in degrees
   clTruthGeoid.SetComponent(0, 0, RadToDeg(clTruthGeoid.GetComponent(0, 0)));
   clTruthGeoid.SetComponent(1, 0, RadToDeg(clTruthGeoid.GetComponent(1, 0)));

   // Fill data arrays
   adErrorValues[0] = pstTruthData_->dGNSSTime;
   adPositionValues[0] = pstTruthData_->dGNSSTime;
   aulErrorPrecision[0] = 3;
   aulPositionPrecision[0] = 3;
   // Position errors
   for (ULONG ulIndex = 1; ulIndex < ulNUM_POSITION_TRENDS; ulIndex++)
   {
      // Error plots
      if (ulIndex < ulNUM_TRENDS_SET1)
      {
         adErrorValues[ulIndex] = clLLFPosError.GetComponent((INT)(ulIndex - 1), 0);
         aulErrorPrecision[ulIndex] = LSCOMMON_ERRORDECIMAL_SIGNIFICANCE;
         adPositionValues[ulIndex] = clCurrPos.GetComponent((INT)(ulIndex - 1), 0);
      }
      // Truth position plots
      else
      {
         adPositionValues[ulIndex] = clTruthGeoid.GetComponent((INT)(ulIndex - ulNUM_TRENDS_SET1), 0);
      }
      aulPositionPrecision[ulIndex] = LSCOMMON_STDEVDECIMAL_SIGNIFICANCE;
   }
   // Position standard deviations estimations
   for (ULONG ulIndex = 0; ulIndex < ulNUM_TRENDS_SET1 - 1; ulIndex++)
   {
      adErrorValues[2 * ulIndex + ulNUM_TRENDS_SET1] = sqrt(clLLFPosCov.GetComponent((INT)ulIndex, (INT)ulIndex));
      adErrorValues[2 * ulIndex + ulNUM_TRENDS_SET1 + 1] = -sqrt(clLLFPosCov.GetComponent((INT)ulIndex, (INT)ulIndex));
      aulErrorPrecision[2 * ulIndex + ulNUM_TRENDS_SET1] = LSCOMMON_STDEVDECIMAL_SIGNIFICANCE;
      aulErrorPrecision[2 * ulIndex + ulNUM_TRENDS_SET1 + 1] = LSCOMMON_STDEVDECIMAL_SIGNIFICANCE;
   }

   // Output files
   PrintLSPlots(pmpFiles_->at(sErrOut), &adErrorValues[0], &aulErrorPrecision[0], ulNUM_ERROR_TRENDS);
   PrintLSPlots(pmpFiles_->at(sPosOut), &adPositionValues[0], &aulPositionPrecision[0], ulNUM_POSITION_TRENDS);
   PrintLSPlots(pmpFiles_->at(sTrajOut), &adPositionValues[1], &aulPositionPrecision[1], ulLATITUDE_VS_LONGITUDE);
   PrintLSPlots(pmpFiles_->at(sTruthOut), &adPositionValues[4], &aulPositionPrecision[4], ulLATITUDE_VS_LONGITUDE);

   return *pclCurrPos_;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeStateCovariance(CMatrix* pclObsCovMatrix_, CMatrix* pclDesignMatrixA_)
{
   // Normal matrix (A^T * Cl^(-1) * A)
   CMatrix clN(pclDesignMatrixA_->Transpose() * pclObsCovMatrix_->NumericInverse2() * (*pclDesignMatrixA_));
   // State covariance matrix (A^T * Cl^(-1) * A) ^ (-1)
   CMatrix clStateCov(clN.NumericInverse2());

   return clStateCov;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeStateCorrection(CMatrix* pclObsCovMatrix_, CMatrix* pclDesignMatrixA_, CMatrix* pclMisclosure_)
{
   // Zero matrix
   CMatrix clZero("Zero Matrix", pclDesignMatrixA_->GetNumCols(), 1);
   clZero.SetZero();
   // Normal matrix
   CMatrix clN(pclDesignMatrixA_->Transpose() * pclObsCovMatrix_->NumericInverse2() * (*pclDesignMatrixA_));
   // Correction vector
   CMatrix clDelta(clN.NumericInverse2() * pclDesignMatrixA_->Transpose() * pclObsCovMatrix_->NumericInverse2() * (*pclMisclosure_));
   clDelta = clZero - clDelta;
   return clDelta;
}
