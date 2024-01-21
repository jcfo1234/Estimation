//----------------------------------------------------------------------
// KalmanCommon.cpp
// Source file of KalmanCommon.h
// Source code for Kalman filter positioning algorithm
//----------------------------------------------------------------------
#include "..\inc\KalmanCommon.h"

//------------------------------------------------------------------------------------------------
CMatrix InitialStateVector(string sKalmanModel_, DOUBLE* pdConstants_)
{
   map<string, INT> mpDimension;
   mpDimension.operator[]("POSRANDWALK") = 4;
   mpDimension.operator[]("VELRANDWALK") = 8;
   mpDimension.operator[]("GAUSSMARKOV") = 8;

   // Search for the Kalman state matrix configuration model. Abort if not found
   if (mpDimension.find(sKalmanModel_) == mpDimension.end())
   {
      throw 9;
   }

   CMatrix clCurrPos("Position Vector", 4, 1);
   CMatrix clECEFPos("ECEF Position Vector", 4, 1);
   CMatrix clStateVector("Position State Vector", mpDimension.at(sKalmanModel_), 1);
   clCurrPos.SetComponent(0, 0, DegToRad(pdConstants_[DATATYPEDEFINITIONS_LATITUDEPOINTEXPANSION_INDEX]));
   clCurrPos.SetComponent(1, 0, DegToRad(pdConstants_[DATATYPEDEFINITIONS_LONGITUDEPOINTEXPANSION_INDEX]));
   clCurrPos.SetComponent(2, 0, pdConstants_[DATATYPEDEFINITIONS_HEIGHTPOINTEXPANSION_INDEX]);
   clCurrPos.SetComponent(3, 0, pdConstants_[DATATYPEDEFINITIONS_TIMESHIFTPOINTEXPANSION_INDEX]);
   GeoidToECEF(&clCurrPos, &clECEFPos);
   // Fill position states (Latitude, Longitude, Height, Clock bias)
   clStateVector.SetComponent(0, 0, clECEFPos.GetComponent(0, 0));
   clStateVector.SetComponent(1, 0, clECEFPos.GetComponent(1, 0));
   clStateVector.SetComponent(2, 0, clECEFPos.GetComponent(2, 0));
   clStateVector.SetComponent(3, 0, clECEFPos.GetComponent(3, 0));
   // Velocity random walk or gauss markov Kalman models
   if (sKalmanModel_ != "POSRANDWALK")
   {
      // ECEF Velocities X, Y, Z and Clock Drift
      clStateVector.SetComponent(4, 0, pdConstants_[DATATYPEDEFINITIONS_VELXECEFPOINTEXPANSION_INDEX]);
      clStateVector.SetComponent(5, 0, pdConstants_[DATATYPEDEFINITIONS_VELYECEFPOINTEXPANSION_INDEX]);
      clStateVector.SetComponent(6, 0, pdConstants_[DATATYPEDEFINITIONS_VELZECEFPOINTEXPANSION_INDEX]);
      clStateVector.SetComponent(7, 0, pdConstants_[DATATYPEDEFINITIONS_CLOCKDRIFTPOINTEXPANSION_INDEX]);
   }
   mpDimension.clear();
   return clStateVector;
}

//------------------------------------------------------------------------------------------------
CMatrix InitialStateCovP(string sKalmanModel_, DOUBLE* pdConstants_)
{
   map<string, INT> mpDimension;
   mpDimension.operator[]("POSRANDWALK") = 4;
   mpDimension.operator[]("VELRANDWALK") = 8;
   mpDimension.operator[]("GAUSSMARKOV") = 8;

   // Search for the Kalman state matrix configuration model. Abort if not found
   if (mpDimension.find(sKalmanModel_) == mpDimension.end())
   {
      throw 10;
   }

   CMatrix clCovP("State Covariance Matrix", mpDimension.at(sKalmanModel_), mpDimension.at(sKalmanModel_));
   clCovP.SetZero();

   clCovP.SetComponent(0, 0, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARX_INDEX], 2));
   clCovP.SetComponent(1, 1, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARY_INDEX], 2));
   clCovP.SetComponent(2, 2, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARZ_INDEX], 2));
   clCovP.SetComponent(3, 3, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARCLKBIAS_INDEX], 2));
   // Velocity states covariance
   if (sKalmanModel_ != "POSRANDWALK")
   {
      clCovP.SetComponent(4, 4, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARVELX_INDEX], 2));
      clCovP.SetComponent(5, 5, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARVELY_INDEX], 2));
      clCovP.SetComponent(6, 6, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARVELZ_INDEX], 2));
      clCovP.SetComponent(7, 7, pow(pdConstants_[DATATYPEDEFINITIONS_INITVARCLKDRIFT_INDEX], 2));
   }
   return clCovP;
}

//------------------------------------------------------------------------------------------------
void TransformStateToLLF(CMatrix* pclECEFState_, CMatrix* pclECEFCovP_, CMatrix* pclGeoidState_, CMatrix* pclLLFCovP_, string sKalmanModel_)
{
   CMatrix clECEFPos(pclECEFState_->SubMatrix(*pclECEFState_, 0, 4, 0, 1));
   CMatrix clECEFVel("ECEF Velocity", 4, 1);
   CMatrix clLLFPos("LLF Position", 4, 1);
   CMatrix clLLFVel("LLF Velocity", 4, 1);
   // Local Level Frame coordinates
   ECEFToGeoid(&clECEFPos, &clLLFPos);
   // Transformation matrix
   CMatrix clRLE(ECEFToLLF(&clLLFPos));
   pclGeoidState_->SetComponent(0, 0, clLLFPos.GetComponent(0, 0));
   pclGeoidState_->SetComponent(1, 0, clLLFPos.GetComponent(1, 0));
   pclGeoidState_->SetComponent(2, 0, clLLFPos.GetComponent(2, 0));
   pclGeoidState_->SetComponent(3, 0, clLLFPos.GetComponent(3, 0));
   // Velocity state augmentation
   if (sKalmanModel_ != "POSRANDWALK")
   {
      clECEFVel = pclECEFState_->SubMatrix(*pclECEFState_, 4, 8, 0, 1);
      clLLFVel = clRLE * clECEFVel;
      clRLE.SetWarningEnable(FALSE);
      clRLE = clRLE.AugmentDim(clRLE, 4, 4);
      // Local Level Frame state vector
      pclGeoidState_->SetComponent(4, 0, clLLFVel.GetComponent(0, 0));
      pclGeoidState_->SetComponent(5, 0, clLLFVel.GetComponent(1, 0));
      pclGeoidState_->SetComponent(6, 0, clLLFVel.GetComponent(2, 0));
      pclGeoidState_->SetComponent(7, 0, clLLFVel.GetComponent(3, 0));
   }
   // Local Level Frame covariance matrix
   *pclLLFCovP_ = clRLE * (*pclECEFCovP_) * clRLE.Transpose();
}

//------------------------------------------------------------------------------------------------
CMatrix StateMatrix(string sKalmanModel_, DOUBLE* pdConstants_)
{
   map<string, CMatrix*> mpStateMatrix;
   // State matrices
   CMatrix clF_PositionRandomWalk("PRW State Matrix", 4, 4);
   CMatrix clF_VelocityRandomWalk("VRW State Matrix", 8, 8);
   CMatrix clF_GaussMarkov("GM State Matrix", 8, 8);
   // Position Random Walk Matrix
   clF_PositionRandomWalk.SetZero();
   // Velocity Random Walk Matrix
   clF_VelocityRandomWalk.SetZero();
   clF_VelocityRandomWalk.SetComponent(0, 4, 1);
   clF_VelocityRandomWalk.SetComponent(1, 5, 1);
   clF_VelocityRandomWalk.SetComponent(2, 6, 1);
   clF_VelocityRandomWalk.SetComponent(3, 7, 1);
   // Gauss Markov Matrix
   clF_GaussMarkov = clF_VelocityRandomWalk;
   clF_GaussMarkov.SetComponent(4, 4, -1 / pdConstants_[DATATYPEDEFINITIONS_CORRVXTIMECONSTANT_INDEX]);
   clF_GaussMarkov.SetComponent(5, 5, -1 / pdConstants_[DATATYPEDEFINITIONS_CORRVYTIMECONSTANT_INDEX]);
   clF_GaussMarkov.SetComponent(6, 6, -1 / pdConstants_[DATATYPEDEFINITIONS_CORRVZTIMECONSTANT_INDEX]);
   clF_GaussMarkov.SetComponent(7, 7, -1 / pdConstants_[DATATYPEDEFINITIONS_CORRCLKDRIFTTIMECONSTANT_INDEX]);
   // Matrix Map
   mpStateMatrix.operator[]("POSRANDWALK") = &clF_PositionRandomWalk;
   mpStateMatrix.operator[]("VELRANDWALK") = &clF_VelocityRandomWalk;
   mpStateMatrix.operator[]("GAUSSMARKOV") = &clF_GaussMarkov;
   return (*mpStateMatrix.at(sKalmanModel_));
}

//------------------------------------------------------------------------------------------------
CMatrix InputMatrix(string sKalmanModel_)
{
   map<string, CMatrix*> mpInputMatrix;
   // Spectral density shape input matrices
   CMatrix clG_PositionRandomWalk("PRW Input Matrix", 4, 4);
   CMatrix clG_VelocityRandomWalk("VRW Input Matrix", 8, 4);
   CMatrix clG_GaussMarkov("GM Input Matrix", 8, 4);
   // Position random walk input matrix
   clG_PositionRandomWalk.SetIdentity();
   // Velocity random walk input matrix;
   clG_VelocityRandomWalk.SetComponent(4, 0, 1);
   clG_VelocityRandomWalk.SetComponent(5, 1, 1);
   clG_VelocityRandomWalk.SetComponent(6, 2, 1);
   clG_VelocityRandomWalk.SetComponent(7, 3, 1);
   // Gauss Markov input matrix
   clG_GaussMarkov = clG_VelocityRandomWalk;
   // Matrix Map
   mpInputMatrix.operator[]("POSRANDWALK") = &clG_PositionRandomWalk;
   mpInputMatrix.operator[]("VELRANDWALK") = &clG_VelocityRandomWalk;
   mpInputMatrix.operator[]("GAUSSMARKOV") = &clG_GaussMarkov;
   return (*mpInputMatrix.at(sKalmanModel_));
}

//----------------------------------------------------------------------------
CMatrix SpectralDensity(string sKalmanModel_, DOUBLE* pdConstants_)
{
   map<string, CMatrix*> mpSpectralMatrix;
   // Spectral density shape input matrices
   CMatrix clQ_PositionRandomWalk("PRW Spectral Densities", 4, 4);
   CMatrix clQ_VelocityRandomWalk("VRW Spectral Densities", 4, 4);
   CMatrix clQ_GaussMarkov("GM Spectral Densities", 4, 4);
   // Position random walk input matrix
   clQ_PositionRandomWalk.SetZero();
   clQ_PositionRandomWalk.SetComponent(0, 0, pow(pdConstants_[DATATYPEDEFINITIONS_PRWXSPECTRALDENSITY_INDEX], 2));
   clQ_PositionRandomWalk.SetComponent(1, 1, pow(pdConstants_[DATATYPEDEFINITIONS_PRWYSPECTRALDENSITY_INDEX], 2));
   clQ_PositionRandomWalk.SetComponent(2, 2, pow(pdConstants_[DATATYPEDEFINITIONS_PRWZSPECTRALDENSITY_INDEX], 2));
   clQ_PositionRandomWalk.SetComponent(3, 3, pdConstants_[DATATYPEDEFINITIONS_PRWCLKBIASSPECTRALDENSITY_INDEX]);
   // Velocity random walk input matrix;
   clQ_VelocityRandomWalk.SetZero();
   clQ_VelocityRandomWalk.SetComponent(0, 0, pow(pdConstants_[DATATYPEDEFINITIONS_VRWXSPECTRALDENSITY_INDEX], 2));
   clQ_VelocityRandomWalk.SetComponent(1, 1, pow(pdConstants_[DATATYPEDEFINITIONS_VRWYSPECTRALDENSITY_INDEX], 2));
   clQ_VelocityRandomWalk.SetComponent(2, 2, pow(pdConstants_[DATATYPEDEFINITIONS_VRWZSPECTRALDENSITY_INDEX], 2));
   clQ_VelocityRandomWalk.SetComponent(3, 3, pdConstants_[DATATYPEDEFINITIONS_VRWCLKDRIFTSPECTRALDENSITY_INDEX]);
   // Gauss Markov input matrix
   clQ_GaussMarkov = clQ_VelocityRandomWalk;
   // Matrix Map
   mpSpectralMatrix.operator[]("POSRANDWALK") = &clQ_PositionRandomWalk;
   mpSpectralMatrix.operator[]("VELRANDWALK") = &clQ_VelocityRandomWalk;
   mpSpectralMatrix.operator[]("GAUSSMARKOV") = &clQ_GaussMarkov;
   return (*mpSpectralMatrix.at(sKalmanModel_));
}

//------------------------------------------------------------------------------------------------
CMatrix ProcessNoiseQ(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_, DOUBLE dStepSize_)
{
   CMatrix clQ("Process Noise", pclF_->GetNumRows(), pclF_->GetNumCols());
   // Scaled versions of the state matrix
   CMatrix clF0("Scaled by Zero", pclF_->GetNumRows(), pclF_->GetNumCols());
   CMatrix clFHalfStep("Scaled by half step size", pclF_->GetNumRows(), pclF_->GetNumCols());
   CMatrix clFStep("Scaled by step size", pclF_->GetNumRows(), pclF_->GetNumCols());
   // Scale factors used in the Simpson integration rule
   const DOUBLE dHalfStepSize = dStepSize_ / 2;
   const DOUBLE dSimpsonFactor = dStepSize_ / 6;
   // Scale the state matrix
   clF0 = *pclF_ * 0;
   clFHalfStep = *pclF_ * dHalfStepSize;
   clFStep = *pclF_ * dStepSize_;
   // Compute the scaled state transition matrices
   CMatrix clPHI0(clF0.MatrixExponential2());
   CMatrix clPHIHalfStep(clFHalfStep.MatrixExponential2());
   CMatrix clPHIStep(clFStep.MatrixExponential2());
   // Compute the scaled process noise
   CMatrix clQ0(clPHI0 * (*pclG_) * (*pclQ_) * pclG_->Transpose() * clPHI0.Transpose());
   CMatrix clQHalfStep(clPHIHalfStep * (*pclG_) * (*pclQ_) * pclG_->Transpose() * clPHIHalfStep.Transpose());
   CMatrix clQStep(clPHIStep * (*pclG_) * (*pclQ_) * pclG_->Transpose() * clPHIStep.Transpose());
   // Compute process Noise by the Simpson integration rule
   CMatrix clTempQ(clQ0 + clQHalfStep * 4 + clQStep);
   clQ = clTempQ * dSimpsonFactor;

   return clQ;
}

//----------------------------------------------------------------------------
CMatrix KalmanGain(CMatrix* pclP_, CMatrix* pclH_, CMatrix* pclR_)
{
   CMatrix clTemp("Temporal", pclR_->GetNumRows(), pclR_->GetNumCols());
   CMatrix clK("Kalman Gain", pclP_->GetNumRows(), pclR_->GetNumCols());

   clTemp = *pclH_ * (*pclP_) * pclH_->Transpose() + (*pclR_);
   clK = *pclP_ * pclH_->Transpose() * clTemp.NumericInverse2(1e-6);

   return clK;
}

//------------------------------------------------------------------------------------------------
void IEKF(CMatrix* pclPredictedState_,
          CMatrix* pclPredStateCovP_,
          CMatrix* pclObsCovR_,
          CMatrix* pclUpdatedState_,
          CMatrix* pclUpdatedCovP_,
          vector<SatelliteData>* pvEpochData_)
{
   CMatrix clIterState(*pclUpdatedState_);
   CMatrix clOldIterState(*pclPredictedState_);
   CMatrix clIterDiff("Iteration Difference", pclPredictedState_->GetNumRows(), pclPredictedState_->GetNumCols());
   CMatrix clIdentity("Identity", pclUpdatedCovP_->GetNumRows(), pclUpdatedCovP_->GetNumCols());
   clIdentity.SetIdentity();
   vector<ULONG> vEpochPRN;
   // Misclosures
   CMatrix clIterMisc(ComputeMisclosure(&clIterState, pvEpochData_, &vEpochPRN));
   // Geometry matrix
   CMatrix clIterA(ComputeDesignMatrix(&clIterState, pvEpochData_));
   // Kalman Gain
   CMatrix clIterK(KalmanGain(pclPredStateCovP_, &clIterA, pclObsCovR_));
   // Old iteration and current iteration difference
   clIterDiff = clIterState - clOldIterState;
   // Iterate for refining the state solution
   while (clIterDiff.MatrixNorm1() > KALMANCOMMON_IEKF_TOLERANCE)
   {
      clOldIterState = clIterState;
      clIterA = ComputeDesignMatrix(&clIterState, pvEpochData_);
      clIterMisc = ComputeMisclosure(&clIterState, pvEpochData_, &vEpochPRN);
      clIterK = KalmanGain(pclPredStateCovP_, &clIterA, pclObsCovR_);
      clIterState = *pclPredictedState_ - clIterK * (clIterMisc + clIterA * (*pclPredictedState_ - clIterState));
      clIterDiff = clIterState - clOldIterState;
   }
   // Updated states and covariances
   *pclUpdatedState_ = clIterState;
   *pclUpdatedCovP_ = (clIdentity - clIterK * clIterA) * (*pclPredStateCovP_);
}

//------------------------------------------------------------------------------------------------
CMatrix KalmanEpochSolution(string sKalmanModel_, map<string, fstream*>* pmpFiles_, 
                            vector<SatelliteData>* pvEpochData_, TruthData* pstTruthData_, 
                            string sOutDirPath_, DOUBLE* pdConstants_, 
                            CMatrix* pclStateVector_, CMatrix* pclStateCovP_,
                            DOUBLE& dTimeStep_)
{
   // Output files full string
   string sErrOut = sOutDirPath_;
   string sPosOut = sOutDirPath_;
   string sTrajOut = sOutDirPath_;
   string sTruthOut = sOutDirPath_;
   sErrOut.append("\\KalmanErrorENU.GIT");
   sPosOut.append("\\KalmanPosition.GIT");
   sTrajOut.append("\\KalmanTrajectory.GIT");
   sTruthOut.append("\\KalmanTruth.GIT");
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

   // [X,Y,Z,c*DeltaT,V_X,V_Y,V_Z,(c*DeltaT)^Dot]
   CMatrix clECEFState("ECEF Point of Expansion", pclStateVector_->GetNumRows(), pclStateVector_->GetNumCols());
   CMatrix clPMinus("Apriori State Covariance", pclStateCovP_->GetNumRows(), pclStateCovP_->GetNumCols());
   // Geoid state space
   CMatrix clGeoidState("Geoid State Space", pclStateVector_->GetNumRows(), pclStateVector_->GetNumCols());
   CMatrix clGeoidP("Geoid State Covariance", pclStateCovP_->GetNumRows(), pclStateCovP_->GetNumCols());

   // State power spectral density matrix
   CMatrix clSpectralDensityQ(SpectralDensity(sKalmanModel_, pdConstants_));
   // State matrix
   CMatrix clF(StateMatrix(sKalmanModel_, pdConstants_));
   // Input noise state matrix
   CMatrix clG(InputMatrix(sKalmanModel_));
   // Identity Matrix
   CMatrix clI("Identity", pclStateCovP_->GetNumRows(), pclStateCovP_->GetNumCols());
   clI.SetIdentity();

   // Prediction loop
   CMatrix clPHI((clF * dTimeStep_).MatrixExponential2());
   CMatrix clQ(ProcessNoiseQ(&clF, &clG, &clSpectralDensityQ, dTimeStep_));
   clECEFState = clPHI * (*pclStateVector_);
   clPMinus = clPHI * (*pclStateCovP_) * clPHI.Transpose() + clQ;

   // Misclosures
   CMatrix clMisc(ComputeMisclosure(&clECEFState, pvEpochData_, &vEpochPRN));
   // Geometry matrix
   CMatrix clA(ComputeDesignMatrix(&clECEFState, pvEpochData_));
   // Observation covariance matrix
   CMatrix clR(ComputeObsCovMatrix(pvEpochData_, pdConstants_));
   // Kalman Gain
   CMatrix clK(KalmanGain(&clPMinus, &clA, &clR));

   // Correction loop
   CMatrix clTemp(clI - clK * clA);
   // Corrected state vector
   *pclStateVector_ = clECEFState - clK * clMisc;
   // Corrected state covariance matrix
   *pclStateCovP_ = clTemp * clPMinus;

   // Perform Iterations if IEKF is enabled
   if (pdConstants_[DATATYPEDEFINITIONS_IEKFENABLE_INDEX] > 0)
   {
      IEKF(&clECEFState, &clPMinus, &clR, pclStateVector_, pclStateCovP_, pvEpochData_);
   }

   // Print to output files
   TransformStateToLLF(pclStateVector_, pclStateCovP_, &clGeoidState, &clGeoidP, sKalmanModel_);
   CMatrix clGeoidPos(clGeoidState.SubMatrix(clGeoidState, 0, 4, 0, 1));
   CMatrix clTruthGeoid("Truth LLF Position", 4, 4);
   CMatrix clTruthECEF("Truth ECEF Position", 4, 4);
   CMatrix clRLE(ECEFToLLF(&clGeoidPos));
   CMatrix clECEFError("ECEF Position Error Vector", 4, 1);
   CMatrix clLLFError("LLF Position Error Vector", 4, 1);
   // Set position error components in ECEF
   clECEFError.SetComponent(0, 0, pstTruthData_->dXCoord - pclStateVector_->GetComponent(0, 0));
   clECEFError.SetComponent(1, 0, pstTruthData_->dYCoord - pclStateVector_->GetComponent(1, 0));
   clECEFError.SetComponent(2, 0, pstTruthData_->dZCoord - pclStateVector_->GetComponent(2, 0));
   clECEFError.SetComponent(3, 0, pclStateVector_->GetComponent(3, 0));
   // Truth Data
   clTruthECEF.SetComponent(0, 0, pstTruthData_->dXCoord);
   clTruthECEF.SetComponent(1, 0, pstTruthData_->dYCoord);
   clTruthECEF.SetComponent(2, 0, pstTruthData_->dZCoord);
   ECEFToGeoid(&clTruthECEF, &clTruthGeoid);
   // Set Geoid curvilinear position vector in degrees
   clGeoidPos.SetComponent(0, 0, RadToDeg(clGeoidState.GetComponent(0, 0)));
   clGeoidPos.SetComponent(1, 0, RadToDeg(clGeoidState.GetComponent(1, 0)));
   clTruthGeoid.SetComponent(0, 0, RadToDeg(clTruthGeoid.GetComponent(0, 0)));
   clTruthGeoid.SetComponent(1, 0, RadToDeg(clTruthGeoid.GetComponent(1, 0)));
   // Set position error components in LLF
   clLLFError = clRLE * clECEFError;

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
         adErrorValues[ulIndex] = clLLFError.GetComponent((INT)(ulIndex - 1), 0);
         aulErrorPrecision[ulIndex] = KALMANCOMMON_ERRORDECIMAL_SIGNIFICANCE;
         adPositionValues[ulIndex] = clGeoidPos.GetComponent((INT)(ulIndex - 1), 0);
      }
      // Truth position plots
      else
      {
         adPositionValues[ulIndex] = clTruthGeoid.GetComponent((INT)(ulIndex - ulNUM_TRENDS_SET1), 0);
      }
      aulPositionPrecision[ulIndex] = KALMANCOMMON_STDEVDECIMAL_SIGNIFICANCE;
   }
   // Position standard deviations estimations
   for (ULONG ulIndex = 0; ulIndex < ulNUM_TRENDS_SET1 - 1; ulIndex++)
   {
      adErrorValues[2 * ulIndex + ulNUM_TRENDS_SET1] = sqrt(clGeoidP.GetComponent((INT)ulIndex, (INT)ulIndex));
      adErrorValues[2 * ulIndex + ulNUM_TRENDS_SET1 + 1] = -sqrt(clGeoidP.GetComponent((INT)ulIndex, (INT)ulIndex));
      aulErrorPrecision[2 * ulIndex + ulNUM_TRENDS_SET1] = KALMANCOMMON_STDEVDECIMAL_SIGNIFICANCE;
      aulErrorPrecision[2 * ulIndex + ulNUM_TRENDS_SET1 + 1] = KALMANCOMMON_STDEVDECIMAL_SIGNIFICANCE;
   }

   // Output files
   PrintLSPlots(pmpFiles_->at(sErrOut), &adErrorValues[0], &aulErrorPrecision[0], ulNUM_ERROR_TRENDS);
   PrintLSPlots(pmpFiles_->at(sPosOut), &adPositionValues[0], &aulPositionPrecision[0], ulNUM_POSITION_TRENDS);
   PrintLSPlots(pmpFiles_->at(sTrajOut), &adPositionValues[1], &aulPositionPrecision[1], ulLATITUDE_VS_LONGITUDE);
   PrintLSPlots(pmpFiles_->at(sTruthOut), &adPositionValues[4], &aulPositionPrecision[4], ulLATITUDE_VS_LONGITUDE);

   return (*pclStateVector_);
}

//------------------------------------------------------------------------------------------------
void KalmanProcess(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, DOUBLE* pdConstants_, string sKalmanModel_)
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
   DOUBLE dTimeStepSize = 0;
   // Flags indicate to continue reading at file
   BOOLEANO bReadMeas = TRUE;
   BOOLEANO bReadTruth = TRUE;
   // Store the current state vector and covariance matrix
   CMatrix clCurrState(InitialStateVector(sKalmanModel_, pdConstants_));
   CMatrix clP(InitialStateCovP(sKalmanModel_, pdConstants_));
   CMatrix clGeoidState("LLF Position", clCurrState.GetNumRows(), clCurrState.GetNumCols());
   CMatrix clGeoidP("LLF Covariance", clP.GetNumRows(), clP.GetNumCols());
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
         // Compute the time step size for updating the state space and computing the state transition matrix
         dTimeStepSize = dNextGNSSTime - dCurrGNSSTime;
         clCurrState = KalmanEpochSolution(sKalmanModel_, pmpFiles_, &vEpochSatelliteData, &stTruthData, pvFilePaths_->at(ulOUTPUTFILES_DIRECTORY_INDEX), pdConstants_, &clCurrState, &clP, dTimeStepSize);
         // Transform to geodetical coordinates
         TransformStateToLLF(&clCurrState, &clP, &clGeoidState, &clGeoidP, sKalmanModel_);
         // Print Position Coordinates to confirm thread is alive
         cout << fixed << setprecision(3) << vEpochSatelliteData.at(0).dGNSSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidState.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidState.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(3) << clGeoidState.GetComponent(2, 0) << endl;
         // Clear epoch data list
         vEpochSatelliteData.clear();
      }
   }
}
