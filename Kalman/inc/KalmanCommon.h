//----------------------------------------------------------------------
// KalmanCommon.h
// Module that handles common least squares functions
//----------------------------------------------------------------------

#ifndef KALMANCOMMON_H
#define KALMANCOMMON_H

#include <stdio.h>
#include <tchar.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "..\..\DataTypeDefinitions.h"
#include "..\..\Fetch\inc\FetchCommon.h"
#include "..\..\Fetch\inc\CommonFunctions.h"
#include "..\..\Matrix\inc\Matrix.h"
using namespace std;

#define KALMANCOMMON_ERRORDECIMAL_SIGNIFICANCE ((ULONG) 9)
#define KALMANCOMMON_STDEVDECIMAL_SIGNIFICANCE ((ULONG) 6)
#define KALMANCOMMON_IEKF_TOLERANCE ((DOUBLE) 1e-3)

//----------------------------------------------------------------------------
// Reads measurements file and truth file to plot least squares position errors
// pvFilePaths_: Pointer to vector with input file paths
// pmpFiles_: Map from file paths to pointer to file for handling input/output
// pdConstatns_: Useful values for iterating on the Least squares algorithm
// sKalmanModel_: Type of dynamical model to use in Kalman state matrix
//----------------------------------------------------------------------------
void KalmanProcess(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, DOUBLE* pdConstants_, string sKalmanModel_);

//----------------------------------------------------------------------------
// Transforms state vector from ECEF to LLF coordinates
// pclECEFState_: Pointer to ECEF state vector
// pclECEFCovP_: Pointer to ECEF state covariance matrix P
// pclGeoidState_: Pointer to Geoid state vector (latitude, longitude, height)
// pclLLFCovP_: Pointer to Geoid state covariance matrix
// sKalmanModel_: Kalman model
//----------------------------------------------------------------------------
void TransformStateToLLF(CMatrix* pclECEFState_, CMatrix* pclECEFCovP_, CMatrix* pclGeoidState_, CMatrix* pclLLFCovP_, string sKalmanModel_);

//----------------------------------------------------------------------------
// Returns the initial state vector for Kalman processing
// sKalmanModel_: Kalman model
// pdConstants_: Constants used to set the initial point of expansion
//----------------------------------------------------------------------------
CMatrix InitialStateVector(string sKalmanModel_, DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Returns the initial state vector for Kalman processing
// sKalmanModel_: Kalman model
// pdConstants_: Constants used to set the initial state covariance matrix
//----------------------------------------------------------------------------
CMatrix InitialStateCovP(string sKalmanModel_, DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Returns the state matrix for Kalman processing
// sKalmanModel_: Kalman model
// pdConstants_: Constants used to set the state matrix
//----------------------------------------------------------------------------
CMatrix StateMatrix(string sKalmanModel_, DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Returns the state input matrix of noise for Kalman processing
// sKalmanModel_: Kalman model
// pdConstants_: Constants used to set the state matrix
//----------------------------------------------------------------------------
CMatrix InputMatrix(string sKalmanModel_);

//----------------------------------------------------------------------------
// Returns the state spectral noise density for Kalman processing
// sKalmanModel_: Kalman model
// pdConstants_: Constants used to set the state matrix
//----------------------------------------------------------------------------
CMatrix SpectralDensity(string sKalmanModel_, DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Returns the discrete time process noise matrix Q
// pclF_: State matrix
// pclG_: Input noise matrix
// pclQ_: Power spectral density matrix
// dStepSize_: In-between epoch differences
//----------------------------------------------------------------------------
CMatrix ProcessNoiseQ(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_, DOUBLE dStepSize_);

//----------------------------------------------------------------------------
// Returns the discrete Kalman gain matrix
// pclP_: Apriori state covariance matrix
// pclH_: Design geometry matrix
// pclR_: Observation covariance matrix
//----------------------------------------------------------------------------
CMatrix KalmanGain(CMatrix* pclP_, CMatrix* pclH_, CMatrix* pclR_);

//----------------------------------------------------------------------------
// Computes Kalman state vector and state covariance matrix and plots in files
// sKalmanModel_: PRW, VRW or Gauss Markov dynamic model
// pmpFiles_: Map to Input/Output files
// pvEpochData_: Measurement data read from file
// pstTruthData_: Truth data read from file
// sOutDirPath_: Output directory path were results are plotted
// pdConstants_: Constants used to control the Kalman filter processor
// pclStateVector_: Last Kalman state vector that is updated after measurement is processed
// pclStateCovP_: Last Kalman state covariance matrix updated after measurement is processed
// dTimeStep_: Difference between last and current epoch
//----------------------------------------------------------------------------
CMatrix KalmanEpochSolution(string sKalmanModel_,
                            map<string, fstream*>* pmpFiles_,
                            vector<SatelliteData>* pvEpochData_, 
                            TruthData* pstTruthData_, 
                            string sOutDirPath_, 
                            DOUBLE* pdConstants_,
                            CMatrix* pclStateVector_, 
                            CMatrix* pclStateCovP_,
                            DOUBLE& dTimeStep_);

//----------------------------------------------------------------------------
// Computes Kalman state vector and state covariance matrix iterating on single epoch
// pclPredictedState_: From the Kalman filter prediction loop state vector
// pclPredStateCovP_: From Kalman filter prediction loop state covariance matrix
// pclObsCovR_: Current epoch observation covariance matrix
// pclUpdatedState_: Kalman filter corrected state vector
// pclUpdatedCovP_: Kalman filter corrected state covariance matrix
// pvEpochData_: Current epoch measurement data vector
//----------------------------------------------------------------------------
void IEKF(CMatrix* pclPredictedState_,
          CMatrix* pclPredStateCovP_, 
          CMatrix* pclObsCovR_,
          CMatrix* pclUpdatedState_,
          CMatrix* pclUpdatedCovP_,
          vector<SatelliteData>* pvEpochData_);

#endif
