//----------------------------------------------------------------------
// LSCommon.h
// Module that handles common least squares functions
//----------------------------------------------------------------------

#ifndef LSCOMMON_H
#define LSCOMMON_H

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

#define LSCOMMON_ERROR_TOLERANCE ((DOUBLE) 1e-6)
#define LSCOMMON_ERRORDECIMAL_SIGNIFICANCE ((ULONG) 9)
#define LSCOMMON_STDEVDECIMAL_SIGNIFICANCE ((ULONG) 6)

//----------------------------------------------------------------------------
// Reads measurements file and truth file to plot least squares position errors
// pvFilePaths_: Pointer to vector with input file paths
// pmpFiles_: Map from file paths to pointer to file for handling input/output
// pdConstatns_: Useful values for iterating on the Least squares algorithm
// sKalmanModel_: Unused in least squares
//----------------------------------------------------------------------------
void LeastSquaresProcess(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, DOUBLE* pdConstants_, string sKalmanModel_);

//----------------------------------------------------------------------------
// Least squares single point solution
// pmpFiles_: Map from file paths to pointer to file for handling input/output
// pclCurrPos_: Current point of expansion, (Latitude, Longitude in degrees, heigth n meters)
// pvEpochData_: Measurements and satellite data vector
// pstTruthData_: Truth data structure
// sOutDirPath_: Output directory
// pdConstants_: Constants used by the least squares processing algorithm
//----------------------------------------------------------------------------
CMatrix SPEpochSolution(map<string, fstream*>* pmpFiles_,
                        CMatrix* pclCurrPos_,
                        vector<SatelliteData>* pvEpochData_, 
                        TruthData* pstTruthData_, 
                        string sOutDirPath_,
                        DOUBLE* pdConstants_);

//----------------------------------------------------------------------------
// Computes the state covariance matrix
// pclObsCovMatrix_: Observation covariance matrix
// pclDesignMatrixA_: Input design matrix
//----------------------------------------------------------------------------
CMatrix ComputeStateCovariance(CMatrix* pclObsCovMatrix_, CMatrix* pclDesignMatrixA_);

//----------------------------------------------------------------------------
// Computes the state correction vector
// pclObsCovMatrix_: Observation covariance matrix
// pclDesignMatrixA_: Input design matrix
// pclMisclusure_: Misclosure vector
//----------------------------------------------------------------------------
CMatrix ComputeStateCorrection(CMatrix* pclObsCovMatrix_, CMatrix* pclDesignMatrixA_, CMatrix* pclMisclosure_);

#endif