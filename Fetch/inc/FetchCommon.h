//----------------------------------------------------------------------
// FetchCommon.h
// Module that handles common file functions
//----------------------------------------------------------------------

#ifndef FETCHCOMMON_H
#define FETCHCOMMON_H

#include <stdio.h>
#include <tchar.h>
#include <Windows.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
using namespace std;

#define FETCH_ZERO_TOLERANCE ((DOUBLE) 1e-3)

// Satellite and measurement data structure
#define FETCH_SATELLITEDATA_NUMMEMBERS ((ULONG) 8)
#define FETCH_SVGNSSTIME_INDEX ((ULONG) 0)
#define FETCH_SVPRN_INDEX ((ULONG) 1)
#define FETCH_SVECEFX_INDEX ((ULONG) 2)
#define FETCH_SVECEFY_INDEX ((ULONG) 3)
#define FETCH_SVECEFZ_INDEX ((ULONG) 4)
#define FETCH_SVPSR_INDEX ((ULONG) 5)
#define FETCH_SVADR_INDEX ((ULONG) 6)
#define FETCH_SVELEVATION_INDEX ((ULONG) 7)

// Truth file position data structure
#define FETCH_TRUTHPOSITION_NUMMEMBERS ((ULONG) 4)
#define FETCH_TRUTHGNSSTIME_INDEX ((ULONG) 0)
#define FETCH_TRUTHECEFX_INDEX ((ULONG) 1)
#define FETCH_TRUTHECEFY_INDEX ((ULONG) 2)
#define FETCH_TRUTHECEFZ_INDEX ((ULONG) 3)

//--------------------------------------------------------------------------------------
// Function Name:
// bOpenFiles
// Arguments:
// pvFilePaths_ (In): Vector with full file paths to input files
// pmpFiles_ (Out): Map from file paths to file pointer handler
// ulFlag_ (In): Type of file, whether for read or write
// sOutDir_ (In): Path of the output directory to save result files
// Function returns TRUE when both input files are opened successfully
BOOLEANO bOpenFiles(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, ULONG ulFlag_, string sOutDir_="");

//--------------------------------------------------------------------------------------
// Function Name:
// bReadMeasEpoch
// Arguments:
// pmpFiles_ (In): Map from file paths to file pointer handler
// pvFileKeys_ (In): Vector with full file paths to input files
// pvEpochSVData_ (Out): Vector with satellite data measurements per epoch
// dCurrGPSTime_ (In/Out): Marks current GPS Epoch
// dNextGPSTime_ (In/Out): Marks next GPS Epoch
// Function returns TRUE when satellite data per epoch is successfully read
BOOLEANO bReadMeasEpoch(map<string, fstream*>* pmpFiles_, vector<string>* pvFileKeys_, vector<SatelliteData>* pvEpochSVData_, DOUBLE& dCurrGPSTime_, DOUBLE& dNextGPSTime_);

//--------------------------------------------------------------------------------------
// Function Name:
// bReadTruthEpoch
// Arguments:
// pmpFiles_ (In): Map from file paths to file pointer handler
// pvFileKeys_ (In): Vector with full file paths to input files
// pstTruthData (Out): Structure with truth position data per epoch
// dCurrGPSTime_ (In/Out): Marks current GPS Epoch
// dMeasGPSTime_ (In): Variable used to line up truth file with measurement file
// Function returns TRUE when truth data per epoch is successfully read
BOOLEANO bReadTruthEpoch(map<string, fstream*>* pmpFiles_, vector<string>* pvFileKeys_, TruthData* pstTruthData, DOUBLE& dCurrGPSTime_, DOUBLE dMeasGPSTime_);

//--------------------------------------------------------------------------------------
// Function Name:
// ParseReadLine
// Arguments:
// pdValues_ (Out): Stores double values in an array
// ulNumStructMembers_ (In): Number of elements in structure
// sData_ (In): Text file line read from measurement or truth file
void ParseReadLine(DOUBLE* pdValues_, ULONG ulNumStructMembers_, string sData_);

void PrintLSErrorHeader(fstream* pfOutFile_);

void PrintLSPositionHeader(fstream* pfOutFile_);

void PrintLSTrajectoryHeader(fstream* pfOutFile_);

void PrintLSPlots(fstream* pfOutFile_, DOUBLE* pdValues_, ULONG* pulPrecision_, ULONG ulSeriesSize_);

#endif
