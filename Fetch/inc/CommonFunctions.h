//----------------------------------------------------------------------
// CommonFunctions.h
// Module that handles common vector operation functions
//----------------------------------------------------------------------

#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

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
#include "..\..\Matrix\inc\Matrix.h"
using namespace std;

#define COORDTRANS_WGS84_SEMIMAJORAXIS ((DOUBLE) 6378137.0)
#define COORDTRANS_WGS84_SEMIMINORAXIS ((DOUBLE) 6356752.31425)
#define COORDTRANS_PI ((DOUBLE) 3.14159265358979323846)

//----------------------------------------------------------------------------
// Calculate N (Prime vertical radius) at a give latitude for a given ellipsoid
// pvLat_: Pointer to vector with latitude values
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
CMatrix CalcN(CMatrix* pvLat_, DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS);

//----------------------------------------------------------------------------
// Calculate dN/dlatitude at a give latitude for a given ellipsoid
// dLat_: Latitude in radians
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
DOUBLE CalcdNdLatitude(DOUBLE dLatitude_, DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS);

//----------------------------------------------------------------------------
// Calculate Eccentricity of Ellipsoid
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
DOUBLE CalcEccentricity(DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS);

//----------------------------------------------------------------------------
// Calculate M (Meridian radius of curvature) at a give latitude 
// for a given ellipsoid
// pvLat_: Pointer to vector with latitude values
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
CMatrix CalcM(CMatrix* pvLat_, DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS);

//----------------------------------------------------------------------------
// Transforms Latitude, Longitude and Height to ECEF coordinates
// for a given ellipsoid
// pvGeoidCoord_: Pointer to vector with Latitude (rad), Longitude (rad) and Height (m)
// pvECEFCoord_: Pointer to vector with ECEF coordinates X (m), Y(m), Z (m)
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
void GeoidToECEF(CMatrix* pvGeoidCoord_, CMatrix* pvECEFCoord_, DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS);

//----------------------------------------------------------------------------
// Transforms ECEF coordinates to Latitude, Longitude and Height
// for a given ellipsoid
// pvECEFCoord_: Pointer to vector with ECEF coordinates X (m), Y(m), Z (m)
// pvGeoidCoord_: Pointer to vector with Latitude (rad), Longitude (rad) and Height (m)
// dSemiMajorAxis: Ellipsoid semi-major axis
// dSemiMinorAxis: Ellipsoid semi-minor axis
//----------------------------------------------------------------------------
void ECEFToGeoid(CMatrix* pvECEFCoord_, CMatrix* pvGeoidCoord_, DOUBLE dSemiMajorAxis_ = COORDTRANS_WGS84_SEMIMAJORAXIS, DOUBLE dSemiMinorAxis_ = COORDTRANS_WGS84_SEMIMINORAXIS, DOUBLE CONVERGENCE_ = 0.000001);

//----------------------------------------------------------------------------
// Returns rotation matrix that rotates vectors from ECEF to Local Level Frame
// pvGeoidCoord_: Pointer to vector with Latitude (rad), Longitude (rad) and Height (m)
//----------------------------------------------------------------------------
CMatrix ECEFToLLF(CMatrix* pvGeoidCoord_);

//----------------------------------------------------------------------------
// Returns rotation matrix that rotates vectors from ECEF to Local Level Frame
// pvGeoidCoord_: Pointer to vector with Latitude (rad), Longitude (rad) and Height (m)
//----------------------------------------------------------------------------
CMatrix ECEFToLLF2(CMatrix* pvGeoidCoord_);

//----------------------------------------------------------------------------
// Transforms angular quantity from degrees to radians
// dDeg_: Angle given in degrees
DOUBLE DegToRad(DOUBLE dDeg_);

//----------------------------------------------------------------------------
// Transforms angular quantity from radians to degrees
// dRad_: Angle given in radians
DOUBLE RadToDeg(DOUBLE dRad_);

//----------------------------------------------------------------------------
// Computes misclosure from measurements and estimated geometric ranges
// pclPointExpansion_ (In): ECEF initial coordinates
// pvEpochSVData_ (In): Vector with satellite measurement data
// pvPRNOrder_ (Out): Vector storing the PRN order in the measurements
// Function outputs the misclosure vector of the measurements
CMatrix ComputeMisclosure(CMatrix* pclPointExpansion_, vector<SatelliteData>* pvEpochSVData_, vector<ULONG>* pvPRNOrder_);

//----------------------------------------------------------------------------
// Computes misclosure from measurements and estimated geometric ranges
// pclPointExpansion_ (In): ECEF initial coordinates
// pvEpochSVData_ (In): Vector with satellite measurement data
// Function outputs the desing matrix
CMatrix ComputeDesignMatrix(CMatrix* pclPointExpansion_, vector<SatelliteData>* pvEpochSVData_);

//----------------------------------------------------------------------------
// Computes misclosure from measurements and estimated geometric ranges
// pvEpochSVData_ (In): Vector with satellite measurement data
// Function outputs the desing matrix
CMatrix ComputeObsCovMatrix(vector<SatelliteData>* pvEpochSVData_, DOUBLE* pdConstants_);

#endif
