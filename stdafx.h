// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"
#include "DataTypeDefinitions.h"
#include "Matrix\inc\Matrix.h"
#include "LS\inc\LSCommon.h"
#include "Kalman\inc\KalmanCommon.h"

#include <stdio.h>
#include <tchar.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#define ESTIMATION_NUM_OPTIONS ((ULONG) 35)
// Estimation options
#define ESTIMATION_POINTEXPANSIONLATITUDE_OPTIONINDEX ((ULONG) 0)
#define ESTIMATION_POINTEXPANSIONLONGITUDE_OPTIONINDEX ((ULONG) 1)
#define ESTIMATION_POINTEXPANSIONHEIGHT_OPTIONINDEX ((ULONG) 2)
#define ESTIMATION_POINTEXPANSIONTIMESHIFT_OPTIONINDEX ((ULONG) 3)
#define ESTIMATION_APRIORIVARIANCEFACTOR_OPTIONINDEX ((ULONG) 4)
#define ESTIMATION_POINTEXPANSIONCLOCKDRIFT_OPTIONINDEX ((ULONG) 5)
#define ESTIMATION_MEASUREMENTFILE_OPTIONINDEX ((ULONG) 6)
#define ESTIMATION_TRUTHFILE_OPTIONINDEX ((ULONG) 7)
#define ESTIMATION_FUNCTIONTYPE_OPTIONINDEX ((ULONG) 8)
#define ESTIMATION_OUTPATHDIRECTORY_OPTIONINDEX ((ULONG) 9)
#define ESTIMATION_VELXECEFPOINTEXPANSION_INDEX ((ULONG) 10)
#define ESTIMATION_VELYECEFPOINTEXPANSION_INDEX ((ULONG) 11)
#define ESTIMATION_VELZECEFPOINTEXPANSION_INDEX ((ULONG) 12)
#define ESTIMATION_PRWXSPECTRALDENSITY_INDEX ((ULONG) 13)
#define ESTIMATION_PRWYSPECTRALDENSITY_INDEX ((ULONG) 14)
#define ESTIMATION_PRWZSPECTRALDENSITY_INDEX ((ULONG) 15)
#define ESTIMATION_PRWCLKBIASSPECTRALDENSITY_INDEX ((ULONG) 16)
#define ESTIMATION_VRWXSPECTRALDENSITY_INDEX ((ULONG) 17)
#define ESTIMATION_VRWYSPECTRALDENSITY_INDEX ((ULONG) 18)
#define ESTIMATION_VRWZSPECTRALDENSITY_INDEX ((ULONG) 19)
#define ESTIMATION_VRWCLKDRIFTSPECTRALDENSITY_INDEX ((ULONG) 20)
#define ESTIMATION_CORRVXTIMECONSTANT_INDEX ((ULONG) 21)
#define ESTIMATION_CORRVYTIMECONSTANT_INDEX ((ULONG) 22)
#define ESTIMATION_CORRVZTIMECONSTANT_INDEX ((ULONG) 23)
#define ESTIMATION_CORRCLKDRIFTTIMECONSTANT_INDEX ((ULONG) 24)
#define ESTIMATION_INITVARX_INDEX ((ULONG) 25)
#define ESTIMATION_INITVARY_INDEX ((ULONG) 26)
#define ESTIMATION_INITVARZ_INDEX ((ULONG) 27)
#define ESTIMATION_INITVARVELX_INDEX ((ULONG) 28)
#define ESTIMATION_INITVARVELY_INDEX ((ULONG) 29)
#define ESTIMATION_INITVARVELZ_INDEX ((ULONG) 30)
#define ESTIMATION_INITVARCLKBIAS_INDEX ((ULONG) 31)
#define ESTIMATION_INITVARCLKDRIFT_INDEX ((ULONG) 32)
#define ESTIMATION_KALMANSTATEMODEL_INDEX ((ULONG) 33)
#define ESTIMATION_IEKFENABLE_INDEX ((ULONG) 34)

void AssignConstants(DOUBLE* pdConstants_);

void DeallocateHeap(map<string, fstream*>* pmpFiles_);
// TODO: reference additional headers your program requires here
