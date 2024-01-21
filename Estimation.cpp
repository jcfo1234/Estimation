// Estimation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
using namespace std;

// Default point of expansion
DOUBLE dTheLatitude = 51.1225;
DOUBLE dTheLongitude = -114.0133;
DOUBLE dTheHeight = 1084;
DOUBLE dTheVelX = 0;
DOUBLE dTheVelY = 0;
DOUBLE dTheVelZ = 0;
DOUBLE dGPSTimeShift = 0;
DOUBLE dGPSClockDrift = 0;

// Apriori variance factor [m^2]
DOUBLE dVarFactor = 1;
// Initial position state standard deviations [m]
DOUBLE dTheVarX = 1000;
DOUBLE dTheVarY = 1000;
DOUBLE dTheVarZ = 1000;
// Initial velocity state standard deviations [m/s]
DOUBLE dTheVarVelX = 10;
DOUBLE dTheVarVelY = 10;
DOUBLE dTheVarVelZ = 10;
// Initial clock bias state standard deviation [m]
DOUBLE dTheVarClkBias = 100;
// Initial clock drift state standard deviation [m/s]
DOUBLE dTheVarClkDrift = 100;
// Position Random Walk spectral densities [m/Hz]
DOUBLE dPRWXSpectralDensity = 10;
DOUBLE dPRWYSpectralDensity = 10;
DOUBLE dPRWZSpectralDensity = 10;
DOUBLE dPRWh0 = 1e-21;
// Velocity Random Walk spectral densities
DOUBLE dVRWXSpectralDensity = 1000;
DOUBLE dVRWYSpectralDensity = 1000;
DOUBLE dVRWZSpectralDensity = 1000;
DOUBLE dPRWh1 = 1e-20;
// Gauss-Markov time constants
DOUBLE dCorrVxTimeConstant = 1e-3;
DOUBLE dCorrVyTimeConstant = 1e-3;
DOUBLE dCorrVzTimeConstant = 1e-3;
DOUBLE dCorrClkDriftTimeConstant = 1e-3;
BOOLEANO bIEKFEnable = FALSE;

// Input files
string sMeasFile("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO629\\Labs\\Lab1\\satpos_meas.txt");
string sPosFile("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO629\\Labs\\Lab1\\refpos.txt");
string sOutPath("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO629\\Labs\\Lab1");
// Data Processing to perform
string sDataProcess("LEASTSQUARES");
string sKalmanModel("POSRANDWALK");

//-----------------------------------------------------------------------------------
// Control switch options
// Control switches are:
// /-latpexp<Latitude in Degrees>             : Change latitude default value
// /-lngpexp<Longitude in Degrees>            : Change longitude default value
// /-hgtpexp<Height in Degrees>               : Change height default value
// /-clkbpexp<Time shift in meters>           : Change time shift default value
// /-clkdpexp<Clock drift in m/s>             : Change the clock drift default value
// /-filerov<Full path to rover file>         : Change rover file path default value
// /-filetruth<Full path to truth file>       : Change truth file path default value
// /-fnctn<LEASTSQUARES or KALMANFILTER>      : Data Processing solution type
// /-outpath<Output files path>               : Directory where results are output
// /-velxpexp<ECEF X Velocity>                : ECEF velocity X component
// /-velypexp<ECEF Y Velocity>                : ECEF velocity Y component
// /-velzpexp<ECEF Z Velocity>                : ECEF velocity Z component
// /-prwqx<ECEF X spectral density>           : Power spectral density of X
// /-prwqy<ECEF Y spectral density>           : Power spectral density of Y
// /-prwqz<ECEF Z spectral density>           : Power spectral density of Z
// /-prwqclkb<Clock h0>                       : Clock oscilator h0 parameter
// /-vrwqx<ECEF Velocity X PSD>               : Power spectral density of velocity X
// /-vrwqy<ECEF Velocity Y PSD>               : Power spectral density of velocity Y
// /-vrwqz<ECEF Velocity Z PSD>               : Power spectral density of velocity Z
// /-vrwqclkd<Clock h_-2>                     : Clock oscillator h_-2 parameter
// /-gmtcx<Time constant>                     : Gauss-Markov velocity X time constant
// /-gmtcy<Time constant>                     : Gauss-Markov velocity Y time constant
// /-gmtcz<Time constant>                     : Gauss-Markov velocity Z time constant
// /-gmtcclkd<Time constant>                  : Gauss-Markov clock drift time constant
// /-varx<standard deviation>                 : Initial X standard deviation
// /-vary<standard deviation>                 : Initial Y standard deviation
// /-varz<standard deviation>                 : Initial Z standard deviation
// /-varvelx<standard deviation>              : Initial velocity X standard deviation
// /-varvely<standard deviation>              : Initial velocity Y standard deviation
// /-varvelz<standard deviation>              : Initial velocity Z standard deviation
// /-varclkb<standard deviation>              : Initial clock bias standard deviation
// /-varclkd<standard deviation>              : Initial clock drift standard deviation
// /-kalmod<string>                           : POSRANDWALK, VELRANDWALK or GAUSSMARKOV
// /-IEKF                                     : Enable iterated extended Kalman filter
//-----------------------------------------------------------------------------------
const CHAR* acOptions[ESTIMATION_NUM_OPTIONS] = { "LATPEXP", "LNGPEXP", "HGTPEXP", "CLKBPEXP",
                                                  "CLKDPEXP", "VARFACTOR", "FILEROV", "FILETRUTH", 
                                                  "FNCTN", "OUTPATH", "VELXPEXP", "VELYPEXP",
                                                  "VELZPEXP", "PRWQX", "PRWQY", "PRWQZ",
                                                  "PRWQCLKB", "VRWQX", "VRWQY", "VRWQZ",
                                                  "VRWQCLKD", "GMTCX", "GMTCY", "GMTCZ",
                                                  "GMTCCLKD", "VARX", "VARY", "VARZ",
                                                  "VARVELX", "VARVELY", "VARVELZ", "VARCLKB", 
                                                  "VARCLKD", "KALMOD", "IEKF"};
vector<string> vOptions(acOptions, acOptions + ESTIMATION_NUM_OPTIONS);
map<string, void*> mpTheSwitchOptions;
// Error codes
map<INT, string> mpTheErrorCodes;

// Function pointer to the appropriate data processing
typedef void(*DataProcessing)(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, DOUBLE* pdConstants_, string sKalmanModel_);

//-----------------------------------------------------------------------------------
INT main(INT iargc_, CHAR** ppcargv_)
{
   // Useful program constants
   DOUBLE adConstants[DATATYPEDEFINITIONS_CONSTANTSARRAYSIZE];
   // Vector with full paths to input and output files
   vector<string> vFilePaths;
   // Map with input/output file handlers
   map<string, fstream*> mpFiles;
   // Map with data processing function pointers
   map<string, DataProcessing> mpFunc;
   mpFunc["LEASTSQUARES"] = &LeastSquaresProcess;
   mpFunc["KALMANFILTER"] = &KalmanProcess;
   // Initialize control switch options
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POINTEXPANSIONLATITUDE_OPTIONINDEX)] = &dTheLatitude;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POINTEXPANSIONLONGITUDE_OPTIONINDEX)] = &dTheLongitude;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POINTEXPANSIONHEIGHT_OPTIONINDEX)] = &dTheHeight;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POINTEXPANSIONTIMESHIFT_OPTIONINDEX)] = &dGPSTimeShift;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POINTEXPANSIONCLOCKDRIFT_OPTIONINDEX)] = &dGPSClockDrift;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_APRIORIVARIANCEFACTOR_OPTIONINDEX)] = &dVarFactor;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_MEASUREMENTFILE_OPTIONINDEX)] = &sMeasFile;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_TRUTHFILE_OPTIONINDEX)] = &sPosFile;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_FUNCTIONTYPE_OPTIONINDEX)] = &sDataProcess;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_OUTPATHDIRECTORY_OPTIONINDEX)] = &sOutPath;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VELXECEFPOINTEXPANSION_INDEX)] = &dTheVelX;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VELYECEFPOINTEXPANSION_INDEX)] = &dTheVelY;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VELZECEFPOINTEXPANSION_INDEX)] = &dTheVelZ;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_PRWXSPECTRALDENSITY_INDEX)] = &dPRWXSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_PRWYSPECTRALDENSITY_INDEX)] = &dPRWYSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_PRWZSPECTRALDENSITY_INDEX)] = &dPRWZSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_PRWCLKBIASSPECTRALDENSITY_INDEX)] = &dPRWh0;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VRWXSPECTRALDENSITY_INDEX)] = &dVRWXSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VRWYSPECTRALDENSITY_INDEX)] = &dVRWYSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VRWZSPECTRALDENSITY_INDEX)] = &dVRWZSpectralDensity;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VRWCLKDRIFTSPECTRALDENSITY_INDEX)] = &dPRWh1;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CORRVXTIMECONSTANT_INDEX)] = &dCorrVxTimeConstant;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CORRVYTIMECONSTANT_INDEX)] = &dCorrVyTimeConstant;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CORRVZTIMECONSTANT_INDEX)] = &dCorrVzTimeConstant;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CORRCLKDRIFTTIMECONSTANT_INDEX)] = &dCorrClkDriftTimeConstant;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARX_INDEX)] = &dTheVarX;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARY_INDEX)] = &dTheVarY;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARZ_INDEX)] = &dTheVarZ;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARVELX_INDEX)] = &dTheVarVelX;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARVELY_INDEX)] = &dTheVarVelY;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARVELZ_INDEX)] = &dTheVarVelZ;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARCLKBIAS_INDEX)] = &dTheVarClkBias;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVARCLKDRIFT_INDEX)] = &dTheVarClkDrift;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_KALMANSTATEMODEL_INDEX)] = &sKalmanModel;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_IEKFENABLE_INDEX)] = &bIEKFEnable;
   // Error codes
   mpTheErrorCodes[1] = "Option used is invalid (main)";
   mpTheErrorCodes[2] = "Switch option used is invalid (main)";
   mpTheErrorCodes[3] = "Number of elements read from input file is greater than the expected (ParseReadLine)";
   mpTheErrorCodes[4] = "String argument could not be transformed to a number (ParseReadLine)";
   mpTheErrorCodes[5] = "Matrix provided must be one dimensional (CalcN)";
   mpTheErrorCodes[6] = "Matrix provided must be one dimensional (CalcM)";
   mpTheErrorCodes[7] = "Vector provided must be dimension 4";
   mpTheErrorCodes[8] = "Vector provided must be dimension 4";
   mpTheErrorCodes[9] = "Unrecognized Kalman model configuration (InitialStateVector)";
   mpTheErrorCodes[10] = "Unrecognized Kalman model configuration (InitialStateCovP)";

   try
   {
      // Parse program command arguments
      for (INT i = 1; i < iargc_; i++)
      {
         // Transform one string to uppercase and keep original string
         string stOption(ppcargv_[i]);
         string stOptionOrig(ppcargv_[i]);
         BOOLEANO bOptionFound = FALSE;
         transform(stOption.begin(), stOption.end(), stOption.begin(), toupper);
         switch (stOption.at(0))
         {
            // Allowed switch delimiters
         case '/':
         case '-':
            for (vector<string>::iterator itOptions = vOptions.begin(); itOptions != vOptions.end(); itOptions++)
            {
               // Option found
               if (stOption.find(*itOptions) != string::npos)
               {
                  // File path assignment changes from default
                  if (*itOptions == vOptions.at(ESTIMATION_MEASUREMENTFILE_OPTIONINDEX) ||
                      *itOptions == vOptions.at(ESTIMATION_TRUTHFILE_OPTIONINDEX) ||
                      *itOptions == vOptions.at(ESTIMATION_OUTPATHDIRECTORY_OPTIONINDEX))
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     *static_cast<string*>(mpTheSwitchOptions.at(*itOptions)) = stChangeValue;
                  }
                  // Data Processing function call
                  else if (*itOptions == vOptions.at(ESTIMATION_FUNCTIONTYPE_OPTIONINDEX) ||
                           *itOptions == vOptions.at(ESTIMATION_KALMANSTATEMODEL_INDEX))
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     transform(stChangeValue.begin(), stChangeValue.end(), stChangeValue.begin(), toupper);
                     *static_cast<string*>(mpTheSwitchOptions.at(*itOptions)) = stChangeValue;
                  }
                  // Enable iterated extended kalman filter
                  else if (*itOptions == vOptions.at(ESTIMATION_IEKFENABLE_INDEX))
                  {
                     *((BOOLEANO*)mpTheSwitchOptions.at(*itOptions)) = TRUE;
                  }
                  // Latitude, Longitude, Height, Time shift or Measurement standard deviation assignment
                  else
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     *((DOUBLE*)mpTheSwitchOptions.at(*itOptions)) = stod(stChangeValue);
                  }
                  bOptionFound = TRUE;
                  break;
               }
            }
            // Option not found
            if (!bOptionFound)
            {
               throw 1;
            }
            break;
            // Invalid switch delimiter
         default:
            throw 2;
         }
      }
      // Store file paths in vector
      vFilePaths.push_back(sMeasFile);
      vFilePaths.push_back(sPosFile);
      vFilePaths.push_back(sOutPath);
      // Assign program constants
      AssignConstants(&adConstants[0]);
      // Execute data processing program
      mpFunc.at(sDataProcess)(&vFilePaths, &mpFiles, &adConstants[0], sKalmanModel);
      // De-allocate heap memory
      DeallocateHeap(&mpFiles);
   }
   // Catch errors discovered in the program
   catch (INT iError)
   {
      DeallocateHeap(&mpFiles);
      cout << mpTheErrorCodes.at(iError) << endl;
      return EXIT_FAILURE;
   }
   // Catch unknown errors discovered in sub-functions and dependencies
   catch (...)
   {
      DeallocateHeap(&mpFiles);
      cout << "Unknown Error" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------------
void AssignConstants(DOUBLE* pdConstants_)
{
   // Position vector in curvilinear geoid coordinates
   pdConstants_[DATATYPEDEFINITIONS_LATITUDEPOINTEXPANSION_INDEX] = dTheLatitude;
   pdConstants_[DATATYPEDEFINITIONS_LONGITUDEPOINTEXPANSION_INDEX] = dTheLongitude;
   pdConstants_[DATATYPEDEFINITIONS_HEIGHTPOINTEXPANSION_INDEX] = dTheHeight;
   pdConstants_[DATATYPEDEFINITIONS_TIMESHIFTPOINTEXPANSION_INDEX] = dGPSTimeShift;
   pdConstants_[DATATYPEDEFINITIONS_CLOCKDRIFTPOINTEXPANSION_INDEX] = dGPSClockDrift;
   // Measurements variance factor
   pdConstants_[DATATYPEDEFINITIONS_APRIORIVARIANCE_INDEX] = dVarFactor;
   // Velocity vector ECEF point of expansion
   pdConstants_[DATATYPEDEFINITIONS_VELXECEFPOINTEXPANSION_INDEX] = dTheVelX;
   pdConstants_[DATATYPEDEFINITIONS_VELYECEFPOINTEXPANSION_INDEX] = dTheVelY;
   pdConstants_[DATATYPEDEFINITIONS_VELZECEFPOINTEXPANSION_INDEX] = dTheVelZ;
   // Position random walk spectral density
   pdConstants_[DATATYPEDEFINITIONS_PRWXSPECTRALDENSITY_INDEX] = dPRWXSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_PRWYSPECTRALDENSITY_INDEX] = dPRWYSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_PRWZSPECTRALDENSITY_INDEX] = dPRWZSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_PRWCLKBIASSPECTRALDENSITY_INDEX] = pow(SPEED_OF_LIGHT, 2) * dPRWh0 / 2;
   // Velocity random walk spectral density
   pdConstants_[DATATYPEDEFINITIONS_VRWXSPECTRALDENSITY_INDEX] = dVRWXSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_VRWYSPECTRALDENSITY_INDEX] = dVRWYSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_VRWZSPECTRALDENSITY_INDEX] = dVRWZSpectralDensity;
   pdConstants_[DATATYPEDEFINITIONS_VRWCLKDRIFTSPECTRALDENSITY_INDEX] = pow(SPEED_OF_LIGHT, 2) * dPRWh1 * 2 * pow(DATATYPEDEFINITIONS_PI, 2);
   // Autocorrelation time constants
   pdConstants_[DATATYPEDEFINITIONS_CORRVXTIMECONSTANT_INDEX] = dCorrVxTimeConstant;
   pdConstants_[DATATYPEDEFINITIONS_CORRVYTIMECONSTANT_INDEX] = dCorrVyTimeConstant;
   pdConstants_[DATATYPEDEFINITIONS_CORRVZTIMECONSTANT_INDEX] = dCorrVzTimeConstant;
   pdConstants_[DATATYPEDEFINITIONS_CORRCLKDRIFTTIMECONSTANT_INDEX] = dCorrClkDriftTimeConstant;
   // Initial position state standard deviations
   pdConstants_[DATATYPEDEFINITIONS_INITVARX_INDEX] = dTheVarX;
   pdConstants_[DATATYPEDEFINITIONS_INITVARY_INDEX] = dTheVarY;
   pdConstants_[DATATYPEDEFINITIONS_INITVARZ_INDEX] = dTheVarZ;
   // Initial velocity state standard deviations
   pdConstants_[DATATYPEDEFINITIONS_INITVARVELX_INDEX] = dTheVarVelX;
   pdConstants_[DATATYPEDEFINITIONS_INITVARVELY_INDEX] = dTheVarVelY;
   pdConstants_[DATATYPEDEFINITIONS_INITVARVELZ_INDEX] = dTheVarVelZ;
   // Initial clock bias standard deviation
   pdConstants_[DATATYPEDEFINITIONS_INITVARCLKBIAS_INDEX] = dTheVarClkBias;
   // Initial clock drift standard deviation
   pdConstants_[DATATYPEDEFINITIONS_INITVARCLKDRIFT_INDEX] = dTheVarClkDrift;
   // Enable/Disable iterated extended kalman filter
   (bIEKFEnable == TRUE) ? pdConstants_[DATATYPEDEFINITIONS_IEKFENABLE_INDEX] = 1.0 : pdConstants_[DATATYPEDEFINITIONS_IEKFENABLE_INDEX] = -1.0;
}

//-----------------------------------------------------------------------------------
void DeallocateHeap(map<string, fstream*>* pmpFiles_)
{
   // Iterate through all the opened file object array
   for (map<string, fstream*>::iterator pIterFiles = pmpFiles_->begin(); pIterFiles != pmpFiles_->end(); pIterFiles++)
   {
      // Delete any allocated heap memory
      if (pIterFiles->second)
      {
         pIterFiles->second->close();
         delete pIterFiles->second;
      }
   }
   pmpFiles_->clear();
}