//----------------------------------------------------------------------
// FetchCommon.cpp
// Source file of FetchCommon.h
// Source code for reading input files
//----------------------------------------------------------------------
#include "..\inc\FetchCommon.h"

//------------------------------------------------------------------------------------------------
BOOLEANO bOpenFiles(vector<string>* pvFilePaths_, map<string, fstream*>* pmpFiles_, ULONG ulFlag_, string sOutDir_)
{
   BOOLEANO bFileOpened = TRUE;
   DWORD ftyp = GetFileAttributes(sOutDir_.c_str());
   // Create the output directory if it does not exist
   if (ulFlag_ == ios::out)
   {
      if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES)
      {
         CreateDirectory(sOutDir_.c_str(), NULL);
      }
   }

   // Full file path names are stored in pvFilePaths_
   // Map will store the full file path name with the file handler for further reading
   for (vector<string>::iterator pIterFilePaths = pvFilePaths_->begin(); pIterFilePaths != pvFilePaths_->end(); pIterFilePaths++)
   {
      // If the file not already opened
      if (pmpFiles_->find(*pIterFilePaths) == pmpFiles_->end())
      {
         pmpFiles_->operator[](*pIterFilePaths) = new fstream;
         pmpFiles_->at(*pIterFilePaths)->open(*pIterFilePaths, ulFlag_);
         bFileOpened = bFileOpened & TRUE;
      }
      // File is already opened
      else
      {
         bFileOpened = bFileOpened & FALSE;
      }
   }

   return bFileOpened;
}

//------------------------------------------------------------------------------------------------
void ParseReadLine(DOUBLE* pdValues_, ULONG ulNumStructMembers_, string sData_)
{
   UINT uiValueIndex = 0;
   istringstream ss(sData_);

   while (getline(ss, sData_, ' '))
   {
      if (!sData_.empty())
      {
         try
         {
            (uiValueIndex < ulNumStructMembers_) ? pdValues_[uiValueIndex] = stod(sData_.data()) : throw 3;
         }
         catch (std::invalid_argument ex)
         {
            throw 4;
         }
         uiValueIndex++;
      }
   }
}

//------------------------------------------------------------------------------------------------
BOOLEANO bReadMeasEpoch(map<string, fstream*>* pmpFiles_, 
                        vector<string>* pvFileKeys_,
                        vector<SatelliteData>* pvEpochSVData_, 
                        DOUBLE& dCurrGPSTime_, 
                        DOUBLE& dNextGPSTime_)
{
   // Input file indexes
   const ULONG ulMEASUREMENT_INPUTFILE_INDEX = 0;
   // Next GPS Time
   DOUBLE dNextGPSTime = 0;
   BOOLEANO bReadNextEntrie = TRUE;
   BOOLEANO bEndOfFile = FALSE;
   // Read text file line by line
   string sMeasLine;
   // File pointer position
   streamoff ulMeasOffset = 0;

   //# Read one epoch at a time
   while (bReadNextEntrie && !pmpFiles_->at(pvFileKeys_->at(ulMEASUREMENT_INPUTFILE_INDEX))->eof())
   {
      SatelliteData stSatelliteData;
      ulMeasOffset = pmpFiles_->at(pvFileKeys_->at(ulMEASUREMENT_INPUTFILE_INDEX))->tellg();
      // Read measurement file line by line
      if (getline(*pmpFiles_->at(pvFileKeys_->at(ulMEASUREMENT_INPUTFILE_INDEX)), sMeasLine, '\n'))
      {
         DOUBLE adValues[FETCH_SATELLITEDATA_NUMMEMBERS];
         ParseReadLine(&adValues[0], FETCH_SATELLITEDATA_NUMMEMBERS, sMeasLine);
         // Fill satellite data structure
         dNextGPSTime = adValues[FETCH_SVGNSSTIME_INDEX];
         stSatelliteData.dGNSSTime = adValues[FETCH_SVGNSSTIME_INDEX];
         stSatelliteData.ulPRN = (ULONG)adValues[FETCH_SVPRN_INDEX];
         stSatelliteData.dPRNX = adValues[FETCH_SVECEFX_INDEX];
         stSatelliteData.dPRNY = adValues[FETCH_SVECEFY_INDEX];
         stSatelliteData.dPRNZ = adValues[FETCH_SVECEFZ_INDEX];
         stSatelliteData.dPSR = adValues[FETCH_SVPSR_INDEX];
         stSatelliteData.dADR = adValues[FETCH_SVADR_INDEX];
         stSatelliteData.dElevation = adValues[FETCH_SVELEVATION_INDEX];
         // Transition to next epoch boundary
         if (abs(dNextGPSTime - dCurrGPSTime_) > FETCH_ZERO_TOLERANCE && dCurrGPSTime_ > FETCH_ZERO_TOLERANCE && dNextGPSTime > dCurrGPSTime_)
         {
            bReadNextEntrie = FALSE;
            dNextGPSTime_ = dNextGPSTime;
         }
         else
         {
            dCurrGPSTime_ = dNextGPSTime;
            pvEpochSVData_->push_back(stSatelliteData);
         }
      }
      // No more data left to read
      else
      {
         bEndOfFile = TRUE;
         break;
      }
   }
   // Seek for old epoch so that program starts reading right at the new epoch character
   if (!bEndOfFile)
   {
      pmpFiles_->at(pvFileKeys_->at(ulMEASUREMENT_INPUTFILE_INDEX))->seekg(ulMeasOffset, ios::beg);
   }

   return !bEndOfFile;
}

//------------------------------------------------------------------------------------------------
BOOLEANO bReadTruthEpoch(map<string, fstream*>* pmpFiles_, 
                         vector<string>* pvFileKeys_, 
                         TruthData* pstTruthData, 
                         DOUBLE& dCurrGPSTime_, 
                         DOUBLE dMeasGPSTime_)
{
   // Input file indexes
   const ULONG ulTRUTH_INPUTFILE_INDEX = 1;
   // Next GPS Time
   DOUBLE dNextGPSTime = 0;
   BOOLEANO bReadNextEntrie = TRUE;
   BOOLEANO bEndOfFile = FALSE;
   // Read text file line by line
   string sTruthLine;
   // File pointer position
   streamoff ulTruthOffset = 0;

   //# Read one epoch at a time
   while (bReadNextEntrie && !pmpFiles_->at(pvFileKeys_->at(ulTRUTH_INPUTFILE_INDEX))->eof())
   {
      ulTruthOffset = pmpFiles_->at(pvFileKeys_->at(ulTRUTH_INPUTFILE_INDEX))->tellg();
      // Read measurement file line by line
      if (getline(*pmpFiles_->at(pvFileKeys_->at(ulTRUTH_INPUTFILE_INDEX)), sTruthLine, '\n'))
      {
         DOUBLE adValues[FETCH_TRUTHPOSITION_NUMMEMBERS];
         ParseReadLine(&adValues[0], FETCH_TRUTHPOSITION_NUMMEMBERS, sTruthLine);
         // Fill satellite data structure
         dNextGPSTime = adValues[FETCH_SVGNSSTIME_INDEX];
         // Measurement and truth file data GNSS time must align
         if (abs(dNextGPSTime - dMeasGPSTime_) >= 0 && abs(dNextGPSTime - dMeasGPSTime_) < FETCH_ZERO_TOLERANCE)
         {
            pstTruthData->dGNSSTime = adValues[FETCH_TRUTHGNSSTIME_INDEX];
            pstTruthData->dXCoord = adValues[FETCH_TRUTHECEFX_INDEX];
            pstTruthData->dYCoord = adValues[FETCH_TRUTHECEFY_INDEX];
            pstTruthData->dZCoord = adValues[FETCH_TRUTHECEFZ_INDEX];
            dCurrGPSTime_ = dNextGPSTime;
         }
         // Transition to next epoch boundary
         else if (abs(dNextGPSTime - dCurrGPSTime_) > FETCH_ZERO_TOLERANCE && dCurrGPSTime_ > FETCH_ZERO_TOLERANCE)
         {
            bReadNextEntrie = FALSE;
         }
      }
      // No more data left to read
      else
      {
         bEndOfFile = TRUE;
         break;
      }
   }
   // Seek for old epoch so that program starts reading right at the new epoch character
   if (!bEndOfFile)
   {
      pmpFiles_->at(pvFileKeys_->at(ulTRUTH_INPUTFILE_INDEX))->seekg(ulTruthOffset, ios::beg);
   }

   return !bEndOfFile;
}

//------------------------------------------------------------------------------------------------
void PrintLSErrorHeader(fstream* pfOutFile_)
{
   // Write header of file
   *pfOutFile_ << "[header]" << endl;
   *pfOutFile_ << "title=East, North, Up Position Errors" << endl;
   *pfOutFile_ << "xlabel=Time(s)" << endl;
   *pfOutFile_ << "ylabel=Position Errors Stats (m)" << endl;
   *pfOutFile_ << "grid=both" << endl;
   *pfOutFile_ << "pointsize=5" << endl;
   *pfOutFile_ << "xcolumn=0" << endl;
   *pfOutFile_ << "series=9" << endl;
   *pfOutFile_ << "miny=-10" << endl;
   *pfOutFile_ << "maxy=10" << endl;
   *pfOutFile_ << "\n[series1]" << endl;
   *pfOutFile_ << "name=East-Error" << endl;
   *pfOutFile_ << "column=1" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series2]" << endl;
   *pfOutFile_ << "name=North-Error" << endl;
   *pfOutFile_ << "column=2" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series3]" << endl;
   *pfOutFile_ << "name=Up-Error" << endl;
   *pfOutFile_ << "column=3" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series4]" << endl;
   *pfOutFile_ << "name=East-StdDev+" << endl;
   *pfOutFile_ << "column=4" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series5]" << endl;
   *pfOutFile_ << "name=East-StdDev-" << endl;
   *pfOutFile_ << "column=5" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series6]" << endl;
   *pfOutFile_ << "name=North-StdDev+" << endl;
   *pfOutFile_ << "column=6" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series7]" << endl;
   *pfOutFile_ << "name=North-StdDev-" << endl;
   *pfOutFile_ << "column=7" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series8]" << endl;
   *pfOutFile_ << "name=Up-StdDev+" << endl;
   *pfOutFile_ << "column=8" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series9]" << endl;
   *pfOutFile_ << "name=Up-StdDev-" << endl;
   *pfOutFile_ << "column=9" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[data]" << endl;
}

//------------------------------------------------------------------------------------------------
void PrintLSPlots(fstream* pfOutFile_, DOUBLE* pdValues_, ULONG* pulPrecision_, ULONG ulSeriesSize_)
{
   for (ULONG ulIndex = 0; ulIndex < ulSeriesSize_; ulIndex++)
   {
      if (ulIndex < ulSeriesSize_ - 1)
         *pfOutFile_ << fixed << setprecision(pulPrecision_[ulIndex]) << pdValues_[ulIndex] << ",";
      else
         *pfOutFile_ << fixed << setprecision(pulPrecision_[ulIndex]) << pdValues_[ulIndex] << "," << endl;
   }
}

//------------------------------------------------------------------------------------------------
void PrintLSPositionHeader(fstream* pfOutFile_)
{
   // Write header of file
   *pfOutFile_ << "[header]" << endl;
   *pfOutFile_ << "title=Position (Latitude, Longitude, Height)" << endl;
   *pfOutFile_ << "xlabel=Time(s)" << endl;
   *pfOutFile_ << "ylabel=Position [Deg or m]" << endl;
   *pfOutFile_ << "grid=both" << endl;
   *pfOutFile_ << "pointsize=5" << endl;
   *pfOutFile_ << "xcolumn=0" << endl;
   *pfOutFile_ << "series=6" << endl;
   *pfOutFile_ << "\n[series1]" << endl;
   *pfOutFile_ << "name=Latitude" << endl;
   *pfOutFile_ << "column=1" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series2]" << endl;
   *pfOutFile_ << "name=Longitude" << endl;
   *pfOutFile_ << "column=2" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series3]" << endl;
   *pfOutFile_ << "name=Height" << endl;
   *pfOutFile_ << "column=3" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series4]" << endl;
   *pfOutFile_ << "name=Truth-Latitude" << endl;
   *pfOutFile_ << "column=4" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series5]" << endl;
   *pfOutFile_ << "name=Truth-Longitude" << endl;
   *pfOutFile_ << "column=5" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[series6]" << endl;
   *pfOutFile_ << "name=Truth-Height" << endl;
   *pfOutFile_ << "column=6" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[data]" << endl;
}

//------------------------------------------------------------------------------------------------
void PrintLSTrajectoryHeader(fstream* pfOutFile_)
{
   // Write header of file
   *pfOutFile_ << "[header]" << endl;
   *pfOutFile_ << "title=Position (2-D Trajectory)" << endl;
   *pfOutFile_ << "xlabel=Time(s)" << endl;
   *pfOutFile_ << "ylabel=Position [Deg or m]" << endl;
   *pfOutFile_ << "grid=both" << endl;
   *pfOutFile_ << "pointsize=5" << endl;
   *pfOutFile_ << "xcolumn=0" << endl;
   *pfOutFile_ << "series=1" << endl;
   *pfOutFile_ << "\n[series1]" << endl;
   *pfOutFile_ << "name=2DTraj" << endl;
   *pfOutFile_ << "column=1" << endl;
   *pfOutFile_ << "point=dot" << endl;
   *pfOutFile_ << "connect=1" << endl;
   *pfOutFile_ << "connectdist=-1" << endl;
   *pfOutFile_ << "linewidth=1" << endl;
   *pfOutFile_ << "\n[data]" << endl;
}
