/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#include "SetHandler.h"
#include <cmath>

SetHandler::SetHandler(unsigned int maxPSMs) : maxPSMs_(maxPSMs) {}

SetHandler::~SetHandler() {
  reset();
}

/**
 * Frees the feature rows of all PSMs, which are only needed for training.
 */
void SetHandler::releaseFeatures() {
  for (const auto& subset : subsets_) {
    subset->clearFeaturePointers();
  }
  featurePool_.destroyPool();
}

void SetHandler::reset() {
  subsets_.clear();
  psmStore_.clear();
  featurePool_.destroyPool();
  DataSet::resetFeatureNames();
}

std::string SetHandler::describeLine(unsigned int lineNr) const {
  std::ostringstream oss;
  oss << "line ";
  // the last input file whose lines start at or before lineNr
  auto it = std::upper_bound(
      inputSegments_.begin(), inputSegments_.end(), lineNr,
      [](unsigned int line, const InputSegment& segment) {
        return line < segment.firstLine;
      });
  if (it == inputSegments_.begin()) {
    oss << lineNr;
  } else {
    --it;
    oss << it->firstFileLine + (lineNr - it->firstLine) << " of "
        << it->fileName;
  }
  return oss.str();
}

void SetHandler::rethrowWithSource(const MyException& e,
                                   unsigned int lineNr) const {
  if (inputSegments_.empty()) throw;
  std::ostringstream oss;
  oss << e.what() << "Line " << lineNr << " of the concatenated input is "
      << describeLine(lineNr) << "." << std::endl;
  throw MyException(oss.str());
}

void SetHandler::incrementLineNr(unsigned int& lineNr) {
  if (lineNr >= static_cast<unsigned int>(INT_MAX)) {
    std::ostringstream oss;
    oss << "ERROR: The input has more than " << INT_MAX << " lines. "
        << "Percolator counts lines and PSMs with 32-bit integers and "
        << "cannot process more than " << INT_MAX << " PSMs in one run."
        << std::endl;
    throw MyException(oss.str());
  }
  ++lineNr;
}

void SetHandler::warnIgnoredLabel(unsigned int lineNr,
                                  unsigned int& numIgnored) const {
  if (++numIgnored <= kMaxLabelWarnings) {
    std::cerr << "Warning: the PSM on " << describeLine(lineNr)
              << " has a label not in {1,-1} and will be ignored." << std::endl;
  }
}

bool SetHandler::addLabel(unsigned char& seen, bool isDecoy, bool sampling) {
  const unsigned char bit = labelBit(isDecoy);
  if (sampling && (seen & bit) && !(seen & kSeenTwice)) {
    seen |= kSeenTwice;
    ++sampledSpectraWithRepeatedLabel_;
  }
  // a label not seen yet means that the spectrum has both labels
  const bool sameLabels = (seen & bit) != 0;
  seen |= bit;
  return sameLabels;
}

void SetHandler::warnIgnoredLabelTotal(unsigned int numIgnored) const {
  if (numIgnored > kMaxLabelWarnings) {
    std::cerr << "Warning: " << numIgnored << " PSMs have a label not in "
              << "{1,-1} and were ignored; only the first "
              << kMaxLabelWarnings << " are listed above." << std::endl;
  }
}
/**
 * Gets the vector index of the DataSet matching the label
 * @param label DataSet label
 * @return index of matching DataSet
 */
unsigned int SetHandler::getSubsetIndexFromLabel(LabelType label) {
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    if (subsets_[ix]->getLabel() == label) return ix;
  }
  ostringstream temp;
  temp << "Error: No DataSet found with label " << label << std::endl;
  throw MyException(temp.str());
}

/**
 * Insert DataSet object into this SetHandler, which takes ownership
 * @param ds DataSet to be inserted
 */
void SetHandler::push_back_dataset(std::unique_ptr<DataSet> ds) {
  subsets_.push_back(std::move(ds));
}

void SetHandler::populateScoresWithPSMs(vector<ScoreHolder> &scores, LabelType label) {
  subsets_[getSubsetIndexFromLabel(label)]->fillScores(scores);
}

void SetHandler::normalizeFeatures(Normalizer& norm) {
  std::vector<double*> featuresV, rtFeaturesV;
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    subsets_[ix]->fillFeatures(featuresV);
  }
  norm.setSet(featuresV,
      rtFeaturesV,
      FeatureNames::getNumFeatures(),
      0);
  norm.normalizeSet(featuresV, rtFeaturesV);
}

LabelType SetHandler::getLabel(int setPos) {
  assert(setPos >= 0 && setPos < (signed int)subsets_.size());
  return subsets_[static_cast<std::size_t>(setPos)]->getLabel();
}

int SetHandler::readTab(istream& dataStream, SanityCheck& check) {
  std::vector<double> noWeights;
  Scores noScores(true);
  return readAndScoreTab(dataStream, noWeights, noScores, check);
}

int SetHandler::getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields) {
  TabReader reader(headerLine);
  reader.skip(2u); // discard id, label
  bool hasScannr = false;
  while (!reader.error()) {
    std::string optionalHeader = reader.readString();
    // transform to lower case for case insensitive matching
    std::transform(optionalHeader.begin(), optionalHeader.end(), 
                   optionalHeader.begin(), ::tolower);
    if (optionalHeader == "scannr") {
      optionalFields.push_back(SCANNR);
      hasScannr = true;
    } else if (optionalHeader == "expmass") {
      optionalFields.push_back(EXPMASS);
    } else if (optionalHeader == "calcmass") {
      optionalFields.push_back(CALCMASS);
    } else if ((optionalHeader == "rt" || optionalHeader == "retentiontime")) {
      optionalFields.push_back(RETTIME);
    } else if ((optionalHeader == "filename" || optionalHeader == "spectrafile")) {
      optionalFields.push_back(FILENAME);
    } else {
      break;
    }
  }
  if (!hasScannr && VERB > 0) {
    cerr << "\nWARNING: Tab delimited input does not contain ScanNr column," <<
            "\n         scan numbers will be assigned automatically.\n" << endl;
  }
  return static_cast<int>(optionalFields.size());
}

bool SetHandler::isDefaultDirectionLine(const std::string& defaultDirectionLine) {
  TabReader reader(defaultDirectionLine);
  
  std::string psmid = reader.readString();

  std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
  return (psmid == "defaultdirection");
}

int SetHandler::getNumFeatures(const std::string& line, int optionalFieldCount) {
  TabReader reader(line);
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount)); // remove id, label and optional fields
  
  reader.readDouble();
  int numFeatures = 0;
  while (!reader.error()){
    ++numFeatures;
    reader.readDouble();
  }

  return numFeatures;
}

void SetHandler::getFeatureNames(const std::string& headerLine, 
    int numFeatures, int optionalFieldCount, FeatureNames& featureNames) {
  TabReader reader(headerLine);
  // removes enumerator, label and if present optional fields
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount));
  int numFeatLeft = numFeatures;
  while (!reader.error() && numFeatLeft > 0) {
    std::string tmp = reader.readString();
    if (!tmp.empty()) { // Ensure that the string is not empty or consider other validity checks
      featureNames.insertFeature(tmp);
      numFeatLeft--;
    }
  }
  // Check if we have read the correct number of features
  if (numFeatLeft != 0) {
    // thrown rather than exit(), so that the caller's destructors still
    // remove the temporary file of a multi-file input
    ostringstream temp;
    temp << "ERROR: Reading tab file, the header line has fewer columns than "
         << "the first PSM: " << numFeatures - numFeatLeft
         << " column names for its " << numFeatures << " features." << std::endl;
    throw MyException(temp.str());
  }
  featureNames.initFeatures();
  assert(numFeatures == DataSet::getNumFeatures());
}

/**
 * Checks that the Peptide column of the header, if it has one, follows the
 * feature columns of the first PSM line. Otherwise the features would get
 * the names of other columns.
 */
void SetHandler::checkPeptideColumn(const std::string& headerLine,
    int numFeatures, int optionalFieldCount,
    unsigned int firstPsmLineNr) const {
  const int firstFeatureColumn = 2 + optionalFieldCount;
  TabReader reader(headerLine);
  for (int column = 0; !reader.error(); ++column) {
    std::string name = reader.readString();
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (column < firstFeatureColumn || name != "peptide") continue;
    if (column == firstFeatureColumn + numFeatures) return;
    ostringstream temp;
    temp << "ERROR: Reading tab file, the header line has "
         << column - firstFeatureColumn << " feature columns before its "
         << "Peptide column, but the first PSM, on "
         << describeLine(firstPsmLineNr) << ", has " << numFeatures
         << " features." << std::endl;
    Globals::failUnlessNoTerminate(temp.str(),
        "No-terminate flag set: ignoring the error, the feature names may be "
        "wrong.");
    return;
  }
}

bool SetHandler::getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values) {
  TabReader reader(defaultDirectionLine);
  // removes enumerator, label and if present optional fields
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount));
  
  bool hasDefaultValues = false;
  unsigned int ix = 0;
  double a = reader.readDouble();
  while (!reader.error()) {
    if (a != 0.0) hasDefaultValues = true;
    if (VERB > 2) {
      std::cerr << "Initial direction for " << 
                   DataSet::getFeatureNames().getFeatureName(ix) << " is " << 
                   a << std::endl;
    }
    init_values.push_back(a);
    a = reader.readDouble();
    ix++;
  }
  return hasDefaultValues;
}

void SetHandler::writeTab(const string& dataFN, const SanityCheck& check) {
  ofstream dataStream(dataFN.data(), ios::out);
  dataStream << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\t";
  if (PSMDescription::hasSpectrumFileName()) {
    dataStream << "filename\t";
  }
  dataStream << DataSet::getFeatureNames().getFeatureNames()
      << "\tPeptide\tProteins" << std::endl;
  const vector<double>& initial_values = check.getDefaultWeights();
  if (initial_values.size() > 0) {
    dataStream << "DefaultDirection\t-\t-\t-\t-";
    for (size_t i = 0; i < initial_values.size(); ++i) {
      dataStream << '\t' << initial_values[i];
    }
    dataStream << std::endl;
  }
  for (const auto& subset : subsets_) {
    subset->writeTabData(dataStream);
  }
  closeOutputFile(dataStream, dataFN);
}

void SetHandler::readPSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, bool& concatenatedSearch,
    std::vector<OptionalField>& optionalFields) {
  // owned here until reading succeeded
  std::unique_ptr<DataSet> targetSet(new DataSet());
  targetSet->setLabel(LabelType::TARGET);
  std::unique_ptr<DataSet> decoySet(new DataSet());
  decoySet->setLabel(LabelType::DECOY);
  
  unsigned int lineNr = (hasInitialValueRow ? 3u : 2u);
  unsigned int numIgnoredLabels = 0u;
  size_t numSampledLines = 0u;
  sampledLines_ = sampledSpectra_ = sampledSpectraWithRepeatedLabel_ = 0u;
  if (psmLine.size() == 0) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, could not find any PSMs." << std::endl;
    Globals::failUnlessNoTerminate(temp.str(),
        "No-terminate flag set: ignoring error and continuing without PSMs.");
  } else if (maxPSMs_ > 0u) { // reservoir sampling to create subset of size maxPSMs_
    std::priority_queue<PSMDescriptionPriority> subsetPSMs;
    // scan id -> (priority of its first PSM, labels seen as in addLabel);
    // all PSMs of a scan id share its priority
    ScanIdPriorityTable scanIdLookUp;
    unsigned int upperLimit = UINT_MAX;
    try {  // not indented: rethrowWithSource adds the file of the line
    do {
      if (lineNr % 1000000 == 0 && VERB > 1) {
        std::cerr << "Processing line " << lineNr << std::endl;
      }
      psmLine = rtrim(psmLine);
      
      int label = 0;
      ScanId scanId;
      TabReader reader(psmLine);
      PsmHeader header;
      const bool simple = readScanIdAndLabel(psmLine, lineNr, optionalFields,
                                             reader, header, scanId, label);
      bool isDecoy = (label == -1);
      const bool sampling = numSampledLines < kRepeatSampleLines;
      // the generator's values are below PseudoRandom::kRandMax < 2^32
      uint32_t priority = 0u;
      bool seen = false;
      uint8_t& labels = scanIdLookUp.findOrInsert(
          scanId.first, scanId.second, PseudoRandom::lcg_rand, priority, seen);
      if (!seen) {
        labels = labelBit(isDecoy);
      } else if (!addLabel(labels, isDecoy, sampling)) {
        concatenatedSearch = false;
      }
      size_t randIdx = priority;
      if (sampling && ++numSampledLines == kRepeatSampleLines) {
        sampledSpectra_ = scanIdLookUp.size();
      }

      if (label != 1 && label != -1) {
        warnIgnoredLabel(lineNr, numIgnoredLabels);
      } else if (subsetPSMs.size() < maxPSMs_ || randIdx < upperLimit) {
        PSMDescriptionPriority psmPriority;
        bool readProteins = false;
        if (!simple ||
            !DataSet::readPsmBody(reader, header, readProteins, psmPriority.psm,
                                  psmPriority.label, featurePool_, psmStore_,
                                  decoyPrefix_)) {
          psmPriority.label = DataSet::readPsm(psmLine, lineNr, optionalFields,
                                   readProteins, psmPriority.psm, featurePool_,
                                   psmStore_, decoyPrefix_);
        }
        psmPriority.priority = randIdx;
        subsetPSMs.push(psmPriority);
        if (subsetPSMs.size() > maxPSMs_) {
          PSMDescriptionPriority del = subsetPSMs.top();
          upperLimit = static_cast<unsigned int>(del.priority);
          featurePool_.deallocate(del.psm->features);
          psmStore_.release(del.psm);
          subsetPSMs.pop();
        }
      }
      incrementLineNr(lineNr);
    } while (getline(dataStream, psmLine));
    } catch (const MyException& e) {
      rethrowWithSource(e, lineNr);
    }
    if (numSampledLines < kRepeatSampleLines) {
      sampledSpectra_ = scanIdLookUp.size();
    }
    sampledLines_ = numSampledLines;
    
    addQueueToSets(subsetPSMs, targetSet.get(), decoySet.get());
  } else { // simply read all PSMs
    // only needed until a separate search is detected
    SeparateSearchDetector separateSearchDetector;
    // scan id -> labels seen as in addLabel, only for the first
    // kRepeatSampleLines lines; the priorities are not used
    std::unique_ptr<ScanIdPriorityTable> sampledScanIds(
        new ScanIdPriorityTable());
    try {  // not indented: rethrowWithSource adds the file of the line
    do {
      if (lineNr % 1000000 == 0 && VERB > 1) {
        std::cerr << "Reading line " << lineNr << std::endl;
      }
      psmLine = rtrim(psmLine);
      int label = 0;
      ScanId scanId;
      TabReader reader(psmLine);
      PsmHeader header;
      const bool simple = readScanIdAndLabel(psmLine, lineNr, optionalFields,
                                             reader, header, scanId, label);
      bool isDecoy = (label == -1);
      if (concatenatedSearch &&
          separateSearchDetector.add(scanId.first, scanId.second, isDecoy)) {
        concatenatedSearch = false;
      }
      if (numSampledLines < kRepeatSampleLines) {
        uint32_t unusedPriority = 0u;
        bool seen = false;
        uint8_t& labels = sampledScanIds->findOrInsert(
            scanId.first, scanId.second, [] { return 0u; }, unusedPriority,
            seen);
        if (!seen) {
          labels = labelBit(isDecoy);
        } else {
          addLabel(labels, isDecoy, true);
        }
        if (++numSampledLines == kRepeatSampleLines) {
          sampledSpectra_ = sampledScanIds->size();
          sampledScanIds.reset();
        }
      }
      if (label == 1 || label == -1) {
        DataSet* dataSet = label == 1 ? targetSet.get() : decoySet.get();
        if (!simple || !dataSet->readPsmBody(reader, header, featurePool_,
                                             psmStore_, decoyPrefix_)) {
          dataSet->readPsm(psmLine, lineNr, optionalFields, featurePool_,
                           psmStore_, decoyPrefix_);
        }
      } else {
        warnIgnoredLabel(lineNr, numIgnoredLabels);
      }
      incrementLineNr(lineNr);
    } while (getline(dataStream, psmLine));
    } catch (const MyException& e) {
      rethrowWithSource(e, lineNr);
    }
    if (concatenatedSearch &&
        separateSearchDetector.foundSeparateSearches()) {
      concatenatedSearch = false;
    }
    if (numSampledLines < kRepeatSampleLines) {
      sampledSpectra_ = sampledScanIds->size();
    }
    sampledLines_ = numSampledLines;
  }
  warnIgnoredLabelTotal(numIgnoredLabels);

  numPsmLines_ = lineNr - (hasInitialValueRow ? 3u : 2u);
  if (VERB > 1) {
    std::cerr << "Found " << numPsmLines_ << " PSMs" << std::endl;
  }
  
  push_back_dataset(std::move(targetSet));
  push_back_dataset(std::move(decoySet));
}

void SetHandler::addQueueToSets(
    std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet) {
  while (!subsetPSMs.empty()) {
    PSMDescriptionPriority psmPriority = subsetPSMs.top();
    // DataSet::readPsm only returns TARGET or DECOY
    if (psmPriority.label == LabelType::TARGET) {
      targetSet->registerPsm(psmPriority.psm);
    } else {
      decoySet->registerPsm(psmPriority.psm);
    }
    subsetPSMs.pop();
  }
}

/**
 * Reads the label and the scan id of a PSM line, if possible with
 * DataSet::readPsmHeader, which leaves reader and header ready for
 * DataSet::readPsmBody, and otherwise with getScanId. Returns whether
 * readPsmHeader succeeded.
 */
bool SetHandler::readScanIdAndLabel(const std::string& psmLine,
    unsigned int lineNr, std::vector<OptionalField>& optionalFields,
    TabReader& reader, PsmHeader& header, ScanId& scanId, int& label) {
  if (DataSet::readPsmHeader(reader, lineNr, optionalFields, header)) {
    label = header.label;
    scanId = ScanId(static_cast<int>(header.psm.scan), header.psm.expMass);
    return true;
  }
  scanId = getScanId(psmLine, label, optionalFields, lineNr);
  return false;
}

// DataSet::readPsmHeader must fail for every line on which this function
// reports an error or changes a value.
ScanId SetHandler::getScanId(const std::string& psmLine, int& label,
    std::vector<OptionalField>& optionalFields, unsigned int lineNr) {
  ScanId scanId;
  TabReader reader(psmLine);
  
  reader.skip();
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read PSMid." << std::endl;
    throw MyException(temp.str());
  }
  
  label = reader.readInt();
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read label." << std::endl;
    throw MyException(temp.str());
  }
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        scanId.first = reader.readInt();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number on line " 
              << lineNr << ". Check if scan number is an integer between "
              << INT_MIN << " and " << INT_MAX << "." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        scanId.second = reader.readDouble();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading experimental mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        if (!std::isfinite(scanId.second)) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, the experimental mass of PSM "
              << psmLine.substr(0, psmLine.find('\t')) << " on line " << lineNr
              << " is " << scanId.second << " instead of a finite number."
              << std::endl;
          Globals::failUnlessNoTerminate(temp.str(),
              "No-terminate flag set: setting it to 0 and ignoring the error.");
          scanId.second = 0.0;
        }
        break;
      } case CALCMASS: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading calculated mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case RETTIME: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading retention time on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case DELTAMASS: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading delta mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case FILENAME: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading spectra file name on line " 
              << lineNr << "." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) scanId.first = static_cast<int>(lineNr);
  
  return scanId;
}

int SetHandler::readAndScoreTab(istream& dataStream,
    std::vector<double>& rawWeights, Scores& allScores, SanityCheck& check,
    const std::vector<std::vector<double> >& foldWeights,
    const std::vector<std::pair<unsigned int, unsigned int> >& lineFolds) {
  if (!dataStream) {
    std::cerr << "ERROR: Cannot open data stream." << std::endl;
    return 0;
  }
  std::string psmLine, headerLine, defaultDirectionLine;
  
  getline(dataStream, headerLine); // line with feature names
  headerLine = rtrim(headerLine);
  if (headerLine.substr(0,5) == "<?xml") {
    std::cerr << "ERROR: Cannot read Tab delimited input from data stream.\n" << 
       "Input file seems to be in XML format, which is deprecated since percolator 3.9." << 
       std::endl;
    return 0;
  }
  
  // Checking for optional headers "ScanNr", "ExpMass", "CalcMass", "Rt"/"RetentionTime" and "dM"/"DeltaMass"
  std::vector<OptionalField> optionalFields;
  int optionalFieldCount = getOptionalFields(headerLine, optionalFields);
  // parse second line for default direction
  getline(dataStream, defaultDirectionLine);
  defaultDirectionLine = rtrim(defaultDirectionLine);
  bool hasInitialValueRow = isDefaultDirectionLine(defaultDirectionLine);

  // count number of features from first PSM
  if (hasInitialValueRow) {
    getline(dataStream, psmLine);
  } else {
    psmLine = defaultDirectionLine;
  }
  psmLine = rtrim(psmLine);
  int numFeatures = getNumFeatures(psmLine, optionalFieldCount);
  FeatureNames& featureNames = DataSet::getFeatureNames();
  // fill in the feature names from the header line
  getFeatureNames(headerLine, numFeatures, optionalFieldCount, featureNames);
  if (!psmLine.empty()) {  // without PSMs the errors below are more apt
    checkPeptideColumn(headerLine, numFeatures, optionalFieldCount,
                       hasInitialValueRow ? 3u : 2u);
  }
  if (numFeatures < 1) {
    ostringstream oss;
    oss << "ERROR: Reading tab file, too few features present." << std::endl;
    Globals::failUnlessNoTerminate(oss.str(),
        "No-terminate flag set: ignoring error and using a pseudo-feature of "
        "zeroes.");
    featurePool_.createPool(1u);
  } else {
    featurePool_.createPool(DataSet::getNumFeatures());
  }  
  
  // fill in the default weights if present
  std::vector<double> init_values;
  bool hasDefaultValues = false;
  if (hasInitialValueRow) {
    hasDefaultValues = getInitValues(defaultDirectionLine, optionalFieldCount, 
                                     init_values);
  }
  if (hasDefaultValues && init_values.size() > static_cast<std::size_t>(numFeatures)) {
    ostringstream oss;
    oss << "ERROR: Reading tab file, too many default values present." << std::endl;
    Globals::failUnlessNoTerminate(oss.str(),
        "No-terminate flag set: ignoring error and trimming default value "
        "vector.");
    init_values.resize(static_cast<std::size_t>(numFeatures));
  }

  // read in the data
  if (rawWeights.size() > 0) {
    readAndScorePSMs(dataStream, psmLine, hasInitialValueRow, optionalFields,
                     rawWeights, allScores, foldWeights, lineFolds);
  } else {
    // detect if the input came from separate target and decoy searches or 
    // from a concatenated search by looking for scan+expMass combinations
    // that have both at least one target and decoy PSM
    bool concatenatedSearch = true;
    
    readPSMs(dataStream, psmLine, hasInitialValueRow, concatenatedSearch, optionalFields);
    
    check.checkAndSetDefaultDir();
    if (hasDefaultValues) check.addDefaultWeights(init_values); 
    check.setConcatenatedSearch(concatenatedSearch);
  }
  return 1;
}

void SetHandler::readAndScorePSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields, 
    std::vector<double>& rawWeights, Scores& allScores,
    const std::vector<std::vector<double> >& foldWeights,
    const std::vector<std::pair<unsigned int, unsigned int> >& lineFolds) {
  unsigned int lineNr = (hasInitialValueRow ? 3u : 2u);
  bool readProteins = true;
  size_t nextFold = 0;  // first entry of lineFolds not before lineNr
  // (index in allScores, fold) of the PSMs scored with a fold's weights
  std::vector<std::pair<size_t, unsigned int> > foldRows;
  foldRows.reserve(lineFolds.size());
  try {  // not indented: rethrowWithSource adds the file of the line
  do {
    if (lineNr % 1000000 == 0 && VERB > 1) {
      std::cerr << "Processing line " << lineNr << std::endl;
    }
    psmLine = rtrim(psmLine);
    // DataSet::readPsm would read any label other than 1 as a decoy, while
    // the training data ignore such lines
    TabReader reader(psmLine);
    PsmHeader header;
    const bool simple =
        DataSet::readPsmHeader(reader, lineNr, optionalFields, header);
    int label = header.label;
    if (!simple) {
      TabReader labelReader(psmLine);
      labelReader.skip(1u);  // id
      label = labelReader.readInt();
    }
    if (label == 1 || label == -1) {
      ScoreHolder sh;
      if (!simple ||
          !DataSet::readPsmBody(reader, header, readProteins, sh.pPSM, sh.label,
                                featurePool_, psmStore_, decoyPrefix_)) {
        sh.label = DataSet::readPsm(psmLine, lineNr, optionalFields, readProteins,
                                    sh.pPSM, featurePool_, psmStore_, decoyPrefix_);
      }
      while (nextFold < lineFolds.size() && lineFolds[nextFold].first < lineNr) {
        ++nextFold;
      }
      const bool inFold = nextFold < lineFolds.size() &&
                          lineFolds[nextFold].first == lineNr;
      if (inFold) {
        foldRows.emplace_back(
            static_cast<size_t>(allScores.end() - allScores.begin()),
            lineFolds[nextFold].second);
      }
      allScores.scoreAndAddPSM(
          sh, inFold ? foldWeights.at(lineFolds[nextFold].second) : rawWeights,
          featurePool_);
    }
    incrementLineNr(lineNr);
  } while (getline(dataStream, psmLine));
  } catch (const MyException& e) {
    rethrowWithSource(e, lineNr);
  }
  rescaleFoldScores(allScores, foldRows, foldWeights.size());

  if (VERB > 1) {
    std::cerr << "Found " << lineNr - (hasInitialValueRow ? 3u : 2u) << " PSMs" << std::endl;
  }
}

/**
 * Puts the scores of the PSMs scored with a fold's weights on the scale of
 * the other PSMs, which the averaged weights scored: for each fold, an
 * affine map gives its PSMs the mean and standard deviation of the scores
 * of all other PSMs. Both groups are random samples of the PSMs, and the
 * map uses no labels. With the merge scheme test, the weights of each fold
 * were rescaled by merge() with its test set's labels (score 0 at q = fdr),
 * which on a small test set favours a fold whose decoys happen to score low;
 * the map undoes that. With the other schemes (--merge-normalization) it
 * puts the fold's scores on the scale of the averaged weights in the same
 * way. Nothing changes if there are no other PSMs.
 * @param rows (index in allScores, fold) of the PSMs scored by a fold
 */
void SetHandler::rescaleFoldScores(
    Scores& allScores, const std::vector<std::pair<size_t, unsigned int> >& rows,
    size_t numFolds) {
  const size_t numScores = static_cast<size_t>(allScores.end() - allScores.begin());
  if (rows.empty() || rows.size() >= numScores || numFolds == 0) return;
  std::vector<double> n(numFolds, 0.0), sum(numFolds, 0.0), sumSq(numFolds, 0.0);
  double restN = 0.0, restSum = 0.0, restSumSq = 0.0;
  size_t next = 0;
  Scores::iterator it = allScores.begin();
  for (size_t ix = 0; ix < numScores; ++ix, ++it) {
    const double s = it->score;
    if (next < rows.size() && rows[next].first == ix) {
      const unsigned int fold = rows[next++].second;
      n[fold] += 1.0;
      sum[fold] += s;
      sumSq[fold] += s * s;
    } else {
      restN += 1.0;
      restSum += s;
      restSumSq += s * s;
    }
  }
  if (restN < 2.0) return;
  const double restMean = restSum / restN;
  const double restSd = std::sqrt(std::max(0.0,
      (restSumSq - restN * restMean * restMean) / (restN - 1.0)));
  std::vector<double> scale(numFolds, 1.0), offset(numFolds, 0.0);
  for (size_t fold = 0; fold < numFolds; ++fold) {
    if (n[fold] < 2.0) continue;
    const double mean = sum[fold] / n[fold];
    const double sd = std::sqrt(std::max(0.0,
        (sumSq[fold] - n[fold] * mean * mean) / (n[fold] - 1.0)));
    if (!(sd > 0.0) || !(restSd > 0.0)) continue;
    scale[fold] = restSd / sd;
    offset[fold] = restMean - mean * scale[fold];
  }
  for (const auto& row : rows) {
    ScoreHolder& sh = *(allScores.begin() + static_cast<std::ptrdiff_t>(row.first));
    sh.score = sh.score * scale[row.second] + offset[row.second];
  }
}

std::string& SetHandler::rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
  return s;
}
