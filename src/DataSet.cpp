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

#include "DataSet.h"
#include <cmath>

FeatureNames DataSet::featureNames_;
bool DataSet::decoyWarningTripped_ = false;
bool DataSet::flankWarningTripped_ = false;

DataSet::DataSet() {}

DataSet::~DataSet() {}

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    if (psm->hasSpectrumFileName()) {
      out << '\t' << psm->getSpectrumFileName();
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << psm->getFullPeptide() << '\t';
    psm->printProteins(out);
    out << '\n';
  }
  return true;
}

// TODO: find a way to make these three functions generic
void DataSet::fillScores(std::vector<ScoreHolder>& scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm));
  }
}

void DataSet::clearFeaturePointers() {
  for (PSMDescription* psm : psms_) {
    psm->features = NULL;
  }
}

void DataSet::fillFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

/**
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
void DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, FeatureMemoryPool& featurePool,
    PsmStore& psmStore, const std::string& decoyPrefix) {
  PSMDescription* myPsm = NULL;
  bool readProteins = true;
  readPsm(line, lineNr, optionalFields, readProteins, myPsm, featurePool,
          psmStore, decoyPrefix);
  registerPsm(myPsm);
}

LabelType DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool,
    PsmStore& psmStore, const std::string& decoyPrefix) {
  TabReader reader(line);
  
  // Parsed into a local object, which only goes into the store when the
  // line is complete. Until then its strings are views into line.
  PSMDescription psm;
  psm.setId(reader.readView());
  LabelType label = reader.readInt() == 1 ? LabelType::TARGET : LabelType::DECOY;
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        psm.scan = static_cast<unsigned int>(reader.readInt());
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number of PSM " 
              << psm.getId() << ". Check if scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        psm.expMass = reader.readDouble();
        // A NaN would break the strict weak order that the comparators in
        // ScoreHolder.h need. SetHandler::getScanId checks every line first and
        // reports such a value also with NO_TERMINATE, so in that case it is
        // set to 0 here without a second message.
        if (!isfinite(psm.expMass)) {
          if (!NO_TERMINATE) {
            ostringstream temp;
            temp << "ERROR: Reading tab file, the experimental mass of PSM "
                << psm.getId() << " on line " << lineNr << " is "
                << psm.expMass << " instead of a finite number." << std::endl;
            throw MyException(temp.str());
          }
          psm.expMass = 0.0;
        }
        break;
      } case CALCMASS: {
        psm.calcMass = reader.readDouble();
        break;
      } case RETTIME: {
        psm.setRetentionTime(reader.readDouble());
        break;
      } case DELTAMASS: {
          ostringstream temp;
          temp << "ERROR: The parser tries to allocate a deprecated delta mass field for the PSM " 
              << psm.getId() << std::endl;
          throw MyException(temp.str());
        break;
      } case FILENAME: {
        psm.setSpectrumFileName(reader.readView());
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) psm.scan = lineNr;
  psm.setLineNr(lineNr);
  
  unsigned int numFeatures = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  double* featureRow = featurePool.allocate();
  psm.features = featureRow;
  for (unsigned int j = 0; j < numFeatures; j++) {
    featureRow[j] = reader.readDouble();
    if (!isfinite(featureRow[j])) {
      ostringstream oss;
      oss << "ERROR: Reached strange feature with val=" << featureRow[j]
           << " col=" << j << " for PSM with id " << psm.getId() << endl;
      if (NO_TERMINATE) {
        std::cerr << oss.str();
        std::cerr << "No-terminate flag set: setting value to 0 and ignoring the error." << std::endl;
        featureRow[j] = 0.0;
      } else {
        throw MyException(oss.str() + "Terminating.\n");
      }
    }
  }
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in feature vector of PSM " 
      << psm.getId() << ". Check if there are enough features on this line and "
      << "if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  psm.setPeptide(reader.readView());
  const std::string_view peptide_seq = psm.getFullPeptide();
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM " << psm.getId() 
      << ". Check if a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (!PSMDescription::hasFlanks(peptide_seq)) {
    checkUnflankedPeptide(peptide_seq, psm.getId(), lineNr);
  }
  
  std::vector<std::string_view> proteins;
  if (readProteins) {
    if (PSMDescription::getProteinNameSeparator() == "\t") {
      while (!reader.error()) {
        std::string_view tmp = reader.readView();
        if (tmp.size() > 0) proteins.push_back(tmp);
      }
    } else {
      const std::string_view names = reader.readView();
      const std::string& separator = PSMDescription::getProteinNameSeparator();
      size_t start = 0, pos = 0;
      while ((pos = names.find(separator, start)) != std::string_view::npos) {
        if (pos > start) proteins.push_back(names.substr(start, pos - start));
        start = pos + separator.length();
      }
      if (start < names.size()) proteins.push_back(names.substr(start));

      // check if end of line was reached after reading the protein string
      // if this is not the case, it means the user probably incorrectly provided proteins as a tab separated input
      if (!reader.error()) {
        ostringstream temp;
        temp << "ERROR: Reading tab file, error reading PSM " << psm.getId() 
          << ". Found tab character after reading proteins even though " << PSMDescription::getProteinNameSeparator() 
          << " was specified as --protein-name-separator." << std::endl;
        if (NO_TERMINATE) {
          std::cerr << temp.str();
          std::cerr << "No-terminate flag set: ignoring error and continuing "
                    << "without PSMs." << std::endl;
        } else {
          throw MyException(temp.str());
        }
      }
    }
  }

  if (label == LabelType::DECOY && VERB > 1 && !decoyWarningTripped_) {
    for (auto const& proteinId: proteins) { 
      bool startsWithDecoyPrefix = (proteinId.rfind(decoyPrefix, 0) == 0);
      if (!startsWithDecoyPrefix && !decoyWarningTripped_) {
        std::cerr << "Warning: protein decoy prefix " << decoyPrefix 
                  << " doesn't match the decoy protein identifier " 
                  << proteinId << "." << std::endl;
        decoyWarningTripped_ = true;
      }
    }
  }
  psm.setId(psmStore.storeString(psm.getId()));
  psm.setPeptide(psmStore.storeString(psm.getFullPeptide()));
  psm.setProteins(psmStore.storeProteins(proteins), proteins.size());
  myPsm = psmStore.add(std::move(psm));
  return label;
}

/**
 * Handles a peptide that is not of the form X.PEPTIDE.Y. Protein inference
 * needs the flanks, and a number in the peptide column means that the line
 * has more feature columns than the first PSM line, so both are errors.
 * Otherwise the peptide is used as it is (see
 * PSMDescription::getPeptideSequenceView), with one warning per run.
 */
void DataSet::checkUnflankedPeptide(std::string_view peptide,
                                    std::string_view psmId,
                                    unsigned int lineNr) {
  // a digit is required, as the peptides NAN and INF would parse as numbers
  std::string field(peptide);
  TabReader numberReader(field);
  numberReader.readDouble();
  if (!numberReader.error() &&
      std::any_of(field.begin(), field.end(),
                  [](unsigned char c) { return std::isdigit(c) != 0; })) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide column of PSM " << psmId
         << " on line " << lineNr << " holds the number " << peptide
         << ". The line has more feature columns than the first PSM line ("
         << FeatureNames::getNumFeatures() << " features)." << std::endl;
    throw MyException(temp.str());
  }
  if (ProteinProbEstimator::getCalcProteinLevelProb()) {
    // the protein inference digests the protein sequences and needs flanks
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide " << peptide << " of PSM "
         << psmId << " on line " << lineNr << " is not of the form "
         << "X.PEPTIDE.Y with one flanking amino acid (or '-') on each side, "
         << "which protein inference needs." << std::endl;
    throw MyException(temp.str());
  }
  if (!flankWarningTripped_) {
    std::cerr << "Warning: the peptide " << peptide << " of PSM " << psmId
              << " on line " << lineNr << " is not of the form X.PEPTIDE.Y "
              << "with one flanking amino acid (or '-') on each side. Such "
              << "peptides are used as they are." << std::endl;
    flankWarningTripped_ = true;
  }
}

void DataSet::registerPsm(PSMDescription* myPsm) {
  switch (label_) {
    case LabelType::TARGET: { break; };
    case LabelType::DECOY: { break; };
    default:  { throw MyException("ERROR : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label\n");}
  }
  psms_.push_back(myPsm);
}
