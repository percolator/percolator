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

bool DataSet::isotopeMass = false;
bool DataSet::calcDOC = false;
const string DataSet::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";

string DataSet::ptmAlphabet = "#*@";
FeatureNames DataSet::featureNames;

DataSet::DataSet() {
  numSpectra = 0;
  sqtFN = "";
  pattern = "";
  doPattern = false;
  matchPattern = true;
}

DataSet::~DataSet() {

  for(unsigned i = 0; i < psms.size(); i++)
  {
    if(psms[i]->features)
    {
      delete[] psms[i]->features;
      psms[i]->features = NULL;
    }

    if(psms[i]->retentionFeatures)
    {
      delete[] psms[i]->retentionFeatures;
      psms[i]->retentionFeatures = NULL;
    }

    delete psms[i];
    psms[i] = NULL;
  }

}

PSMDescription* DataSet::getNext(int& pos) {
  pos++;
  if (pos < 0) {
    pos = 0;
  }
  if (pos >= getSize()) {
    return NULL;
  }
  return psms[pos];
}

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out, const string& lab) {
  int pos = -1;
  PSMDescription* pPSM = NULL;
  unsigned int nf = FeatureNames::getNumFeatures();
  if (calcDOC) {
    nf -= DescriptionOfCorrect::numDOCFeatures();
  }
  while ((pPSM = getNext(pos)) != NULL) {
    double* frow = pPSM->features;
    out << pPSM->id << '\t' << lab;
    if (calcDOC) {
      out << '\t' << pPSM->getUnnormalizedRetentionTime() << '\t'
          << pPSM->massDiff;
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << frow[ix];
    }
    out << '\t' << pPSM->peptide;
    for (const auto &proteinId : pPSM->proteinIds) {
      out << '\t' << proteinId;
    }
    out << endl;
  }
  return true;
}

void DataSet::print_features() {
  for (int i = 0; i < getSize(); i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << psms[i]->features[j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features()
{
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (int i = 0; i < 10; i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

void DataSet::print(Scores& test, vector<ResultHolder> &outList)
{
  ostringstream out;
  size_t ix = 0;
  vector<PSMDescription*>::const_iterator psm = psms.begin();
  for (; psm != psms.end(); psm++, ix++) {
    ScoreHolder* pSH = test.getScoreHolder((*psm)->features);
    if (pSH == NULL) {
      continue;
    }
    double score = pSH->score;
    for (const auto &proteinId : (*psm)->proteinIds) {
      out << "\t" << proteinId;
    }
    ResultHolder rh(score, (*psm)->q, (*psm)->pep, (*psm)->id, (*psm)->peptide, out.str());
    outList.push_back(rh);
    out.str("");
  }
}

/*
double DataSet::isPngasef(const string& peptide) {
  bool isDecoy;
  switch (label) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("ERROR : class DataSet has not been initiated\
		to neither target nor decoy label\n");}
  }
  return isPngasef( peptide, isDecoy);
}
//
double DataSet::isPngasef(const string& peptide, bool isDecoy ) {
  size_t next_pos = 0, pos;
  while ((pos = peptide.find("N*", next_pos)) != string::npos) {
    next_pos = pos + 1;
    if (! isDecoy) {
      pos += 3;
      if (peptide[pos] == '#') {
        pos += 1;
      }
    } else {
      pos -= 2;
      if (peptide[pos] == '#') {
        pos -= 1;
      }
    }
    if (peptide[pos] == 'T' || peptide[pos] == 'S') {
      return 1.0;
    }
  }
  return 0.0;
}
*/

void DataSet::readTabData(ifstream& is, const vector<unsigned int>& ixs) {
  // count number of features
  string tmp, line;
  is.clear();
  is.seekg(0, ios::beg);
  getline(is, line); // skip over line with column names
  getline(is, line);
  unsigned int m = 0, n = ixs.size();
  istringstream buff(line);
  double a;
  buff >> tmp >> tmp; // remove id and label
  while (true) {
    buff >> a;
    if (buff.good()) {
      m++;
    } else {
      buff >> tmp;
    }
    if (!buff) {
      break;
    }
    buff.clear();
  }

  if (m < 1) {
    is.close();
    throw MyException("Error : Reading tab file, too few features present.");
  }
  if (calcDOC) {
    m -= 2;
  }

  initFeatureTables((calcDOC ? m + DescriptionOfCorrect::numDOCFeatures(): m),calcDOC);
  if (calcDOC) {
    getFeatureNames().setDocFeatNum(m);
  }
  string seq;
  
  // read in the data
  is.clear();
  is.seekg(0, ios::beg);
  getline(is, line); // id line
  unsigned int ix = 0;
  getline(is, line);
  for (unsigned int i = 0; i < n; i++) {
    while (ix < ixs[i]) {
      getline(is, line);
      ix++;
    }
    PSMDescription  *myPsm = new PSMDescription();
    buff.str(line);
    buff.clear();
    buff >> myPsm->id;
    buff >> tmp; // get rid of label
    double *featureRow = new double[m];
    myPsm->features = featureRow;
    if (calcDOC) {
      buff >> myPsm->retentionTime;
      buff >> myPsm->massDiff;
    }
    for (register unsigned int j = 0; j < m; j++) {
      buff >> featureRow[j];
    }
    std::string peptide_seq = "";
    buff >> peptide_seq;
    cerr << peptide_seq << endl;
    //NOTE to check if the peptide sequence contains flanks or not
    if(peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.')
    {
      is.close();
      ostringstream temp;
      temp << "Error : Reading tab file, the peptide sequence " << peptide_seq << " \
      does not contain one or two of its flaking amino acids." << std::endl;
      throw MyException(temp.str());
    }
    myPsm->peptide = peptide_seq;

    while (!!buff) {
      buff >> tmp;
      if (tmp.size() > 0) {
        myPsm->proteinIds.insert(tmp);
      }
    }
    if (calcDOC) {
      DescriptionOfCorrect::calcRegressionFeature(*psms[i]);
      featureRow[m] = abs(psms[i]->pI - 6.5);
      featureRow[m + 1] = abs(psms[i]->massDiff);
      featureRow[m + 2] = 0;
    }
    psms.push_back(myPsm);
    ++numSpectra;
  }
}

unsigned int DataSet::peptideLength(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (aaAlphabet.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

unsigned int DataSet::cntPTMs(const string& pep) {
  unsigned int len = 0;
  for (string::size_type pos = 2; (pos + 2) < pep.size(); pos++) {
    if (ptmAlphabet.find(pep.at(pos)) != string::npos) {
      len++;
    }
  }
  return len;
}

void DataSet::initFeatureTables(const unsigned int numFeat, bool __regressionTable)
{
  FeatureNames::setNumFeatures(numFeat);
  regresionTable = __regressionTable;
  psms.clear();
}

#ifdef XML_SUPPORT
// Convert a peptide with or without modifications into a string
std::string DataSet::decoratePeptide(const ::percolatorInNs::peptideType& peptide) {
  std::list<std::pair<int,std::string> > mods;
  std::string peptideSeq = peptide.peptideSequence();
  for(const auto &mod_ref : peptide.modification()){
    std::stringstream ss;
    if (mod_ref.uniMod().present()) {
      ss << "[UNIMOD:" << mod_ref.uniMod().get().accession() << "]";
      mods.push_back(std::pair<int,std::string>(mod_ref.location(),ss.str()));
    }
    if (mod_ref.freeMod().present()) {
      ss << "[" << mod_ref.freeMod().get().moniker() << "]";
      mods.push_back(std::pair<int,std::string>(mod_ref.location(),ss.str()));
    }
  }
  mods.sort(greater<std::pair<int,std::string> >());
  std::list<std::pair<int,std::string> >::const_iterator it;
  for(it=mods.begin();it!=mods.end();++it) {
    peptideSeq.insert(it->first,it->second);
  }
  return peptideSeq;
}

void DataSet::readPsm(const percolatorInNs::peptideSpectrumMatch& psm, unsigned scanNumber)
{
  bool isDecoy;
  switch (label) {
    case 1: { isDecoy = false; break; };
    case -1: { isDecoy = true; break; };
    default:  { throw MyException("Error : Reading PSM, class DataSet has not been initiated\
		to neither target nor decoy label\n");}
  }

  if(psm.isDecoy() != isDecoy)
  {
    ostringstream temp;
    temp << "Error : adding PSM " << psm.id() << " to the dataset.\n\
    The label isDecoy of the PSM is not the same in the dataset." << std::endl;
    throw MyException(temp.str());
  }
  else
  {
      PSMDescription  *myPsm = new PSMDescription();
      string mypept = decoratePeptide(psm.peptide());

      if(psm.occurence().size() <= 0)
      {
	ostringstream temp;
	temp << "Error: adding PSM " << psm.id() << " to the dataset.\n\
	The PSM does not contain protein occurences." << std::endl;
	throw MyException(temp.str());
      }

      for( const auto & oc : psm.occurence() )
      {
        myPsm->proteinIds.insert( oc.proteinId() );
        // adding n-term and c-term residues to peptide
	//NOTE the residues for the peptide in the PSMs are always the same for every protein
        myPsm->peptide = oc.flankN() + "." + mypept + "." + oc.flankC();
      }

      myPsm->id = psm.id();
      myPsm->charge = psm.chargeState();
      myPsm->scan = scanNumber;
      myPsm->expMass = psm.experimentalMass();
      myPsm->calcMass = psm.calculatedMass();
      if ( psm.observedTime().present() ) {
        myPsm->retentionTime = psm.observedTime().get();
      }

      const ::percolatorInNs::features::feature_sequence & featureS = psm.features().feature();
      int featureNum = 0;

      myPsm->features = new double[FeatureNames::getNumFeatures()];
      if (regresionTable)
      {
	myPsm->retentionFeatures = new double[RTModel::totalNumRTFeatures()];
      }

      for ( ::percolatorInNs::features::feature_const_iterator featureIter = featureS.begin(); featureIter != featureS.end(); featureIter++ ) {
        myPsm->features[featureNum]=*featureIter;
        featureNum++;
      }

      // myPsm.peptide = psmIter->peptide().peptideSequence();
      myPsm->massDiff = MassHandler::massDiff(psm.experimentalMass() ,psm.calculatedMass(),psm.chargeState());

      if (calcDOC)
      {
        DescriptionOfCorrect::calcRegressionFeature(*myPsm);
        myPsm->features[featureNum++] = abs( myPsm->pI - 6.5);
        myPsm->features[featureNum++] = abs( myPsm->massDiff);
        myPsm->features[featureNum++] = 0;
        myPsm->features[featureNum++] = 0;
      }

      psms.push_back(myPsm);
      ++numSpectra;
  }
}
#endif // XML_SUPPORT
