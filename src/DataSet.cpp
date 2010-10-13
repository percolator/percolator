/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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

int DataSet::hitsPerSpectrum = 1;
bool DataSet::calcQuadraticFeatures = false;
bool DataSet::calcAAFrequencies = false;
bool DataSet::calcPTMs = false;
bool DataSet::isotopeMass = false;
bool DataSet::calcDOC = false;
bool DataSet::pngasef = false;
const string DataSet::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";

string DataSet::ptmAlphabet = "#*@";
FeatureNames DataSet::featureNames;

static char buf[4096];

DataSet::DataSet() {
  feature = NULL;
  regressionFeature = NULL;
  numSpectra = 0;
  sqtFN = "";
  pattern = "";
  doPattern = false;
  matchPattern = true;
  psmNum = 0;
}

DataSet::~DataSet() {
  if (feature) {
    delete[] feature;
    feature = NULL;
  }
  if (regressionFeature) {
    delete[] regressionFeature;
    regressionFeature = NULL;
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
  return &psms[pos];
}

const double* DataSet::getFeatures(const int pos) const {
  return &feature[DataSet::rowIx(pos)];
}

bool DataSet::writeTabData(ofstream& out, const string& lab) {
  int pos = -1;
  PSMDescription* pPSM = NULL;
  unsigned int nf = FeatureNames::getNumFeatures();
  if (calcDOC) {
    nf -= DescriptionOfCorrect::numDOCFeatures();
  }
  while ((pPSM = getNext(pos)) != NULL) {
    double* frow = pPSM->features;
    out << psms[pos].id << '\t' << lab;
    if (calcDOC) {
      out << '\t' << psms[pos].getUnnormalizedRetentionTime() << '\t'
          << psms[pos].massDiff;
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << frow[ix];
    }
    out << "\t" << pPSM->peptide;
    set<string>::const_iterator it = pPSM->proteinIds.begin();
    for (; it != pPSM->proteinIds.end(); it++) {
      out << "\t" << *it;
    }
    out << endl;
  }
  return true;
}

void DataSet::print_features() {
  for (int i = 0; i < getSize(); i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << feature[DataSet::rowIx(i) + j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features() {
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (int i = 0; i < 10; i++) {
    for (unsigned int j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << feature[DataSet::rowIx(i) + j] << "\t";
    }
    cerr << endl;
  }
}

void DataSet::print(Scores& test, vector<ResultHolder> &outList) {
  ostringstream out;
  size_t ix = 0;
  vector<PSMDescription>::const_iterator psm = psms.begin();
  for (; psm != psms.end(); psm++, ix++) {
    ScoreHolder* pSH = test.getScoreHolder(psm->features);
    if (pSH == NULL) {
      continue;
    }
    double score = pSH->score;
    double q = psm->q;
    double pep = psm->pep;
    set<string> prots = psm->proteinIds;
    set<string>::const_iterator it = psm->proteinIds.begin();
    for (; it != psm->proteinIds.end(); it++) {
      out << "\t" << *it;
    }
    ResultHolder rh(score, q, pep, psm->id, psm->peptide, out.str());
    outList.push_back(rh);
    out.str("");
  }
}

double DataSet::isPngasef(const string& peptide) {
  bool isDecoy;
  switch (label) {
  case 1: { isDecoy = false; break; };
  case -1: { isDecoy = true; break; };
  default:  { fprintf(stderr,"programming error 123234123\n"); exit(EXIT_FAILURE); } 
  }
  return isPngasef( peptide, isDecoy);
}

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

void DataSet::readTabData(ifstream& is, const vector<unsigned int>& ixs) {
  string tmp, line;
  is.clear();
  is.seekg(0, ios::beg);
  getline(is, line);
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
    cerr << "To few features in Tab data file";
    exit(-1);
  }
  if (calcDOC) {
    m -= 2;
  }
  initFeatureTables((calcDOC ? m + DescriptionOfCorrect::numDOCFeatures()
  : m), n, calcDOC);
  if (calcDOC) {
    getFeatureNames().setDocFeatNum(m);
  }
  string seq;
  is.clear();
  is.seekg(0, ios::beg);
  getline(is, line); // id line
  //  getFeatureNames().setFeatures(line,2,m);
  unsigned int ix = 0;
  getline(is, line);
  for (unsigned int i = 0; i < n; i++) {
    while (ix < ixs[i]) {
      getline(is, line);
      ix++;
    }
    buff.str(line);
    buff.clear();
    buff >> psms[i].id;
    buff >> tmp; // get rid of label
    double* featureRow = &feature[rowIx(i)];
    psms[i].features = featureRow;
    if (calcDOC) {
      buff >> psms[i].retentionTime;
      buff >> psms[i].massDiff;
    }
    for (register unsigned int j = 0; j < m; j++) {
      buff >> featureRow[j];
    }
    buff >> psms[i].peptide;
    while (!!buff) {
      buff >> tmp;
      if (tmp.size() > 0) {
        psms[i].proteinIds.insert(tmp);
      }
    }
    if (calcDOC) {
      DescriptionOfCorrect::calcRegressionFeature(psms[i]);
      featureRow[m] = abs(psms[i].pI - 6.5);
      featureRow[m + 1] = abs(psms[i].massDiff);
      featureRow[m + 2] = 0;
    }
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

void DataSet::readFragSpectrumScans( const ::percolatorInNs::fragSpectrumScan & fss) {

  bool isDecoy;
  switch (label) {
  case 1: { isDecoy = false; break; };
  case -1: { isDecoy = true; break; };
  default:  { fprintf(stderr,"programming error 123234123\n"); exit(EXIT_FAILURE); } 
  }
  const ::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_sequence & psmSeq = fss.peptideSpectrumMatch();
  for ( ::percolatorInNs::fragSpectrumScan::peptideSpectrumMatch_const_iterator psmIter = psmSeq.begin(); psmIter != psmSeq.end(); ++psmIter) {

    if ( psmIter->isDecoy() == isDecoy ) {

      assert( psms.size() > psmNum );
      PSMDescription & myPsm = psms[psmNum];

      // rng:oneOrMore so the assert should always be true
      assert( psmIter->occurence().size() > 0 );


      BOOST_FOREACH( const percolatorInNs::occurence & oc,  psmIter->occurence() )  {
        myPsm.proteinIds.insert( oc.proteinId() );
        // adding n-term and c-term residues to peptide
        myPsm.peptide = oc.flankN() + "." + psmIter->peptide().peptideSequence() + "." + oc.flankC();
      }
      myPsm.id = psmIter->id();
      myPsm.scan = fss.scanNumber();
      if(MassHandler::monoisotopic == true){
        if(! psmIter->experimentalMassToCharge().present()){
          cerr << "\nYou have selected the -M option, but no experimental mass "
              << "information is available in the pin file you have given in input "
              << "for peptideSpectrumMatch elements.\n";
          exit(-1);
        }
        myPsm.expMass = psmIter->experimentalMassToCharge().get();
      }
      myPsm.calcMass = psmIter->calculatedMassToCharge();
      if ( psmIter->observedTime().present() ) {
        myPsm.retentionTime = psmIter->observedTime().get();
      }

      const ::percolatorInNs::features::feature_sequence & featureS = psmIter->features().feature();
      int featureNum = 0;

      for ( ::percolatorInNs::features::feature_const_iterator featureIter = featureS.begin(); featureIter != featureS.end(); featureIter++ ) {
        myPsm.features[featureNum]=*featureIter;
        featureNum++;
      }

      // myPsm.peptide = psmIter->peptide().peptideSequence();

      myPsm.massDiff =
          MassHandler::massDiff(fss.experimentalMassToCharge() ,
              psmIter->calculatedMassToCharge(),
              psmIter->chargeState(),
              myPsm.peptide.substr(2, myPsm.peptide.size()
                  - 4));

      if (calcDOC) {
        // These features will be set before each iteration
        DescriptionOfCorrect::calcRegressionFeature(myPsm);
        myPsm.features[featureNum++] = abs( myPsm.pI - 6.5);
        myPsm.features[featureNum++] = abs( myPsm.massDiff);
        // myPsm.features[featureNum++]=abs(psm.retentionTime);
        myPsm.features[featureNum++] = 0;
        myPsm.features[featureNum++] = 0;
      }
      ++psmNum;
    } 
  }
  return;
}


void DataSet::initFeatureTables(const unsigned int numFeat,
    const unsigned int numSpec,
    bool regressionTable) {
  FeatureNames::setNumFeatures(numFeat);
  numSpectra = numSpec;

  feature = new double[numFeat * numSpec];
  if (regressionTable) {
    regressionFeature
    = new double[RTModel::totalNumRTFeatures() * numSpec];
  }
  psms.resize(numSpectra);
  for (int ix = 0; ix < numSpectra; ++ix) {
    psms[ix].features = &feature[rowIx(ix)];
  }
  if (regressionTable) {
    double* ptr = regressionFeature;
    size_t nf = RTModel::totalNumRTFeatures();
    for (int ix = 0; ix < numSpectra; ++ix, ptr += nf) {
      psms[ix].retentionFeatures = ptr;
    }
  }
}
