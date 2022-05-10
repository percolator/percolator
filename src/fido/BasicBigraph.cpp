// Written by Oliver Serang 2009
// see license for more information

#include "BasicBigraph.h"

BasicBigraph::BasicBigraph(): PsmThreshold(0.0), PeptideThreshold(1e-3),
  ProteinThreshold(1e-3) {}

BasicBigraph::BasicBigraph(double __psmthreshold, double __peptidethreshold, 
    double __proteinthreshold):
  PsmThreshold(__psmthreshold), PeptideThreshold(__peptidethreshold),
  ProteinThreshold(__proteinthreshold) {}

BasicBigraph::~BasicBigraph() {}

void BasicBigraph::read(Scores* fullset, bool multiple_labeled_peptides) {
  string pepName, protName;
  double value =  -10;
  int pepIndex = -1;
  StringTable PSMNames, proteinNames;

  vector<ScoreHolder>::iterator psm = fullset->begin();
  for (; psm!= fullset->end(); ++psm) {
    // e peptide_string
    pepName = psm->pPSM->getFullPeptide();
    
    if ( pepName[1] == '.' ) {
      // trim off the cleavage events
      pepName = pepName.substr(2, pepName.size() - 4);
    }
    //NOTE fido will keep only one peptide in the case that a target and a decoy peptide
    //      contain the same sequence
    //     this is a bit of a tricky situation because for separate target-decoy searches
    //      two unique peptides ( target and decoy ) might have the same sequence, therefore
    //     fido assumes that peptide is connected to the correspondent target and decoy proteins
    //     however, this is not the real situation. This scenario is not common but as a very rough
    //    and quick wordaround I am appending a * to all the decoy peptides to distinguish them
    //    from the target peptides. 
    if (psm->isDecoy() && multiple_labeled_peptides) {
      pepName += "*";
    }
    
    if ( PSMNames.lookup(pepName) == -1 ){
      add(PSMsToProteins, PSMNames, pepName);
    }
    pepIndex = PSMNames.lookup(pepName);

    // r proteins
    std::vector<string>::const_iterator pid = psm->pPSM->proteinIds.begin();
    for (; pid!= psm->pPSM->proteinIds.end(); ++pid) {
      protName = getRidOfUnprintablesAndUnicode(*pid);
      if (proteinNames.lookup(protName) == -1) {
        add(proteinsToPSMs, proteinNames, protName);
      }
      connect(PSMNames, pepName, proteinNames, protName);
    }
    // p probability of the peptide match to the spectrum
    value = 1 - psm->pep;
    PSMsToProteins.weights[ pepIndex ] = max(PSMsToProteins.weights[pepIndex], value);
 }

  PSMsToProteins.names = PSMNames.getItemsByNumber();
  proteinsToPSMs.names = proteinNames.getItemsByNumber();
  
  //NOTE this function is assigning PeptideThreshold probablity to all the PSMs with a prob below PeptideThreshold
  /**pseudoCountPSMs();**/
}


void BasicBigraph::read(istream & is, bool multiple_labeled_peptides) {
  char instr;
  string pepName, protName;
  double value =  -10;

  int pepIndex = -1;
  int state = 'e';

  StringTable PSMNames, proteinNames;

  while (is >> instr) {
    if (instr == 'e' && (state == 'e' || state == 'p')) {
      if (state == 'p') {
        cerr << "Warning: no peptide score for peptide entry " << pepName << ", using last score (" << value << ")" << endl;
        value = -1;
        if ( value == -10 ) {
          throw MyException("Error: No previous peptide entry to use");
        }
        PSMsToProteins.weights[ pepIndex ] = max(value, PSMsToProteins.weights[ pepIndex ]);
      }
      is >> pepName;
      //pepName = cleanPeptideSequence(pepName);
      
      if (PSMNames.lookup(pepName) == -1) {
        add(PSMsToProteins, PSMNames, pepName);
      }
      
      pepIndex = PSMNames.lookup(pepName);
      state = 'c';
    } else if (instr == 'c' && state == 'c') {
      state = 'r';
    } else if ( instr == 'r' && ( state == 'c' || state == 'r' || state == 'p' ) ) {
      is >> protName;

      if ( proteinNames.lookup(protName) == -1 )
        add(proteinsToPSMs, proteinNames, protName);

      connect(PSMNames, pepName, proteinNames, protName);
      state = 'p';
    } else if ( instr == 'p' && state == 'p' ) {
      is >> value;
      // this option scores peptides using only the best match
      PSMsToProteins.weights[ pepIndex ] = max(PSMsToProteins.weights[pepIndex], value);
      state = 'e';
    } else if ( instr == '#' ) {
      // comment line
      string garbage;
      getline(is, garbage);
    } else {
      throw MyException("");
    }
  }

  PSMsToProteins.names = PSMNames.getItemsByNumber();
  proteinsToPSMs.names = proteinNames.getItemsByNumber();

  //NOTE this function is assigning PeptideThreshold probablity to all the PSMs with a prob below PeptideThreshold
  /**pseudoCountPSMs();**/
}

void BasicBigraph::printGraphStats() const
{
  cout << "There are \t" << PSMsToProteins.size() << " PSMs" << endl;
  cout << "      and \t" << proteinsToPSMs.size() << " proteins" << endl;

  int edgeCount = 0;
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      edgeCount += static_cast<int>(PSMsToProteins.associations[k].size());
    }

  cout << "      and \t" << edgeCount << " edges" << endl;
}


void BasicBigraph::saveSeveredProteins() {
  severedProteins = Array<string>();
  for (int k = 0; k < proteinsToPSMs.size(); k++) {
    if (proteinsToPSMs.associations[k].size() == 0) {
      severedProteins.add( proteinsToPSMs.names[k] );
    }
  }
}

void BasicBigraph::prune() {
  removePoorPSMs();
  removePoorProteins();
  saveSeveredProteins();
  reindex();
  markSectionPartitions();
  cloneMultipleMarkedPSMs();
  reindex();
}


void BasicBigraph::cloneMultipleMarkedPSMs() {
  numberClones = 0;
  // compute this once, since it will change as you add clones
  int N = PSMsToProteins.size();

  for (int k = 0; k < N; k++) {
    const GraphNode & psm = PSMsToProteins[k];
    // the way the marking procedure works, it will only multiple
    // mark PSMs with a score <= PeptideThreshold
    if (psm.sectionMark.size() > 1) {
      clonePSM(k);
    }
  }
}

void BasicBigraph::clonePSM(int pepIndex) {
  // add a new copy of this PSM for every section
  // first find the sections that this PSM associates with
  Set sections;
  const Set & s = PSMsToProteins.associations[pepIndex];
  for (int k = 0; k < s.size(); k++) {
    int sect = proteinsToPSMs[ s[k] ].section;
    sections |= Set::SingletonSet(sect);
  }
  // index the proteins by the sections they belong to
  Array<Set> associatedProteinsBySection(sections.size());
  for (int k = 0; k < s.size(); k++) {
    int sect = proteinsToPSMs.sections[ s[k] ];
    int ind = sections.find(sect);
    associatedProteinsBySection[ ind ] |= Set::SingletonSet(s[k]);
  }
  // add a clone for each of the sections-- include the first one,
  // since it will be easier to build them all than to have a special
  // case. 

  for (int k = 0; k < associatedProteinsBySection.size(); k++) {
    int sect = sections[k];
    ostringstream ost;
    ost << PSMsToProteins[ pepIndex ].name << "_clone_" << sect;
    
    PSMsToProteins.names.add( ost.str() );
    PSMsToProteins.associations.add( associatedProteinsBySection[k] );
    PSMsToProteins.weights.add( PSMsToProteins[ pepIndex ].weight );
    PSMsToProteins.sections.add(sect);

    // add the association from this section's proteins to the new peptide
    for (int j = 0; j<associatedProteinsBySection[k].size(); j++)
      proteinsToPSMs.associations[ associatedProteinsBySection[k][j] ] |= Set::SingletonSet( PSMsToProteins.size()-1 );
  }

  // afterward, erase the original
  // (remove the associations from proteins)
  for (int k = 0; k < associatedProteinsBySection.size(); k++) {
    for (int j = 0; j < associatedProteinsBySection[k].size(); j++) {
      int prot = associatedProteinsBySection[k][j];
      proteinsToPSMs[prot].association = proteinsToPSMs[prot].association.without( Set::SingletonSet(pepIndex) );
    }
  }

  // (remove the associations to proteins)
  PSMsToProteins[ pepIndex ].association = Set();
  numberClones += static_cast<int>(sections.size()-1);
}

void BasicBigraph::reindex() {
  int k;

  // note: later check to see if you can do such that in Array using a
  // pointer to a member function
  Set connectedPSMs;
  for (k = 0; k < PSMsToProteins.size(); k++) {
    if (!PSMsToProteins.associations[k].isEmpty())
      connectedPSMs.add(k);
  }

  Set connectedProteins;
  for (k = 0; k < proteinsToPSMs.size(); k++) {
    if (!proteinsToPSMs.associations[k].isEmpty())
      connectedProteins.add(k);
  }

  // note: ahh! this is bad design. You need to remake this code
  int backupNumberClones = numberClones;
  Array<string> backupSeveredProteins = severedProteins;
  double backupPeptideThreshold = PeptideThreshold;

  *this = buildSubgraph(connectedProteins, connectedPSMs);
  numberClones = backupNumberClones;
  severedProteins = backupSeveredProteins;
  PeptideThreshold = backupPeptideThreshold;
}

BasicBigraph BasicBigraph::buildSubgraph(const Set & connectedProteins, const Set & connectedPSMs) {
  BasicBigraph result;

  result.PSMsToProteins.names = PSMsToProteins.names[ connectedPSMs ];
  result.PSMsToProteins.associations = PSMsToProteins.associations[ connectedPSMs ];
  result.PSMsToProteins.weights = PSMsToProteins.weights[ connectedPSMs ];
  result.PSMsToProteins.sections = PSMsToProteins.sections[ connectedPSMs ];

  for (int k = 0; k < result.PSMsToProteins.associations.size(); k++) {
    result.PSMsToProteins.associations[k] = result.PSMsToProteins.associations[k].reindexToFind( connectedProteins );
  }

  result.proteinsToPSMs.names = proteinsToPSMs.names[ connectedProteins ];
  result.proteinsToPSMs.associations = proteinsToPSMs.associations[ connectedProteins ];
  result.proteinsToPSMs.weights = proteinsToPSMs.weights[ connectedProteins ];
  result.proteinsToPSMs.sections = proteinsToPSMs.sections[ connectedProteins ];

  for (int k = 0; k < result.proteinsToPSMs.associations.size(); k++) {
    result.proteinsToPSMs.associations[k] = result.proteinsToPSMs.associations[k].reindexToFind( connectedPSMs );
  }

  return result;
}

void BasicBigraph::removePoorPSMs() {
  for (int k = 0; k < PSMsToProteins.size(); k++) {
    if (PSMsToProteins.weights[k] < PsmThreshold) {
      disconnectPSM(k);
    }
  }
}

void BasicBigraph::removePoorProteins() {
  for (int k = 0; k < proteinsToPSMs.size(); k++) {
    if (Vector(PSMsToProteins.weights[ proteinsToPSMs.associations[k] ]).max() < ProteinThreshold) {
      disconnectProtein(k);
    }
  }
}

void BasicBigraph::disconnectPSM(int k) {
  Set& as = PSMsToProteins.associations[k];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++) {
    Set& setRef = proteinsToPSMs.associations[ *iter ];
    setRef = setRef.without( Set::SingletonSet( k ) );
  }
  as = Set();
}

void BasicBigraph::disconnectProtein(int k) {
  Set & as = proteinsToPSMs.associations[k];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++) {
    Set & setRef = PSMsToProteins.associations[ *iter ];
    setRef = setRef.without( Set::SingletonSet( k ) );
  }
  as = Set();
}

void BasicBigraph::add(GraphLayer & gl, StringTable & st, const string & item) {
  if ( st.lookup(item) == -1 ) {
    // if the string is not already known, then add a new node for it
    st.add(item);
    gl.associations.add( Set() );
    gl.weights.add( -1.0 );
    gl.sections.add(-1);
  }
}

void BasicBigraph::connect(const StringTable & pepNames, const string & pepName, const StringTable & proteinNames, const string & protName)
{
  int pepIndex = pepNames.lookup(pepName);
  int protIndex = proteinNames.lookup(protName);

  // performance note: currently O(n^2) worstcase. Later use a bitset
  // and then after the graph is read, pack it into a set. 

  PSMsToProteins.associations[ pepIndex ] |= Set::SingletonSet(protIndex);
  proteinsToPSMs.associations[ protIndex ] |= Set::SingletonSet(pepIndex);
}

void BasicBigraph::printProteinWeights() const
{
  const Array<string> & protNames = proteinsToPSMs.names;
  Array<double> sorted = proteinsToPSMs.weights;
  Array<int> indices = sorted.sort();
  
  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      cout << sorted[k] << Array<string>( 1, protNames[ indices[k] ] ) << endl;
    }
}

void BasicBigraph::traceConnected(int index, GraphLayer & gl, int sectionNumber) {
  if ( gl.sections[index] == sectionNumber )
    return;
  
  // if it has not already been marked by this section, do so
  gl.sections[index] = sectionNumber;
  // add mark to set of marks for this node, set can only be larger than 1 for 
  // peptides below PeptideThreshold
  gl.sectionMarks[index] |= Set::SingletonSet(sectionNumber);
  
  // do not follow edges with PSM probability below PeptideThreshold
  if (&gl == &PSMsToProteins && gl.weights[index] <= PeptideThreshold)
    return;

  const Set& as = gl.associations[index];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++) {
    if ( &gl == &proteinsToPSMs )
      traceConnected(*iter, PSMsToProteins, sectionNumber);
    else
      traceConnected(*iter, proteinsToPSMs, sectionNumber);
  }
}

int BasicBigraph::markSectionPartitions() {
  // returns the number of sections that are found
  proteinsToPSMs.sectionMarks = Array<Set>(proteinsToPSMs.size());
  PSMsToProteins.sectionMarks = Array<Set>(PSMsToProteins.size());

  PSMsToProteins.sections = Array<int>(PSMsToProteins.size(), -1);
  proteinsToPSMs.sections = Array<int>(proteinsToPSMs.size(), -1);
  
  // MT: make sure proteins with equal peptide evidence end up in the same section
  Array<Set> groups = ReplicateIndexer<Set>::replicates(Set::sumSetElements, proteinsToPSMs.associations );
  std::map<int, int> mapToFirstInGroup;
  
  for (int k = 0; k < groups.size(); k++) {
    Set reps = groups[k].without( Set::SingletonSet(groups[k][0]) );
    for (Set::Iterator iter = reps.begin(); iter != reps.end(); iter++) {
      mapToFirstInGroup[*iter] = groups[k][0];
    }
  }
  
  int section = 0;
  for (int k = 0; k < proteinsToPSMs.size(); k++) {
    if (mapToFirstInGroup.find(k) != mapToFirstInGroup.end()) {
      int sectionNumber = proteinsToPSMs.sections[mapToFirstInGroup[k]];
      proteinsToPSMs.sections[k] = sectionNumber;
      // add mark to set of marks for this node, set can only be larger than 1 for 
      // peptides below PeptideThreshold
      proteinsToPSMs.sectionMarks[k] |= Set::SingletonSet(sectionNumber);
    } else if (proteinsToPSMs[k].section == -1) {
      traceConnected(k, proteinsToPSMs, section);
      section++;
    }
  }

  return section;
}

Array<BasicBigraph> BasicBigraph::partitionSections() {
  int numSections = markSectionPartitions();

  Array<Set> proteinSubsets(numSections), PSMSubsets(numSections);
  for (int k = 0; k < proteinsToPSMs.size(); k++) {
    proteinSubsets[ proteinsToPSMs[k].section ].add(k);
  }

  for (int k = 0; k < PSMsToProteins.size(); k++) {
    PSMSubsets[ PSMsToProteins[k].section ].add(k);
  }

  // now reindex them to their proper sets

  Array<BasicBigraph> result;
  for (int k = 0; k < numSections; k++) {
    result.add( buildSubgraph( proteinSubsets[k], PSMSubsets[k] ) );
  }

  return result;
}

/* this function is currently not in use */
void BasicBigraph::pseudoCountPSMs() {
  for (int k = 0; k < PSMsToProteins.size(); k++) {
    if (PSMsToProteins.weights[k] < PeptideThreshold) {
      PSMsToProteins.weights[k] = PeptideThreshold;
    }
  }  
}

void BasicBigraph::setPsmThreshold(double __psm_threshold)
{
  PsmThreshold = __psm_threshold;
}

void BasicBigraph::setPeptideThreshold(double __peptide_threshold)
{
  PeptideThreshold = __peptide_threshold;
}

void BasicBigraph::setProteinThreshold(double __protein_threshold)
{
  ProteinThreshold = __protein_threshold;
}

double BasicBigraph::getPeptideThreshold()
{
  return PeptideThreshold;
}

double BasicBigraph::getProteinThreshold()
{
  return ProteinThreshold;
}

double BasicBigraph::getPsmThreshold()
{
  return PsmThreshold; 
}
