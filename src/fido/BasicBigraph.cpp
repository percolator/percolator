// Written by Oliver Serang 2009
// see license for more information

#include "BasicBigraph.h"

BasicBigraph::BasicBigraph() :
  ProteinIdentifier()
{
  //
}

// Mattia Tomasoni
inline string getRidOfUnprintablesAndUnicode(string inpString) {
  string outputs = "";
  for (int jj = 0; jj < inpString.size(); jj++) {
    char ch = inpString[jj];
    if (((int)ch) >= 32 && ((int)ch) <= 128) {
      outputs += ch;
    }
  }
  return outputs;
}

//Mattia Tomasoni
#include "Globals.h"
/**
 * read graph protein-psm graph from percolator's Scores class intead of file
 */
void BasicBigraph::read(Scores* fullset, bool debug){
	string pepName, protName;
	double value =  -10;
	int pepIndex = -1;
	StringTable PSMNames, proteinNames;

	ofstream of;
	if(debug) {
	  string s = string(TEMP_DIR) + "1percolator_psm_graph_file.txt";
	  of.open(s.c_str());
	}

	vector<ScoreHolder>::iterator psm = fullset->begin();
	for (; psm!= fullset->end(); ++psm) {

	  // e peptide_string
	  pepName = psm->pPSM->peptide;
	  if ( pepName[1] == '.' ) {
	    // trim off the cleavage events
	    pepName = pepName.substr(2, pepName.size()-4 );
	  }
	  if ( PSMNames.lookup(pepName) == -1 ){
	    //cout << "Adding e " << pepName << endl;
	    add(PSMsToProteins, PSMNames, pepName);
	  }
	  pepIndex = PSMNames.lookup(pepName);
	  if(debug) of << "e " << pepName << endl;

	  // r proteins
	  set<string>::const_iterator pid = psm->pPSM->proteinIds.begin();
	  for (; pid!= psm->pPSM->proteinIds.end(); ++pid) {
	    protName = getRidOfUnprintablesAndUnicode(*pid);
	    if ( proteinNames.lookup(protName) == -1 ){
	      //cout << "Adding r " << pepName << endl;
	      add(proteinsToPSMs, proteinNames, protName);
	    }
	    connect(PSMNames, pepName, proteinNames, protName);
	    if(debug) of << "r " << protName << endl;
	  }

	  // p probability of the peptide match to the spectrum
	  value = 1- psm->pPSM->pep;
	  PSMsToProteins.weights[ pepIndex ] = max(PSMsToProteins.weights[pepIndex], value);
	  if(debug) of << "p " << value << endl;
	}

	if(debug) of.close();

	PSMsToProteins.names = PSMNames.getItemsByNumber();
	proteinsToPSMs.names = proteinNames.getItemsByNumber();
}

void BasicBigraph::read(istream & is)
{
	char instr;
	string pepName, protName;
	double value =  -10;

	int pepIndex = -1;

	int state = 'e';

	StringTable PSMNames, proteinNames;

	while ( is >> instr )
	{
		//      cout << "Processing line: " << instr << endl;
		if ( instr == 'e' && (state == 'e' || state == 'p') )
		{
			is >> pepName;

			if ( pepName[1] == '.' )
			{
				// trim off the cleavage events
				pepName = pepName.substr(2, pepName.size()-4 );
			}

			if ( PSMNames.lookup(pepName) == -1 )
			{
				//	      cout << "Adding e " << pepName << endl;
				add(PSMsToProteins, PSMNames, pepName);
			}

			pepIndex = PSMNames.lookup(pepName);

			if ( state == 'p' )
			{
				cerr << "Warning: no peptide score for peptide entry " << pepName << ", using last score (" << value << ")" << endl;

				// using -10 since -1 may be used by peptide prophet
				if ( value == -10 )
				{
					cerr << "Error: No previous peptide entry to use" << endl;
					throw FormatException();
				}
				PSMsToProteins.weights[ pepIndex ] = max(value, PSMsToProteins.weights[ pepIndex ]);
			}

			state = 'r';

			//	  cout << "Read peptide: " << pepName << endl;
		}
		else if ( instr == 'r' && ( state == 'r' || state == 'p' ) )
		{
			is >> protName;

			if ( proteinNames.lookup(protName) == -1 )
				add(proteinsToPSMs, proteinNames, protName);

			// note: this is a hack to only get the nodes connected the
			// first time the peptide is seen; (it's initial weight is
			// -1.0 and will always be changed to >= 0.0 when it is
			// first finished being seen. Might want to make this better
			// later.

			//	  if ( PSMsToProteins.weights[ pepIndex ] == -1.0 )
			connect(PSMNames, pepName, proteinNames, protName);

			//	  cout << "Read protein: " << protName << endl;
			state = 'p';
		}
		else if ( instr == 'p' && state == 'p' )
		{
			is >> value;
			//	  cout << "Read value " << value << endl;

			// this option scores peptides using only the best match
			PSMsToProteins.weights[ pepIndex ] = max(PSMsToProteins.weights[pepIndex], value);

			//	  cout << "Gets value " << PSMsToProteins.weights[ pepIndex ] << endl;

			// this option scores peptides using all matches (and will
			// be at least as good as the best match)
			/***
	  if ( PSMsToProteins.weights[ pepIndex ] == -1 )
	    {
	      PSMsToProteins.weights[ pepIndex ] = value;
	    }
	  else
	    {
	      PSMsToProteins.weights[ pepIndex ] = 1 - ( 1-PSMsToProteins.weights[ pepIndex] )*(1 - value);
	    }
			 ***/

			state = 'e';
		}
		else if ( instr == '#' )
		{
			// comment line

			string garbage;
			getline(is, garbage);
		}
		else
		{
			throw FormatException();
		}
	}

	PSMsToProteins.names = PSMNames.getItemsByNumber();
	proteinsToPSMs.names = proteinNames.getItemsByNumber();

	//  printGraphStats();

	//  exit(0);
}

void BasicBigraph::printGraphStats() const
{
  cout << "There are \t" << PSMsToProteins.size() << " PSMs" << endl;
  cout << "      and \t" << proteinsToPSMs.size() << " proteins" << endl;

  int edgeCount = 0;
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      edgeCount += PSMsToProteins.associations[k].size();
    }

  cout << "      and \t" << edgeCount << " edges" << endl;
}

void BasicBigraph::readFromMCMC(istream & graph, istream & pepProph)
{
  string pepName, protName;
  double garbage;
  double value;

  int pepIndex, protIndex;

  StringTable PSMNames, proteinNames;

  while ( graph >> pepName >> protName >> garbage )
    {
      //      cerr << "\t" << pepName << " <-- " << protName << endl;

      //      cout << "\t\tPepstuffs" << endl;
      if ( PSMNames.lookup(pepName) == -1 )
	add(PSMsToProteins, PSMNames, pepName);
	  
      pepIndex = PSMNames.lookup(pepName);

      //      cout << "\t\tProtstuffs" << endl;
      if ( proteinNames.lookup(protName) == -1 )
	add(proteinsToPSMs, proteinNames, protName);

      protIndex = proteinNames.lookup(protName);

      //      cout << "\t\tConnecting" << endl;
      if ( PSMsToProteins.associations[ pepIndex ].find(protIndex) == -1 )
	connect(PSMNames, pepName, proteinNames, protName);
    }

  /***
  // set all peptide prophet scores to 0.0
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      int len = PSMNames.getItemsByNumber()[k].length();
      
      if ( len >= 10 && len <= 28 )
      	PSMsToProteins.weights[k] = 0.0;
    }
  ***/

  while ( pepProph >> pepName >> value )
    {
      //      cout << "Peptide " << pepName << " " << value << endl;

      pepIndex = PSMNames.lookup(pepName);

      if ( pepIndex != -1 )
	PSMsToProteins.weights[ pepIndex ] = max( PSMsToProteins.weights[ pepIndex ] , value );
    }

  PSMsToProteins.names = PSMNames.getItemsByNumber();
  proteinsToPSMs.names = proteinNames.getItemsByNumber();

  //  print();

  printGraphStats();
}

void BasicBigraph::saveSeveredProteins()
{
  severedProteins = Array<string>();
  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( proteinsToPSMs.associations[k].size() == 0 )
	{
	  // this protein will be excised out by reindex
	  //	  cout << "Protein with no assocs " << proteinsToPSMs.names[k] << endl;
	  severedProteins.add( proteinsToPSMs.names[k] );
	}
    }

  //  cout << "Severed protein list: " << severedProteins << endl;
}

void BasicBigraph::prune()
{
  //  cout << "Traditional prune and reindex" << endl;
  //  cout << "Pruning with threshold " << PeptideThreshold << endl;

  //  displayDotty("Before");
  
  removePoorPSMs();
  removePoorProteins();
  
  saveSeveredProteins();

  reindex();

  //  floorLowPSMs();
  markSectionPartitions();

  //  cout << "Cloning low-scoring PSMs" << endl;
  cloneMultipleMarkedPSMs();

  //  cout << "Number clones is (after cloning) " << numberClones << endl;

  //  cout << "Reindexing" << endl;
  reindex();

  //  displayDotty("After");

  //  cout << "Number clones is (after reindex) " << numberClones << endl;
}

void BasicBigraph::floorLowPSMs()
{
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      if ( PSMsToProteins.weights[k] <= PeptideThreshold )
	{
	  PSMsToProteins.weights[k] = 0.0;
	}
    }
}

void BasicBigraph::cloneMultipleMarkedPSMs()
{
  numberClones = 0;

  // compute this once, since it will change as you add clones
  int N = PSMsToProteins.size();

  for (int k=0; k<N; k++)
    {
      const GraphNode & psm = PSMsToProteins[k];

      // the way the marking procedure works, it will only multiple
      // mark PSMs with a score of 0.0
      if ( psm.sectionMark.size() > 1 )
	{
	  clonePSM( k );
	}
    }
}

void BasicBigraph::clonePSM(int pepIndex)
{
  //  cout << "\tcloning one PSM: " << PSMsToProteins.names[ pepIndex ] << endl;
  // add a new copy of this PSM for every section

  // first find the sections that this PSM associates with
  Set sections;
  const Set & s = PSMsToProteins.associations[pepIndex];
  int k, j;
  for (k=0; k<s.size(); k++)
    {
      int sect = proteinsToPSMs[ s[k] ].section;
      sections |= Set::SingletonSet(sect);
    }

  // index the proteins by the sections they belong to
  Array<Set> associatedProteinsBySection(sections.size());
  for (k=0; k<s.size(); k++)
    {
      int sect = proteinsToPSMs.sections[ s[k] ];
      int ind = sections.find(sect);
      associatedProteinsBySection[ ind ] |= Set::SingletonSet(s[k]);
    }

  // add a clone for each of the sections-- include the first one,
  // since it will be easier to build them all than to have a special
  // case. 

  for (k=0; k<associatedProteinsBySection.size(); k++)
    {
      int sect = sections[k];
      ostringstream ost;
      ost << PSMsToProteins[ pepIndex ].name << "_clone_" << sect;
      
      PSMsToProteins.names.add( ost.str() );
      PSMsToProteins.associations.add( associatedProteinsBySection[k] );
      PSMsToProteins.weights.add( PSMsToProteins[ pepIndex ].weight );
      PSMsToProteins.sections.add(sect);

      // add the association from this section's proteins to the new peptide
      for (int j=0; j<associatedProteinsBySection[k].size(); j++)
	proteinsToPSMs.associations[ associatedProteinsBySection[k][j] ] |= Set::SingletonSet( PSMsToProteins.size()-1 );
    }

  // afterward, erase the original
  // (remove the associations from proteins)
  for (k=0; k<associatedProteinsBySection.size(); k++)
    {
      for (j=0; j<associatedProteinsBySection[k].size(); j++)
	{
	  int prot = associatedProteinsBySection[k][j];

	  proteinsToPSMs[prot].association = proteinsToPSMs[prot].association.without( Set::SingletonSet(pepIndex) );
	}
    }

  // (remove the associations to proteins)
  PSMsToProteins[ pepIndex ].association = Set();

  // count the number of clones added (not including the original)

  //  cout << "Adding " << numberClones << " clones" << endl;
  numberClones += sections.size()-1;
}

void BasicBigraph::reindex()
{
  int k;

  // note: later check to see if you can do suchThat in Array using a
  // pointer to a member function
  Set connectedPSMs;
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      if ( ! PSMsToProteins.associations[k].isEmpty() )
	connectedPSMs.add(k);
    }

  Set connectedProteins;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( ! proteinsToPSMs.associations[k].isEmpty() )
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

BasicBigraph BasicBigraph::buildSubgraph(const Set & connectedProteins, const Set & connectedPSMs)
{
  int k;
  BasicBigraph result;

  result.PSMsToProteins.names = PSMsToProteins.names[ connectedPSMs ];
  result.PSMsToProteins.associations = PSMsToProteins.associations[ connectedPSMs ];
  result.PSMsToProteins.weights = PSMsToProteins.weights[ connectedPSMs ];
  result.PSMsToProteins.sections = PSMsToProteins.sections[ connectedPSMs ];

  for (k=0; k<result.PSMsToProteins.associations.size(); k++)
    {
      result.PSMsToProteins.associations[k] = result.PSMsToProteins.associations[k].reindexToFind( connectedProteins );
    }

  result.proteinsToPSMs.names = proteinsToPSMs.names[ connectedProteins ];
  result.proteinsToPSMs.associations = proteinsToPSMs.associations[ connectedProteins ];
  result.proteinsToPSMs.weights = proteinsToPSMs.weights[ connectedProteins ];
  result.proteinsToPSMs.sections = proteinsToPSMs.sections[ connectedProteins ];

  for (k=0; k<result.proteinsToPSMs.associations.size(); k++)
    {
      result.proteinsToPSMs.associations[k] = result.proteinsToPSMs.associations[k].reindexToFind( connectedPSMs );
    }

  return result;
}

void BasicBigraph::removePoorPSMs()
{
  int k;
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      //      if ( PSMsToProteins.weights[k] < BasicBigraph::PeptideThreshold )

      // note: asms hack
      //      if ( PSMsToProteins.weights[k] < 0.2 )
      if ( PSMsToProteins.weights[k] < 0.0 )
	{
	  //	  cout << "\tRemoving bad PSM " << PSMsToProteins.names[k] << " which had weight " << PSMsToProteins.weights[k] << endl;
	  disconnectPSM(k);
	}
    }
}

void BasicBigraph::removePoorProteins()
{
  //  cout << "Pruning with ProteinThreshold = " << ProteinThreshold << endl;
  int k;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( Vector(PSMsToProteins.weights[ proteinsToPSMs.associations[k] ]).max() < ProteinThreshold )
	{
	  //	  cout << "\tDisconnecting protein " << proteinsToPSMs.names[k] << endl;
	  //	  cout << "\t" << PSMsToProteins.weights[ proteinsToPSMs.associations[k] ] << endl;
	  disconnectProtein(k);
	}
    }
}

void BasicBigraph::disconnectPSM(int k)
{
  Set & as = PSMsToProteins.associations[k];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++)
    {
      Set & setRef = proteinsToPSMs.associations[ *iter ];
      setRef = setRef.without( Set::SingletonSet( k ) );
    }
  as = Set();
}

void BasicBigraph::disconnectProtein(int k)
{
  Set & as = proteinsToPSMs.associations[k];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++)
    {
      Set & setRef = PSMsToProteins.associations[ *iter ];
      setRef = setRef.without( Set::SingletonSet( k ) );
    }
  as = Set();

  //  cout << "Just discon protein " << proteinsToPSMs.names[k] << endl;
  //  cout << "\tnow has associations.size(): " << proteinsToPSMs.associations[k].size() << endl;
}

void BasicBigraph::add(GraphLayer & gl, StringTable & st, const string & item)
{
  if ( st.lookup(item) == -1 )
    {
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

void BasicBigraph::traceConnected(int index, GraphLayer & gl, int sectionNumber)
{
  if ( gl.sections[index] == sectionNumber )
    return;

  gl.sections[index] = sectionNumber;

  // if it has not already been marked by this section, do so
  gl.sectionMarks[index] |= Set::SingletonSet(sectionNumber);

  if ( & gl == & PSMsToProteins && gl.weights[index] <= PeptideThreshold )
    {
      //      cerr << "Bounding trace through PSM " << PSMsToProteins.names[index] << " with score " << gl.weights[index] << endl;
      return;
    }

  const Set & as = gl.associations[index];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++)
    {
      if ( & gl == & proteinsToPSMs )
	traceConnected( *iter, PSMsToProteins, sectionNumber );
      else
	{
	  traceConnected( *iter, proteinsToPSMs, sectionNumber );
	}
    }
}

int BasicBigraph::markSectionPartitions()
{
  // returns the number of sections that are found

  //  cout << "Marking with threshold " << PeptideThreshold << endl;

  // clear & allocate the sectionMarks
  proteinsToPSMs.sectionMarks = Array<Set>(proteinsToPSMs.size());
  PSMsToProteins.sectionMarks = Array<Set>(PSMsToProteins.size());

  PSMsToProteins.sections = Array<int>(PSMsToProteins.size(), -1);
  proteinsToPSMs.sections = Array<int>(proteinsToPSMs.size(), -1);

  int section = 0;
  int k;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( proteinsToPSMs[k].section == -1 )
	{
	  traceConnected( k, proteinsToPSMs, section );
	  //	  cout << "Just started a new section: " << section << endl;
	  section++;
	}
    }
  
  //  cout << "There are " << section << " partitions" << endl;

  return section;
}

Array<BasicBigraph> BasicBigraph::partitionSections()
{
  int numSections = markSectionPartitions();

  Array<Set> proteinSubsets(numSections), PSMSubsets(numSections);
  int k;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      proteinSubsets[ proteinsToPSMs[k].section ].add(k);
    }

  for (k=0; k<PSMsToProteins.size(); k++)
    {
      PSMSubsets[ PSMsToProteins[k].section ].add(k);
    }

  // now reindex them to their proper sets

  Array<BasicBigraph> result;
  for (k=0; k<numSections; k++)
    {
      /***
      if ( proteinSubsets[k].size() > 10 )
	{
	  cout << "\tProt Section: " << proteinSubsets[k].size();
	}
      ***/

      result.add( buildSubgraph( proteinSubsets[k], PSMSubsets[k] ) );
    }

  return result;
}

void BasicBigraph::outputDotty(ofstream & fout, const string & name) const
{
  fout << "graph " << name << " {" << endl;
  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      const Set & s = proteinsToPSMs.associations[k];

      string protName = proteinsToPSMs.names[k];

      for (int i=0; i<s.size(); i++)
	{
	  string pepName =  PSMsToProteins.names[ s[i] ];

	  ostringstream ost;
	  ost << "\\n" << PSMsToProteins.weights[ s[i] ];
	    
	  pepName += " " + ost.str();

	  fout << "R" << k << "[label=\"" << protName << "\"]" << endl;
	  fout << "E" << s[i] << "[label=\"" << pepName << "\"]" << endl;
	  fout << "R" << k << " -- " << "E" << s[i] << ";" << endl;
	}

    }
  fout << "}" << endl;
}

void BasicBigraph::displayDotty(const string & name) const
{
  string graphName = "/tmp/displayGraph_" + name + ".dot";
    
  ofstream fout(graphName.c_str());
  outputDotty(fout, name);

  string cmd = "dotty " + graphName + " ";
  if(system( cmd.c_str() ) != 0)
  {
    cerr << "Error: doing system call : " << cmd << endl;
  }

  //    getchar();
}
