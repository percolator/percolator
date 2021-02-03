#ifndef _BasicBigraph_H
#define _BasicBigraph_H

#include <fstream>
#include "Scores.h"
#include "StringTable.h"
#include "Array.h"
#include "Vector.h"
#include "ReplicateIndexer.h"

/*
* GraphNode represents a PSM or protein with associated probability and
* keeps a list of nodes it is connected to in the bigraph
*
*/
struct GraphNode {
  const string & name;
  Set & association;
  double & weight;
  int & section;
  Set & sectionMark;

  GraphNode(const string & n, Set & a, double & w, int & s, Set & sM):
    name(n), association(a), weight(w), section(s), sectionMark(sM) { }
    
  ~GraphNode() { }
};

/*
* GraphLayer represents a collection of PSMs or proteins
*
*/
struct GraphLayer {
  Array<string> names;
  Array<Set> associations;
  Array<double> weights;
  Array<int> sections;
  Array<Set> sectionMarks;

  ~GraphLayer() { }
    
  friend ostream & operator <<(ostream & os, const GraphLayer & gl) {
    os << "\t" << gl.associations << endl << "\t" << gl.weights << endl << "\t" << gl.names << endl;
    return os;
  }
  
  GraphNode operator [] (int k) {
    return GraphNode( names[k], associations[k], weights[k], sections[k], sectionMarks[k]);
  }
  
  int size() const {
    return static_cast<int>(associations.size());
  }
};

/*
* BasicBigraph represents a (sub)bigraph of PSMs and proteins
*
*/
class BasicBigraph {
 public:
  BasicBigraph();
  BasicBigraph(double __psmthreshold, double __peptidethreshold, 
	       double __proteinthreshold);
  virtual ~BasicBigraph();
  
  void read(Scores* fullset, bool multiple_labeled_peptides = false);
  void read(istream & is, bool multiple_labeled_peptides = false);
  void prune();
  void printGraph();
  void printProteinWeights() const;
  void printGraphStats() const;
  void print() const {
    cout << "PSM graph layer: " << endl;
    cout << PSMsToProteins << endl;
    cout << "Protein graph layer: " << endl;
    cout << proteinsToPSMs << endl;
    cout << proteinsToPSMs.names << endl << endl;
  }
  
  Array<BasicBigraph> partitionSections();
  void setPsmThreshold(double psm_threshold);
  double getPsmThreshold();
  void setPeptideThreshold(double peptide_threshold);
  double getPeptideThreshold();
  void setProteinThreshold(double protein_threshold);
  double getProteinThreshold();
  
protected:
  
  void add(GraphLayer & gl, StringTable & st, const string & item);
  void connect(const StringTable & PSMNames, const string & pepStr, 
	       const StringTable & proteinNames, const string & protStr);
  void disconnectProtein(int k);
  void disconnectPSM(int k);
  void pseudoCountPSMs();
  void floorLowPSMs();
  int markSectionPartitions();
  void removeDegeneratePSMs();
  void cloneDegeneratePSMs();
  void removePoorPSMs();
  void removePoorProteins();
  void reindex();
  void cloneMultipleMarkedPSMs();
  void clonePSM(int pepIndex);
  void saveSeveredProteins();
  
  BasicBigraph buildSubgraph(const Set & connectedProteins, const Set & connectedPSMs);
  void traceConnected(int index, GraphLayer & gl, int sectionNumber);

  double PsmThreshold;
  double PeptideThreshold;
  double ProteinThreshold;
  
public:
  //NOTE these should be protected and access them with functions
  int numberClones;
  Array<string> severedProteins;
  GraphLayer proteinsToPSMs, PSMsToProteins;
  
};

#endif

