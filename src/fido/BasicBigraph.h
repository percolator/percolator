// Written by Oliver Serang 2009
// see license for more information

#ifndef _BasicBigraph_H
#define _BasicBigraph_H

#include "ProteinIdentifier.h"
#include <fstream>

class BasicBigraph : public ProteinIdentifier
{
  // note: hack; should be protected later
public:

  struct GraphNode
  {
    const string & name;
    Set & association;
    double & weight;
    int & section;
    Set & sectionMark;

  GraphNode(const string & n, Set & a, double & w, int & s, Set & sM):
    name(n),
      association(a),
      weight(w),
      section(s),
      sectionMark(sM)
    {
    }
  };

  struct GraphLayer
  {
    Array<string> names;
    Array<Set> associations;
    Array<double> weights;
    Array<int> sections;
    Array<Set> sectionMarks;

    friend ostream & operator <<(ostream & os, const GraphLayer & gl)
    {
      os << "\t" << gl.associations << endl << "\t" << gl.weights << endl << "\t" << gl.names << endl;
      return os;
    }
    GraphNode operator [] (int k)
    {
      return GraphNode( names[k], associations[k], weights[k], sections[k], sectionMarks[k] );
    }
    int size() const
    {
      // all should be the same
      return associations.size();
    }
  };

  void read(istream & is);
  void readFromMCMC(istream & graph, istream & pepProph);
  void add(GraphLayer & gl, StringTable & st, const string & item);
  void connect(const StringTable & PSMNames, const string & pepStr, const StringTable & proteinNames, const string & protStr);

  void disconnectProtein(int k);
  void disconnectPSM(int k);
  void floorLowPSMs();
  int markSectionPartitions();
  void removePoorPSMs();
  void removePoorProteins();
  void reindex();
  void prune();

  void cloneMultipleMarkedPSMs();
  void clonePSM(int pepIndex);

  int numberClones;

  Array<string> severedProteins;

  void saveSeveredProteins();

public:
  void outputDotty(ofstream & fout, const string & name) const;
  void displayDotty(const string & name) const;

  GraphLayer proteinsToPSMs, PSMsToProteins;

  void printGraphStats() const;
  void print() const
  {
    cout << "PSM graph layer: " << endl;
    cout << PSMsToProteins << endl;
    cout << "Protein graph layer: " << endl;
    cout << proteinsToPSMs << endl;
    cout << proteinsToPSMs.names << endl << endl;
  }

  BasicBigraph();

  void printProteinWeights() const;

  Array<BasicBigraph> partitionSections();
  BasicBigraph buildSubgraph(const Set & connectedProteins, const Set & connectedPSMs);

  void traceConnected(int index, GraphLayer & gl, int sectionNumber);
};

#endif

