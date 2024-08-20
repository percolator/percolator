#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// based on the djb2 algorithm
struct DJB2Hash {
  size_t operator()(const std::string& s) const {
    unsigned long hash = 5381;
    for (char c : s) {
      hash = ((hash << 5) + hash) + c;
    }
    return hash;
  }
};

class ScoreHolder;
class Scores;

class CompositionSorter {
 public:
  // int addPSMs(Scores& psms, bool useTDC=true);
  int addPSMs(const Scores& psms, bool useTDC = false);
  std::string generateCompositionSignature(const std::string& peptide);
  int sortScorePerPeptide();
  int inCompositionCompetition(Scores& bestScoreHolders,
                               unsigned int decoysPerTarget = 1);
  int psmAndPeptide(const Scores& scores,
                    Scores& winnerPeptides,
                    unsigned int decoysPerTarget = 1);
  static void psmsOnly(const Scores& scores, Scores& winnerPeptides);
  static void retainRepresentatives(const Scores& psms,
                                    Scores& winnerPeptides,
                                    double selectionFDR,
                                    unsigned int decoysPerTarget,
                                    bool useCompositionMatch);

 protected:
  std::unordered_map<std::string,
                     std::map<std::string, std::vector<ScoreHolder*>>,
                     DJB2Hash>
      compositionToPeptidesToScore_;
};