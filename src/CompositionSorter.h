#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

class ScoreHolder;
class Scores;

class CompositionSorter {
 public:
  // int addPSMs(Scores& psms, bool useTDC=true);
  int addPSMs(const Scores& psms, bool useTDC = false);
  // the residues and modifications of a peptide with their counts, e.g.
  // "A2C1K1[+57.02]1": in the order of their strings, each followed by
  // its count
  std::string generateCompositionSignature(std::string_view peptide);
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
  // composition signature -> number, in the order in which the signatures
  // first appear. The compositions are processed in the order of their
  // numbers; the map is only used for lookups, so its hash function and
  // iteration order do not affect the result.
  std::unordered_map<std::string, unsigned int> compositionIndex_;
  // the PSMs added, with the number of their composition and their peptide,
  // which views the PSM's own sequence; sortScorePerPeptide orders them by
  // composition, peptide and decreasing score
  struct CompositionPsm {
    unsigned int composition;
    std::string_view peptide;
    std::size_t order;  // position in which the PSM was added
    ScoreHolder* psm;
  };
  std::vector<CompositionPsm> compositionPsms_;
};