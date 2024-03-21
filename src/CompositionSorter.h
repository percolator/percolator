#include <iostream>
#include <string>
#include <map>
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

class ScoreHolder ;
class Scores;

class CompositionSorter {
    public:
        int addPSMs(Scores& psms); 
        std::string generateCompositionSignature(const std::string& peptide);
        int sortScorePerPeptide();
        int inCompositionCompetition(Scores& winnerPeptides, unsigned int decoysPerTarget=1);
        int psmAndPeptide(Scores& scores, Scores& winnerPeptides, unsigned int decoysPerTarget=1);
    protected:
        std::unordered_map<std::string, std::map<std::string,std::vector<const ScoreHolder*>>, DJB2Hash> compositionToPeptidesToScore_;
};