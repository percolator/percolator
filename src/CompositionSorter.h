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

class PeptidepairScore;
class ScoreBase;
class Scores;

class CompositionSorter {
    public:
        int addPSMs(std::vector<ScoreBase>& psms); 
        std::string generateCompositionSignature(const std::string& peptide);
        int fillPeptidepairs(Scores& scores); 
    protected:
        std::unordered_map<std::string, std::map<std::string,std::vector<ScoreBase&>>, DJB2Hash> targetCompositionToPeptides_;
        std::unordered_map<std::string, std::map<std::string,std::vector<ScoreBase&>>, DJB2Hash> decoyCompositionToPeptides_;
};