#include <vector>
#include <algorithm> // For std::max_element

class ScoreBase {
public:
    virtual double getScore() const = 0;
};

class SpectrumGroup {
public:
    ScoreBase* target, decoy1, decoy2;
    ScoreBase* highestScoring;

    SpectrumGroup() : target(nullptr), decoy1(nullptr), decoy2(nullptr), highestScoring(nullptr) {}

    void addTarget(ScoreBase* t) {
        target = t;
        updateHighestScoring();
    }

    void addDecoy1(ScoreBase* d) {
        decoy1=d;
        updateHighestScoring();
    }

    void addDecoy2(ScoreBase* d) {
        decoy2=d;
        updateHighestScoring();
    }


private:
    void updateHighestScoring() {
        highestScoring = target; // Assume target is the highest at the start
        if (decoy1 && (decoy1->getScore() > highestScoring->getScore())) {
            highestScoring = decoy1;
        }
        if (decoy2 && (decoy2->getScore() > highestScoring->getScore())) {
            highestScoring = decoy2;
        }
    }
};
