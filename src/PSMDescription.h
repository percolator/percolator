/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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
#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "Enzyme.h"

/*
 * PSMDescription
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 *
 */
class PSMDescription {
   public:
    PSMDescription() = default;
    // the characters of peptide are not copied and must outlive the PSM
    explicit PSMDescription(std::string_view peptide) : peptide(peptide) {}

    // view of the protein names of a PSM, in input order
    class ProteinList {
     public:
        class const_iterator {
         public:
            explicit const_iterator(const std::string* const* p) : p_(p) {}
            const std::string& operator*() const { return **p_; }
            const std::string* operator->() const { return *p_; }
            const_iterator& operator++() { ++p_; return *this; }
            bool operator==(const const_iterator& other) const { return p_ == other.p_; }
            bool operator!=(const const_iterator& other) const { return p_ != other.p_; }
         private:
            const std::string* const* p_;
        };
        ProteinList(const std::string* const* first, size_t n) : first_(first), n_(n) {}
        const_iterator begin() const { return const_iterator(first_); }
        const_iterator end() const { return const_iterator(first_ + n_); }
        size_t size() const { return n_; }
        bool empty() const { return n_ == 0; }
        const std::string& operator[](size_t i) const { return *first_[i]; }
     private:
        const std::string* const* first_;
        size_t n_;
    };


    // TODO: move these static functions somewhere else
    static std::string removePTMs(const std::string& peptideSeq);
    static std::string removeFlanks(const std::string& peptideSeq) {
        return peptideSeq.substr(2, peptideSeq.size() - 4);
    }
    static bool ptrLess(PSMDescription* one, PSMDescription* other) {
        return *one < *other;
    }
    static bool ptrEqual(PSMDescription* one, PSMDescription* other) {
        return *one == *other;
    }

    std::string getPeptideSequence() const { return std::string(getPeptideSequenceView()); }
    // same as getPeptideSequence, without copying
    std::string_view getPeptideSequenceView() const {
        return peptide.substr(2, peptide.size() - 4);
    }
    std::string_view getFullPeptideSequence() const { return peptide; }
    std::string getFlankN() const { return std::string(peptide.substr(0, 1)); }
    std::string getFlankC() const {
        return std::string(peptide.substr(peptide.size() - 1, peptide.size()));
    }

    friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
    void printProteins(std::ostream& out) const;

    bool operator<(const PSMDescription& other) const {
        return (peptide < other.peptide) ||
               (peptide == other.peptide && getRetentionTime() < other.getRetentionTime());
    }

    bool operator==(const PSMDescription& other) const {
        return (peptide == other.peptide);
    }

    // the characters of id and pep_seq are not copied and must outlive the
    // PSM, e.g. by being stored in a PsmStore
    inline void setId(std::string_view id) { id_ = id; }
    inline std::string_view getId() const { return id_; }

    std::string_view getFullPeptide() const { return peptide; }
    void setPeptide(std::string_view pep_seq) { peptide = pep_seq; }

    ProteinList proteins() const { return ProteinList(proteins_, numProteins_); }
    // the list is not copied and must outlive the PSM
    void setProteins(const std::string* const* list, size_t n) {
        proteins_ = list;
        numProteins_ = static_cast<unsigned int>(n);
    }
    static inline void setProteinNameSeparator(const std::string sep) {
        proteinNameSeparator_ = sep;
    }
    static inline const std::string& getProteinNameSeparator() { return proteinNameSeparator_; }

    void setSpectrumFileName(std::string_view fileName) {
        // consecutive PSMs usually come from the same file
        size_t& lastIndex = lastSpectrumFileIndex_;
        if (lastIndex >= spectraFileNames_.size() || spectraFileNames_[lastIndex] != fileName) {
            auto specFilePos = std::find(spectraFileNames_.begin(), spectraFileNames_.end(), fileName);
            lastIndex = specFilePos - spectraFileNames_.begin();
            if (specFilePos == spectraFileNames_.end()) {
                spectraFileNames_.push_back(std::string(fileName));
            }
        }
        specFileNr = static_cast<unsigned int>(lastIndex);
    }
    // forgets the spectrum file names of an earlier run in the same process
    static void resetSpectrumFileNames() {
        spectraFileNames_.clear();
        lastSpectrumFileIndex_ = 0;
    }
    inline const std::string& getSpectrumFileName() {
        static const std::string noFileName;
        return hasSpectrumFileName() ? spectraFileNames_.at(specFileNr) : noFileName;
    }
    inline bool static hasSpectrumFileName() { return !spectraFileNames_.empty(); }

    inline void setRetentionTime(const double retentionTime) {
        retentionTime_ = retentionTime;
    }
    inline double getRetentionTime() const { return retentionTime_; }

    // line of the PSM in the input, used as the last tie-break when sorting
    inline void setLineNr(unsigned int lineNr) { lineNr_ = lineNr; }
    inline unsigned int getLineNr() const { return lineNr_; }

    double* features = NULL;  // owned by a FeatureMemoryPool instance, no need to delete
    double expMass = 0.0, calcMass = 0.0;
    double retentionTime_ = std::numeric_limits<double>::quiet_NaN();
    unsigned int scan = 0u;
    unsigned int specFileNr = 0u;

   protected:
    std::string_view id_;
    std::string_view peptide;
    const std::string* const* proteins_ = nullptr;
    unsigned int numProteins_ = 0u;
    unsigned int lineNr_ = 0u;  // fills the padding after numProteins_
    static std::string proteinNameSeparator_;
    static vector<std::string> spectraFileNames_;
    static size_t lastSpectrumFileIndex_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
    out << "Peptide: " << psm.peptide << endl;
    out << "Spectrum scan number: " << psm.scan << endl;
    out << endl;
    return out;
}

#endif /*PSMDESCRIPTION_H_*/
