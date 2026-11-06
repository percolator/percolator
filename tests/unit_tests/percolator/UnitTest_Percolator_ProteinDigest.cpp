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
/*
 * Unit tests for the FASTA reading, digestion and peptide index used by the
 * detection of fragment and duplicate proteins (src/picked_protein).
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Enzyme.h"
#include "ProteinDigest.h"
// the PercolatorCrux digestion that ProteinDigest replaces
#include "Peptide.h"
#include "PeptideConstraint.h"
#include "Protein.h"
#include "ProteinPeptideIterator.h"

using PickedProtein::DigestConstraint;
using PickedProtein::DigestType;
using PickedProtein::FastaProtein;
using PickedProtein::PeptideIndex;
using PickedProtein::PeptideSpan;

namespace {

typedef std::vector<std::pair<int, int> > Spans;  // (0-based start, length)

/*
 * Reference: the peptide selection of the PercolatorCrux
 * ProteinPeptideIterator (adapted from crux), which ProteinDigest replaces,
 * transcribed with the same arrays of cleavage positions, the same pair
 * selection with its skipping of C-terminal positions, and the same array
 * bounds. Only the enzyme rule is taken from Enzyme.h. The residue after the
 * C terminus is '\0', as in the null-terminated sequence of the original.
 */
Spans legacyDigest(const std::string& seq, const Enzyme& enzyme,
                   DigestType digest, int minLength, int maxLength,
                   int numMisCleavage) {
  const int len = static_cast<int>(seq.size());
  std::vector<int> cleavage(static_cast<size_t>(len) + 1);
  std::vector<int> nonCleavage(static_cast<size_t>(len) + 1);
  std::vector<int> all(static_cast<size_t>(len) + 1);
  std::vector<int> cumulative;  // cleavage positions before residue i, plus 1
  int numCleavage = 0, numNonCleavage = 0;
  cleavage[static_cast<size_t>(numCleavage++)] = 0;
  for (int i = 0; i < len; ++i) {
    cumulative.push_back(numCleavage);
    if (enzyme.isEnzymatic(seq[static_cast<size_t>(i)],
                           seq.c_str()[i + 1])) {
      cleavage[static_cast<size_t>(numCleavage++)] = i + 1;
    } else {
      nonCleavage[static_cast<size_t>(numNonCleavage++)] = i + 1;
    }
    all[static_cast<size_t>(i)] = i;
  }
  if (cleavage[static_cast<size_t>(numCleavage - 1)] != len) {
    cleavage[static_cast<size_t>(numCleavage++)] = len;
  }
  all[static_cast<size_t>(len)] = len;

  Spans out;
  auto select = [&](const int* nterm, int numNterm, const int* cterm,
                    int numCterm, int skip) {
    int previousStart = 0;
    for (int n = 0; n < numNterm; ++n) {
      int nextStart = previousStart;
      bool noNewStart = true;
      for (int c = previousStart; c < numCterm; ++c) {
        if (cumulative[static_cast<size_t>(cterm[c] - 1)] -
                cumulative[static_cast<size_t>(nterm[n])] >
            skip) {
          break;
        }
        if (cterm[c] <= nterm[n]) continue;
        const int length = cterm[c] - nterm[n];
        if (length < minLength) {
          continue;
        } else if (length > maxLength) {
          break;
        } else if (noNewStart) {
          nextStart = c;
          noNewStart = false;
        }
        out.push_back(std::make_pair(nterm[n], length));
      }
      previousStart = nextStart;
    }
  };
  switch (digest) {
    case DigestType::Full:
      select(cleavage.data(), numCleavage - 1, cleavage.data() + 1,
             numCleavage - 1, numMisCleavage);
      break;
    case DigestType::Partial:
      select(all.data(), len, cleavage.data() + 1, numCleavage - 1,
             numMisCleavage);
      select(cleavage.data(), numCleavage - 1, nonCleavage.data(),
             numNonCleavage - 1, numMisCleavage);
      break;
    case DigestType::NonSpecific:
      select(all.data(), len, all.data() + 1, len, 500);
      break;
  }
  return out;
}

Spans newDigest(const std::string& seq, const Enzyme& enzyme,
                DigestType digest, int minLength, int maxLength,
                int numMisCleavage) {
  std::vector<PeptideSpan> spans;
  PickedProtein::digestProtein(
      seq, enzyme, DigestConstraint{digest, minLength, maxLength,
                                    numMisCleavage},
      spans);
  Spans out;
  for (const PeptideSpan& s : spans) {
    out.push_back(std::make_pair(static_cast<int>(s.start),
                                 static_cast<int>(s.length)));
  }
  return out;
}

// Percolator is compiled with _GLIBCXX_PARALLEL, whose std::sort and vector
// comparison start a team of threads for each of the many small calls of
// these tests; they sort and compare sequentially instead. (libc++, as on
// macOS, ignores _GLIBCXX_PARALLEL and has no __gnu_parallel.)
Spans sorted(Spans s) {
#if defined(_GLIBCXX_PARALLEL) && defined(__GLIBCXX__)
  std::sort(s.begin(), s.end(), __gnu_parallel::sequential_tag());
#else
  std::sort(s.begin(), s.end());
#endif
  return s;
}

template <typename T>
bool same(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!(a[i] == b[i])) return false;
  }
  return true;
}

// The old class itself, for the enzymes it supports.
PercolatorCrux::ENZYME_T cruxEnzyme(Enzyme::EnzymeType type) {
  switch (type) {
    case Enzyme::NO_ENZYME: return PercolatorCrux::NO_ENZYME;
    case Enzyme::TRYPSIN: return PercolatorCrux::TRYPSIN;
    case Enzyme::TRYPSINP: return PercolatorCrux::TRYPSINP;
    case Enzyme::CHYMOTRYPSIN: return PercolatorCrux::CHYMOTRYPSIN;
    case Enzyme::ELASTASE: return PercolatorCrux::ELASTASE;
    case Enzyme::LYSN: return PercolatorCrux::LYSN;
    case Enzyme::LYSC: return PercolatorCrux::LYSC;
    case Enzyme::ARGC: return PercolatorCrux::ARGC;
    case Enzyme::ASPN: return PercolatorCrux::ASPN;
    case Enzyme::GLUC: return PercolatorCrux::GLUC;
    default: return PercolatorCrux::INVALID_ENZYME;
  }
}

Spans oldIteratorDigest(const std::string& seq, Enzyme::EnzymeType type,
                        DigestType digest, int minLength, int maxLength,
                        int numMisCleavage) {
  PercolatorCrux::DIGEST_T cruxDigest =
      digest == DigestType::Full      ? PercolatorCrux::FULL_DIGEST
      : digest == DigestType::Partial ? PercolatorCrux::PARTIAL_DIGEST
                                      : PercolatorCrux::NON_SPECIFIC_DIGEST;
  PercolatorCrux::Protein protein;
  protein.setId("p");
  protein.setSequence(seq.c_str());
  protein.setLength(static_cast<unsigned int>(seq.size()));
  PercolatorCrux::PeptideConstraint constraint(
      cruxEnzyme(type), cruxDigest, minLength, maxLength, numMisCleavage);
  PercolatorCrux::ProteinPeptideIterator it(&protein, &constraint);
  Spans out;
  while (it.hasNext()) {
    PercolatorCrux::Peptide* peptide = it.next();
    out.push_back(std::make_pair(
        static_cast<int>(peptide->getSequencePointer() -
                         protein.getSequencePointer()),
        static_cast<int>(peptide->getLength())));
    PercolatorCrux::Peptide::free(peptide);
  }
  return out;
}

const std::vector<std::string>& realProteins() {
  // canonical sequences from UniProtKB/Swiss-Prot (release 2026_04)
  static const std::vector<std::string> proteins = {
    // sp|P0CG48|UBC_HUMAN, 685 residues
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKE"
    "STLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLED"
    "GRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQ"
    "QRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKA"
    "KIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLE"
    "VEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIF"
    "VKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLH"
    "LVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTL"
    "SDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLI"
    "FAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQD"
    "KEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGV",
    // sp|P62805|H4_HUMAN, 103 residues
    "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLE"
    "NVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG",
    // sp|P04637|P53_HUMAN, 393 residues
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAP"
    "RMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSP"
    "ALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQ"
    "HLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIIT"
    "LEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKK"
    "KPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLM"
    "FKTEGPDSD",
    // sp|P69905|HBA_HUMAN, 142 residues
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVA"
    "DALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDK"
    "FLASVSTVLTSKYR",
    // sp|P02768|ALBU_HUMAN, 609 residues
    "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHV"
    "KLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQ"
    "HKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTEC"
    "CQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVS"
    "KLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVEND"
    "EMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKC"
    "CAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLV"
    "EVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCF"
    "SALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFA"
    "AFVEKCCKADDKETCFAEEGKKLVAASQAALGL",
    // sp|P00924|ENO1_YEAST, 437 residues
    "MAVSKVYARSVYDSRGNPTVEVELTTEKGVFRSIVPSGASTGVHEALEMRDGDKSKWMGKGVLH"
    "AVKNVNDVIAPAFVKANIDVKDQKAVDDFLISLDGTANKSKLGANAILGVSLAASRAAAAEKNV"
    "PLYKHLADLSKSKTSPYVLPVPFLNVLNGGSHAGGALALQEFMIAPTGAKTFAEALRIGSEVYH"
    "NLKSLTKKRYGASAGNVGDEGGVAPNIQTAEEALDLIVDAIKAAGHDGKIKIGLDCASSEFFKD"
    "GKYDLDFKNPNSDKSKWLTGPQLADLYHSLMKRYPIVSIEDPFAEDDWEAWSHFFKTAGIQIVA"
    "DDLTVTNPKRIATAIEKKAADALLLKVNQIGTLSESIKAAQDSFAAGWGVMVSHRSGETEDTFI"
    "ADLVVGLRTGQIKTGAPARSERLAKLNQLLRIEEELGDNAVFAGENFHHGDKL",
    // sp|P00330|ADH1_YEAST, 348 residues
    "MSIPETQKGVIFYESHGKLEYKDIPVPKPKANELLINVKYSGVCHTDLHAWHGDWPLPVKLPLV"
    "GGHEGAGVVVGMGENVKGWKIGDYAGIKWLNGSCMACEYCELGNESNCPHADLSGYTHDGSFQQ"
    "YATADAVQAAHIPQGTDLAQVAPILCAGITVYKALKSANLMAGHWVAISGAAGGLGSLAVQYAK"
    "AMGYRVLGIDGGEGKEELFRSIGGEVFIDFTKEKDIVGAVLKATDGGAHGVINVSVSEAAIEAS"
    "TRYVRANGTTVLVGMPAGAKCCSDVFNQVVKSISIVGSYVGNRADTREALDFFARGLVKSPIKV"
    "VGLSTLPEIYEKMEKGQIVGRYVVDTSK",
    // sp|P02994|EF1A_YEAST, 458 residues
    "MGKEKSHINVVVIGHVDSGKSTTTGHLIYKCGGIDKRTIEKFEKEAAELGKGSFKYAWVLDKLK"
    "AERERGITIDIALWKFETPKYQVTVIDAPGHRDFIKNMITGTSQADCAILIIAGGVGEFEAGIS"
    "KDGQTREHALLAFTLGVRQLIVAVNKMDSVKWDESRFQEIVKETSNFIKKVGYNPKTVPFVPIS"
    "GWNGDNMIEATTNAPWYKGWEKETKAGVVKGKTLLEAIDAIEQPSRPTDKPLRLPLQDVYKIGG"
    "IGTVPVGRVETGVIKPGMVVTFAPAGVTTEVKSVEMHHEQLEQGVPGDNVGFNVKNVSVKEIRR"
    "GNVCGDAKNDPPKGCASFNATVIVLNHPGQISAGYSPVLDCHTAHIACRFDELLEKNDRRSGKK"
    "LEDHPKFLKSGDAALVKFVPSKPMCVEAFSEYPPLGRFAVRDMRQTVAVGVIKSVDKTEKAAKV"
    "TKAAQKAAKK",
  };
  return proteins;
}

// Random proteins over an alphabet rich in the residues of the cleavage
// rules, of lengths 0 to 300, and some fixed corner cases.
std::vector<std::string> randomProteins(unsigned seed, size_t count) {
  std::vector<std::string> proteins = {"", "K", "M", "P", "KP", "MK", "KK",
                                       "MKR", "RPK", "DDEE", "PPPP",
                                       "MKKKKKKKKKKKKKKKKK", "KRKRKRKRKRKRK",
                                       "AAAAAAAAAAAAAAK", "MAAAAKAAAAAAAR"};
  const std::string alphabet = "KRPKRPDEMFLWYAGIVSTNQCHX";
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> letter(0, alphabet.size() - 1);
  std::uniform_int_distribution<size_t> length(0, 300);
  while (proteins.size() < count) {
    std::string s(length(rng), 'A');
    for (char& c : s) c = alphabet[letter(rng)];
    if (!s.empty() && rng() % 3 == 0) s[0] = 'M';
    proteins.push_back(s);
  }
  return proteins;
}

const Enzyme::EnzymeType kAllEnzymes[] = {
    Enzyme::NO_ENZYME, Enzyme::TRYPSIN, Enzyme::TRYPSINP,
    Enzyme::CHYMOTRYPSIN, Enzyme::THERMOLYSIN, Enzyme::PROTEINASEK,
    Enzyme::PEPSIN, Enzyme::ELASTASE, Enzyme::LYSN, Enzyme::LYSC,
    Enzyme::ARGC, Enzyme::ASPN, Enzyme::GLUC};
const DigestType kAllDigests[] = {DigestType::Full, DigestType::Partial,
                                  DigestType::NonSpecific};

struct Limits {
  int minLength, maxLength, misCleavages;
};
const Limits kLimits[] = {{1, 1, 0},  {1, 4, 0},  {2, 6, 1},  {6, 50, 2},
                          {7, 30, 0}, {7, 30, 3}, {4, 60, 8}, {7, 255, 1},
                          {0, 20, 2}, {10, 5, 2}, {3, 12, -1}};

// Non-specific digests give about length x maxLength peptides and depend on
// the enzyme only through the limit of 500 missed cleavages, so they are
// tested with two enzymes and short peptides.
bool skipCase(Enzyme::EnzymeType type, DigestType digest, const Limits& lim) {
  return digest == DigestType::NonSpecific &&
         (lim.maxLength > 30 ||
          (type != Enzyme::NO_ENZYME && type != Enzyme::TRYPSIN));
}

}  // namespace

// The digester gives the same peptides as the transcribed PercolatorCrux
// selection for every enzyme, digestion type and length/missed-cleavage
// limit, on real and random proteins.
TEST(ProteinDigestTest, SameAsLegacyDigestion) {
  std::vector<std::string> proteins = realProteins();
  std::vector<std::string> random = randomProteins(11, 200);
  proteins.insert(proteins.end(), random.begin(), random.end());
  size_t numSpans = 0;
  for (Enzyme::EnzymeType type : kAllEnzymes) {
    std::unique_ptr<Enzyme> enzyme = Enzyme::createEnzyme(type);
    for (DigestType digest : kAllDigests) {
      for (const Limits& lim : kLimits) {
        if (skipCase(type, digest, lim)) continue;
        for (const std::string& seq : proteins) {
          Spans expected = sorted(legacyDigest(seq, *enzyme, digest,
              lim.minLength, lim.maxLength, lim.misCleavages));
          Spans got = sorted(newDigest(seq, *enzyme, digest, lim.minLength,
              lim.maxLength, lim.misCleavages));
          ASSERT_TRUE(same(expected, got))
              << enzyme->getStringEnzyme() << " digest "
              << static_cast<int>(digest) << " limits " << lim.minLength
              << "," << lim.maxLength << "," << lim.misCleavages
              << " protein " << seq;
          numSpans += got.size();
        }
      }
    }
  }
  EXPECT_GT(numSpans, 1000000u);
}

// The transcription above reproduces the old ProteinPeptideIterator, peptide
// for peptide and in the same order, for all enzymes the iterator supports.
TEST(ProteinDigestTest, LegacyDigestionMatchesProteinPeptideIterator) {
  std::vector<std::string> proteins = realProteins();
  std::vector<std::string> random = randomProteins(12, 200);
  proteins.insert(proteins.end(), random.begin(), random.end());
  for (Enzyme::EnzymeType type : kAllEnzymes) {
    if (cruxEnzyme(type) == PercolatorCrux::INVALID_ENZYME) continue;
    std::unique_ptr<Enzyme> enzyme = Enzyme::createEnzyme(type);
    for (DigestType digest : kAllDigests) {
      for (const Limits& lim : kLimits) {
        if (skipCase(type, digest, lim)) continue;
        for (const std::string& seq : proteins) {
          Spans legacy = legacyDigest(seq, *enzyme, digest, lim.minLength,
                                      lim.maxLength, lim.misCleavages);
          Spans old = oldIteratorDigest(seq, type, digest, lim.minLength,
                                        lim.maxLength, lim.misCleavages);
          ASSERT_TRUE(same(legacy, old))
              << enzyme->getStringEnzyme() << " digest "
              << static_cast<int>(digest) << " limits " << lim.minLength
              << "," << lim.maxLength << "," << lim.misCleavages
              << " protein " << seq;
        }
      }
    }
  }
}

// Spot checks of the cleavage rules on a protein that ends in a cleavage
// residue: for trypsin, peptides end after K or R unless P follows.
TEST(ProteinDigestTest, TrypticPeptides) {
  std::unique_ptr<Enzyme> trypsin = Enzyme::createEnzyme(Enzyme::TRYPSIN);
  const std::string seq = "MAKPGRAAKDDR";
  Spans full = sorted(newDigest(seq, *trypsin, DigestType::Full, 1, 50, 0));
  // MAKPGR, AAK, DDR
  EXPECT_TRUE(same(Spans({{0, 6}, {6, 3}, {9, 3}}), full));
  Spans missed = sorted(newDigest(seq, *trypsin, DigestType::Full, 1, 50, 1));
  EXPECT_TRUE(same(Spans({{0, 6}, {0, 9}, {6, 3}, {6, 6}, {9, 3}}),
                   missed));
  // a limit of 0 allows no peptides
  EXPECT_TRUE(newDigest(seq, *trypsin, DigestType::Full, 1, 0, 0).empty());
}

namespace {

// The protein lists of each peptide as PickedProteinCaller built them before:
// a map from the sequence to the protein indices, with repeats.
std::map<std::string, std::vector<size_t> > mapIndex(
    const std::vector<FastaProtein>& proteins, const Enzyme& enzyme,
    const DigestConstraint& c, std::vector<std::vector<std::string> >& peps,
    std::vector<size_t>& numPeptides) {
  std::map<std::string, std::vector<size_t> > index;
  peps.assign(proteins.size(), std::vector<std::string>());
  numPeptides.assign(proteins.size(), 0);
  for (size_t p = 0; p < proteins.size(); ++p) {
    const std::string& s = proteins[p].sequence;
    Spans spans = legacyDigest(s, enzyme, c.digest, c.minLength, c.maxLength,
                               c.maxMissedCleavages);
    numPeptides[p] = spans.size();
    for (const auto& span : spans) {
      std::string pep = s.substr(static_cast<size_t>(span.first),
                                 static_cast<size_t>(span.second));
      index[pep].push_back(p);
      peps[p].push_back(pep);
      if (pep[0] == 'M' && span.first == 0 &&
          span.second - 1 >= c.minLength) {
        index[pep.substr(1)].push_back(p);
        peps[p].push_back(pep.substr(1));
      }
    }
  }
  return index;
}

}  // namespace

// The flat index gives the same protein intersections (with repeats) and
// peptide counts as the map of strings it replaces.
TEST(ProteinDigestTest, PeptideIndexMatchesMap) {
  // proteins with shared peptides: random proteins, copies, fragments and
  // variants without the N-terminal methionine
  std::vector<std::string> seqs = randomProteins(13, 80);
  std::vector<std::string> real = realProteins();
  seqs.insert(seqs.end(), real.begin(), real.end());
  std::mt19937 rng(5);
  const size_t numBase = seqs.size();
  for (size_t i = 0; i < numBase; ++i) {
    const std::string s = seqs[i];
    if (s.size() < 4) continue;
    switch (rng() % 5) {
      case 0: seqs.push_back(s); break;
      case 1: seqs.push_back(s.substr(rng() % (s.size() / 2))); break;
      case 2: seqs.push_back(s.substr(0, s.size() / 2 + rng() % (s.size() / 2))); break;
      case 3: if (s[0] == 'M') seqs.push_back(s.substr(1)); break;
      default: break;
    }
  }
  std::shuffle(seqs.begin(), seqs.end(), rng);
  std::vector<FastaProtein> proteins;
  for (size_t i = 0; i < seqs.size(); ++i) {
    proteins.push_back(FastaProtein{"p" + std::to_string(i), seqs[i]});
  }
  std::unique_ptr<Enzyme> trypsin = Enzyme::createEnzyme(Enzyme::TRYPSIN);
  size_t numMulti = 0;
  for (const DigestConstraint& c :
       {DigestConstraint{DigestType::Full, 6, 50, 2},
        DigestConstraint{DigestType::Full, 4, 20, 0},
        DigestConstraint{DigestType::Partial, 5, 15, 1}}) {
    std::vector<std::vector<std::string> > peps;
    std::vector<size_t> numPeptides;
    std::map<std::string, std::vector<size_t> > map =
        mapIndex(proteins, *trypsin, c, peps, numPeptides);
    PeptideIndex index(proteins, std::vector<size_t>(), *trypsin, c);
    std::vector<size_t> got;
    for (size_t p = 0; p < proteins.size(); ++p) {
      ASSERT_EQ(numPeptides[p], index.numPeptides(p));
      std::vector<size_t> expected;
      for (size_t k = 0; k < peps[p].size(); ++k) {
        const std::vector<size_t>& list = map[peps[p][k]];
        if (k == 0) {
          expected = list;
        } else {
          std::vector<size_t> next;
          std::set_intersection(expected.begin(), expected.end(),
                                list.begin(), list.end(),
                                std::back_inserter(next));
          expected.swap(next);
        }
      }
      index.proteinsContainingAllPeptides(p, got);
      if (expected.size() > 1) {
        ASSERT_TRUE(same(expected, got)) << "protein " << p;
        ++numMulti;
      } else {
        ASSERT_LT(got.size(), 2u) << "protein " << p;
      }
    }
    // the same on a subset of the proteins, as for the candidate groups
    std::vector<size_t> subset;
    for (size_t p = 0; p < proteins.size(); p += 3) subset.push_back(p);
    PeptideIndex subIndex(proteins, subset, *trypsin, c);
    std::map<std::string, std::vector<size_t> > subMap;
    for (size_t p : subset) {
      for (const std::string& pep : peps[p]) subMap[pep].push_back(p);
    }
    for (size_t p : subset) {
      ASSERT_EQ(numPeptides[p], subIndex.numPeptides(p));
      std::vector<size_t> expected;
      for (size_t k = 0; k < peps[p].size(); ++k) {
        const std::vector<size_t>& list = subMap[peps[p][k]];
        if (k == 0) {
          expected = list;
        } else {
          std::vector<size_t> next;
          std::set_intersection(expected.begin(), expected.end(),
                                list.begin(), list.end(),
                                std::back_inserter(next));
          expected.swap(next);
        }
      }
      subIndex.proteinsContainingAllPeptides(p, got);
      if (expected.size() > 1) {
        ASSERT_TRUE(same(expected, got)) << "protein " << p;
      } else {
        ASSERT_LT(got.size(), 2u) << "protein " << p;
      }
    }
  }
  EXPECT_GT(numMulti, 50u);
}

// FASTA parsing: headers, sequence formatting, '>' inside a line and
// proteins longer than the 39,998 residues the old reader allowed.
TEST(ProteinDigestTest, ReadFasta) {
  std::string longSeq(100000, 'A');
  for (size_t i = 7; i < longSeq.size(); i += 13) longSeq[i] = 'K';
  std::ostringstream fasta;
  fasta << "junk before the first protein >notaprotein\n"
        << ">sp|P1|ONE_HUMAN description > with a gt\r\n"
        << "MKWV tfis\r\n12 LLF*\n\n"
        << ">  P2 leading blanks\n"
        << "ACDE>P3 starts inside a sequence line\nGHIK\n"
        << ">P4\n";
  for (size_t i = 0; i < longSeq.size(); i += 60) {
    fasta << longSeq.substr(i, 60) << "\n";
  }
  fasta << ">P5 empty\n>P6 no newline at the end\nMMM";
  std::istringstream in(fasta.str());
  std::vector<FastaProtein> proteins;
  ASSERT_TRUE(PickedProtein::readFasta(in, proteins));
  ASSERT_EQ(6u, proteins.size());
  EXPECT_EQ("sp|P1|ONE_HUMAN", proteins[0].id);
  EXPECT_EQ("MKWVTFISLLF", proteins[0].sequence);
  EXPECT_EQ("P2", proteins[1].id);
  EXPECT_EQ("ACDE", proteins[1].sequence);
  EXPECT_EQ("P3", proteins[2].id);
  EXPECT_EQ("GHIK", proteins[2].sequence);
  EXPECT_EQ("P4", proteins[3].id);
  EXPECT_EQ(longSeq, proteins[3].sequence);
  EXPECT_EQ("P5", proteins[4].id);
  EXPECT_EQ("", proteins[4].sequence);
  EXPECT_EQ("P6", proteins[5].id);
  EXPECT_EQ("MMM", proteins[5].sequence);

  std::istringstream none("no header line\n");
  proteins.clear();
  EXPECT_TRUE(PickedProtein::readFasta(none, proteins));
  EXPECT_TRUE(proteins.empty());
}
