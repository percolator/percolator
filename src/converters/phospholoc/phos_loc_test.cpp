#include <cstdio>
#include "phos_loc_spectrum.h"
#include "phos_loc_peptide.h"
#include "phos_loc_parameters.h"
#include "phos_loc_ion_series.h"
#include "phos_loc_match.h"
#include "phos_loc_ascore.h"
#include "phos_loc_phosphors.h"

using namespace phos_loc;

int main()
{
  Parameters paras;
  paras.Print();

  // Spectrum spec("s-pxsp-1.2068.2068.3.dta"); // AEHVAEADK.2.23153.dta
  Spectrum spec("s-pxsp-1.1771.1771.2.dta", paras.activation_type_);
  spec.Preprocess(TOPN_WINDOW, paras.preproc_parameters_,
                  Tolerance(0.5, 0, DA_TOL, MONO));
  spec.Print();


  std::vector<LocationMod> var_mods;
  // var_mods.push_back(LocationMod(5,21)); // unimod_id("Phospho") = 21
  // var_mods.push_back(LocationMod(13,21));
  // std::string pep_seq = "AGEPNSPDAEEANSPDVTAGCDPAGVHPPR"; // AEHVAEADK
  std::string pep_seq = "SQDSYPGSPSLSPR";
  var_mods.push_back(LocationMod(7,21));
  var_mods.push_back(LocationMod(11,21));

  IonSeries ion_series(pep_seq, var_mods, paras, spec.precursor().charge);
  int n = ion_series.GetNumPredictedIons();
  fprintf(stdout, "Number of all predicted ions: %d\n", n);
  std::map<IonType, std::vector<Ion> >::const_iterator it;
  for (it = ion_series.ion_matrix().begin();
       it != ion_series.ion_matrix().end(); ++it) {
    fprintf(stdout, "%s (%d+):\t", it->first.ion_name().c_str(), it->first.charge());
    for (std::vector<Ion>::const_iterator it_j = it->second.begin();
         it_j != it->second.end(); ++it_j) {
      fprintf(stdout, "%f ", it_j->m_over_z());
    }
    fprintf(stdout, "\n");
  }

  std::vector<Ion>::const_iterator it_k;
  for (it_k = ion_series.sorted_ions().begin();
       it_k != ion_series.sorted_ions().end(); ++it_k) {
    fprintf(stdout, "%s.%d: %f; ", it_k->ion_type().ion_name().c_str(),
            it_k->cleavage_site(), it_k->m_over_z());
  }
  fprintf(stdout, "\n");

  std::vector<std::vector<LocationMod> > var_mod_combs;
  var_mod_combs.push_back(var_mods);
  Ascore ascore(spec, pep_seq, var_mod_combs, paras);
  std::vector<int> phospho_locs;
  phospho_locs.push_back(var_mods[0].aa_idx);
  phospho_locs.push_back(var_mods[1].aa_idx);
  double score = ascore.GetPeptideScore(phospho_locs);
  double lscore = ascore.GetPhosphoSiteScore(phospho_locs);
  fprintf(stdout, "%.6f %.6f\n", score, lscore);

  PhosphoRS phosphors(spec, pep_seq, var_mod_combs, paras);
  double s1 = phosphors.GetPeptideScore(phospho_locs);
  double s2 = phosphors.GetPhosphoSiteScore(phospho_locs);
  fprintf(stdout, "%.6f %.6f\n", s1, s2);
}
