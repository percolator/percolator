#include <cstdio>
#include "spectrum.h"
#include "peptide.h"
#include "parameters.h"
#include "ion_series.h"
#include "match.h"
#include "ascore.h"

int main()
{
  Parameters paras;
  paras.Print();

  Spectrum spec("s-pxsp-1.2068.2068.3.dta"); // AEHVAEADK.2.23153.dta
  spec.Preprocess(TOPN_WINDOW, paras.preproc_parameters_,
                  Tolerance(0.5, 0, DA_TOL, MONO));
  spec.Print();


  std::vector<LocationMod> var_mods;
  var_mods.push_back(LocationMod(5,0));
  var_mods.push_back(LocationMod(13,0));
  std::string pep_seq = "AGEPNSPDAEEANSPDVTAGCDPAGVHPPR"; // AEHVAEADK
  IonSeries ion_series(pep_seq, var_mods, paras, spec.precursor().charge);

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
  double score = ascore.GetPeptideScore(var_mods, paras);
  fprintf(stdout, "PeptideScore=%.6f\n", score);
}
