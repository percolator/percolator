Spectra from:

Moruz, Luminita, et al. "Mass fingerprinting of complex mixtures: protein inference from high-resolution peptide masses and predicted retention times." Journal of proteome research 12.12 (2013): 5730-5741.

Input files:

- 103111-Yeast-2hr-01.RAW: original data set, not uploaded to github
- small.ms2: random subsample of 107 spectra corresponding to 292 mass-charge states assigned by bullseye
- Searched with Tide, MS-GF+ and X! Tandem, against the swissprot_yeast database concatenated with a mimic database of 9x the original size to calibrate the statistics, 10ppm, no miscleavages, full digestion, up to 2 methionine oxidations and variable modification of an acetylation of the peptide N terminal.

Test cases sqt2pin:

- separate target decoy searches (target.sqt and decoy.sqt)
- concatenated fasta database search (combined.sqt)
- PTMs (oxidation of methionine, K.QISSIM\[16\]SK.R, SpecId target_3664_2_1)
- Peptide with both a target and decoy protein (SpecId combined_8892_2_1 in combined.sqt)

Test cases msgf2pin:

- separate target decoy searches (target.mzid and decoy.mzid)
- concatenated fasta database search (combined.mzid)
- PTMs (oxidation of methionine, R.SM\[UNIMOD:35\]NGISIC\[UNIMOD:4\]GK.N, SpecId target_SII_10_1_288002_2_1)

Test cases tandem2pin:

- separate target decoy searches (target.t.xml and decoy.t.xml)
- concatenated fasta database search (combined.t.xml)
- PTMs (oxidation of methionine, R.M\[15.99492\]PFSHEVAM\[15.99492\]NGGIIVK.L, SpecId target_8_4_1)
- Peptide with both a target and decoy protein (SpecId combined_84_2_1 in combined.t.xml)
