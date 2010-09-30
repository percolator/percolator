/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@cbr.su.se>

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
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/* This files stores the implementations of the methods for the EludeCaller class */
#include <iostream>
#include <sstream>
#include <algorithm>

#include "EludeCaller.h"
#include "Normalizer.h"
#include "Option.h"
#include "Enzyme.h"
#include "Globals.h"
#include "DataManager.h"

string EludeCaller::library_path_ = "models/";
double EludeCaller::lts_coverage_ = 0.95;
double EludeCaller::hydrophobicity_diff_ = 5.0;

EludeCaller::EludeCaller():automatic_model_sel_(false), append_model_(false),
                           linear_calibration_(true), remove_duplicates_(false),
                           remove_in_source_(false), remove_non_enzymatic_(false),
                           context_format_(false), test_includes_rt_(false),
                           remove_common_peptides_(false), train_features_table_(NULL),
                           test_features_table_(NULL), processed_test_(false) {
  Normalizer::setType(Normalizer::UNI);
  Enzyme::setEnzyme(Enzyme::TRYPSIN);
  string basic_aa[] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};
  train_aa_alphabet_.insert(basic_aa, basic_aa + 20);
  test_aa_alphabet_.insert(basic_aa, basic_aa + 20);
}

EludeCaller::~EludeCaller() {
  if (train_features_table_) {
    DataManager::CleanUpTable(train_psms_, train_features_table_);
    train_features_table_ = NULL;
  }
  if (test_features_table_) {
    DataManager::CleanUpTable(test_psms_, test_features_table_);
    test_features_table_ = NULL;
  }
}

/* introductory message */
string EludeCaller::Greeter() const {
  ostringstream oss;
  oss << "Elude version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under MIT License" << endl;
  oss << "Written by Lukas Käll (lukas.kall@cbr.su.se) "
         "and Luminita Moruz (lumi@sbc.su.se)" << endl;
  oss << "Center for Biomembrane Research." << endl;
  oss << "Dept. of Biochemistry, Stockholm University, Stockholm." << endl;
  oss << "Usage:" << endl;
  oss << "   elude [options]" << endl << endl;
  return oss.str();
}

/* parse the command line arguments */
bool EludeCaller::ParseOptions(int argc, char** argv) {
  ostringstream intro;
  intro << Greeter() << endl;
  CommandLineParser cmd(intro.str());
  // define available options
  cmd.defineOption("v",
                   "verbose",
                   "Set verbosity of output: 0 = no processing info, 5 = all, default is 2.",
                   "level");
  cmd.defineOption("t",
                   "train",
                   "Specifies the file including the training data.",
                   "filename");
  cmd.defineOption("e",
                   "evaluate",
                   "Specifies the file including the test data.",
                   "filename");
  cmd.defineOption("s",
                   "save-model",
                   "Specifies the file in which the model will be saved.",
                   "filename");
  cmd.defineOption("l",
                   "load-model",
                   "Specifies a file including a SVR model to be loaded.",
                   "filename");
  cmd.defineOption("o",
                   "out",
                   "File to save the predictions.",
                   "filename");
  cmd.defineOption("a",
                   "auto",
                   "The SVR model is selected automatically from the library.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("b",
                   "lib-path",
                   "Specifies the path to the library.",
                   "filename");
  cmd.defineOption("d",
                   "append",
                   "Append current model to library.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("j",
                   "no_linear_adjust",
                   "The model will not be linearly adjusted.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("c",
                   "lts-coverage",
                   "Specifies the fraction of data used in calibrating a model via LTS. "
                   "This option is not valid when the -j option is used.",
                   "value");
  cmd.defineOption("u",
                   "unique",
                   "Remove all redundant peptides from the test set",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("k",
                   "common",
                   "Remove the peptides from the train that are also in the test set",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("y",
                   "no-in-source",
                   "Specifies that in source fragments should be removed from the test set."
                   "This option can be used only when the test set includes retention time.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("i",
                   "save-in-source",
                   "The file where the detected in source fragments are stored",
                   "filename");
  cmd.defineOption("z",
                   "enzyme",
                   "The enzyme used for digestion. Possible values {NO_ENZYME,TRYPSIN,CHYMOTRYPSIN,ELASTASE}."
                   "By default: TRYPSIN",
                   "value");
  cmd.defineOption("x",
                   "remove-non-enzymatic",
                   "All non enzymatic peptides will be removed from both train and test."
                   "The option is available only when the sequence is given as A.XXX.B",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("r",
                   "retention-index",
                   "File to save the trained retention index.",
                   "filename");
  cmd.defineOption("f",
                   "context-format",
                   "The peptides are given in the format A.XXX.B, where XXX is the"
                   "peptide sequence and A and B are the previous and next amino acid"
                   "in the protein sequence. If the peptide is at the beginning (end)"
                   "of the protein, then A(B) will be -.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("g",
                   "test-rt",
                   "The test file includes rt. In this case the in source-fragments in"
                   "the test data can be detected and various performance measures for"
                   "the test data will be displayed.",
                   "",
                   TRUE_IF_SET);

  cmd.parseArgs(argc, argv);

  // process options
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (cmd.optionSet("t")) {
    train_file_ = cmd.options["t"];
  }
  if (cmd.optionSet("e")) {
    test_file_ = cmd.options["e"];
  }
  if (cmd.optionSet("s")) {
    save_model_file_ = cmd.options["s"];
  }
  if (cmd.optionSet("l")) {
    load_model_file_ = cmd.options["l"];
  }
  if (cmd.optionSet("o")) {
    output_file_ = cmd.options["o"];
  }
  if (cmd.optionSet("a")) {
    automatic_model_sel_ = true;
  }
  if (cmd.optionSet("b")) {
    library_path_ = cmd.options["b"];
  }
  if (cmd.optionSet("d")) {
    append_model_ = true;
  }
  if (cmd.optionSet("j")) {
    linear_calibration_ = false;
  }
  if (cmd.optionSet("c")) {
    lts_coverage_ = cmd.getDouble("r", 0.0, 1.0);
  }
  if (cmd.optionSet("u")) {
    remove_duplicates_ = true;
  }
  if (cmd.optionSet("k")) {
    remove_common_peptides_ = true;
  }
  if (cmd.optionSet("y")) {
    remove_in_source_ = true;
  }
  if (cmd.optionSet("i")) {
    in_source_file_ = cmd.options["i"];
  }
  if (cmd.optionSet("z")) {
    SetEnzyme(cmd.options["z"]);
  }
  if (cmd.optionSet("x")) {
    remove_non_enzymatic_ = true;
  }
  if (cmd.optionSet("r")) {
    index_file_ = cmd.options["r"];
  }
  if (cmd.optionSet("f")) {
    context_format_ = true;
  }
  if (cmd.optionSet("g")) {
    test_includes_rt_ = true;
  }
  return true;
}

void EludeCaller::SetEnzyme(const string &enzyme) {
  if ((enzyme.compare("CHYMOTRYPSIN") == 0) ||
      (enzyme.compare("chymotrypsin") == 0)) {
    Enzyme::setEnzyme(Enzyme::CHYMOTRYPSIN);
  } else if ((enzyme.compare("ELASTASE") == 0) ||
             (enzyme.compare("elastase") == 0)) {
    Enzyme::setEnzyme(Enzyme::ELASTASE);
  } else if ((enzyme.compare("TRYPSIN") == 0) ||
             (enzyme.compare("trypsin") == 0)) {
    Enzyme::setEnzyme(Enzyme::TRYPSIN);
  } else {
    Enzyme::setEnzyme(Enzyme::NO_ENZYME);
    if (VERB >= 3) {
      cerr << "Warning: Enzyme " + enzyme + " not recognized. No enzyme set. Please use"
           << "one of the values {NO_ENZYME, TRYPSIN, CHYMOTRYPSIN, ELASTASE}." << endl;
    }
  }
}

/* process train data when a model is trained*/
int EludeCaller::ProcessTrainData() {
  // load the training peptides
  DataManager::LoadPeptides(train_file_, true, context_format_, train_psms_, train_aa_alphabet_);
  // remove duplicates
  DataManager::RemoveDuplicates(train_psms_);
  // remove in source fragments
  if (!test_file_.empty()) {
    // load the test peptides
    if (test_includes_rt_) {
      DataManager::LoadPeptides(test_file_, true, context_format_, test_psms_, test_aa_alphabet_);
    } else {
      DataManager::LoadPeptides(test_file_, false, context_format_, test_psms_, test_aa_alphabet_);
    }
    if (remove_duplicates_) {
      DataManager::RemoveDuplicates(test_psms_);
    }
    // remove from the train the peptides in the test
    if (remove_common_peptides_) {
      DataManager::RemoveCommonPeptides(test_psms_, train_psms_);
    }
    processed_test_= true;
  }
  // remove in source fragmentation
  vector< pair<PSMDescription, string> > in_source_fragments;
  if (test_includes_rt_) {
    in_source_fragments = DataManager::RemoveInSourceFragments(hydrophobicity_diff_,
        RetentionFeatures::kKyteDoolittle, remove_in_source_, train_psms_, test_psms_);
  } else {
    vector<PSMDescription> tmp;
    in_source_fragments = DataManager::RemoveInSourceFragments(hydrophobicity_diff_,
        RetentionFeatures::kKyteDoolittle, false, train_psms_, tmp);
  }
  if (!in_source_file_.empty()) {
    DataManager::WriteInSourceToFile(in_source_file_, in_source_fragments);
  }
  // remove non enzymatic
  if (remove_non_enzymatic_) {
    if (context_format_) {
      DataManager::RemoveNonEnzymatic(train_psms_);
      DataManager::RemoveNonEnzymatic(test_psms_);
    } else {
      if (VERB >= 4) {
        cerr << "Warning: non-enzymatic peptides cannot be detected unless the peptides"
             << " are give in the format A.XXX.Y. All peptides will be included in the"
             << " subsequent analyses" << endl;
      }
    }
  }
}

int EludeCaller::NormalizeRetentionTimes(vector<PSMDescription> &psms) {
  PSMDescription::setPSMSet(psms);
  PSMDescription::normalizeRetentionTimes(psms);
}


/* build the retention index by training a linear model */
/*
int EludeCaller::BuildRetentionIndex() {
  RetentionModel model;

  // set the amino acids alphabet
  model.SetAlphabet(train_aa_alphabet_);
  // set the aa as the active features
  model.SetAAFeatures();
  // compute the aa features
  model.ComputeAAFeatures(train_psms_);
  // normalize the features

  // initialize a linear svr
  model.InitSVR(true);


}*/
