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
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <ctime>

#include "EludeCaller.h"
#include "Normalizer.h"
#include "Option.h"
#include "Enzyme.h"
#include "Globals.h"
#include "DataManager.h"

const double EludeCaller::kFractionPeptides = 0.95;
string EludeCaller::library_path_ = "models/";
double EludeCaller::lts_coverage_ = 0.95;
double EludeCaller::hydrophobicity_diff_ = 5.0;

EludeCaller::EludeCaller():automatic_model_sel_(false), append_model_(false),
                           linear_calibration_(true), remove_duplicates_(false),
                           remove_in_source_(false), remove_non_enzymatic_(false),
                           context_format_(false), test_includes_rt_(false),
                           remove_common_peptides_(false), train_features_table_(NULL),
                           test_features_table_(NULL), processed_test_(false),
                           rt_model_(NULL), ignore_ptms_(false), lts(NULL) {
  Normalizer::setType(Normalizer::UNI);
  Enzyme::setEnzyme(Enzyme::TRYPSIN);
  RetentionFeatures::set_ignore_ptms(false);
  the_normalizer_ = Normalizer::getNormalizer();
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
  if (rt_model_) {
    delete rt_model_;
    rt_model_ = NULL;
  }
  if (lts != NULL) {
    delete lts;
  }
  DeleteRTModels();
}

void EludeCaller::DeleteRTModels() {
  if (rt_models_.size() <= 0) {
    return;
  }
  vector<RetentionModel*>::iterator it = rt_models_.begin();
  for( ; it != rt_models_.end(); ++it) {
    delete (*it);
  }
  rt_models_.clear();
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
                   "File to save the ions.",
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
                   "The model will NOT be linearly adjusted.",
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
  cmd.defineOption("p",
                   "ignore-ptms",
                   "If there are ptms in the test set that were not present when "
                   "training the model, they will be ignored and the index value of "
                   "the unmodified amino acid is used ",
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
    int n = library_path_.length();
    if (library_path_[n - 1] != '\\' && library_path_[n - 1] != '/') {
      library_path_ += "/";
    }
  }
  if (cmd.optionSet("d")) {
    append_model_ = true;
  }
  if (cmd.optionSet("j")) {
    linear_calibration_ = false;
  }
  if (cmd.optionSet("c")) {
    double coverage = cmd.getDouble("r", 0.0, 1.0);
    LTSRegression::setCoverage(coverage);
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
  if (cmd.optionSet("p")) {
    RetentionFeatures::set_ignore_ptms(true);
    ignore_ptms_ = true;
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
  DataManager::LoadPeptides(train_file_, true, context_format_,
      train_psms_, train_aa_alphabet_);
  // remove duplicates
  DataManager::RemoveDuplicates(train_psms_);
  // remove in source fragments
  if (!test_file_.empty()) {
    // load the test peptides
    DataManager::LoadPeptides(test_file_, test_includes_rt_, context_format_,
        test_psms_, test_aa_alphabet_);
    // remove duplicates
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

  // shuffle the training data
  random_shuffle(train_psms_.begin(), train_psms_.end());
}

int EludeCaller::NormalizeRetentionTimes(vector<PSMDescription> &psms) {
  PSMDescription::setPSMSet(psms);
  PSMDescription::normalizeRetentionTimes(psms);
}

/* train model */
int EludeCaller::TrainRetentionModel() {
  /* initialize model */
  rt_model_ = new RetentionModel(the_normalizer_);
  /* build the retention index */
  retention_index_ = rt_model_->BuildRetentionIndex(train_aa_alphabet_, false, train_psms_);
  /* train the model */
  rt_model_->TrainRetentionModel(train_aa_alphabet_, retention_index_, true, train_psms_);
}

/* process test data */
int EludeCaller::ProcessTestData() {
  // load the peptides
  DataManager::LoadPeptides(test_file_, test_includes_rt_, context_format_, test_psms_, test_aa_alphabet_);
  // remove duplicates
  if (remove_duplicates_) {
    DataManager::RemoveDuplicates(test_psms_);
  }
  if (remove_common_peptides_)  {
    DataManager::RemoveCommonPeptides(test_psms_, train_psms_);
  }
  // remove in source
  if (!in_source_file_.empty() || remove_in_source_) {
    if (test_includes_rt_) {
      vector< pair<PSMDescription, string> > in_source_fragments;
      in_source_fragments = DataManager::RemoveInSourceFragments(hydrophobicity_diff_,
          RetentionFeatures::kKyteDoolittle, remove_in_source_, train_psms_, test_psms_);
      if (!in_source_file_.empty()) {
        DataManager::WriteInSourceToFile(in_source_file_, in_source_fragments);
      }
    } else {
      if (VERB >= 4) {
        cerr << "Warning: in-source fragments in the test data cannot be identified "
             << "since the retention time is not included. All peptides are included "
             << "in the subsequent analyses" << endl;
      }
    }
  }
  // remove non enzymatic
  if (remove_non_enzymatic_) {
    if (context_format_) {
      DataManager::RemoveNonEnzymatic(test_psms_);
    } else {
      if (VERB >= 4) {
        cerr << "Warning: non-enzymatic peptides cannot be detected unless the peptides"
             << " are give in the format A.XXX.Y. All peptides are included in the"
             << " subsequent analyses" << endl;
      }
    }
  }
  processed_test_ = true;
  return 0;
}

// get retention times as vectors
pair<vector<double> , vector<double> > GetRTs(const vector<PSMDescription> &psms) {
  vector<double> rts, prts;
  vector<PSMDescription>::const_iterator it = psms.begin();
  for ( ; it != psms.end(); ++it) {
    rts.push_back(it->retentionTime);
    prts.push_back(it->predictedTime);
  }
  return make_pair(prts, rts);
}

/* add a model to the library */
int EludeCaller::AddModelLibrary() const {
  string file_name;
  if (!load_model_file_.empty()) {
    file_name = load_model_file_;
  } else if (!save_model_file_.empty()) {
    file_name = save_model_file_;
  } else if (!train_file_.empty()) {
    file_name = GetFileName(train_file_);
  } else {
    time_t current_time = time(NULL);
    struct tm* time_info = localtime(&current_time);
    file_name = asctime(time_info);
    file_name = file_name.substr(4, (file_name.length() - 10));
    int found = file_name.find_first_of(" ");
      while (found != string::npos) {
        file_name.replace(found, 1, "_");
        found = file_name.find_first_of(" ");
      }
  }
  file_name = library_path_ + file_name;
  rt_model_->SaveModelToFile(file_name);
  return 0;
}

/* main function in Elude */
int EludeCaller::Run() {
  pair<int, double> best_model(-1, -1.0);
  // build a rt model, either by training one or by loading one
  if (!train_file_.empty() && !load_model_file_.empty() && VERB >= 4
      && !linear_calibration_) {
    cerr << "Warning: a model can be either trained or loaded from a file. "
         << "The two options should not be used together, unless linear calibration "
         << "should be carried out. In such a case please use the -j option. "
         << "The model will be trained using the peptides in " << train_file_ << endl;
  }
  if (!train_file_.empty()) {
    ProcessTrainData();
    // initialize the feature table
    train_features_table_ = DataManager::InitFeatureTable(
         RetentionFeatures::kMaxNumberFeatures, train_psms_);
    if (automatic_model_sel_) {
      best_model = AutomaticModelSelection();
    } else if (load_model_file_.empty()){
      TrainRetentionModel();
    }
  }
  // load a model from a file
  if (!load_model_file_.empty() && !automatic_model_sel_) {
    rt_model_ = new RetentionModel(the_normalizer_);
    rt_model_->LoadModelFromFile(load_model_file_);
  }
  // save the model
  if (!save_model_file_.empty()) {
    if (rt_model_ != NULL && !rt_model_->IsModelNull()) {
      rt_model_->SaveModelToFile(save_model_file_);
    } else if (VERB >= 2) {
      cerr << "Warning: No trained model available. Nothing to save to "
           << save_model_file_ << endl;
    }
  }
  // append a file to the library
  if (append_model_) {
    if (automatic_model_sel_) {
      if (VERB >= 3) {
        cerr << "Warning: The model should already be in the library if "
             << "the automatic model selection option is employed. No model "
             << "will be added to the library"<< endl;
      }
    } else if (rt_model_ == NULL) {
      if (VERB >= 3) {
        cerr << "Warning: No model available, nothing to append to the library."
             << endl;
      }
    } else {
      AddModelLibrary();
    }
  }

  // test a model
  if (!test_file_.empty()) {
    // process the test data
    if (!processed_test_) {
      ProcessTestData();
    }
    // initialize the feature table
    test_features_table_ = DataManager::InitFeatureTable(
            RetentionFeatures::kMaxNumberFeatures, test_psms_);
    if (automatic_model_sel_) {
      int index = best_model.first;
      if (index < 0) {
        if (VERB >= 2) {
          cerr << "Error: No model available to predict rt. Execution aborted." << endl;
        }
        exit(1);
      }
      rt_models_[index]->PredictRT(test_aa_alphabet_, ignore_ptms_, test_psms_);
      if (linear_calibration_ && train_psms_.size() > 1) {
        rt_models_[index]->PredictRT(train_aa_alphabet_, ignore_ptms_, train_psms_);
      }
    } else {
      rt_model_->PredictRT(test_aa_alphabet_, ignore_ptms_, test_psms_);
      if (linear_calibration_ && train_psms_.size() > 1) {
        rt_model_->PredictRT(train_aa_alphabet_, ignore_ptms_, train_psms_);
      }
    }
    // linear calibration is performed only for automatic model selection or when
    // loading a model from a file
    if (linear_calibration_ && (automatic_model_sel_ || !load_model_file_.empty())) {
      if (train_psms_.size() <= 1) {
        if (VERB >= 3) {
          cerr << "Warning: at least 2 training psms are needed to calibrate the model. "
               << "No calibration performed. " << endl;
         }
       } else {
         // get the a and b coefficients
         pair<vector<double> , vector<double> > rts = GetRTs(train_psms_);
         lts = new LTSRegression();
         lts->setData(rts.first, rts.second);
         lts->runLTS();
         AdjustLinearly(test_psms_);
       }
    }
    // compute performance measures
    if (test_includes_rt_) {
      double rank_correl = ComputeRankCorrelation(test_psms_);
      double pearson_correl = ComputePearsonCorrelation(test_psms_);
      double win = ComputeWindow(test_psms_);
      if (VERB >= 3) {
        cerr << "Performance measures for the test data: " << endl;
        cerr << "  Pearson's correlation r = " << pearson_correl << endl;
        cerr << "  Spearman's rank correlation rho = " << rank_correl << endl;
        cerr << "  Delta_t 95% = " << win << endl;
      }
    }
    // write the predictions to file
    if (!output_file_.empty()) {
      DataManager::WriteOutFile(output_file_, test_psms_, test_includes_rt_);
    }
  }
}

int EludeCaller::AdjustLinearly(vector<PSMDescription> &psms) {
  vector<PSMDescription>::iterator it = psms.begin();
  for( ; it != psms.end(); ++it) {
    it->predictedTime = lts->predict(it->predictedTime);
  }
}

/* Load the best model from the library; the function returns a pair consisting of
 * the index of this model in the vector of models and the rank correlation
 * obtained on the calibration peptides using this model */
pair<int, double> EludeCaller::AutomaticModelSelection() {
  vector<string> model_files = ListDirFiles(library_path_);
  if (model_files.size() == 0) {
    if (VERB >= 2) {
      cerr << "Error: No model available in " << library_path_ << endl;
    }
    return make_pair(-1, -1.0);
  } else {
    RetentionModel *m = new RetentionModel(the_normalizer_);
    m->LoadModelFromFile(model_files[0]);
    rt_models_.push_back(m);
    // no calibration psms available
    if (!train_psms_.size() > 0) {
      if (VERB >= 3) {
        cerr << "Warning: no calibration psms available. First model available "
             << "in the library is selected (" << model_files[0] << ")." << endl;
      }
      return make_pair(0, -1.0);
    } else {
      m->PredictRT(train_aa_alphabet_, ignore_ptms_, train_psms_);
      double rank_correl = ComputeRankCorrelation(train_psms_);
      double best_correl = rank_correl;
      int best_index = 0;
      //cout << model_files[0] << ", correl = " << rank_correl << endl;
      for(int i = 1; i < model_files.size(); ++i) {
        m = new RetentionModel(the_normalizer_);
        m->LoadModelFromFile(model_files[i]);
        rt_models_.push_back(m);
        m->PredictRT(train_aa_alphabet_, ignore_ptms_, train_psms_);
        rank_correl = ComputeRankCorrelation(train_psms_);
        //cout << model_files[i] << ", correl = " << rank_correl << endl;
        if (rank_correl > best_correl) {
          best_correl = rank_correl;
          best_index = i;
        }
      }
      //cout << "Select model: " << model_files[best_index] << endl;
      return make_pair(best_index, best_correl);
    }
  }
}

/* Return a list of files in a directory dir_name*/
vector<string> EludeCaller::ListDirFiles(const std::string &dir_name) {
  string dir = dir_name;
  char last_char = dir[dir.length() - 1];
  if (last_char != '/' && last_char != '\\') {
    dir += "/";
  }
  DIR* dp = opendir(dir.c_str());
  if (dp == NULL) {
    if (VERB >= 1) {
      cerr << "Error: Unable to open " << dir << endl;
      cerr << "Execution aborted. " << endl;
    }
    exit(1);
  }
  vector<string> files;
  string file_name;
  struct stat info;
  struct dirent* dirp = readdir(dp);
  while (dirp != NULL) {
    file_name = string(dirp->d_name);
    if (!(file_name == ".") && !(file_name == "..")) {
      stat((dir + file_name).c_str(), &info);
      if (!S_ISDIR(info.st_mode)) {
        files.push_back(dir + file_name);
      }
    }
    dirp = readdir(dp);
  }
  closedir(dp);
  return files;
}

bool EludeCaller::FileExists(const string &file) {
  struct stat file_info;
  int int_stat;
  if (stat(file.c_str(), &file_info) != 0) {
    return false;
  } else {
    return true;
  }
}

/* find the best line that fits the data (basically the coefficients a, b)
 * The predicted time is the explanatory variable */
void EludeCaller::FindLeastSquaresSolution(const vector<PSMDescription> &psms,
    double& a, double& b) {
  double sum_x = 0.0, sum_y = 0.0;
  double sum_x_squared = 0.0, sum_xy = 0.0;
  double x, y;
  vector<PSMDescription>::const_iterator it = psms.begin();
  for ( ; it != psms.end(); ++it) {
    y = it->retentionTime;
    x = it->predictedTime;
    sum_x += x;
    sum_y += y;
    sum_x_squared += x * x;
    sum_xy += x * y;
  }
  double n = (double)psms.size();
  double avg_x = sum_x / n;
  double avg_y = sum_y / n;
  double ssxx = sum_x_squared - (n * avg_x * avg_x);
  double ssxy = sum_xy - (n * avg_x * avg_y);
  a = ssxy / ssxx;
  b = avg_y - (a * avg_x);
}

bool ComparePsmsDeltaRT(PSMDescription psm1, PSMDescription psm2) {
  double delta_rt1 = psm1.getPredictedRetentionTime() - psm1.getRetentionTime();
  double delta_rt2 = psm2.getPredictedRetentionTime() - psm2.getRetentionTime();
  return delta_rt1 < delta_rt2;
}

/*compute Delta t(95%) window */
double EludeCaller::ComputeWindow(vector<PSMDescription> &psms) {
  double win, diff;
  int nr = (int) round(kFractionPeptides * (double) psms.size());
  int k = psms.size() - nr, i = 1;

  sort(psms.begin(), psms.end(), ComparePsmsDeltaRT);
  win = (psms[nr - 1].getPredictedRetentionTime() - psms[nr - 1].getRetentionTime())
      - (psms[0].getPredictedRetentionTime() - psms[0].getRetentionTime());
  while (i <= k) {
    diff = (psms[i + nr - 1].getPredictedRetentionTime()
        - psms[i + nr - 1].getRetentionTime()) - (psms[i].getPredictedRetentionTime()
        - psms[i].getRetentionTime());
    if (diff < win) {
      win = diff;
    }
    i++;
  }

  return win;
}

/* Compare 2 psms according to retention time */
bool ComparePsmsRT(PSMDescription psm1, PSMDescription psm2) {
  return psm1.getRetentionTime() < psm2.getRetentionTime();
}

/* Compare 2 psms according to predicted retention time */
bool ComparePsmsPRT(pair<PSMDescription, int> psm1, pair<PSMDescription,
    int> psm2) {
  double prt1, prt2;
  prt1 = psm1.first.getPredictedRetentionTime();
  prt2 = psm2.first.getPredictedRetentionTime();
  return prt1 < prt2;
}

/* calculate Spearman's rank correlation */
double EludeCaller::ComputeRankCorrelation(vector<PSMDescription> &psms) {
  double corr = 0.0, d = 0.0, avg_rank, rank_p, rank_o;
  int i, j;
  int n = psms.size();
  vector<pair<PSMDescription, double> > rankedPsms;

  // sort peptides according to observed retention time
  sort(psms.begin(), psms.end(), ComparePsmsRT);
  // record ranks
  i = 0;
  while (i < n) {
    avg_rank = j = i + 1;
    while ((j < n) && (psms[i].getRetentionTime()
        == psms[j].getRetentionTime())) {
      avg_rank += ++j;
    }
    avg_rank = avg_rank / (double)(j - i);
    for (int k = i; k < j; ++k) {
      rankedPsms.push_back(make_pair(psms[k], avg_rank));
    }
    i = j;
  }
  // sort peptides according to predicted rt
  sort(rankedPsms.begin(), rankedPsms.end(), ComparePsmsPRT);
  // calculate sum of squared differences btw ranks
  i = 0;
  while (i < n) {
    // calculate rank of predicted rt
    rank_p = j = i + 1;
    while ((j < n) && (rankedPsms[i].first.getPredictedRetentionTime()
        == rankedPsms[j].first.getPredictedRetentionTime())) {
      rank_p += ++j;
    }
    rank_p = rank_p / (double)(j - i);
    // calculate and add squared difference
    for (int k = i; k < j; ++k) {
      d += pow(rankedPsms[k].second - rank_p, 2);
    }
    // increase i
    i = j;
  }

  return 1.0 - ((6.0 * d) / (double)(n * (pow(n, 2.) - 1)));
}

/* Compute Pearson's correlation coefficient */
double EludeCaller::ComputePearsonCorrelation(vector<PSMDescription> & psms) {
  int no_psms = psms.size();
  // calculate means
  double sum_obs = 0.0, sum_pred = 0.0;
  for (int i = 0; i < no_psms; i++) {
    sum_obs += psms[i].retentionTime;
    sum_pred += psms[i].predictedTime;
  }
  double avg_obs = sum_obs / no_psms;
  double avg_pred = sum_pred / no_psms;
  // calculate stdevs
  sum_obs = 0.0;
  sum_pred = 0.0;
  double dev_obs, dev_pred;
  double numerator = 0.0;
  for (int i = 0; i < no_psms; i++) {
    dev_obs = psms[i].retentionTime - avg_obs;
    dev_pred = psms[i].predictedTime - avg_pred;
    numerator += dev_obs * dev_pred;
    sum_obs += pow(dev_obs, 2);
    sum_pred += pow(dev_pred, 2);
  }
  double sd_obs = 1, sd_pred = 1;
  sd_obs = sqrt(sum_obs / (no_psms - 1));
  sd_pred = sqrt(sum_pred / (no_psms - 1));

  return numerator / ((no_psms - 1) * sd_obs * sd_pred);
}

/* get the filename from a path (excluding the extension) */
string EludeCaller::GetFileName(const string &path) {
  int pos1, pos2;
  if ((pos1 = path.find_last_of("/")) == string::npos) {
    if ((pos1 = path.find_last_of("\\"))== string::npos) {
      pos1 = -1;
    }
  }
  string file_name = path.substr(pos1 + 1, path.length());
  if ((pos2 = file_name.find_last_of(".")) != string::npos);
    file_name = file_name.substr(0, file_name.length() - pos2);
  return file_name;
}

/**********************************************************************/
/*** Additional functions from the provious implementation of Elude ***/
/**********************************************************************/
/*


// write peptides and retention times to a file
void RTPredictor::writeRetentionTimeFile(const char* fileName, vector<
    PSMDescription> & psms) {
  ofstream out(fileName);
  vector<PSMDescription>::iterator it;
  for (it = psms.begin(); it != psms.end(); ++it) {
    out << (it->peptide) << " " << (it->retentionTime) << "\n";
  }
  out.close();
}

// write the value of the features to a file
void RTPredictor::writeFeaturesFile(const char* file, vector<
    PSMDescription> & psms, bool unnormalized) {
  ofstream out(file);
  vector<PSMDescription>::iterator it;
  double predictedRT, observedRT;
  out << "# File generated by Elude " << endl;
  out << "# " << __DATE__ << " , " << __TIME__ << endl;
  out << "Peptide\tpredicted_retention_time\tobserved_retention_time"
      << endl;
  if (unnormalized) {
    vector<double*> tmp = PSMDescription::getRetFeatures(psms);
    theNormalizer->unNormalizeSet(tmp);
  }
  for (it = psms.begin(); it != psms.end(); ++it) {
    predictedRT = it->predictedTime;
    observedRT = it->retentionTime;
    out << (it->peptide) << "\t" << predictedRT << "\t" << observedRT;
    for (size_t i = 0; i < *(theNormalizer->getNumRetFeatures()); ++i) {
      out << "\t" << it->retentionFeatures[i];
    }
    out << "\n";
  }
  out.close();
}

// print a vector of PSMs
void RTPredictor::printPsms(vector<PSMDescription> & psms) {
  cerr << endl << "Printing..." << endl;
  vector<PSMDescription>::iterator it;
  for (it = psms.begin(); it != psms.end(); it++) {
    cerr << it->peptide << it->retentionTime << endl;
  }
  cerr << "Done." << endl << endl;
}

// calculate ms between observed and predicted retention times
double RTPredictor::computeMS(vector<PSMDescription> & psms) {
  double ms = 0.0, diff = 0.0;
  int noPsms;
  if (VERB > 2) {
    cerr << endl << "Computing ms..." << endl;
  }
  noPsms = psms.size();
  for (int i = 0; i < noPsms; i++) {
    diff = psms[i].predictedTime - psms[i].retentionTime;
    ms += diff * diff;
  }
  ms = ms / noPsms;
  if (VERB > 2) {
    cerr << "Done." << endl;
  }
  return ms;
}
/**********************************************************************/
