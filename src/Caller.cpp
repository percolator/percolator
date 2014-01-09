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

#include "Caller.h"
#include "unistd.h"
#include <iomanip>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>

using namespace std;
#ifdef XML_SUPPORT
using namespace xercesc;
#endif //XML_SUPPORT

const unsigned int Caller::xval_fold = 3; /* number of folds for cross validation*/
const double requiredIncreaseOver2Iterations = 0.01; /* checks cross validation convergence */

#ifdef XML_SUPPORT
/** some constants to be used to compare xml strings **/

//databases
static const XMLCh databasesStr[] = {
      chLatin_d, chLatin_a, chLatin_t, chLatin_a, chLatin_b, chLatin_a,
      chLatin_s, chLatin_e, chLatin_s, chNull };
     
      
//calibration
static const XMLCh calibrationStr[] = { chLatin_c, chLatin_a,
      chLatin_l, chLatin_i, chLatin_b, chLatin_r, chLatin_a,
      chLatin_t, chLatin_i, chLatin_o, chLatin_n, chNull };

//proteins
static const XMLCh proteinsStr[] = { chLatin_p, chLatin_r,
      chLatin_o, chLatin_t, chLatin_e, chLatin_i, chLatin_n,
      chLatin_s, chNull };

//protein      
static const XMLCh proteinStr[] = { chLatin_p, chLatin_r,
      chLatin_o, chLatin_t, chLatin_e, chLatin_i, chLatin_n, chNull };
      
//fragSpectrumScan 
static const XMLCh fragSpectrumScanStr[] = { chLatin_f, chLatin_r,
      chLatin_a, chLatin_g, chLatin_S, chLatin_p, chLatin_e,
      chLatin_c, chLatin_t, chLatin_r, chLatin_u, chLatin_m,
      chLatin_S, chLatin_c, chLatin_a, chLatin_n, chNull };   
#endif //XML_SUPPORT      
      
/** some constants to be used to compare xml strings **/
      
Caller::Caller() :
        xmlOutputFN(""), pNorm(NULL), pCheck(NULL), svmInput(NULL), protEstimator(NULL),
        forwardTabInputFN(""), decoyWC(""), resultFN(""), tabFN(""),
        xmlInputFN(""), weightFN(""), tabInput(false), readStdIn(false),
        docFeatures(false), quickValidation(false), reportPerformanceEachIteration(false),
        reportUniquePeptides(true), calculateProteinLevelProb(false),
        schemaValidation(true), hasProteins(false), target_decoy_competition(false),
        test_fdr(0.01), selectionfdr(0.01), selectedCpos(0), selectedCneg(0),
        threshTestRatio(0.3), trainRatio(0.6), niter(10) {

  /*fido parameters*/
  fido_alpha = -1;
  fido_beta = -1;
  fido_gamma = -1;
  fido_nogrouProteins = false; 
  fido_noprune = false;
  fido_noseparate = false;
  fido_reduceTree = false;
  fido_truncate = false;
  fido_trivialGrouping = false;
  fido_depth = 0;
  fido_mse_threshold = 0.1;
  /* general protein probabilities options */
  tiesAsOneProtein = false;
  usePi0 = false;
  outputEmpirQVal = false;  
  decoy_prefix = "random";
}

Caller::~Caller() {
  if (pNorm) {
    delete pNorm;
  }
  pNorm = NULL;
  if (pCheck) {
    delete pCheck;
  }
  pCheck = NULL;
  if (svmInput) {
    delete svmInput;
  }
  svmInput = NULL;
  if (protEstimator) {
    delete protEstimator;
  }
  protEstimator = NULL;
  if(readStdIn) {
    rmdir(xmlInputDir);
    delete xmlInputDir;
  }
  xmlInputDir = NULL;
}

string Caller::extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss << "Started " << ctime(&startTime) << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << " on " << host << endl;
  oss << "Hyperparameters fdr=" << selectionfdr;
  oss << ", Cpos=" << selectedCpos << ", Cneg=" << selectedCneg
      << ", maxNiter=" << niter << endl;
  return oss.str();
}

string Caller::greeter() {
  ostringstream oss;
  oss << "Percolator version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2006-9 University of Washington. All rights reserved.\n"
      << "Written by Lukas Käll (lukall@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n";
  return oss.str();
}

bool Caller::parseOptions(int argc, char **argv) {
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  call = call.substr(0,call.length()-1); // trim ending carriage return
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   percolator [-X pout.xml] [other options] pin.xml" << endl;
  intro << "Where pin.xml is the output file generated by sqt2pin; pout.xml is where" << endl;
  intro << "the output will be written (ensure to have read and write access on the file)." << endl;
  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("X",
      "xmloutput",
      "path to file in xml-output format (pout)",
      "filename");
  cmd.defineOption("e",
      "stdinput",
      "read xml-input format (pin) from standard input",
      "",
      TRUE_IF_SET);
  cmd.defineOption("Z",
      "decoy-xml-output",
      "Include decoys (PSMs, peptides and/or proteins) in the xml-output. Only available if -X is used.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("p",
      "Cpos",
      "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.",
      "value");
  cmd.defineOption("n",
      "Cneg",
      "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified or -p not specified.",
      "value");
  cmd.defineOption("F",
      "trainFDR",
      "False discovery rate threshold to define positive examples in training. Set by cross validation if 0. Default is 0.01.",
      "value");
  cmd.defineOption("t",
      "testFDR",
      "False discovery rate threshold for evaluating best cross validation result and the reported end result. Default is 0.01.",
      "value"); 
  cmd.defineOption("i",
      "maxiter",
      "Maximal number of iterations",
      "number");
  cmd.defineOption("x",
      "quick-validation",
      "Quicker execution by reduced internal cross-validation.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("f",
      "train-ratio",
      "Fraction of the negative data set to be used as train set when only providing one negative set, \
      remaining examples will be used as test set. Set to 0.6 by default.",
      "value");
  cmd.defineOption("J",
      "tab-out",
      "Output the computed features to the given file in tab-delimited format. A file with the features with the given file name will be created",
      "file name");
  cmd.defineOption("j",
      "tab-in",
      "Input files are given as a tab delimited file. In this case the only argument should be a file name \
      of the data file. The tab delimited fields should be id <tab> label <tab> feature1 \
      <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM \
      Labels are interpreted as 1 -- positive set \
      and test set, -1 -- negative set.\
      When the --doc option the first and second feature (third and fourth column) should contain \
      the retention time and difference between observed and calculated mass",
      "filename");
  cmd.defineOption("w",
      "weights",
      "Output final weights to the given file",
      "filename");
  cmd.defineOption("W",
      "init-weights",
      "Read initial weights from the given file (one per line)",
      "filename");
  cmd.defineOption("V",
      "default-direction",
      "The most informative feature given as feature number, can be negated to indicate that a lower value is better.",
      "featureNum");
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("u",
      "unitnorm",
      "Use unit normalization [0-1] instead of standard deviation normalization",
      "",
      TRUE_IF_SET);
  cmd.defineOption("R",
      "test-each-iteration",
      "Measure performance on test set each iteration",
      "",
      TRUE_IF_SET);
  cmd.defineOption("O",
      "override",
      "Override error check and do not fall back on default score vector in case of suspect score vector",
      "",
      TRUE_IF_SET);
  cmd.defineOption("S",
      "seed",
      "Setting seed of the random number generator. Default value is 1",
      "value");
  cmd.defineOption("K",
      "klammer",
      "Retention time features calculated as in Klammer et al.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("D",
      "doc",
      "Include description of correct features.",
      "",
      MAYBE,
      "15");
  cmd.defineOption("r",
      "results",
      "Output tab delimited results to a file instead of stdout",
      "filename");
  cmd.defineOption("B",
      "decoy-results",
      "Output tab delimited results for decoys into a file",
      "filename");
  cmd.defineOption("U",
      "only-psms",
      "Do not remove redundant peptides, keep all PSMS and exclude peptide level probabilities.",
      "",
      FALSE_IF_SET);
  cmd.defineOption("s",
      "no-schema-validation",
      "skip validation of input file against xml schema.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("A",
      "protein",
      "output protein level probabilities",
      "",
      TRUE_IF_SET);
  cmd.defineOption("a",
      "fido-alpha",
      "Probability with which a present protein emits an associated peptide (to be used jointly with the -A option) \
       Set by grid search if not specified.",
      "value");
  cmd.defineOption("b",
      "fido-beta",
      "Probability of the creation of a peptide from noise (to be used jointly with the -A option). Set by grid search if not specified",
      "value");
  cmd.defineOption("G",
      "fido-gamma",
      "Prior probability of that a protein is present in the sample ( to be used with the -A option). Set by grid search if not specified",
      "value");
  cmd.defineOption("g",
      "allow-protein-group",
      "treat ties as if it were one protein (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("I",
      "protein-level-pi0",
      "use pi_0 value when calculating empirical q-values (no effect if option Q is activated) (Only valid if option -A is active).",
      "", 
      TRUE_IF_SET);
  cmd.defineOption("q",
      "empirical-protein-q", 		   
      "output empirical q-values and p-values (from target-decoy analysis) (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("N",
      "fido-no-group-proteins", 		   
      "disactivates the grouping of proteins with similar connectivity, \
       for example if proteins P1 and P2 have the same peptides matching both of them, P1 and P2 will not be grouped as one protein \
       (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("E",
      "fido-no-separate-proteins", 		   
      "Proteins graph will not be separated in sub-graphs (Only valid if option -A is active).",
      "",
      TRUE_IF_SET); 
  cmd.defineOption("C",
      "fido-no-prune-proteins", 		   
      "it does not prune peptides with a very low score (~0.0) which means that if a peptide with a very low score is matching two proteins,\
       when we prune the peptide,it will be duplicated to generate two new protein groups (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("d",
      "fido-gridsearch-depth",
      "Setting depth 0 or 1 or 2 from low depth to high depth(less computational time) \
       of the grid search for the estimation Alpha,Beta and Gamma parameters for fido(Only valid if option -A is active). Default value is 0",
      "value");
  cmd.defineOption("P",
      "pattern",
      "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that idenfifies the decoys in the database \
       is not the default (by default : ramdom) (Only valid if option -A  is active).",
      "value");
  cmd.defineOption("T",
      "fido-reduce-tree-in-gridsearch",
      "Reduce the tree of proteins (removing low scored proteins) in order to estimate alpha,beta and gamma faster.(Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("Y",
      "post-processing-tdcn",
      "Use target decoy competition to compute peptide probabilities.(recommended when using -A).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("H",
      "grid-search-mse-threshold",
      "Q-value threshold that will be used in the computation of the MSE and ROC AUC score in the grid search (recommended 0.05 for normal size datasets and 0.1 for big size datasets).(Only valid if option -A is active).",
      "",
      "value");
  cmd.defineOption("W",
      "fido-truncation",
      "Proteins with a very low score (< 0.001) will be truncated (assigned 0.0 probability).(Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("Q",
      "fido-protein-group-level-inference",
      "Uses protein group level inference, each cluster of proteins is either present or not, therefore when grouping proteins discard all possible combinations for each group.(Only valid if option -A is active and -N is inactive).",
      "",
      TRUE_IF_SET);
  
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  // now query the parsing results
  if (cmd.optionSet("X")) xmlOutputFN = cmd.options["X"];
  
  if (cmd.optionSet("B")) {
    decoyOut = cmd.options["B"];
  }
  
  if (cmd.optionSet("r")) {
    resultFN = cmd.options["r"];
  }
  
  if (cmd.optionSet("U")) {
    if (cmd.optionSet("A")){
      cerr
      << "The -U option cannot be used in conjunction with -A: peptide level statistics\n"
      << "are needed to calculate protein level ones.";
      return 0;
    }
    reportUniquePeptides = false;
  }

  if (cmd.optionSet("A")) {
  
    calculateProteinLevelProb = true;
    tiesAsOneProtein = cmd.optionSet("g");
    usePi0 = cmd.optionSet("I");
    outputEmpirQVal = cmd.optionSet("q");
    fido_nogrouProteins = cmd.optionSet("N"); 
    fido_noprune = cmd.optionSet("C");
    fido_noseparate = cmd.optionSet("E");
    fido_reduceTree = cmd.optionSet("T");
    fido_truncate = cmd.optionSet("W");
    fido_trivialGrouping = cmd.optionSet("Q");
    if (cmd.optionSet("P"))  decoy_prefix = cmd.options["P"];
    if (cmd.optionSet("d"))  fido_depth = cmd.getInt("d", 0, 2);
    if (cmd.optionSet("a"))  fido_alpha = cmd.getDouble("a", 0.00, 1.0);
    if (cmd.optionSet("b"))  fido_beta = cmd.getDouble("b", 0.00, 1.0);
    if (cmd.optionSet("G"))  fido_gamma = cmd.getDouble("G", 0.00, 1.0);
    if (cmd.optionSet("H"))  fido_mse_threshold = cmd.getDouble("H",0.001,1.0);

  }
  
  if (cmd.optionSet("e")) {

    readStdIn = true;
    string str = "";
    try
    {
      boost::filesystem::path ph = boost::filesystem::unique_path();
      boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
      boost::filesystem::path file("pin-tmp.xml");
      xmlInputFN = std::string((dir / file).string()); 
      str =  dir.string();
      xmlInputDir = new char[str.size() + 1];
      std::copy(str.begin(), str.end(), xmlInputDir);
      xmlInputDir[str.size()] = '\0';
      if (boost::filesystem::is_directory(dir)) {
	      boost::filesystem::remove_all(dir);
      }
      boost::filesystem::create_directory(dir);
    } 
    catch (boost::filesystem::filesystem_error &e)
    {
      std::cerr << e.what() << std::endl;
      return 0;
    }
  }
  
  if (cmd.optionSet("p")) {
    selectedCpos = cmd.getDouble("p", 0.0, 1e127);
  }
  if (cmd.optionSet("n")) {
    selectedCneg = cmd.getDouble("n", 0.0, 1e127);
    if(selectedCpos == 0)
    {
      std::cerr << "Warning : the positive penalty(cpos) is 0, therefore both the "  
		 << "positive and negative penalties are going "
		 << "to be cros-validated. The option --Cneg has to be used together "
		 << "with the option --Cpos" << std::endl;
    }
  }
  if (cmd.optionSet("J")) {
    tabFN = cmd.options["J"];
  }
  if (cmd.optionSet("j")) {
    tabInput = true;
    forwardTabInputFN = cmd.options["j"];
  }
  if (cmd.optionSet("w")) {
    weightFN = cmd.options["w"];
  }
  if (cmd.optionSet("W")) {
    SanityCheck::setInitWeightFN(cmd.options["W"]);
  }
  if (cmd.optionSet("V")) {
    SanityCheck::setInitDefaultDir(cmd.getInt("V", -100, 100));
  }
  if (cmd.optionSet("f")) {
    double frac = cmd.getDouble("f", 0.0, 1.0);
    trainRatio = frac;
  }
  if (cmd.optionSet("u")) {
    Normalizer::setType(Normalizer::UNI);
  }
  if (cmd.optionSet("O")) {
    SanityCheck::setOverrule(true);
  }
  if (cmd.optionSet("R")) {
    reportPerformanceEachIteration = true;
  }
  if (cmd.optionSet("x")) {
    quickValidation=true;
  }
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (cmd.optionSet("F")) {
    selectionfdr = cmd.getDouble("F", 0.0, 1.0);
  }
  if (cmd.optionSet("t")) {
    test_fdr = cmd.getDouble("t", 0.0, 1.0);
  }
  if (cmd.optionSet("S")) {
    Scores::setSeed(cmd.getInt("S", 1, 20000));
  }
  if (cmd.optionSet("K")) {
    DescriptionOfCorrect::setKlammer(true);
  }
  if (cmd.optionSet("D")) {
    docFeatures = true;
    DataSet::setCalcDoc(true);
    DescriptionOfCorrect::setDocType(cmd.getInt("D", 0, 15));
  }
  if (cmd.optionSet("Z")) {
    Scores::setOutXmlDecoys(true);
  }
  if (cmd.optionSet("s")) {
    schemaValidation = false;
  }
  showExpMass = true;
  Scores::setShowExpMass(showExpMass);
  if (cmd.optionSet("Y")) {
    target_decoy_competition = true; 
  }
  // if there are no arguments left...
  if (cmd.arguments.size() == 0) {
    if(!cmd.optionSet("j") && !cmd.optionSet("e") ){ // unless the input comes from -j option or -e option
      cerr << "Error: too few arguments.";
      cerr << "\nInvoke with -h option for help\n";
      return 0; // ...error
    }
  }
  // if there is one argument left...
  if (cmd.arguments.size() == 1) {
#ifdef XML_SUPPORT 
    xmlInputFN = cmd.arguments[0]; // then it's the pin input
    if(cmd.optionSet("j")){ // and if the tab input is also present
      cerr << "Error: use one of either pin or tab-delimited input format.";
      cerr << "\nInvoke with -h option for help.\n";
      return 0; // ...error
    }
    if(cmd.optionSet("e")){ // if stdin pin file is present
      cerr << "Error: the pin file has already been given as stdinput argument.";
      cerr << "\nInvoke with -h option for help.\n";
      return 0; // ...error
    }
#else //XML_SUPPORT
    tabInput = true;
    forwardTabInputFN = cmd.arguments[0]; // then it's the tab input
#endif //XML_SUPPORT
  }
  // if there is more then one argument left...
  if (cmd.arguments.size() > 1) {
    cerr << "Error: too many arguments.";
    cerr << "\nInvoke with -h option for help\n";
    return 0; // ...error
  }

  return true;
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w normal vector
 */
void Caller::printWeights(ostream & weightStream, vector<double>& w) {
  weightStream
  << "# first line contains normalized weights, second line the raw weights"
  << endl;
  weightStream << DataSet::getFeatureNames().getFeatureNames() << "\tm0"
      << endl;
  weightStream.precision(3);
  weightStream << w[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << w[ix];
  }
  weightStream << endl;
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(w, ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

/**
 * Reads in the files from XML (must be enabled at compile time) or tab format
 */
int Caller::readFiles() { 
  if (xmlInputFN.size() != 0) {    
#ifdef XML_SUPPORT 
    xercesc::XMLPlatformUtils::Initialize();
    
    DataSet * targetSet = new DataSet();
    assert(targetSet);
    targetSet->setLabel(1);
    DataSet * decoySet = new DataSet();
    assert(decoySet);
    decoySet->setLabel(-1);
    
    try {
      
      namespace xml = xsd::cxx::xml;
      std::ifstream xmlInStream;
      xmlInStream.exceptions(ifstream::badbit | ifstream::failbit);
      xmlInStream.open(xmlInputFN.c_str());

      string schemaDefinition= Globals::getInstance()->getXMLDir()+PIN_SCHEMA_LOCATION+string("percolator_in.xsd");
      string schema_major = PIN_VERSION_MAJOR;
      string schema_minor = PIN_VERSION_MINOR;
      parser p;
      xml_schema::dom::auto_ptr<xercesc::DOMDocument> doc(p.start(
          xmlInStream, xmlInputFN.c_str(), Caller::schemaValidation,
          schemaDefinition, schema_major, schema_minor));

      doc = p.next();
      // read enzyme element
      // the enzyme element is a subelement but CodeSynthesis Xsd does not
      // generate a class for it. (I am trying to find a command line option
      // that overrides this decision). As for now special treatment is needed
      char* value = XMLString::transcode(doc->getDocumentElement()->getTextContent());
      
      if(VERB > 1) std::cerr << "enzyme=" << value << std::endl;
      
      Enzyme::setEnzyme(value);
      XMLString::release(&value);
      doc = p.next();

      //checking if database is present to jump it
      if(XMLString::equals(databasesStr, doc->getDocumentElement()->getTagName()))
      {
        //NOTE I dont really need this info, do I? good to have it though
        /*
          std::auto_ptr< ::percolatorInNs::databases > 
	      databases( new ::percolatorInNs::databases(*doc->getDocumentElement()));
        */
        doc = p.next();
        Caller::hasProteins = true;
      }
      
      // read process_info element
      percolatorInNs::process_info
      processInfo(*doc->getDocumentElement());
      otherCall = processInfo.command_line();
      doc = p.next();


      if (XMLString::equals(calibrationStr,doc->getDocumentElement()->getTagName())) 
      {
	//NOTE the calibration should define the initial direction
        //percolatorInNs::calibration calibration(*doc->getDocumentElement());
        doc = p.next();
      };

      percolatorInNs::featureDescriptions featureDescriptions(*doc->getDocumentElement());

      //I want to get the initial values that are present in feature descriptions
      vector<double> init_values;
      for( const auto & descr : featureDescriptions.featureDescription() ) {
          if(descr.initialValue().present()){
              if(VERB >2){
                  std::cerr << "Initial direction for " << descr.name() << " is " << descr.initialValue().get() << std::endl;
              }
              init_values.push_back(descr.initialValue().get());
          }
      }
      
      FeatureNames& feNames = DataSet::getFeatureNames();
      feNames.setFromXml(featureDescriptions, docFeatures);
      targetSet->initFeatureTables(feNames.getNumFeatures(), docFeatures);
      decoySet->initFeatureTables(feNames.getNumFeatures(), docFeatures);

      // import info from xml: read Fragment Spectrum Scans
      for (doc = p.next(); doc.get()!= 0 && 
	XMLString::equals(fragSpectrumScanStr, doc->getDocumentElement()->getTagName()); doc = p.next()) 
      {
        percolatorInNs::fragSpectrumScan fragSpectrumScan(*doc->getDocumentElement());
	      for (const auto &psm : fragSpectrumScan.peptideSpectrumMatch())
	      {
	        if(psm.isDecoy())
	        {
	          decoySet->readPsm(psm,fragSpectrumScan.scanNumber());
	        }
	        else
	        {
	          targetSet->readPsm(psm,fragSpectrumScan.scanNumber());
	        }
	      }
      }

      // import info from xml: read database proteins
      // only read them if they are present and the option of using mayusfdr is activated
      unsigned readProteins = 0;
      for (doc = p.next(); doc.get()!= 0 
	&& Caller::hasProteins && Caller::calculateProteinLevelProb /*&& Caller::protEstimator->getMayuFdr()*/
	&& XMLString::equals(proteinStr, doc->getDocumentElement()->getTagName()); doc = p.next()) 
      {
        std::auto_ptr< ::percolatorInNs::protein > protein( new ::percolatorInNs::protein(*doc->getDocumentElement()));
        protEstimator->addProteinDb(*protein);
        ++readProteins;
      }
      
      /*if(Caller::calculateProteinLevelProb && Caller::protEstimator->getMayuFdr() && readProteins <= 0)
      {
	std::cerr << "Warning : options -Q and -A are activated but the number of proteins found in the input file is zero.\n\
		       Did you run converters with the flag -F ?\n" << std::endl;
	Caller::protEstimator->setMayusFDR(false);
      }*/
      
      //maybe better to do :
      //SanityCheck::addDefaultWeights(init_values);
      pCheck = SanityCheck::initialize(otherCall);
      assert(pCheck);
      pCheck->addDefaultWeights(init_values);
      setHandler.push_back_dataset(targetSet);
      setHandler.push_back_dataset(decoySet);
      xmlInStream.close();
    }

    catch (const xml_schema::exception& e) {
      std::cerr << e << endl;
      return 0;
    } catch (const std::ios_base::failure&) {
      std::cerr << "unable to open or read failure" << std::endl;
      return 0;
    } catch (const xercesc::DOMException& e) {
      char * tmpStr = XMLString::transcode(e.getMessage());
      std::cerr << "catch  xercesc::DOMException=" << tmpStr << std::endl;
      XMLString::release(&tmpStr);
      return 0;
    }
#else //XML_SUPPORT
    std::cerr << "Warning: Compiler flag XML_SUPPORT was off, trying to process input as tab delimited file" << std::endl;
    pCheck = new SanityCheck();
    setHandler.readTab(forwardTabInputFN);
    std::cerr << "Features:\n" << DataSet::getFeatureNames().getFeatureNames() << std::endl;
#endif //XML_SUPPORT
  } else if (tabInput) {
    pCheck = new SanityCheck();
    setHandler.readTab(forwardTabInputFN);
    std::cerr << "Features:\n" << DataSet::getFeatureNames().getFeatureNames() << std::endl;
  } 
  return true;
}

/** 
 * Train one of the crossvalidation bins 
 * @param set identification number of the bin that is processed
 * @param w list of normal vectors (in the linear algebra sense) of the hyperplane from SVM, one for each bin
 * @param updateDOC boolean deciding to calculate retention features @see DescriptionOfCorrect
 * @param cpos_vec vector with soft margin parameter for positives
 * @param cfrac_vec vector with soft margin parameter for fraction negatives / positives
 * @param best_cpos best soft margin parameter for positives
 * @param best_cfrac best soft margin parameter for fraction negatives / positives
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
*/
int Caller::xv_process_one_bin(unsigned int set, vector<vector<double> >& w, bool updateDOC, vector<double>& cpos_vec, 
                               vector<double>& cfrac_vec, double &best_cpos, double &best_cfrac, vector_double* pWeights,
                               options * pOptions) {
  int bestTP = 0;
  if (VERB > 2) {
    cerr << "cross validation - fold " << set + 1 << " out of "
         << xval_fold << endl;
  }
  
  vector<double> ww = w[set]; // normal vector initial guess and result holder
  vector<double> bestW = w[set]; // normal vector with highest true positive estimate
  xv_train[set].calcScores(ww, selectionfdr);
  if (docFeatures && updateDOC) {
    xv_train[set].recalculateDescriptionOfGood(selectionfdr);
  }
  xv_train[set].generateNegativeTrainingSet(*svmInput, 1.0);
  xv_train[set].generatePositiveTrainingSet(*svmInput, selectionfdr, 1.0);
  if (VERB > 2) {
    cerr << "Calling with " << svmInput->positives << " positives and "
         << svmInput->negatives << " negatives\n";
  }
  
  // Create storage vector for SVM algorithm
  struct vector_double* Outputs = new vector_double;
  Outputs->vec = new double[svmInput->positives + svmInput->negatives];
  Outputs->d = svmInput->positives + svmInput->negatives;
  
  // Find combination of soft margin parameters with highest estimate of true positives
  for (const auto cpos : cpos_vec) {
    for (const auto cfrac : cfrac_vec) {
      if (VERB > 2) cerr << "-cross validation with cpos=" << cpos
          << ", cfrac=" << cfrac << endl;
      int tp = 0;
      for (int ix = 0; ix < pWeights->d; ix++) {
        pWeights->vec[ix] = 0;
      }
      for (int ix = 0; ix < Outputs->d; ix++) {
        Outputs->vec[ix] = 0;
      }
      svmInput->setCost(cpos, (cpos) * (cfrac));
      
      // Call SVM algorithm (see ssl.cpp)
      L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs);
      
      for (int i = FeatureNames::getNumFeatures() + 1; i--;) {
        ww[i] = pWeights->vec[i];
      }
      tp = xv_train[set].calcScores(ww, test_fdr);
      if (VERB > 2) {
        cerr << "- cross validation estimates " << tp
             << " target PSMs over " << test_fdr * 100 << "% FDR level"
             << endl;
      }
      if (tp >= bestTP) {
        if (VERB > 2) {
          cerr << "Better than previous result, store this" << endl;
        }
        bestTP = tp;
        bestW = ww;
        best_cpos = cpos;
        best_cfrac = cfrac;
      }
    }
    if (VERB > 2) cerr << "cross validation estimates " << bestTP
        / (xval_fold - 1) << " target PSMs with q<" << test_fdr
        << " for hyperparameters Cpos=" << best_cpos << ", Cneg="
        << best_cfrac * best_cpos << endl;
  }
  w[set]=bestW;
  delete[] Outputs->vec;
  delete Outputs;
  return bestTP;
}

/** 
 * Executes a cross validation step
 * @param w list of the bins' normal vectors (in linear algebra sense) of the hyperplane from SVM
 * @param updateDOC boolean deciding to calculate retention features @see DescriptionOfCorrect
 * @return Estimation of number of true positives
 */
int Caller::xv_step(vector<vector<double> >& w, bool updateDOC) {
  // Setup
  struct options* pOptions = new options;
  pOptions->lambda = 1.0;
  pOptions->lambda_u = 1.0;
  pOptions->epsilon = EPSILON;
  pOptions->cgitermax = CGITERMAX;
  pOptions->mfnitermax = MFNITERMAX;
  struct vector_double* pWeights = new vector_double;
  pWeights->d = FeatureNames::getNumFeatures() + 1;
  pWeights->vec = new double[pWeights->d];
  int estTP = 0;
  double best_cpos = 1, best_cfrac = 1;
  if (!quickValidation) {
    for (unsigned int set = 0; set < xval_fold; ++set) {
      estTP += xv_process_one_bin(set,w,updateDOC, xv_cposs, xv_cfracs, best_cpos, best_cfrac, pWeights, pOptions);   
    }
  } else {
    // Use limited internal cross validation, i.e take the cpos and cfrac values of the first bin 
    // and use it for the subsequent bins 
    estTP += xv_process_one_bin(0,w,updateDOC, xv_cposs, xv_cfracs, best_cpos, best_cfrac, pWeights, pOptions);
    vector<double> cp(1),cf(1);
    cp[0]=best_cpos; cf[0]= best_cfrac;
    for (unsigned int set = 1; set < xval_fold; ++set) {
      estTP += xv_process_one_bin(set,w,updateDOC, cp, cf, best_cpos, best_cfrac, pWeights, pOptions);   
    }
  }
  delete[] pWeights->vec;
  delete pWeights;
  delete pOptions;
  return estTP / (xval_fold - 1);
}

/** 
 * Train the SVM using several cross validation iterations
 * @param w list of normal vectors
 */
void Caller::train(vector<vector<double> >& w) {
  // iterate
  int foundPositivesOldOld=0, foundPositivesOld=0, foundPositives=0; 
  for (unsigned int i = 0; i < niter; i++) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << " :\t";
    }
    foundPositives = xv_step(w, true);
    if (VERB > 1) {
      cerr << "After the iteration step, " << foundPositives
          << " target PSMs with q<" << selectionfdr
          << " were estimated by cross validation" << endl;
    }
    if (VERB > 2) {
      cerr << "Obtained weights" << endl;
      for (size_t set = 0; set < xval_fold; ++set) {
        printWeights(cerr, w[set]);
      }
    }
    if (foundPositives>0 && foundPositivesOldOld>0 && quickValidation) {
      if ((double)(foundPositives-foundPositivesOldOld)<=(foundPositivesOldOld*requiredIncreaseOver2Iterations)) {
        if (VERB > 1) {
          cerr << "Performance increase over the last two iterations indicate that the algorithm has converged" << endl;
          cerr << "(" << foundPositives << " vs " << foundPositivesOldOld << ")" << endl;
        }
        break;
      }
    }    
    foundPositivesOldOld=foundPositivesOld;    
    foundPositivesOld=foundPositives;
  }
  if (VERB == 2) {
    cerr
    << "Obtained weights (only showing weights of first cross validation set)"
    << endl;
    printWeights(cerr, w[0]);
  }
  foundPositives = 0;
  for (size_t set = 0; set < xval_fold; ++set) {
    if (docFeatures) {
      xv_test[set].getDOC().copyDOCparameters(xv_train[set].getDOC());
      xv_test[set].setDOCFeatures();
    }
    foundPositives += xv_test[set].calcScores(w[set], test_fdr);
  }
  if (VERB > 0) {
    cerr << "After all training done, " << foundPositives << " target PSMs with q<"
        << test_fdr << " were found when measuring on the test set"
        << endl;
  }  
}

/** 
 * Fills in the features previously read from file and normalizes them
 */
void Caller::fillFeatureSets() {
  fullset.fillFeatures(setHandler, reportUniquePeptides);
  if (VERB > 1) {
    cerr << "Train/test set contains " << fullset.posSize()
        << " positives and " << fullset.negSize()
        << " negatives, size ratio=" << fullset.getTargetDecoySizeRatio()
        << " and pi0=" << fullset.getPi0() << endl;
  }
  
  //check for the minimum recommended number of positive and negative hits
  if(fullset.posSize() <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of positive samples read is too small to perform a correct clasification.\n" << std::endl;
  }
  if(fullset.negSize() <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of negative samples read is too small to perform a correct clasification.\n" << std::endl;
  }
  
  //Normalize features
  if (docFeatures) {
    for (auto &subset : setHandler.getSubsets()) {
      subset->setRetentionTime(scan2rt);
    }
  }
  if (tabFN.length() > 0) {
    setHandler.writeTab(tabFN);
  }
  vector<double*> featuresV, rtFeaturesV;
  PSMDescription* pPSM;
  for (auto &subset : setHandler.getSubsets()) {
    subset->fillFeatures(featuresV);
    subset->fillRtFeatures(rtFeaturesV);
  }
  pNorm = Normalizer::getNormalizer();

  pNorm->setSet(featuresV,
      rtFeaturesV,
      FeatureNames::getNumFeatures(),
      docFeatures ? RTModel::totalNumRTFeatures() : 0);
  pNorm->normalizeSet(featuresV, rtFeaturesV);
}

/** 
 * Sets up the SVM classifier: 
 * - divide dataset into training and test sets for each fold
 * - set parameters (fdr, soft margin)
 * @param w list of normal vectors
 * @return number of positives for initial setup
 */
int Caller::preIterationSetup(vector<vector<double> >& w) {
  
  svmInput = new AlgIn(fullset.size(), FeatureNames::getNumFeatures() + 1); // One input set, to be reused multiple times
  assert( svmInput );

  if (selectedCpos >= 0 && selectedCneg >= 0) {
    xv_train.resize(xval_fold);
    xv_test.resize(xval_fold);
    if(xmlInputFN.size() > 0){
    	// take advantage of spectrum information in xml input
    	fullset.createXvalSetsBySpectrum(xv_train, xv_test, xval_fold);
    } else {
    	fullset.createXvalSets(xv_train, xv_test, xval_fold);
    }

    if (selectionfdr <= 0.0) {
      selectionfdr = test_fdr;
    }
    if (selectedCpos > 0) {
      xv_cposs.push_back(selectedCpos);
    } else {
      xv_cposs.push_back(10);
      xv_cposs.push_back(1);
      xv_cposs.push_back(0.1);
      if (VERB > 0) {
        cerr << "selecting cpos by cross validation" << endl;
      }
    }
    if (selectedCpos > 0 && selectedCneg > 0) {
      xv_cfracs.push_back(selectedCneg / selectedCpos);
    } else {
      xv_cfracs.push_back(1.0 * fullset.getTargetDecoySizeRatio());
      xv_cfracs.push_back(3.0 * fullset.getTargetDecoySizeRatio());
      xv_cfracs.push_back(10.0 * fullset.getTargetDecoySizeRatio());
      if (VERB > 0) {
        cerr << "selecting cneg by cross validation" << endl;
      }
    }
    return pCheck->getInitDirection(xv_test, xv_train, pNorm, w, test_fdr);
  } else {
    vector<Scores> myset(1, fullset);
    cerr << "B" << endl;
    return pCheck->getInitDirection(myset, myset, pNorm, w, test_fdr);
  }
}

/** 
 * Subroutine of @see Caller::writeXML() for PSM output
 */
void Caller::writeXML_PSMs() {
  ofstream os;
  xmlOutputFN_PSMs = xmlOutputFN;
  xmlOutputFN_PSMs.append("writeXML_PSMs");
  os.open(xmlOutputFN_PSMs.c_str(), ios::out);

  os << "  <psms>" << endl;
  for (vector<ScoreHolder>::iterator psm = fullset.begin();
      psm != fullset.end(); ++psm) {
      os << *psm;
  }
  os << "  </psms>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see Caller::writeXML() for peptide output
 */
void Caller::writeXML_Peptides() {
  ofstream os;
  xmlOutputFN_Peptides = xmlOutputFN;
  xmlOutputFN_Peptides.append("writeXML_Peptides");
  os.open(xmlOutputFN_Peptides.c_str(), ios::out);
  // append PEPTIDEs
  os << "  <peptides>" << endl;
  for (vector<ScoreHolder>::iterator psm = fullset.begin(); psm
  != fullset.end(); ++psm) {
    os << (ScoreHolderPeptide)*psm;
  }
  os << "  </peptides>" << endl << endl;
  os.close();
}

/** 
 * Subroutine of @see Caller::writeXML() for protein output
 */
void Caller::writeXML_Proteins() {
  xmlOutputFN_Proteins = xmlOutputFN;
  xmlOutputFN_Proteins.append("writeXML_Proteins");
  protEstimator->writeOutputToXML(xmlOutputFN_Proteins, Scores::isOutXmlDecoys());
}

/** 
 * Writes the output of percolator to an pout XML file
 */
void Caller::writeXML(){
  ofstream os;
  const string space = PERCOLATOR_OUT_NAMESPACE;
  string schema_major = POUT_VERSION_MAJOR;
  string schema_minor = POUT_VERSION_MINOR;
  const string schema = space +
      " https://github.com/percolator/percolator/raw/pout-" + schema_major +
      "-" + schema_minor + "/src/xml/percolator_out.xsd";
  os.open(xmlOutputFN.data(), ios::out | ios::binary);
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  os << "<percolator_output "
      << endl << "xmlns=\""<< space << "\" "
      << endl << "xmlns:p=\""<< space << "\" "
      << endl << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
      << endl << "xsi:schemaLocation=\""<< schema <<"\" "
      << endl << "p:majorVersion=\"" << VERSION_MAJOR << "\" p:minorVersion=\""
      << VERSION_MINOR << "\" p:percolator_version=\"Percolator version "
      << VERSION << "\">\n"<< endl;
  os << "  <process_info>" << endl;
  os << "    <command_line>" << call << "</command_line>" << endl;

  os << "    <other_command_line>" << otherCall << "</other_command_line>\n";
  os << "    <pi_0_psms>" << pi_0_psms << "</pi_0_psms>" << endl;
  if(reportUniquePeptides)
    os << "    <pi_0_peptides>" << pi_0_peptides << "</pi_0_peptides>" << endl;
  if(calculateProteinLevelProb) {  
    if(usePi0)
      os << "    <pi_0_proteins>" << protEstimator->getPi0() << "</pi_0_proteins>" << endl;
    /*if(protEstimator->getMayuFdr())
      os << "    <fdr_proteins>" << protEstimator->getFDR() << "</fdr_proteins>" << endl;*/
    os << "    <alpha>" << protEstimator->getAlpha() <<"</alpha>" << endl;
    os << "    <beta>"  << protEstimator->getBeta() <<"</beta>" << endl;
    os << "    <gamma>" << protEstimator->getGamma() <<"</gamma>" << endl;
  }
  os << "    <psms_qlevel>" <<  numberQpsms <<"</psms_qlevel>" << endl;
  if(reportUniquePeptides)
    os << "    <peptides_qlevel>" << fullset.getQvaluesBelowLevel(0.01) << "</peptides_qlevel>" << endl;
  if(calculateProteinLevelProb)
    os << "    <proteins_qlevel>" << protEstimator->getQvaluesBelowLevel(0.01) << "</proteins_qlevel>" << endl;  
  if (docFeatures) {
    os << "    <average_delta_mass>" << fullset.getDOC().getAvgDeltaMass()
                   << "</average_delta_mass>" << endl;
    os << "    <average_pi>" << fullset.getDOC().getAvgPI()
                   << "</average_pi>" << endl;
  }
  os << "  </process_info>" << endl << endl;

  // apppend PSMs
  ifstream ifs_psms(xmlOutputFN_PSMs.data(), ios::in | ios::binary);
  os << ifs_psms.rdbuf();
  ifs_psms.close();
  remove(xmlOutputFN_PSMs.c_str());
  // append Peptides
  if(reportUniquePeptides){
    ifstream ifs_peptides(xmlOutputFN_Peptides.data(), ios::in | ios::binary);
    os << ifs_peptides.rdbuf();
    ifs_peptides.close();
    remove(xmlOutputFN_Peptides.c_str());
  }
  // append Proteins
  if(calculateProteinLevelProb){
    ifstream ifs_proteins(xmlOutputFN_Proteins.data(), ios::in | ios::binary);
    os << ifs_proteins.rdbuf();
    ifs_proteins.close();
    remove(xmlOutputFN_Proteins.c_str());
  }

  os << "</percolator_output>" << endl;
  os.close();
}

/** Calculates the PSM and/or peptide probabilities
 * @param isUniquePeptideRun boolean indicating if we want peptide or PSM probabilities
 * @param procStart clock time when process started
 * @param procStartClock clock associated with procStart
 * @param w list of normal vectors
 * @param diff runtime of the calculations
 * @param TDC boolean for target decoy competition
 */
void Caller::calculatePSMProb(bool isUniquePeptideRun,Scores *fullset, time_t& procStart,
    clock_t& procStartClock, vector<vector<double> >& w, double& diff, bool TDC){
  // write output (cerr or xml) if this is the unique peptide run and the
  // reportUniquePeptides option was switched on OR if this is not the unique
  // peptide run and the option was switched off
  bool writeOutput = (isUniquePeptideRun == reportUniquePeptides);
  
  if (reportUniquePeptides && VERB > 0 && writeOutput) {
    cerr << "Tossing out \"redundant\" PSMs keeping only the best scoring PSM "
        "for each unique peptide." << endl;
  }
  
  if (isUniquePeptideRun) {
    fullset->weedOutRedundant();
  } else {
    fullset->merge(xv_test, selectionfdr);
    if (TDC) {
      fullset->weedOutRedundantTDC();
	    if(VERB > 0) {
	      std::cerr << "Target Decoy Competition yielded " << fullset->posSize() 
	        << " target PSMs and " << fullset->negSize() << " decoy PSMs" << std::endl;
	    }
    }
  }
  
  if (VERB > 0 && writeOutput) {
    std::cerr << "Selecting pi_0=" << fullset->getPi0() << endl;
  }
  if (VERB > 0 && writeOutput) {
    cerr << "Calibrating statistics - calculating q values" << endl;
  }
  int foundPSMs = fullset->calcQ(test_fdr);
  fullset->calcPep();
  if (VERB > 0 && docFeatures && writeOutput) {
    cerr << "For the cross validation sets the average deltaMass are ";
    for (size_t ix = 0; ix < xv_test.size(); ix++) {
      cerr << xv_test[ix].getDOC().getAvgDeltaMass() << " ";
    }
    cerr << "and average pI are ";
    for (size_t ix = 0; ix < xv_test.size(); ix++) {
      cerr << xv_test[ix].getDOC().getAvgPI() << " ";
    }
    cerr << endl;
  }
  if (VERB > 0 && writeOutput) {
    cerr << "New pi_0 estimate on merged list gives " << foundPSMs
        << (reportUniquePeptides ? " peptides" : " PSMs") << " over q="
        << test_fdr << endl;
  }
  if (VERB > 0 && writeOutput) {
    cerr
    << "Calibrating statistics - calculating Posterior error probabilities (PEPs)"
    << endl;
  }
  time_t end;
  time(&end);
  diff = difftime(end, procStart);
  ostringstream timerValues;
  timerValues.precision(4);
  timerValues << "Processing took "
      << ((double)(clock() - procStartClock)) / (double)CLOCKS_PER_SEC;
  timerValues << " cpu seconds or " << diff << " seconds wall time"
      << endl;
  if (VERB > 1 && writeOutput) {
    cerr << timerValues.str();
  }
  if (weightFN.size() > 0) {
    ofstream weightStream(weightFN.data(), ios::out);
    for (unsigned int ix = 0; ix < xval_fold; ++ix) {
      printWeights(weightStream, w[ix]);
    }
    weightStream.close();
  }
  if (resultFN.empty() && writeOutput) {
    setHandler.print(*fullset, NORMAL);
  } else if (!resultFN.empty()) {
    if (writeOutput) {
      ofstream targetStream((resultFN+(reportUniquePeptides ? ".peptides" : ".psms")).data(), ios::out);
      setHandler.print(*fullset, NORMAL, targetStream);
      targetStream.close();
    } else {
      ofstream targetStream((resultFN+".psms").data(), ios::out);
      setHandler.print(*fullset, NORMAL, targetStream);
      targetStream.close();
    }
  }
  if (!decoyOut.empty() && writeOutput) {
    ofstream decoyStream((decoyOut+(reportUniquePeptides ? ".peptides" : ".psms")).data(), ios::out);
    setHandler.print(*fullset, SHUFFLED, decoyStream);
    decoyStream.close();
  } else if(!decoyOut.empty()) {
    ofstream decoyStream((decoyOut+".psms").data(), ios::out);
    setHandler.print(*fullset, SHUFFLED, decoyStream);
    decoyStream.close();
  }
  // set pi_0 value (to be outputted)
  if(isUniquePeptideRun) {
    pi_0_peptides = fullset->getPi0();
  } else {
    pi_0_psms = fullset->getPi0();
    numberQpsms = fullset->getQvaluesBelowLevel(0.01);
  }
}

/** Calculates the protein probabilites by calling Fido and directly writes the results to XML
 */
void Caller::calculateProteinProbabilitiesFido() {
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();  


  protEstimator = new FidoInterface(fido_alpha,fido_beta,fido_gamma,fido_nogrouProteins,fido_noseparate,
				      fido_noprune,fido_depth,fido_reduceTree,fido_truncate,fido_mse_threshold,
				      tiesAsOneProtein,usePi0,outputEmpirQVal,decoy_prefix,fido_trivialGrouping);
  
  if (VERB > 0) {
    cerr << "\nCalculating protein level probabilities with Fido\n";
    cerr << protEstimator->printCopyright();
  }
  
  protEstimator->initialize(&fullset);
  protEstimator->run();
  protEstimator->computeProbabilities();
  protEstimator->computeStatistics();
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff_time = difftime(procStart, startTime);
  
  if (VERB > 1) {  
    cerr << "Estimating Protein Probabilities took : "
    << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
    << " cpu seconds or " << diff_time << " seconds wall time" << endl;
  }
  
  protEstimator->printOut(resultFN,decoyOut);
  if (xmlOutputFN.size() > 0) {
      writeXML_Proteins();
  }
}

/** 
 * Executes the flow of the percolator process:
 * 1. reads in the input file
 * 2. trains the SVM
 * 3. calculate PSM probabilities
 * 4. (optional) calculate peptide probabilities
 * 5. (optional) calculate protein probabilities
 */
int Caller::run() {  

  time(&startTime);
  startClock = clock();
  if (VERB > 0) {
    cerr << extendedGreeter();
  }
  // populate tmp input file with cin information if option is enabled
  if(readStdIn){
    ofstream tmpInputFile;
    tmpInputFile.open(xmlInputFN.c_str());
    while(cin) {
      char buffer[1000];
      cin.getline(buffer, 1000);
      tmpInputFile << buffer << endl;
    }
    tmpInputFile.close();
  }
  
  // Reading input files (pin or temporary file)
  if(!readFiles()) return 0;
  // Copy feature data to Scores object
  fillFeatureSets();

#ifdef XML_SUPPORT
  // terminate xercesc
  if(xmlInputFN.size() != 0){
    xercesc::XMLPlatformUtils::Terminate();
  }
#endif //XML_SUPPORT
  // delete temporary file if reading form stdin
  if(readStdIn){
    remove(xmlInputFN.c_str());
  }
  if(VERB > 2){
    std::cerr << "FeatureNames::getNumFeatures(): "<< FeatureNames::getNumFeatures() << endl;
  }
  vector<vector<double> > w(xval_fold,vector<double> (FeatureNames::getNumFeatures()+ 1)), ww;
  int firstNumberOfPositives = preIterationSetup(w);
  if (VERB > 0) {
    cerr << "Estimating " << firstNumberOfPositives << " over q="
        << test_fdr << " in initial direction" << endl;
  }
  // Set up a first guess of w
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (VERB > 1) cerr << "Reading in data and feature calculation took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  
  // Do the SVM training
  if (VERB > 0) {
    cerr << "---Training with Cpos";
    if (selectedCpos > 0) {
      cerr << "=" << selectedCpos;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", Cneg";
    if (selectedCneg > 0) {
      cerr << "=" << selectedCneg;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", fdr=" << selectionfdr << endl;
  }
  train(w);
  if (!pCheck->validateDirection(w)) {
    fullset.calcScores(w[0]);
  }
  if (VERB > 0) {
    cerr << "Merging results from " << xv_test.size() << " datasets"
        << endl;
  }

  // calculate psms level probabilities
  
  //PSM probabilities TDA or TDC
  calculatePSMProb(false, &fullset, procStart, procStartClock, w, diff, target_decoy_competition);
  if (xmlOutputFN.size() > 0){
    writeXML_PSMs();
  }
  
  // calculate unique peptides level probabilities WOTE
  if(reportUniquePeptides){
    calculatePSMProb(true, &fullset, procStart, procStartClock, w, diff, target_decoy_competition);
    if (xmlOutputFN.size() > 0){
      writeXML_Peptides();
    }
  }
  // calculate protein level probabilities with FIDO
  if(calculateProteinLevelProb){
    calculateProteinProbabilitiesFido();
  }
  // write output to file
  writeXML();  
  return 0;
}
