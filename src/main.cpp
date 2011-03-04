/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"
#include "Caller.h"
using namespace std;

#include <vector>
#include "Scores.h"
#include "Scores.cpp"

int main(int argc, char** argv) {
/*
  vector<ScoreHolder> scores = vector<ScoreHolder>();
  ScoreHolder s1 = ScoreHolder();
  s1.pPSM = new PSMDescription("X.bbb.X", 1);
  s1.pPSM->id = "psm1";
  s1.pPSM->proteinIds.insert("pid1");
  s1.label=1;

  ScoreHolder s2 = ScoreHolder();
  s2.pPSM = new PSMDescription("X.aaa.X", 2);
  s2.pPSM->id= "psm2";
  s2.pPSM->proteinIds.insert("pid2");
  s2.label=1;

  ScoreHolder s3 = ScoreHolder();
  s3.pPSM = new PSMDescription("X.aaa.X", 3);
  s3.pPSM->id ="psm2";
  s3.pPSM->proteinIds.insert("pid3");
  s3.label=-1;
  scores.push_back(s1);
  scores.push_back(s2);
  scores.push_back(s3);
  vector<ScoreHolder>::iterator it = scores.begin();
  for(;it!=scores.end();it++){
    cout << *it;
  }
  cout << endl;

  sort(scores.begin(), scores.end(), lexicOrder());
   // lexicographically order the scores (based on peptides names)
   // new list of scores
   vector<ScoreHolder> uniquePeptideScores = vector<ScoreHolder>();
   string previousPeptide;
   int previousLabel;
   // run a pointer down the scores list
   vector<ScoreHolder>::iterator current = scores.begin();
   for(;current!=scores.end(); current++){
     // compare pointer's peptide with previousPeptide
     string currentPeptide = current->pPSM->getPeptideNoResidues();
     if(previousPeptide.compare(currentPeptide)==0) {
       if((previousLabel == current->label)){
         // if the peptide is a duplicate, append to previously inserted
         vector<ScoreHolder>::iterator last = --uniquePeptideScores.end();
         // the duplicate psm_id...
         last->psms_list.push_back(current->pPSM->id);
         // ... and all its proteins
         BOOST_FOREACH(string pId, current->pPSM->proteinIds){
           last->pPSM->proteinIds.insert(pId);
         }
       }
     } else {
       // otherwise insert as a new score
       current->psms_list.push_back(current->pPSM->id);
       uniquePeptideScores.push_back(*current);

       // update previousPeptide
       previousPeptide = currentPeptide;
       previousLabel = current->label;
     }
   }
   scores = uniquePeptideScores;
   sort(scores.begin(), scores.end(), greater<ScoreHolder> ());


   it = scores.begin();
  for(;it!=scores.end();it++){
    cout << (ScoreHolderPeptide)*it;
  }

  exit(0);
*/

  Caller* pCaller = new Caller();
  int retVal = -1;
  if (pCaller->parseOptions(argc, argv)) {
    retVal = pCaller->run();
  }
  delete pCaller;
  Globals::clean();
  return retVal;
}
