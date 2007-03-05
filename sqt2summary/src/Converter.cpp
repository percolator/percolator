#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;
#include "SequestOut.h"
#include "Converter.h"

Converter::Converter() :
szPeptideLink("PEPTIDELINK=http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=")
{

}

Converter::~Converter()
{
}


void Converter::readFeatures(const string &in, SequestOut& feat,int match) {
  istringstream instr(in),linestr;
  string line,tmp,pep;
  int charge;
  double mass,deltCn,otherXcorr=0,xcorr=0;
  bool gotL=true,gotDeltCn=(match==0);
  int ms=0;
  
  while (getline(instr,line)) {
    if (line[0]=='S') {
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> tmp >> charge >> tmp >> tmp >> mass;
    }
    if (line[0]=='M' && !gotDeltCn) {
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> tmp >> tmp >> deltCn >> otherXcorr;
      feat.dDeltCn=deltCn;
      gotDeltCn = true;
    }
    if ((line[0]=='M') && (match==ms++)) {
      int rSp,matched,expected;
      double cMass,sp;
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> rSp >> cMass >> tmp >> xcorr >> sp >> matched >> expected >> pep;
      
      feat.iRankSp=rSp;      // rank by Sp
      feat.dMass=mass;
      feat.dAMass=cMass; // calc mass
      feat.dXC=xcorr;      // Xcorr
      feat.dSp=sp;      // Sp
      feat.iIon=matched;
      feat.iTot=expected;
      feat.cAA1=pep.at(0);
      feat.cAA2=pep.at(pep.length()-1);
      feat.szSubPep=pep.substr(2,pep.length()-4);
      feat.szPlainPep=feat.szSubPep;
      gotDeltCn = (match!=0);
      gotL = false;
    }
    if (line[0]=='L' && !gotL) {
      if (instr.peek() != 'L') gotL=true;
      string p;
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> p;
      feat.szProt=p;      
    }
  }
}

void Converter::read_sqt(const string fname) {
  int n = 0,charge=0,ms=0;
  string line,tmp,prot;
  istringstream lineParse;  
  ifstream sqtIn;
  sqtIn.open(fname.c_str(),ios::in);
  if (!sqtIn) {
  	cerr << "Could not open file " << fname << endl;
  	exit(-1);
  }
  bool look = false;
  while (getline(sqtIn,line)) {
    if (line[0]=='S' && sqtIn.peek() != 'S') {
         lineParse.str(line);  
         lineParse >> tmp >> tmp >> tmp >> charge;       
         if (charge <= 3) look=true;
         ms=0;
    }
    if (look & line[0]=='L' && ms < hitsPerSpectrum) {
         lineParse.str(line);  
         lineParse >> tmp >> prot;
         ++ms;
         ++n;
    }
  }
  cerr << n << " records in file " << fname << endl;
  sqtIn.clear();
  sqtIn.seekg(0,ios::beg);
  
  int ix=data.size(),lines=0;
  data.resize(ix+n);
  string seq;
  
  string fileId(fname);
  unsigned int spos = fileId.rfind('/');
  if (spos!=string::npos)
    fileId.erase(0,spos+1);
  spos = fileId.find('.');
  if (spos!=string::npos)
    fileId.erase(spos);
  
  ostringstream buff,id;
  
  string scan;
  set<int> theMs;
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
      if(lines>1 && charge<=3) {
        string record=buff.str();
        string idstr = id.str();
        set<int>::const_iterator it;
        for(it=theMs.begin();it!=theMs.end();it++) {
          readFeatures(record,data[ix],*it);
          data[ix].szBaseFileName=idstr;
          data[ix].szFileName=idstr+".dta";
          ix++;
        }
      }
      buff.str("");
      buff.clear();
      id.str("");
      lines=1;
      buff << line << endl;
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scan >> charge;
      id << fileId << '.' << scan << '.'<< scan << '.' << charge ;
      ms=0;
      theMs.clear();
    }
    if (line[0]=='M') {
      ++ms;
      ++lines;
      buff << line << endl;    
    }
    if (line[0]=='L') {
      ++lines;
      buff << line << endl;
      if((int)theMs.size()<hitsPerSpectrum) {
        theMs.insert(ms-1);
      }
    }
  }
  if(lines>1 && charge<=3) {
    string record=buff.str();
    string idstr = id.str();
    set<int>::const_iterator it;
    for(it=theMs.begin();it!=theMs.end();it++) {
      readFeatures(record,data[ix],*it);
      data[ix].szBaseFileName=idstr;
      data[ix].szFileName=idstr+".dta";
      ix++;
    }
  }
  sqtIn.close();
//  cout << "Read File" << endl;
}



 /** printSummary writes to stdout the contents
   *  of the sequestDataStruct and headerStruct in
   *  INTERACT format.
   */
void Converter::printSummary(Header& hdr,string& szCWD)
{
   unsigned int i = 0,
       j = 0,
       iMaxFileNameWidth = 0,
       iMaxRefWidth = 0,
       iMaxSeqWidth = 0,
       lines=data.size();

   // Get max column widths
   for (i = 0; i < lines; i++)
   {
      // FileName
      if (data[i].szBaseFileName.length() > iMaxFileNameWidth)
         iMaxFileNameWidth = data[i].szBaseFileName.length();

      // Ref
      if (data[i].szProt.length()+data[i].szDup.length() > iMaxRefWidth)
         iMaxRefWidth = data[i].szProt.length()+data[i].szDup.length();

      // Sequence
      if (data[i].szSubPep.length() > iMaxSeqWidth)
         iMaxSeqWidth = data[i].szSubPep.length();
   }
   iMaxRefWidth += 1;

   // Print out easy part of the header
   printf("<HTML>\n");
   printf("<HEAD><TITLE>HTML-SUMMARY</TITLE></HEAD>\n");
   printf("<BODY BGCOLOR=\"#FFFFFF\">\n");
   printf("<PRE><FONT COLOR=\"green\">");
   printf("HTML-SUMMARY</FONT>\n");
   printf("Institute for Systems Biology\n");
   printf("Seattle, WA\n");
   printf("%s %s %s %s, %s\n", hdr.szDate.c_str(), hdr.szTime.c_str(), hdr.szTimeSuffix.c_str(), data[0].szDatabase.c_str(), hdr.szMassType.c_str());
   printf("\n");
   printf("<FONT COLOR=\"green\">");

   // Print complex portion of the header
   printf("   #    File");
   for (i = 0; i < iMaxFileNameWidth - 4; i++)
      printf(" ");

   printf("    MH+          XCorr ");
   printf("   dCn      Sp    RSp    Ions");

   printf("   Ref");
   for (i = 0; i < iMaxRefWidth; i++)
      printf(" ");

   printf("  Sequence</FONT>\n");
   printf("<FONT COLOR=\"green\">");
   printf("  ----  ");
   for (i = 0; i < iMaxFileNameWidth; i++)
      printf("-");
   printf("  ------         ------");
   printf("  -----   ------  ---   ------  ");
   for (i = 0; i < iMaxRefWidth; i++)
      printf("-");
   printf("  ");
   for (i = 0; i < iMaxSeqWidth; i++)
      printf("-");
   printf("</FONT>\n");

   for (i = 0; i < lines; i++)
   {

      // Index
      printf(" %5d", i + 1);
      // File
/*
      if (data[i].szFileName[0] != '/')
         printf ("  <A TARGET=\"Win1\" HREF=\"/%s/%s?OutFile=%s/%s\">%s</A>",
             CGI_DIR, OUT_CGI, szCWD, data[i].szFileName, data[i].szBaseFileName);
*/
      if (data[i].szFileName[0] != '/')     /* add extra /./ in path before filename for Pep3D compatibility */
         printf ("  <A TARGET=\"Win1\" HREF=\"/%s/%s?OutFile=%s/./%s\">./%s</A>",    
             CGI_DIR, OUT_CGI, szCWD.c_str(), data[i].szFileName.c_str(), data[i].szBaseFileName.c_str());
      else
         printf ("  <A TARGET=\"Win1\" HREF=\"/%s/%s?OutFile=%s\">%s</A>",
             CGI_DIR, OUT_CGI, data[i].szFileName.c_str(), data[i].szBaseFileName.c_str());
      for (j = 0; j < iMaxFileNameWidth - data[i].szBaseFileName.length(); j++)
         printf(" ");

      if (data[i].szSubPep[0]!='\0')
      {
         // MH+
         printf("  %6.1lf", data[i].dMass);
         // MHError
         printf(" (%+1.1lf)", data[i].dAMass - data[i].dMass);
         // Xcorr
         printf("  %0.4lf", data[i].dXC);
         // dCn 
         if (data[i].dDeltCn >= 0.2)
         {
            if (data[i].bSpecialDeltCn)
               printf("  <FONT COLOR=\"#DD00DD\">%0.3lf*</FONT>", data[i].dDeltCn);
            else
               printf("  <FONT COLOR=\"#DD00DD\">%0.3lf </FONT>", data[i].dDeltCn);
         }
         else
         {
            if (data[i].bSpecialDeltCn)
               printf("  %0.3lf*", data[i].dDeltCn);
            else
               printf("  %0.3lf ", data[i].dDeltCn);
         }
         // Sp
         printf("  %6.1lf", data[i].dSp);
         // Rsp
         if (data[i].iRankSp == 1)
         {
            printf("  <FONT COLOR=\"#DD00DD\">%3d</FONT>", data[i].iRankSp);
         }
         else
         {
            printf("  %3d", data[i].iRankSp);
         }
         // Ions TODO: This needs special parameters based on static/variable 
         // mods
         if (data[i].szBaseFileName[0] != '/')
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s/%s.dta&amp;MassType=%d&amp;NumAxis=1",
                CGI_DIR, PLOT_CGI, szCWD.c_str(), data[i].szBaseFileName.c_str(), data[i].iMassType);
         else
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s.dta&amp;MassType=%d&amp;NumAxis=1",
                CGI_DIR, PLOT_CGI, data[i].szBaseFileName.c_str(), data[i].iMassType);
         if (data[i].dMass1 != 0.0)
            printf("&amp;DMass1=%f", data[i].dMass1);
         if (data[i].dMass2 != 0.0)
            printf("&amp;DMass2=%f", data[i].dMass2);
         if (data[i].dMass3 != 0.0)
            printf("&amp;DMass3=%f", data[i].dMass3);
         if (data[i].dMass4 != 0.0)
            printf("&amp;DMass4=%f", data[i].dMass4);
         if (data[i].dMass5 != 0.0)
            printf("&amp;DMass5=%f", data[i].dMass5);
         if (data[i].dMass6 != 0.0)
            printf("&amp;DMass6=%f", data[i].dMass6);
         if (data[i].szSubPep.length() != data[i].szPlainPep.length())
            printf("&amp;DSite=%s", data[i].szDSite.c_str());
         printf("%s&amp;Pep=%s\">%3d/%3d</A>", data[i].szMod.c_str(), data[i].szPlainPep.c_str(),
                data[i].iIon, data[i].iTot);
   
         // Ref
         printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Ref=%s&amp;Db=%s&amp;NucDb=%d&amp;Pep=%s&amp;MassType=%d\">%s</A>",
             CGI_DIR, COMETDB_CGI, data[i].szProt.c_str(), data[i].szDatabase.c_str(), data[i].bNucDb,
             data[i].szPlainPep.c_str(), data[i].iMassType, data[i].szProt.c_str());

         for (j = 0; j < iMaxRefWidth - data[i].szProt.length() - data[i].szDup.length(); j++)
            printf(" ");
   
         // optional +field
         if (data[i].szDup.length() > 1)
         {
            printf("<A TARGET=\"Win1\" HREF=\"/%s/%s?Db=%s&amp;NucDb=%d&amp;Pep=%s&amp;MassType=%d\">%s</A>",
                CGI_DIR, COMETDB_CGI, data[i].szDatabase.c_str(), data[i].bNucDb, data[i].szPlainPep.c_str(),
                data[i].iMassType, data[i].szDup.c_str());
         }
         printf("  ");
   
         // Seq TODO: Breakup the peptide
         printf("%c.<A TARGET=\"Win1\" HREF=\"%s%s\">%s</A>.%c",
                data[i].cAA1,
                szPeptideLink.c_str(),
                data[i].szPlainPep.c_str(),
                data[i].szSubPep.c_str(),
                data[i].cAA2);

         for (j = 0; j < iMaxSeqWidth - data[i].szSubPep.length() + 1; j++)
            printf(" ");
         printf("\n");

      }
      else
      {
         if (data[i].szBaseFileName[0] != '/')
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s/%s.dta&amp;NumAxis=1",
               CGI_DIR, PLOT_CGI, szCWD.c_str(), data[i].szBaseFileName.c_str());
         else
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s.dta&amp;NumAxis=1",
               CGI_DIR, PLOT_CGI, data[i].szBaseFileName.c_str());

         printf("\">spectrum</A>\n");
      }
   }

   printf("\n</BODY></HTML>\n");
   return;
}

int main(int argc, char *argv[]) {
  int retVal=0;
  Converter *pCaller=new Converter();
  string fName,cwd("");
  for (int i=1;i<argc;i++) {
    fName=argv[i];
    pCaller->read_sqt(fName);
  }
  Header hdr;
  pCaller->printSummary(hdr,cwd);
  delete pCaller;
  return retVal;
}	

