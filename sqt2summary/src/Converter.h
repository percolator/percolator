#ifndef CONVERTER_H_
#define CONVERTER_H_

#define COMETDB_CGI   "comet-fastadb.cgi"

#ifdef  __CYGWIN__
#define CGI_DIR   "isb-bin"
#define PLOT_CGI  "sequest-plot.cgi"
#define OUT_CGI   "sequest-out.cgi"
#else
#define CGI_DIR   "cgi-bin"
#define PLOT_CGI  "sequest-tgz-plot.cgi"
#define OUT_CGI   "sequest-tgz-out.cgi"
#endif

#include <string>
#include <vector>

#include "SequestOut.h"

using namespace std;

class Converter
{
    public:
        Converter();
        virtual ~Converter();
        void printSummary(Header& hdr, string& szCWD);
        void read_sqt(const string fname);
        void readFeatures(const string &in, SequestOut& feat,int match);
    private:
        string szPeptideLink;
        const static int hitsPerSpectrum = 1;
        vector<SequestOut> data;


};
#endif                                            /*CONVERTER_H_*/
