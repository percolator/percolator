#ifndef GLOBALS_H_
#define GLOBALS_H_

#define VERB (Globals::getInstance()->getVerbose())

class Globals
{
public:
	virtual ~Globals();
    static Globals * getInstance();
    int getVerbose() {return verbose;}
    void setVerbose(int verb) {verbose=verb;}
    void decVerbose() {verbose--;}
    void incVerbose() {verbose++;}
private:
	Globals();
    int verbose;
    static Globals * glob; 
};

#endif /*GLOBALS_H_*/
