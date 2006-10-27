#include "Globals.h"

Globals * Globals::glob = 0;

Globals::~Globals()
{
}

Globals::Globals()
{
//    assert(!glob);
    glob=this;
    verbose =2;
}

Globals * Globals::getInstance()
{
    if (!glob) new Globals();
    return glob;
}
