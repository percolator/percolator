/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Globals.cpp,v 1.5 2007/07/16 22:57:21 lukall Exp $
 *******************************************************************************/
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

void Globals::clean()
{
    if (glob)
      delete glob;
    glob = 0;
}
