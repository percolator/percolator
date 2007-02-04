/*******************************************************************************
 * Percolator v 1.00
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas KŠll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Globals.cpp,v 1.2 2007/02/04 04:33:21 lukall Exp $
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
