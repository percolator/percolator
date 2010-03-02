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

#include "Globals.h"

Globals * Globals::glob = 0;

Globals::~Globals() {
}

Globals::Globals() {
  //    assert(!glob);
  glob = this;
  verbose = 2;
}

Globals * Globals::getInstance() {
  if (!glob) new Globals();
  return glob;
}

void Globals::clean() {
  if (glob) delete glob;
  glob = 0;
}
