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

#ifndef CONFIG_H_
#define CONFIG_H_


#ifndef MZIDENTML_VERSION_MAJOR
  #define MZIDENTML_VERSION_MAJOR "@MZIDENTML_VERSION_MAJOR@"
#endif
#ifndef MZIDENTML_VERSION_MINOR
  #define MZIDENTML_VERSION_MINOR "@MZIDENTML_VERSION_MINOR@"
#endif
#ifndef MZIDENTML_NAMESPACE
  #define MZIDENTML_NAMESPACE "@mzIdentML-namespace@"
#endif
#ifndef MZIDENTML_SCHEMA_LOCATION
  #define MZIDENTML_SCHEMA_LOCATION "@MZIDENTML_SCHEMA_LOCATION@"
#endif

#ifndef GAML_TANDEM_VERSION_MAJOR
  #define GAML_TANDEM_VERSION_MAJOR "@GAML_TANDEM_VERSION_MAJOR@"
#endif
#ifndef GAML_TANDEM_VERSION_MINOR
  #define GAML_TANDEM_VERSION_MINOR "@GAML_TANDEM_VERSION_MINOR@"
#endif
#ifndef GAML_TANDEM_NAMESPACE
  #define GAML_TANDEM_NAMESPACE "@gaml_tandem-namespace@"
#endif
#ifndef GAML_TANDEM_SCHEMA_LOCATION
  #define GAML_TANDEM_SCHEMA_LOCATION "@GAML_TANDEM_SCHEMA_LOCATION@"
#endif

#ifndef TANDEM_VERSION
  #define TANDEM_VERSION "@TANDEM_VERSION@"
#endif
#ifndef TANDEM_NAMESPACE
  #define TANDEM_NAMESPACE "@tandem-namespace@"
#endif
#ifndef TANDEM_SCHEMA_LOCATION
  #define TANDEM_SCHEMA_LOCATION "@TANDEM_SCHEMA_LOCATION@"
#endif

#endif /*CONFIG_H_*/
