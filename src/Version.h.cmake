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

#ifndef VERSION_H_
#define VERSION_H_

#ifndef VERSION_MAJOR
  #define VERSION_MAJOR "@CPACK_PACKAGE_VERSION_MAJOR@"
#endif

#ifndef VERSION_MINOR
   #define VERSION_MINOR "@CPACK_PACKAGE_VERSION_MINOR@"
#endif

#ifndef VERSION_PATCH
   #define VERSION_PATCH "@CPACK_PACKAGE_VERSION_PATCH@"
#endif

#ifndef VERSION
  #define VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@.@CPACK_PACKAGE_VERSION_PATCH@"
#endif

#ifndef VERSION_NAME
  #define VERSION_NAME "v@CPACK_PACKAGE_VERSION_MAJOR@-@CPACK_PACKAGE_VERSION_MINOR@"
#endif

#ifndef LVERSION_NAME
  #define LVERSION_NAME L"v@CPACK_PACKAGE_VERSION_MAJOR@-@CPACK_PACKAGE_VERSION_MINOR@"
#endif

#endif /*VERSION_H_*/

