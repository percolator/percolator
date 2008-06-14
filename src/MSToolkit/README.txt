Source Code: MSToolkit
Purpose: The basic code and object structures to facilitate loading, manipulating, and saving MS1 and MS2 files. Also provides
conversion from text based .ms1 files to binary based .bms1 files. Includes ability to read peaks from (not write to) mzXML format.

Version History
------------------------------------------
1.0  - Apr 20 2007 - Stable. Works with Hardklor 1.0 (public). Jesse has older version
1.1  - May 17 2007 - Several bugfixes and new MSReader function (appendFile with MSObject).
1.2  - Sep 18 2007 - Minor bugfixes/improvements.

2.0  - Sep 20 2007 - Added mzXML support - now requires zlib, base64, and ramp libraries and toolkits.
2.1  - Sep 29 2007 - Added file compression. Recommended extensions are .cms1 .cms2 etc.
2.2  - Jan 16 2008 - Fixed bugs in 64-bit file support. Now runs on 64-bit Linux.
2.3  - Apr 03 2008 - Added latest RAMP (2950).
                     Fixed bug in compressed file reader.
                     Added file version numbers to binary and compressed formats (helps if formats change in the future).
                     Changed spectrum header info for binary and compressed formats.
                     * WARNING * New binary and compressed formats are not backwards compatible. All old files must be 
                     re-formatted.
2.31 - Apr 09 2008 - Fixed missing precursor m/z when reading MS/MS from mzXML files.
2.4  - Apr 16 2008 - Added additional changes and documentation as per lab requests.
2.41 - Apr 18 2008 - Bugfixes for mzXML.
2.42 - May 01 2008 - More bugfixes.