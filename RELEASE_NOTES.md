### v3.0.0 -- January 22, 2020
* Add support for netcore SDK 3.1
* BREAKING: Requires .NET Framework 4.7 or greater to run.

### v2.1.1 -- July 18, 2019
* Add a file copy with checksums utility function.

### v2.1.0 -- July 18, 2018
* Add support for generating oligos with ambiguous bases.
* Update FSharp.Core version and use paket management.

### v2.0.22 -- August 3, 2017
* Fix index error in SW alignment.

### v2.0.21 -- February 15, 2017
* Fix error type bug in overlap stitching function.
* Simplify stitching and loopout result types to options.

### v2.0.20 -- February 15, 2017
* Elevate some exception printing methods to ErrorHandling.

### v2.0.19 -- February 15, 2017
* Add OptionExtensions namespace for providing useful extensions to Option.

### v2.0.18 -- February 9, 2017
* Fix Dna hash to use structural hashing of underlying array.

### v2.0.17 -- February 8, 2017
* primerCore gets oligoDesignWithComprise b/c sometimes you just have to break the rules

### v2.0.16 -- January 24, 2016
* Performance improvements to protein translation algorithm,  added test

### v2.0.15 -- December 5, 2016
* Improvements to Dna type.
* Add Length member to SGD feature type.

### v2.0.14 -- December 2, 2016
* Various improvements and extensions to the Dna type.

### v2.0.13 -- December 2, 2016
* Add hashing and append to Dna domain type.

### v2.0.12 -- November 16, 2016
* ENH: SuffixTreeDisk allows parallel access under windows

### v2.0.11 -- November 10, 2016
* ENH: promote a bunch of error handling helpers to Amyris.ErrorHandling.

### v2.0.10 -- October 18, 2016
* BUG: revert behavior of format60 with new performance boost

### v2.0.9 -- October 14, 2016
* BUG: MerHash not wrapping correctly for end of array searches

### v2.0.8 -- October 13, 2016
* utils.format60 reverting return type to string

### v2.0.7 -- October 12, 2016
* merhash functions added, pileup format parser added

### v2.0.6 -- September 27, 2016
* Small bug fix in ePCR.

### v2.0.5 -- September 20, 2016
* Small improvements to stitching module.

### v2.0.4 -- September 12, 2016
* Added stitching module.

### v2.0.3 -- July 17, 2016
* Migrated to Project Scaffold for final release

### v2.0.3-beta -- July 14, 2015
* Fixed warnings, see revision 38 in Mercurial.

### v2.0.2-beta -- July 14, 2015
* Fixed Tm calculation bug in primercore2, see revision 36 in Mercurial.

### v2.0.1-beta -- June 1, 2015
* Fixed Tm calculation bug in primercore2, see revision 33 in Mercurial.

### v2.0.0-beta -- May 13, 2015
* Improved primercore module
* avoids any mono- and dinucleotide repeat
* transferred all hard-coded parameters to the PrimerParam structure
* added ability to specify nucleotide specific penalties.

### v1.0.0 -- April 15, 2015
* Initial NuGet release.
