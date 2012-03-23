neopz_old: this is an old version of neopz. It has been ported to work with cmake. 
This pz lib achieved good scalability results with substruct. It will be 
hosted here while we investigate the performance differences between the old 
and the new version of the neopz.
  NOTES: make sure you use metis-4.0 when building the old version.

When compiling both neopz libs for performance measurements, make sure you:

- use g++ on new MacOS systems
- use the Release option
- use the appropriate Metis lib.
