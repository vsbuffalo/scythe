# Scythe Testing

## File and Operation Overview

 - `Makefile` dispatches cleaning and `run.sh`
 - `run.sh` is what runs other trimmers under a variety of parameters,
   for each of the simulated read files (with all replicates).
 - `analysis.R` contains R methods for matching simulated reads and
   their trimmed corresponding reads (for all trimmers tested) and
   comparing their length.
  
## Trimmers

`run.sh` is designed to use the most recent Scythe version; it
references it in the directory above the testing directory. Other
trimmers must be in the `$PATH` environmental variable.

### Trimmer Versions

 - btrim, not versioned but the 2011-09-16 version of `btrim_mac` was
   downloaded from http://graphics.med.yale.edu/trim/.  The Linux
   version `btrim32` is also not versioned, but the 2011-02-31 version
   was used.
 - cutadapt, version 1.0 was tested. This version is dated 2011-11-04
   and downloaded from
   http://code.google.com/p/cutadapt/downloads/list.
 - scythe version 0.981 was used.


