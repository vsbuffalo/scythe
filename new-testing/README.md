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
