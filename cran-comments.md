## Test environments
* local Windows 10 21H1 install, R 4.1.0
* win-builder (release, oldrelease, devel)
* R-hub (configurations: "with sanitizers", "with valgrind" and "check for CRAN")


## R CMD check results

### local
0 errors | 0 warnings | 1 note

* checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Note: information on .o files for x64 is not available
  File 'D:/Nextcloud/R/PoissonBinomial.Rcheck/PoissonBinomial/libs/i386/PoissonBinomial.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'D:/Nextcloud/R/PoissonBinomial.Rcheck/PoissonBinomial/libs/x64/PoissonBinomial.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)

The cause is unknown; none of the mentioned functions is used.

### R-hub
0 errors | 0 warnings | 3 notes

* NOTE: Package unavailable to check Rd xrefs: 'discreteMTP'

This note seems to be related to the check configuration on R-hub.

### win-builder

0 errors | 0 warnings | 0 notes