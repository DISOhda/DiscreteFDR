## Test environments

* local Manjaro Linux 26.0.2 install, R 4.5.2
* win-builder (release, oldrelease, devel)
* R-hub (configurations: linux, macos-arm64, windows, lto, valgrind), see
  https://github.com/DISOhda/DiscreteFDR/actions/runs/21981039987


## R CMD check results

### local

0 errors | 0 warnings | 0 notes


### win-builder

0 errors | 0 warnings | 1 note

* checking DESCRIPTION meta-information ... NOTE
  Author field differs from that derived from Authors@R
  
  - must be false alarm, there is no additional Authors field in DESCRIPTION
    file
  - happens only for oldrelease checks, not for release or devel


### R-hub

0 errors | 0 warnings | 1 note

* checking sizes of PDF files under ‘inst/doc’ ... NOTE
Unable to find GhostScript executable to run checks on size reduction

  - only happens on some configurations, not on CRAN, so it must be related
    to the affected containers configurations