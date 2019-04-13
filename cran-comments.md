This release removes the use of the unexported S3 method `ggplot2:::predictdf`.
This fixes the bad practice of relying on finding S3 methods on the search path.

## Test environments

* local OS X 10.14.3, R 3.5.3
* win-builder (devel and release)
* Rhub
  * Debian Linux, R-release, GCC
  * Fedora Linux, R-devel, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs.
