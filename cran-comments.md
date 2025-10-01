## Test environments
* Local R installation, windows 11 (R-4.4.2)
* R-hub:
  - [VM] linux (ubuntu-latest, R-devel)
  - [VM] macos (macOS 13, R-release)
  - [VM] windows (Windows-latest, R-release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Notes for CRAN
* This submission fixes issues with Rd cross-references reported in the
  previous CRAN checks (now updated to use `\link[PKG]{FOO}` style
  anchors).
  
* It fixes the issues "Lost braces in \itemize; \value handles \item{}{} directly"

* It fixes the issued from wrongly using the `size` aesthetic for lines in ggplot

* It fixes the error arising from using the function `expression()` as a text label in 
  some ggplot objects
