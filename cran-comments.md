## Resubmission
This is a resubmission. In this version I have:
* Changed one url. https://www.jstor.org/stable/2335637/ --> https://doi.org/10.2307/2335637. So sorry for the issues with this! With 
jstor-urls, I always get the note below which masks any other problems with the
url.

## Test environments
* local R installation, R 4.0.5
* ubuntu 16.04 (on travis-ci), R 4.0.5
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* Found the following (possibly) invalid URLs:
     URL: https://doi.org/10.2307/2335637
       From: man/cor_cyl.Rd
       Status: 403
       Message: Forbidden
