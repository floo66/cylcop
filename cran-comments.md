## Test environments
- local R installation, Windows 10 x64, R 4.2.1
- win-builder (devel and release),
- R-hub windows-x86_64-devel (r-devel)
- R-hub macos-highsierra-release-cran (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results

There were no ERRORs or WARNINGs. There are 2 NOTEs:

### On windows-x86_64-devel (r-devel)

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

As noted in R-hub issue #503, this could be due to a bug/crash in MiKTeX and can likely be ignored.

### On fedora-clang-devel (r-devel)

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable
  
I cannot change that Tidy is not on the path, or update Tidy on the external 
Fedora Linux server.
Same for pacakge V8.
  

## Downstream dependencies
There are no downstream dependencies.
