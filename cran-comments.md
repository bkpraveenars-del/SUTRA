## Test environments

* local Ubuntu 24.04, R 4.3.3 — 0 errors | 0 warnings | 0 notes
* win-builder (devel and release)      — TO BE RUN BEFORE SUBMISSION
* macOS builder (r-release)            — TO BE RUN BEFORE SUBMISSION
* R-hub: linux, macos, windows          — TO BE RUN BEFORE SUBMISSION

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* This is a new submission.
* The DOI cited in the Description field
  (<doi:10.1002/bimj.201100219>) is the methodological reference the package
  implements; it resolves.
* Examples use the package's own simulation engine and run in under 1.5 seconds
  in total, so no example is wrapped in \donttest{} or \dontrun{}.
* The package writes no files and sets no options outside the R session.
* SpATS is used through its exported interface only.
