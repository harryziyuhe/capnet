## CRAN readiness tracker
1. Package structure and metadata
- [x] Valid DESCRIPTION, NAMESPACE, license, and maintainer information (R CMD
      check validated metadata, namespace, and dependencies)
- [x] Unique package name; the maintainer confirmed `capnet` is available on
      CRAN.
- [x] Correct dependencies; `R CMD check` validated package dependencies
- [x] Clean source tarball; inspected its contents after `R CMD build`
- [x] No development-only dependencies or unnecessary files in the source
      tarball; inspected its contents after `R CMD build`

2. Documentation and examples
- [x] Every exported function documented with working examples; `R CMD check`
      reported no missing documentation and all examples passed
- [x] README shows installation and a complete modeling workflow
- [x] The included guide vignette explains contribution caps and passes the
      build and vignette checks
- [x] Examples and vignette run without network data; the guide guards its
      optional `glmnet` and `xtable` dependencies. The guide build takes about
      82 seconds locally.
- [x] Updated the `capnet()` help page to describe projected proximal gradient
      for finite coefficient bounds, L1 soft-thresholding, and `tol` behavior.
      The local roxygen2 version is 8.0.0 versus the recorded 8.1.0, so it
      reported a version mismatch while regenerating documentation.

3. Testing and statistical correctness
- [x] Automated tests cover fitting, prediction, gradients/optimization,
      cross-validation, walk-forward behavior, coefficient bounds, and the
      bounded L1 proximal step; the latest `devtools::test()` passed 164 tests
- [x] Numeric fitted results are compared against trusted baselines: the
      unpenalized Gaussian fit is checked against `stats::lm`, and the uncapped
      elastic-net fit is checked against `glmnet` with explicit tolerances
- [x] Tests cover invalid inputs, convergence retries, standardization, and
      numerical gradient behavior at contribution-cap boundaries
- [x] Reviewed fold construction and walk-forward fitting: CV excludes each
      held-out fold from that fold's training data, and walk-forward fits use
      training data before the evaluation block

4. CRAN checks and portability
- [x] The revised 0.4.0 tarball passed `R CMD check --as-cran` with 0 errors
      and 0 warnings. This local run reported 3 environment notes: CRAN/GitHub
      URLs could not be reached from this machine, the system time could not be
      verified, and R's HTML manual validator flags generated `<main>` tags.
- [x] GitHub Actions checks passed on Windows, macOS, and Linux using release R
- [x] Checks passed with current R and R-devel; the workflow also checks
      oldrel-1 on Linux
- [x] Added a GitHub Actions matrix for release R on Windows/macOS/Linux,
      R-devel on Linux, and oldrel-1 on Linux; all configured jobs passed.
- [x] Cross-platform checks passed in GitHub Actions on Windows, macOS, and
      Linux, including R-devel on Linux. The optional WinBuilder email result
      was not received; the passing workflow matrix provides the R-devel and
      cross-platform checks needed here.
- [x] The updated R-hub Apple-silicon macOS R-devel check passed:
      https://github.com/r-hub2/gentle-falcon-capnet/actions/runs/38800675092
- [x] Source review found no file-writing or network behavior in package
      examples/tests/guide; parallel examples/tests start no workers, and the
      cluster helper cleans up with `on.exit()`

5. Submission readiness
- [x] Verify package hyperlinks and guide citations. The GitHub repository and
      issue URLs and README badges were reviewed; the guide's cited bibliography
      entries were checked, including the OWL-QN/L-BFGS reference.
- [x] Copyright and attribution are consistent: `DESCRIPTION` identifies
      Harry He as the sole author and maintainer and now carries `Copyright:
      2026 Harry He`; `LICENSE.md` and the README identify the same holder.
      No other copyright notices occur in the included package files.
- [x] Built and checked the 0.4.0 source tarball after the baseline,
      input-validation, and bounded-fit tests, the bounded optimizer, new logo
      assets, updated help page, and reviewer-requested DESCRIPTION references
      and vignette option cleanup. The latest manual check completed with
      0 errors and 0 warnings; examples, tests, and vignette checks passed.
- [x] Set the `Date` field in `DESCRIPTION` to the planned submission date
      (`2026-10-28`); update it again if submission happens on a later date.
- [ ] Rebuild and check the 0.4.1 source tarball after the version hike.
- [ ] Submit the final 0.4.1 source tarball through CRAN's submission form and
      respond to the confirmation email.
