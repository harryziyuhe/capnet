## CRAN readiness tracker
1. Package structure and metadata
- [ ] Valid DESCRIPTION, NAMESPACE, license, and maintainer information
- [ ] Unique package name, correct dependencies, and clean source tarball
- [ ] No development-only dependencies or unnecessary files

2. Documentation and examples
- [ ] Every exported function documented with working examples
- [ ] README showing installation and a complete modeling workflow
- [ ] At least one useful vignette explaining contribution caps
- [ ] All examples run quickly and without external data

3. Testing and statistical correctness
- [ ] Comprehensive automated unit tests for fitting, prediction, and optimization
- [ ] Numeric results independently checked against trusted baseline
- [ ] Edge cases, invalid inputs, and numerical stability tested
- [ ] Cross-validation and walk-forward procedures checked for leakage

4. CRAN checks and portability
- [ ] R CMD check --as-cran passes with 0 ERRORs and 0 WARNINGs
- [ ] Every NOTE is resolved or justified
- [ ] Checks pass on Windows, macOS, and Linux
- [ ] Checks pass with current R and R-devel
- [ ] No file-system, parallelization, or runtime policy violations

5. Submission readiness
- [ ] All hyperlinks, copyright notices, and citations verified
- [ ] Final source tarball built and checked
- [ ] cran-comments.md explains the initial release and check results
