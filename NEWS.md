# capnet 0.3.3
* Update function documentations
* Fix parallel processing bug when running on machines with more than 128 cores

# capnet 0.3.2

* Added `cv_capnet()` support for parallel execution via PSOCK cluster.
* Added `family` argument to `capnet()`, `cv_capnet()`, and `walk_capnet()` supporting Gaussian, binomial, Poisson, and Gamma families.
* Added `plot.cv_capnet()` heatmap and slice visualizations for cross-validation error surfaces.
* Added `plot.capnet_path()` for coefficient path plots along a hyperparameter dimension.
* Fixed contribution cap gradient scale-invariance under standardization.

# capnet 0.3.0

* Added `walk_capnet()` for walk-forward evaluation with per-step contribution cap enforcement.
* Added `coef_path()` for fitting models along a hyperparameter grid.
* Added `capnet_violations()` to diagnose which (row, feature) pairs exceed their cap.
* Added `plot_redistribution()` to visualize coefficient and contribution shifts across models.

# capnet 0.2.0

* Added `cv_capnet()` for K-fold cross-validation over `alpha` and `lambda`.
* Added `simulate_capnet_data()` for generating synthetic data with known contribution structure.

# capnet 0.1.0

* Initial release. Core `capnet()` estimator with OWL-QN optimization for elastic net plus contribution cap penalty.
