# mnlfa

**Moderated Nonlinear Factor Analysis (MNLFA) in R**

`mnlfa` implements the stepwise MNLFA procedure for assessing measurement invariance across observed moderator variables (e.g. age, sex, group membership). A key advantage over existing implementations is full support for **multidimensional models** — multiple latent variables (facets) per domain — as well as an optional acquiescence factor for questionnaire data with reversed items.

## Installation

```r
# install.packages("devtools")
devtools::install_github("zieglema/mnlfa_R")
```

## Dependencies

The package uses [OpenMx](https://openmx.ssri.psu.edu/) as its fitting backend via the [mxsem](https://github.com/jhorzek/mxsem) interface. Both are installed automatically.

```r
install.packages("OpenMx", repos = "https://openmx.ssri.psu.edu/OpenMx3/")
install.packages("mxsem")
```

## Pipeline overview

| Step | Function | Purpose |
|------|----------|---------|
| 1 | — | Define per-facet baseline model strings |
| 2 | `run_mxsem()` | Fit baseline CFAs |
| 3 | `add_latent_mean_regression()` + `significant_mean_moderators()` | Screen latent mean moderation |
| 4 | `add_latent_variance_moderation()` + `significant_var_moderators()` | Screen latent variance moderation |
| 5 | `moderate_loadings_and_intercepts()` + `significant_item_moderators()` | Screen item-level DIF |
| 5b | `summarise_anchor_candidates()` + `plot_anchor_candidates()` | Review IC table & plot for anchor selection |
| 6 | — | Select anchor items |
| 7 | `combine_sig_moderators()` | Combine all significant moderations |
| 8 | `create_acquiescence_lines()` | Build acquiescence factor (optional) |
| 9 | `build_full_domain_models()` | Assemble full domain model string |
| 10 | `run_mxsem()` | Fit domain model |
| 11 | `add_covariance_moderation()` | Screen covariance moderation (optional) |
| 12 | `build_all_individual_parameters()` | Compute person-specific parameters |
| 13 | `extract_moderation_table()` | Summary table of all moderation effects |
| 14 | `compute_conditional_parameters()` + `plot_conditional_parameters()` | Interaction plots |

## Usage

See `mnlfa_sim_3factor.R` for a fully annotated simulation script (3 facets × 10 items, 2 moderators), and `mnlfa_worked_example.Rmd` for a narrative worked example with explanations.

## References

- Bauer, D. J. (2017). A more general model for testing measurement invariance and differential item functioning. *Psychological Methods*, 22(3), 507–526. <https://doi.org/10.1080/00273171.2016.1256260>
- Gottfredson, N. C., et al. (2019). Simplifying the implementation of modern scale scoring methods with an application to measure functional impairment. *Psychological Methods*, 24(4), 400–416. <https://doi.org/10.1037/met0000193>
- Kolbe, L., et al. (2022). Investigating measurement invariance with many groups: Representativeness versus power. *Structural Equation Modeling*, 29(3), 380–395. <https://doi.org/10.1080/10705511.2021.1984772>

## Author

Matthias Ziegler — Humboldt-Universität zu Berlin
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--9373--0365-green)](https://orcid.org/0000-0002-9373-0365)
