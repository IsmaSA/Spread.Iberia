**Spread dynamics of invasive species in the Iberian Peninsula**


This repository accompanies the manuscript:

> Soto, I., Novoa, A., Hermoso, V., & Oficialdegui, F. J. *Spread dynamics of
> invasive species in the Iberian Peninsula.* (in review, *Ecological Informatics*).

---

## Overview

For each species, the pipeline aggregates cleaned occurrence records into a
10 × 10 km grid and five-year time windows (1950–2024) and derives a set of
complementary range-expansion metrics:

- **Occupied area** (km²) from per-cluster concave hulls (DBSCAN clustering with a species-specific, adaptive ε);
- **Number of occupied cells** (10 × 10 km);
- **Maximum extent** — greatest pairwise distance between occurrences (invasion "leading edge");
- **Spread rate** (km · year⁻¹) between consecutive windows;
- **Centroid shift** and **dominant direction** of spread;
- A **target-group background correction** to down-/up-weight cells by sampling effort;
- **Sigmoidal (logistic) growth** fits to cumulative occupied cells, yielding lag phase, inflection and saturation points.


### Sources
- **GBIF** occurrence download — DOI/key `0008621-260519110011954`
  (https://doi.org/10.15468/dl.&lt;insert DOI&gt;), retrieved 07 October 2025.
- ~28,000 **expert-validated unpublished records** for the Iberian Peninsula
  (Katsanevakis et al., under review).
- EU **Union List** of invasive alien species (Commission Implementing
  Regulation 2016/1141 and updates to 2025/1422).


## License

Code is released under the MIT License (see `LICENSE`). Occurrence data are
subject to the terms of their original providers (GBIF / data publishers).

## Contact

Ismael Soto — isma-sa@hotmail.com · ismael.soto@ebd.csic.es
