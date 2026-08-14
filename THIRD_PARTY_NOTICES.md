# Third-Party Notices and Research Attribution

This repository contains independently written MATLAB implementations of
selected published equations, constraints, numerical conditions, and methods.
The BSD-3-Clause license in `LICENSE` applies to the original source code and
documentation in this repository unless a file explicitly states otherwise.
It does not grant any rights to the cited journal articles, publisher-formatted
PDFs, article figures, substantial article text, third-party software, or
unpublished/proprietary data.

The Git repository does not redistribute the source-paper PDFs. Consult the
original articles through their DOI links and comply with the terms shown by
the relevant publisher.

## Zhang et al. blackout-zone communication study

- Rouhe Zhang, Kang Wang, Xinran Duan, and Zheng Chen, "Entry trajectory
  optimization considering blackout zone communication constraint,"
  *Advances in Space Research*, vol. 77, pp. 11407-11417, 2026.
- DOI: https://doi.org/10.1016/j.asr.2025.11.062
- Article notice: Copyright 2025, published by Elsevier B.V. on behalf of
  COSPAR. No Creative Commons reuse license is stated in the article copy used
  for this implementation.
- Repository scope: independent implementations and audits of published
  antenna/RAAP geometry, blackout-zone constraints, angle-of-attack and
  aerodynamic expressions, published entry conditions, and selected
  dimensionless equations. Relevant code is principally under
  `+paperstudies/+zhang`, with reduced mission-level adaptations elsewhere.
- Excluded material: the publisher PDF, article figures, substantial text,
  GPOPS-II results, and unpublished vehicle or optimization data.

## Saito et al. controlled small-spacecraft re-entry study

- Takumi Saito, Toshinori Kuwahara, Yuji Saito, and Yuji Sato, "Guidance
  strategies for controlled Earth reentry of small spacecraft in low Earth
  orbit," *Acta Astronautica*, vol. 229, pp. 684-697, 2025.
- DOI: https://doi.org/10.1016/j.actaastro.2024.12.054
- Article notice: Copyright 2024, The Authors; published by Elsevier Ltd on
  behalf of IAA under CC BY-NC-ND 4.0:
  https://creativecommons.org/licenses/by-nc-nd/4.0/
- Repository scope: independent implementations and audits of published
  entry-state tables, de-orbit case grids, spherical entry equations,
  aeroassist guidance expressions, estimator algebra, and uncertainty grids.
  Relevant code is principally under `+paperstudies/+saito`, with reduced
  capsule adaptations elsewhere.
- Excluded material: the article PDF and figures, the unavailable proprietary
  HSRC aerodynamic database, unpublished controller parameters, and raw
  reference trajectories.

## Interpretation and non-endorsement

Source-derived values are classified in the implementation and documentation
as published-exact, convention-required, inferred/project-assumed, surrogate,
or unavailable. These labels are part of the scientific provenance and must
not be removed when results are reported.

This project is an independent research implementation. It is not affiliated
with, sponsored by, endorsed by, or validated by the cited authors, their
institutions, COSPAR, IAA, or Elsevier. Reproducing a deterministic check from
this repository is not, by itself, evidence that every result in either paper
has been numerically reproduced.
