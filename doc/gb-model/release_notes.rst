
..
  SPDX-FileCopyrightText: Open Energy Transition gGmbH and contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
  SPDX-FileCopyrightText: Contributors to gb-open-market-model <https://github.com/open-energy-transition/gb-open-market-model>

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################

Upcoming Release
================

* Tabulated monthly GB powerplant fractional availability profiles (#71).
* Remove unnecessary output in `compose_networks` rule that causes error (#2)
* Tabulated GSP wise powerplant capacities for GB (#4).
* Tabulated EU country level aggregated powerplant capacities (#33)
* Add rule 'retrieve_unavailability_data' to Snakemake workflow for fetching unavailability data from ENTSO-E. (#43)
* Increase number of HTTP download retries to mitigate against Zenodo file retrieval timeouts.
* Keep all retrieved data locally by default to reduce time spent re-downloading data on every run.
* Add FES workbook data download and sheet extraction rule (#50).
* Restructured documentation (#27).
* Added modelling methodology documentation (#20).
* Added GB custom geographic boundary rule and script (#13).
