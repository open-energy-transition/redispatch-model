..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
  SPDX-FileCopyrightText: Contributors to gb-open-market-model <https://github.com/open-energy-transition/gb-open-market-model>

  SPDX-License-Identifier: CC-BY-4.0

#############
Data Sources
#############

gb-open-market-model is compiled from a variety of data sources.
The following table provides an overview of the data sources used exclusively in gb-open-market-model.
For data sources used in PyPSA-Eur, see `this page <../data_sources.html>`_.
Different licenses apply to the data sources.

---------------------------------
The Future Energy Scenarios (FES)
---------------------------------

`The FES <https://www.neso.energy/publications/future-energy-scenarios-fes>`_ is the primary data source for defining the model, both for GB and other European countries.
Here, we use the 2021 FES data workbook.
Tables from the workbook we use are:

- BB1: Building Block Data
- BB2: Building Block Metadata
- SV.34: Installed BECCS generation capacity (GW)


-----------------
GSP coordinates
-----------------
The GSP coordinates <https://api.neso.energy/dataset/963525d6-5d83-4448-a99c-663f1c76330a/resource/41fb4ca1-7b59-4fce-b480-b46682f346c9/download/fes2021_regional_breakdown_gsp_info.csv> is obtained from the NESO website. This is used to assign lat, lon to powerplants extracted from the FES workbook

---------------
EU Supply data
---------------
The EU supply data <https://api.neso.energy/dataset/bd83ce0b-7b1e-4ff2-89e8-12d524c34d99/resource/6563801b-6da4-46e7-b147-3d81c0237779/download/fes2023_es2_v001.csv> is used to retrieve powerplant data of neighbouring countries to GB

---------------
Interconnectors
---------------
Electricity transmission interconnectors between GB regions and neighbouring countries are based on distinct projects considered in the FES (table 9, `FES modelling methods <https://www.neso.energy/document/199916/download>`_).
We combine those projects manually to create a total GB interconnector capacity curve from 2021-2041 that matches the curves given in the FES workbook, sheet SV.37.
The GB region to which those projects connect is based on geolocating the connecting transformer as defined in the NESO `interconnector register <https://www.neso.energy/data-portal/interconnector-register>`_.
For projects not in the register (since some outdated projects are no longer considered), we have used the respective `TYNDP <https://tyndp.entsoe.eu/>`_ project data sheet to estimate their GB onshoring coordinates.
Project definitions and our manually defined start dates for them are user-configurable.

.. note::
  No reasonable combination of projects perfectly matches the FES results.
  However, the combination culminating in the FES results is not publicly available, so the projects we choose is an opinionated assumption.