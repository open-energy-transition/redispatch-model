
..
  SPDX-FileCopyrightText: Open Energy Transition gGmbH and contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
  SPDX-FileCopyrightText: Contributors to gb-open-market-model <https://github.com/open-energy-transition/gb-open-market-model>

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Methodology
##########################################

Here, we detail the existing Network Options Assessment (NOA) Refresh 2021-22 methodology and our plan to reproduce it using entirely open data and tools.

.. note::
    The *NOA Refresh* differs from the *NOA* as it is published 6 months later and includes updates responding to the Holistic Network Design (HND) report.
    The HND maps out the future offshore wind infrastructure requirements to meet defined targets and therefore impacts the NOA recommendations.
    Unlike the initial NOA report, The NOA Refresh is *not* concerned with recommending the future capacity of interconnectors between Great Britain and its neighbours.

.. _existing-method:

--------------------
Existing methodology
--------------------

.. raw:: html

   <iframe src="../_static/existing_methodology_high_level.drawio.html" style="width: 100%; height: 550px; border: none;"></iframe>

   Graphical representation of our understanding of the NOA Refresh 2021-22 methodology.
   Clickable links to methodology / results reports are given for each data source box.

====
Data
====

The existing methodology relies on the Future Energy Scenarios (FES) *Leading the Way* scenario for the majority of its input data.
This sets the basis for the 2030 Great Britain (and neighbours) power system that is then built upon with other datasets and an energy system model that brings it all together.
The FES is enhanced with more detailed baseline transmission network data from the Electricity Ten Year Statement (ETYS) and the Network Options Assessment Interconnector analysis.

Where data is not published by the FES, such as sub-annual supply/demand profiles, a reference database from the modelling team is used.
In the existing methodology, this reference database is proprietary.
The NOA Refresh 2021-22 utilised AFRY databases.
More recent modelling efforts have utilised Energy Exemplar databases.

=====
Model
=====

The model has two distinct phases: dispatch and re-dispatch.
In the unconstrained dispatch phase, a rolling horizon dispatch model is solved while ignoring GB inter-nodal transmission capacities.
In the re-dispatch phase, inter-nodal transmission capacities (a.k.a. "boundary capabilities") are introduced and the rolling horizon optimisation is re-run with the opportunity for generators to increase (offer) or decrease (bid) generation at a pre-defined cost in response to network congestion.
In 2017, the optimisation was specified as running at a 4hr resolution (8 periods per day) for 20 years.
The purpose of the re-dispatch phase is to calculate the system cost incurred by congestion (a.k.a. "constraint costs").
These capacities are varied across scenarios to represent different transmission reinforcement proposals (a.k.a. "network options"), and the resulting constraint costs compared to find the most cost effective options to proceed with.

This section is a summary of the 2017 `Long Term Market and Network Constraint Modelling report <https://www.nationalgrid.com/sites/default/files/documents/Long-term%20Market%20and%20Network%20Constraint%20Modelling.pdf>`_ and the `NOA 2021-22 Methodology <https://www.neso.energy/document/204196/download>`_.


.. _method-existing-boundary-capabilities:

Boundary capabilities
=====================

There are three Transmission Operators (TOs) in Great Britain.
Each of these manages the network infrastructure in their jurisdiction and responds to the ETYS, which highlights the need for future network reinforcements, with proposals for reinforcements across their networks.
These proposals can come in various forms, including as upgrades to existing infrastructure (lines, transformers, etc.) or as new lines on the network.
They can be both onshore and offshore and AC/DC.

The impact of these proposals on boundary capabilities, with boundaries defined by ETYS, is simulated by the TOs and submitted to NESO following a template (the "Systems Requirements Form" - SRF).
These submissions are considered commercially sensitive and are not publicly available.
The SRF includes the expected change in boundary capabilities, the seasonal variation of that change, and the expected financial investment.
Where seasonal variations are not provided, default adjustment factors will be used:

- Spring and autumn thermal:  85%
- Summer thermal 80%
- Summer outage thermal 70%
- Summer outage voltage 90%

The NOA's purpose is to assess the cost-benefit ratio of these proposals.
There are 136 proposals considered in the NOA Refresh 2021-22.
It is evident that all possible combinations of these proposals cannot be simulated and compared (approx. 18,500 scenarios).
Instead, regions of the country are considered in isolation, with only the options in those regions being iterated over.
Boundary capabilities of other regions are held constant, usually using the decisions from the previous NOA.

Bids/Offers
===========

Generator bids and offers are based on their short-run marginal costs (SRMC) scaled by multipliers.

Thermal plants
--------------

SRMC is :math:`\frac{\text{fuel cost} + \text{CO}_{2}\text{ cost}}{\text{plant efficiency}}`.
Multipliers are calculated from historic data (the 5 years preceding the modelling year) as :math:`\frac{\text{annual average bid/offer}}{\text{annual average SRMC}}`.
Values are given in a 2017 report for each type of thermal generator, which *may* be the values used in all subsequent studies:

=====  ====  =====
Asset  Bid   Offer
=====  ====  =====
Gas    0.80  1.63
Coal   0.70  1.56
Oil    0.38  4.16
OCGT   1.37  1.31
Ave.   0.75  1.59
=====  ====  =====

Renewables & new nuclear
------------------------

Since SRMC is low/negligible for renewables and nuclear, a different approach is taken.
Bids/offers are based on the subsidy available for those generators.
The methodology only refers to bids, implying that only *reducing* the output of these generators is possible.
This is a reasonable assumption as the initial dispatch run will prioritise these technologies due to their low/zero SRMC and therefore they are likely to be operating at full capacity when entering the re-dispatch phase.

Storage
-------

No mention is made of storage bid/offer prices in the methodology.

Interconnectors
---------------

There are three possible pricing calculations for bids and offers of interconnectors, based on whether they are importing/exporting/floating (unutilised) in a given time period.
To calculate bid/offer quantities for each of these states, the following is undertaken:

1. Calculate wholesale electricity price in GB and Europe in an unconstrained run.
2. Calculate the interconnector fee (:math:`Fee`) as the market spread between GB and Europe (difference in wholesale electricity price)
3. In each hour of the unconstrained run, calculate the the marginal plant in Europe to assign a bid/offer multiplier to that hour.
   The marginal price combined with the multiplier gives the interconnector :math:`Bid` / :math:`Offer`.
   It is assumed that the import/export quantity in re-dispatch is not sufficient to change the marginal plant in Europe.
4. Apply a cost to different changes in state:
  1. Interconnector importing & reduce imports / float / start exporting: :math:`P_{GB} - P_{foreign} \times (1 + Loss) - Bid \times (1 + Loss)`
  2. Interconnector importing & increase imports: :math:`Fee + Offer \times (1 + Loss)`
  3. Interconnector at float & start exporting: :math:`Fee - Bid \times (1 - Loss)`
  4. Interconnector at float & start importing: :math:`Fee + Offer \times (1 + Loss)`
  5. Interconnector exporting & increase exports: :math:`Fee - Bid \times (1 - Loss)`
  6. Interconnector exporting & reduce exports / float / start importing: :math:`P_{foreign} \times (1 - Loss) - P_{GB} + Offer \times (1 - Loss)`

.. note::

  It is unclear from the documentation what :math:`P` refers to in the above math.


The case of importing/exporting in the unconstrained run and then exporting/importing in the re-dispatch run is poorly represented by these costs (as it switches between two cost formulations).

Plant availability
==================

Data is usually available per month.
Some is supposedly available from the FES, although this is augmented with proprietary datasets.
Availability reflects planned and unplanned outages and is applied to generators as a de-rating of their peak capacity.
This is applied per month, possibly even at a higher resolution if appropriate data is available.
It is unclear how granular plant availability has been applied, beyond monthly de-rating.

Interconnectors have a constant 95% availability applied to them.

Plant representation
====================

Thermal plants
--------------

Thermal plants are grouped if they have similar variable costs.
Characteristics that are represent include "technology type, fuel type, efficiency, start-up cost,
part load efficiency, operating cost, and availability" as well as must-run restrictions (if there are contractual obligations).
Fuel costs do not fluctuate within the year but do vary by European country.
The FES provides fuel and carbon price forecasts.

Combined Heat and Power
-----------------------

The model can distinguish between extraction CHP and backpressure CHP.
For each CHP plant, heat load profiles are specified.

Hydropower
----------

Hydro is split in the model into Reservoir and Run-Of-River.
Inflows are modelled as a cascade with multiple levels, with inflow expectation, the ability of generators to forecast inflows ahead of time, and actual inflow levels.
Mismatches between forecast and actual inflows can then affect the model and are reflected in the constraints and dispatch approach.

Inflows are defined as "unregulated" and "regulated".
Unregulated flows must pass through the plant immediately while regulated ones can be stored in a reservoir.

Hydro reservoirs are grouped into single units per region with water volume represented as stored energy.
Energy spillage is possible (if inflow exceeds capacity) and a maximum release level can be set in given periods.

Since modelled operators account for the uncertainty of future inflows, a medium-term dispatch is likely required.

More description on hydropower modelling can be found `here <https://www.svk.se/4a87f2/siteassets/5.jobba-har/dokument-exjobb/implementation-of-hydropower-modeling-for-electricity-market-simulations-using-bid3.pdf>`_.

Wind and Solar
--------------

Wind and solar generation is based on historical profiles.
The technologies are classed as "must-run" which probably doesn't need to be explicitly forced in the optimisation problem since their short-run marginal cost is negligible so they will be given priority in the optimisation by default.
Curtailment factors and maximum levels of non-synchronous generation *can* be applied but there is no indication that they *are* in this model.

Wind generation profiles are split between existing and future sites.
This is to ensure that capacity factors are uprated from the regional average and instead account for actual turbine siting (which will likely gravitate towards higher capacity factor sites).
Solar generation profiles are not split as they assume the regional average is approximately the same as at individual sites.
Capacity factors are calculated for regions that are not the same as the transmission network regions and are not the same between wind and solar.
Offshore wind regions are not only based on their nearest onshore region.
Some offshore wind regions exist that have no onshore boundaries.

Other considerations
--------------------

A fixed fee `balancing services use of system (BSUoS) charge <https://www.neso.energy/industry-information/charging/balancing-services-use-system-bsuos-charges>`_ is used.

Scarcity rent, Start-up and No load costs, Ramp rates, Temperature dependent start cost, and Minimum on- and off-times are all available but do not look to be used in the model.

European neighbours
===================

Europe is modelled according to a set of assumptions on the state of its energy infrastructure in each modelled year.
The original modelling effort (2016/17) used a proprietary model with in-house expertise to define future infrastructure, considering national policies, emissions targets, GDP growth, stated plans, etc..
The 2016/17 report states that the plan from 2017 was to use the system operator's European model.
This implies that the NOA 2021-22 uses the FES European assumptions.

The European model is run together with the GB model in the unconstrained run, to calculate wholesale electricity prices and marginal plants in neighbouring countries.

Clean Energy Package Constraints
================================

EU/2019/943 Article 13 paragraph 5 of the Clean Energy Package requires that more than 50% of total energy volumes must be renewable (including high-efficiency cogeneration) or that less than 5% of renewables energy volumes are re-dispatched (excluding high-efficiency cogeneration).
That is, so long as the 50% threshold is reached, the 5% re-dispatch threshold does not need to be adhered to.
If the 50% threshold is not reached, the 5% threshold is checked and reinforcement options are iterated upon until it is reached.
It is possible to still be compliant whilst not reaching either threshold (there are other ways to achieve compliance).

See paragraph 2.117 of the NOA 2021-22 methodology for more detail.

================
Decision process
================

A reinforcement option is deemed necessary if it generates a positive Net Present Value (NPV) over its lifetime when introduced into the model.

This section is a summary combining the 2017 `Long Term Market and Network Constraint Modelling report <https://www.nationalgrid.com/sites/default/files/documents/Long-term%20Market%20and%20Network%20Constraint%20Modelling.pdf>`_, the `NOA 2021-22 Methodology <https://www.neso.energy/document/204196/download>`_ and the more recent `tCSNP Refresh Methodology <https://www.neso.energy/document/357916/download>`_ (the successor to the NOA).

Reinforcement option costs
==========================

Net Present Value (NPV) of reinforcement options are calculated for each re-dispatch run by combining constraint costs, investment costs, and societal cost of carbon.

Constraint costs
----------------

The "constraint cost" is calculated as the sum of bids/offers following re-dispatch.
This is calculated for 20 years of re-dispatch modelling (2023 - 2043) plus the final year of modelling repeated 20 times, to account for a 40-year project lifetime.
It is then distributed to individual boundaries in the network by:

1. calculating a "congestion charge" using the shadow price associated with each boundary constraint, i.e. the marginal cost change associated with increasing that boundary's capacity.
2. calculating a "congestion rent delta" as :math:`\text{congestion charge} \times (\text{unconstrained flow} - \text{constrained flow})`.
3. distributing the total constraint cost using the relative congestion rent delta of each boundary.

This distribution can be applied hourly or weekly.
The methodology prefers weekly as hourly includes fluctuations that are within the (re-)dispatch rolling horizon window, in which generators can bid/offer with perfect foresight on congestion for the whole day.

Investment cost
---------------

To calculate the NPV of a system with a reinforcement option, the constraint costs is combined with the option's capital cost, assuming a 40-year project lifetime.
To compare it with an annual constraint cost, it is amortised over its life by applying the `Spackman approach <https://www.ofgem.gov.uk/sites/default/files/docs/2011/10/discounting-for-cost-benefit-analysis-involving-private-investment-but-public-benefit.pdf>`_, using the TO's WACC to annualise and HM Treasury's Social Time Preferential Rate to discount all costs and benefits.
It is an approach that is stated as suitable for projects with private financing but public benefits.

Delay costs
-----------

It is possible for options to be delayed from their Earliest In Service Date (EISD).
This usually comes at a cost, defined by the TO, which is added into the option total cost.

Societal Cost of Carbon
-----------------------

The societal cost of carbon comes into play if operating carbon emissions reduce as a result of a reinforcement option.
In such cases, an economic benefit equal to the quantity of reduced emissions is added to the NPV calculation (assumption: compared to the base case run).
`Traded carbon values <https://www.gov.uk/government/collections/carbon-valuation--2>`_ are likely used for this calculation.

.. note::

  This does not appear in the NOA 2021-22 methodology, only the tCSNP, so we probably do not need to account for it to reproduce the NOA 2021-22.

Modelling steps
===============

.. raw:: html

   <iframe src="../_static/existing_methodology_modelling_steps.drawio.html" style="width: 100%; height: 300px; border: none;"></iframe>

   Graphical representation of our understanding of the NOA Refresh 2021-22 modelling steps.


The modelling process involves adding candidates one a time in a pre-defined priority order for each specific GB region.
As :ref:`aforementioned <method-existing-boundary-capabilities>`, reinforcements in other GB regions are fixed based on the most recent NOA report (here: 2020-2021).
Regions are arbitrarily defined.
The only criterion is that no option may appear in more than one region.

The method for a region goes as follows:

1. Re-dispatch run is undertaken with all proposed options disabled, to generate a *base case* which defines the reinforcement requirements on each boundary in the region for the first year (2023).
2. All reinforcement options are ordered by priority to satisfy the requirements.
3. In order of priority, the options are added one-at-a-time and the 20-year re-dispatch model is re-run to get a new set of constraint costs.
4. If multiple options are equal in priority, they are each modelled independently.
5. For every option modelled at the same priority level, the NPV (constraint cost relative to base case combined with investment and delay costs) is calculated for introducing it at its EISD and for increasing delay years.
   The lowest NPV for each option is chosen as its "optimum" year.
6. If multiple options have been considered, the one with the earliest optimum year is chosen.
   If it is unclear which should be chosen, the modelling is branched at this point.
7. The chosen option is added to the base case and aforementioned steps are repeated with a *new* base case and the next option in the priority order.
8. When the lowest NPV option gives a negative NPV, no more options are considered.

The final set of considered options are referred to as the "reinforcement profile", each having an optimum year that it is in service.
If this in-service year is equal to its EISD, the option is considered "critical".

See paragraphs 2.78 - 2.92 of the NOA 2021-22 methodology report for more detail.