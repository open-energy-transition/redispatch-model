# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import os
import subprocess
from zipfile import ZipFile
from pathlib import Path


# Rule to download and extract ETYS boundary data
rule retrieve_etys_boundary_data:
    output:
        boundary_shp="data/gb-model/etys-boundary-gis-data.zip",
    params:
        url=config["urls"]["gb-etys-boundaries"]
    log:
        logs("retrieve_etys_boundary_data.log")
    resources:
        mem_mb=1000,
    conda: "../envs/shell.yaml"  # This is required to install `curl` into a conda env on Windows
    shell: "curl -sSLvo {output} {params.url}"


# Rule to create region shapes using create_region_shapes.py
rule create_region_shapes:
    input:
        country_shapes=resources("country_shapes.geojson"),
        etys_boundary_lines="data/gb-model/etys-boundary-gis-data.zip",
    output:
        raw_region_shapes=resources("raw_region_shapes.geojson")
    log:
        logs("raw_region_shapes.log")
    resources:
        mem_mb=1000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/create_region_shapes.py"


# Rule to manually merge raw_region_shapes
rule manual_region_merger:
    input:
        raw_region_shapes=resources("raw_region_shapes.geojson"),
    output:
        merged_shapes=resources("merged_shapes.geojson"),
    log:
        logs("manual_region_merger.log")
    resources:
        mem_mb=1000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/manual_region_merger.py"


rule compose_networks:
    input:
        expand(
            resources("networks/composed_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"]
        )


rule compose_network:
    input:
        unpack(input_profile_tech),
        network=resources("networks/base_s_{clusters}.nc"),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
    output:
        network=resources("networks/composed_{clusters}.nc"),
    params:
        countries=config["countries"],
        costs_config=config["costs"],
        electricity=config["electricity"],
        clustering=config["clustering"],
        renewable=config["renewable"],
        lines=config["lines"],
    log:
        logs("compose_network_{clusters}.log")
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/compose_network.py"
