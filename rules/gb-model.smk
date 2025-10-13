# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import os
import subprocess
from zipfile import ZipFile
from pathlib import Path


# Rule to download and extract ETYS boundary data
rule download_data:
    message:
        "Download {wildcards.gb_data} GB model data."
    output:
        downloaded="data/gb-model/downloaded/{gb_data}",
    params:
        url=lambda wildcards: config["urls"][Path(wildcards.gb_data).stem],
    log:
        logs("download_{gb_data}.log"),
    localrule: True
    conda:
        "../envs/gb-model/workflow.yaml"
    shell:
        "curl -sSLvo {output} {params.url}"


# Rule to create region shapes using create_region_shapes.py
rule create_region_shapes:
    input:
        country_shapes=resources("country_shapes.geojson"),
        etys_boundary_lines="data/gb-model/downloaded/gb-etys-boundaries.zip",
    output:
        raw_region_shapes=resources("raw_region_shapes.geojson"),
    log:
        logs("raw_region_shapes.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/create_region_shapes.py"


# Rule to manually merge raw_region_shapes
rule manual_region_merger:
    input:
        raw_region_shapes=rules.create_region_shapes.output.raw_region_shapes,
        country_shapes=resources("country_shapes.geojson"),
    output:
        merged_shapes=resources("merged_shapes.geojson"),
    log:
        logs("manual_region_merger.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/manual_region_merger.py"


# Rule to retrieve generation unit unavailability data from ENTSO-E
rule retrieve_entsoe_unavailability_data:
    message: "Retrieve data from ENTSOE API for generator {wildcards.business_type} unavailability in {wildcards.zone} bidding zone"
    output:
        xml_base_dir = directory("data/gb-model/entsoe_api/{zone}/{business_type}")
    params:
        unavailability = config["entsoe_unavailability"]
    log:
        logs("retrieve_entsoe_unavailability_data_{zone}_{business_type}.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/retrieve_entsoe_unavailability_data.py"

rule process_entsoe_unavailability_data:
    input:
        xml_base_dir = "data/gb-model/entsoe_api/{zone}/{business_type}"
    output:
        unavailability=resources("{zone}_{business_type}_generator_unavailability.csv"),
    log:
        logs("process_entsoe_unavailability_data_{zone}_{business_type}.log"),
    params:
        business_type_codes = config["entsoe_unavailability"]["api_params"]["business_types"]
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/process_entsoe_unavailability_data.py"

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
        logs("compose_network_{clusters}.log"),
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/compose_network.py"


rule compose_networks:
    input:
        expand(
            resources("networks/composed_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        intermediate_data=[
            resources("transmission_availability.csv"),
            expand(resources("fes/{fes_sheet}.csv"), fes_sheet=config["fes-sheet-config"].keys()),
            expand(
                resources("{zone}_{business_type}_generator_unavailability.csv"),
                zone=config["entsoe_unavailability"]["bidding_zones"],
                business_type=config["entsoe_unavailability"]["business_types"]
            ),
            resources("merged_shapes.geojson")
        ]


rule extract_transmission_availability:
    input:
        pdf_report="data/gb-model/downloaded/transmission-availability.pdf",
    output:
        csv=resources("transmission_availability.csv"),
    log:
        logs("extract_transmission_availability.log"),
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/extract_transmission_availability.py"


rule extract_fes_workbook_sheet:
    message:
        "Extract FES workbook sheet {wildcards.fes_sheet} and process into machine-readable, 'tidy' dataframe format according to defined configuration."
    input:
        workbook="data/gb-model/downloaded/fes-workbook.xlsx",
    output:
        csv=resources("fes/{fes_sheet}.csv"),
    params:
        sheet_extract_config=lambda wildcards: config["fes-sheet-config"][
            wildcards.fes_sheet
        ],
    log:
        logs("extract_fes_{fes_sheet}.log"),
    script:
        "../scripts/gb-model/extract_fes_sheet.py"
