# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import os
import subprocess
from zipfile import ZipFile
from pathlib import Path


configfile: "config/gb-model/config.common.yaml"


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


rule extract_naturalearth_data:
    message:
        "Extract natural earth shape data to get lat, lon information of EU countries"
    input:
        zip="data/gb-model/downloaded/world-shapes.zip",
    output:
        shape_files=expand(resources("naturalearth/ne_110m_admin_0_countries.{ext}"),ext=['shp','shx','dbf','prj']),
        country_coordinates=resources("country_coordinates.csv")
    log:
        logs("extract_naturalearth_data.log"),
    script:
        "../scripts/gb-model/extract_naturalearth_data.py"

rule create_powerplants_table:
    message:
        "Tabulate powerplant data GSP-wise from FES workbook sheet BB1 and EU supply data",
    params:
        scenario=config["fes"]["scenario"],
        year=config["fes"]["year"],
    input:
        bb1_sheet=resources("fes/BB1.csv"),
        bb2_sheet=resources("fes/BB2.csv"),
        gsp_coordinates="data/gb-model/downloaded/gsp-coordinates.csv",
        eu_supply="data/gb-model/downloaded/eu-supply-table.csv",
    output:
        csv=resources("fes_powerplants.csv")
    log:
        logs("create_powerplants_table.log")
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/create_powerplants_table.py"