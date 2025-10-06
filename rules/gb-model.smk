# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import os
import subprocess
from zipfile import ZipFile
from pathlib import Path


configfile: "config/gb-model/config.common.yaml"


# Rule to download and extract ETYS boundary data
rule retrieve_etys_boundary_data:
    output:
        boundary_shp="data/gb-model/etys-boundary-gis-data.zip",
    params:
        url=config["urls"]["gb-etys-boundaries"],
    log:
        logs("retrieve_etys_boundary_data.log"),
    localrule: True
    conda:
        "../envs/gb-model/workflow.yaml"
    shell:
        "curl -sSLvo {output} {params.url}"


# Rule to create region shapes using create_region_shapes.py
rule create_region_shapes:
    input:
        country_shapes=resources("country_shapes.geojson"),
        etys_boundary_lines=rules.retrieve_etys_boundary_data.output.boundary_shp,
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


rule download_transmission_availability_pdf:
    output:
        pdf_report="data/gb-model/transmission-availability.pdf",
    params:
        url=config["urls"]["transmission-availability"],
    log:
        logs("transmission_availability.log"),
    localrule: True
    conda:
        "../envs/gb-model/workflow.yaml"
    shell:
        "curl -sSLvo {output} {params.url}"


rule extract_transmission_availability:
    input:
        pdf_report=rules.download_transmission_availability_pdf.output.pdf_report,
    output:
        csv=resources("transmission_availability.csv"),
    log:
        logs("extract_transmission_availability.log"),
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb-model/extract_transmission_availability.py"
