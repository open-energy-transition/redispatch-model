# SPDX-FileCopyrightText:  gb-open-market-model contributors
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
    resources:
        mem_mb=1000,
    conda:
        "../envs/shell.yaml"  # This is required to install `curl` into a conda env on Windows 
    shell:
        "curl -sSLvo {output} {params.url}"


# Rule to create region shapes using create_region_shapes.py
rule create_region_shapes:
    input:
        country_shapes=resources("country_shapes.geojson"),
        etys_boundary_lines="data/gb-model/etys-boundary-gis-data.zip",
    output:
        raw_region_shapes=resources("raw_region_shapes.geojson"),
    log:
        logs("raw_region_shapes.log"),
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
        logs("manual_region_merger.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/manual_region_merger.py"


# Rule to retrieve generation unit unavailability data from ENTSO-E
rule retrieve_unavailability_data:
    output:
        planned_unavailability=resources("planned_unavailability.csv"),
        forced_unavailability=resources("forced_unavailability.csv"),
    log:
        logs("retrieve_unavailability_data.log"),
    resources:
        mem_mb=2000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/gb-model/retrieve_unavailability_data.py"
