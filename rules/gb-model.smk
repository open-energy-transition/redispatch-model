# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


configfile: "config/config.common.yaml"


# Rule to create region shapes using create_region_shapes.py
rule create_region_shapes:
    input:
        country_shapes=resources("country_shapes.geojson"),
        etys_boundary_lines="data/gb-model/etys-boundary-gis-data-mar25/ETYS boundary GIS data Mar25.shp"
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
