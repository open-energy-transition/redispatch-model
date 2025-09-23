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
        boundary_shp="data/gb-model/etys-boundary-gis-data-mar25/ETYS boundary GIS data Mar25.shp",
        boundary_dir=directory("data/gb-model/etys-boundary-gis-data-mar25")
    params:
        url="https://api.neso.energy/dataset/997f4820-1ad4-499b-b1fe-4b8d3d7fbc72/resource/e914fcec-1dc9-4f1f-97e7-59c0d9521bea/download/etys-boundary-gis-data-mar25.zip",
        zip_file="etys-boundary-gis-data-mar25.zip",
    log:
        logs("retrieve_etys_boundary_data.log")
    resources:
        mem_mb=1000,
    run:
        print("Creating output directory...")
        os.makedirs("data/gb-model", exist_ok=True)
        
        print(f"Downloading from {params.url}...")
        result = subprocess.run([
            "curl", "-L", "-o", params.zip_file, params.url
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Download failed: {result.stderr}")
            raise Exception("Download failed")
        
        file_size = os.path.getsize(params.zip_file)
        print(f"Download completed. File size: {file_size} bytes")
        
        print("Extracting zip file...")
        with ZipFile(params.zip_file, 'r') as zip_ref:
            zip_ref.extractall(output.boundary_dir)
        
        if not os.path.exists(output.boundary_shp):
            print("ERROR: Expected shapefile not found!")
            raise Exception("Extraction failed")
        
        os.remove(params.zip_file)
        print("Download and extraction successful!")


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
