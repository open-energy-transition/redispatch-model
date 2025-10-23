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
        raw_region_shapes=resources("gb-model/raw_region_shapes.geojson"),
    log:
        logs("raw_region_shapes.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/create_region_shapes.py"


# Rule to manually merge raw_region_shapes
rule manual_region_merger:
    input:
        raw_region_shapes=rules.create_region_shapes.output.raw_region_shapes,
        country_shapes=resources("country_shapes.geojson"),
    output:
        merged_shapes=resources("gb-model/merged_shapes.geojson"),
    log:
        logs("manual_region_merger.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/manual_region_merger.py"


# Rule to retrieve generation unit unavailability data from ENTSO-E
rule retrieve_entsoe_unavailability_data:
    message:
        "Retrieve data from ENTSOE API for generator {wildcards.business_type} unavailability in {wildcards.zone} bidding zone"
    output:
        xml_base_dir=directory("data/gb-model/entsoe_api/{zone}/{business_type}"),
    params:
        start_date=config["entsoe_unavailability"]["start_date"],
        end_date=config["entsoe_unavailability"]["end_date"],
        bidding_zones=config["entsoe_unavailability"]["bidding_zones"],
        business_types=config["entsoe_unavailability"]["business_types"],
        max_request_days=config["entsoe_unavailability"]["max_request_days"],
        api_params=config["entsoe_unavailability"]["api_params"],
    log:
        logs("retrieve_entsoe_unavailability_data_{zone}_{business_type}.log"),
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/retrieve_entsoe_unavailability_data.py"


rule process_entsoe_unavailability_data:
    input:
        xml_base_dir="data/gb-model/entsoe_api/{zone}/{business_type}",
    output:
        unavailability=resources(
            "gb-model/{zone}_{business_type}_generator_unavailability.csv"
        ),
    log:
        logs("process_entsoe_unavailability_data_{zone}_{business_type}.log"),
    params:
        business_type_codes=config["entsoe_unavailability"]["api_params"][
            "business_types"
        ],
    resources:
        mem_mb=1000,
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/process_entsoe_unavailability_data.py"


rule generator_monthly_unavailability:
    input:
        planned=resources("gb-model/{zone}_planned_generator_unavailability.csv"),
        forced=resources("gb-model/{zone}_forced_generator_unavailability.csv"),
        powerplants=resources("powerplants_s_all.csv"),
    params:
        carrier_mapping=config["entsoe_unavailability"]["carrier_mapping"],
        resource_type_mapping=config["entsoe_unavailability"]["resource_type_mapping"],
        start_date=config["entsoe_unavailability"]["start_date"],
        end_date=config["entsoe_unavailability"]["end_date"],
        max_unavailable_days=config["entsoe_unavailability"]["max_unavailable_days"],
    output:
        csv=resources("gb-model/{zone}_generator_monthly_unavailability.csv"),
    log:
        logs("{zone}_generator_monthly_unavailability.log"),
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/generator_monthly_unavailability.py"


rule extract_transmission_availability:
    input:
        pdf_report="data/gb-model/downloaded/transmission-availability.pdf",
    output:
        csv=resources("gb-model/transmission_availability.csv"),
    log:
        logs("extract_transmission_availability.log"),
    conda:
        "../envs/gb-model/workflow.yaml"
    script:
        "../scripts/gb_model/extract_transmission_availability.py"


rule extract_fes_workbook_sheet:
    message:
        "Extract FES workbook sheet {wildcards.fes_sheet} for FES-{wildcards.fes_year} and process into machine-readable, 'tidy' dataframe format according to defined configuration."
    input:
        workbook="data/gb-model/downloaded/fes-{fes_year}-workbook.xlsx",
    output:
        csv=resources("gb-model/fes/{fes_year}/{fes_sheet}.csv"),
    params:
        sheet_extract_config=lambda wildcards: config["fes-sheet-config"][
            int(wildcards.fes_year)
        ][wildcards.fes_sheet],
    log:
        logs("extract_fes-{fes_year}_{fes_sheet}.log"),
    script:
        "../scripts/gb_model/extract_fes_sheet.py"


rule process_fes_eur_data:
    message:
        "Process FES-compatible European scenario workbook."
    params:
        scenario=config["fes"]["eur"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        countries=config["countries"],
    input:
        eur_supply="data/gb-model/downloaded/eur-supply-table.csv",
    output:
        csv=resources("gb-model/national_eur_data.csv"),
    log:
        logs("process_fes_eur_data.log"),
    script:
        "../scripts/gb_model/process_fes_eur_data.py"


rule process_fes_gsp_data:
    message:
        "Process FES workbook sheet BB1 together with metadata from sheet BB2."
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
    input:
        bb1_sheet=resources("gb-model/fes/2021/BB1.csv"),
        bb2_sheet=resources("gb-model/fes/2021/BB2.csv"),
        gsp_coordinates="data/gb-model/downloaded/gsp-coordinates.csv",
        regions=resources("gb-model/merged_shapes.geojson"),
    output:
        csv=resources("gb-model/regional_gb_data.csv"),
    log:
        logs("process_fes_gsp_data.log"),
    script:
        "../scripts/gb_model/process_fes_gsp_data.py"


rule create_powerplants_table:
    message:
        "Tabulate powerplant data GSP-wise from FES workbook sheet BB1 and EU supply data"
    params:
        gb_config=config["fes"]["gb"],
        eur_config=config["fes"]["eur"],
        default_set=config["fes"]["default_set"],
    input:
        gsp_data=resources("gb-model/regional_gb_data.csv"),
        eur_data=resources("gb-model/national_eur_data.csv"),
    output:
        csv=resources("gb-model/fes_p_nom.csv"),
    log:
        logs("create_powerplants_table.log"),
    script:
        "../scripts/gb_model/create_powerplants_table.py"


rule create_interconnectors_table:
    input:
        regions=resources("gb-model/merged_shapes.geojson"),
    output:
        gsp_data=resources("gb-model/interconnectors_p_nom.csv"),
    params:
        interconnector_config=config["interconnectors"],
        year_range=config["fes"]["year_range_incl"],
    log:
        logs("create_interconnectors_table.log"),
    script:
        "../scripts/gb_model/create_interconnectors_table.py"


rule create_hydrogen_demand_table:
    message:
        "Process hydrogen demand data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_demand_sheets=config["fes"]["hydrogen"]["demand"]["annual_demand_sheets"],
        other_sectors_list=config["fes"]["hydrogen"]["demand"]["other_sectors_list"],
    input:
        demand_sheets=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["demand"][
                "annual_demand_sheets"
            ].items()
            for sheet in sheets.values()
        ],
    output:
        hydrogen_demand=resources("gb-model/fes_hydrogen_demand.csv"),
    log:
        logs("create_hydrogen_demand_table.log"),
    script:
        "../scripts/gb_model/create_hydrogen_demand_table.py"


rule create_grid_electrolysis_table:
    message:
        "Process hydrogen electrolysis data from FES workbook into CSV format"
    input:
        regional_gb_data=resources("gb-model/regional_gb_data.csv"),
    output:
        grid_electrolysis_capacities=resources(
            "gb-model/fes_grid_electrolysis_capacities.csv"
        ),
    log:
        logs("create_grid_electrolysis_table.log"),
    script:
        "../scripts/gb_model/create_grid_electrolysis_table.py"


rule create_hydrogen_supply_table:
    message:
        "Process hydrogen supply data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_supply_sheets=config["fes"]["hydrogen"]["supply"]["supply_sheets"],
        exogeneous_supply_list=config["fes"]["hydrogen"]["supply"][
            "exogeneous_supply_list"
        ],
    input:
        supply_sheets=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["supply"][
                "supply_sheets"
            ].items()
            for sheet in sheets.values()
        ],
    output:
        hydrogen_supply=resources("gb-model/fes_hydrogen_supply.csv"),
    log:
        logs("create_hydrogen_supply_table.log"),
    script:
        "../scripts/gb_model/create_hydrogen_supply_table.py"


rule create_off_grid_electrolysis_demand:
    message:
        "Process electricity demand of off-grid electrolysis from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_supply_sheets=config["fes"]["hydrogen"]["supply"]["supply_sheets"],
    input:
        supply_sheets=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["supply"][
                "supply_sheets"
            ].items()
            for sheet in sheets.values()
        ],
        grid_electrolysis_capacities=resources(
            "gb-model/fes_grid_electrolysis_capacities.csv"
        ),
    output:
        electricity_demand=resources(
            "gb-model/fes_off_grid_electrolysis_electricity_demand.csv"
        ),
    log:
        logs("create_off_grid_electrolysis_demand.log"),
    script:
        "../scripts/gb_model/create_off_grid_electrolysis_demand.py"


rule create_hydrogen_storage_table:
    message:
        "Process hydrogen storage data from FES workbook into CSV format"
    params:
        scenario=config["fes"]["gb"]["scenario"],
        year_range=config["fes"]["year_range_incl"],
        fes_storage_sheets=config["fes"]["hydrogen"]["storage"]["storage_sheets"],
        interpolation_method=config["fes"]["hydrogen"]["storage"][
            "interpolation_method"
        ],
    input:
        storage_sheet=lambda wildcards: [
            resources(f"gb-model/fes/{year}/{sheet}.csv")
            for year, sheets in config["fes"]["hydrogen"]["storage"][
                "storage_sheets"
            ].items()
            for sheet in sheets.values()
        ],
    output:
        hydrogen_storage=resources("gb-model/fes_hydrogen_storage.csv"),
    log:
        logs("create_hydrogen_storage_table.log"),
    script:
        "../scripts/gb_model/create_hydrogen_storage_table.py"


rule compose_network:
    input:
        unpack(input_profile_tech),
        network=resources("networks/base_s_{clusters}.nc"),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
        intermediate_data=[
            resources("gb-model/transmission_availability.csv"),
            expand(
                resources(
                    "gb-model/{zone}_{business_type}_generator_unavailability.csv"
                ),
                zone=config["entsoe_unavailability"]["bidding_zones"],
                business_type=config["entsoe_unavailability"]["business_types"],
            ),
            resources("gb-model/merged_shapes.geojson"),
            resources("gb-model/fes_p_nom.csv"),
            resources("gb-model/interconnectors_p_nom.csv"),
            resources("gb-model/GB_generator_monthly_unavailability.csv"),
            resources("gb-model/fes_hydrogen_demand.csv"),
            resources("gb-model/fes_grid_electrolysis_capacities.csv"),
            resources("gb-model/fes_hydrogen_supply.csv"),
            resources("gb-model/fes_off_grid_electrolysis_electricity_demand.csv"),
            resources("gb-model/fes_hydrogen_storage.csv"),
        ],
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
        "../scripts/gb_model/compose_network.py"


rule compose_networks:
    input:
        expand(
            resources("networks/composed_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
