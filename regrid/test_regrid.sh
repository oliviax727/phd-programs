#!/bin/bash
# Run a test of the OSKAR regridder
cd ./oskar_run_stage

# Source oskar.bashrc
source ../oskar.bashrc

# Set up run stage
rm -rf output/*
rm -rf .nv
rm -rf oskar_*.log
rm -rf test_*.ini
rm -rf test_*_osm

# Define the presets
test_presets=("yuxiang1_zenith")

# Interferometer and image
for preset in ${test_presets[@]}; do
    # Move OSM folder to directory
    cp -r "../regrid/osm_output/${preset}_osm" "./test_${preset}_osm"

    # Get first OSM file name
    cd "test_${preset}_osm"
    files=(*)
    cd ..

    # Create output fits directory
    rm -r "../regrid/fits_output/${preset}_fits"
    mkdir -p "../regrid/fits_output/${preset}_fits"

    for file in ${files[@]}; do
        # Generate the INI files
        cp ../regrid/test_intif_inis/test_intif_gen_bash.ini "test_intif_${preset}.ini"
        cp ../regrid/test_intif_inis/test_img_gen_bash.ini "test_img_${preset}.ini"

        # Replace sky model location in INI file
        ofname="oskar_sky_model\/file=test_${preset}_osm\/${file}"
        sed -i "s/^fileset.*/${ofname}/" "test_intif_${preset}.ini"

        # Replace frequency bin
        word=$(sed '5!d' "test_${preset}_osm/${file}")
        ofname="start_frequency_hz=${word:45:9}"
        sed -i "s/^freqset.*/${ofname}/" "test_intif_${preset}.ini"

        # Then, run the interferometry and imager simulations
        oskar_bash -g -i -f "test_intif_${preset}.ini" -c
        oskar_bash -g -I -f "test_img_${preset}.ini" -c

        # Copy output data to regrid folder
        cp output/sim_image_I.fits "../regrid/fits_output/${preset}_fits/${file}.fits"

        # Clear OSM folder and INI from directory
        rm "test_intif_${preset}.ini"
        rm "test_img_${preset}.ini"
    done

    rm -r "test_${preset}_osm"
done

cd ..