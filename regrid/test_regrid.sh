#!/bin/bash
# Run a test of the OSKAR regridder
cd ./oskar_run_stage

# OSKAR Basic Command - add to .bashrc
function oskar() {

    gflag=0
    prog=""
    file=""
    outf="./"
    prevd=$PWD
    
    while [ $# -gt 0 ]; do
        case $1 in
            -g | --global | -s | --sample)
                gflag=1
            ;;
            -l | --local)
                gflag=0
            ;;
            -i | --intf)
                prog="oskar_sim_interferometer"
            ;;
            -b | --beam)
                prog="oskar_sim_beam_pattern"
            ;;
            -I | --image)
                prog="oskar_imager"
            ;;
            -f | --file)
                file=$2
                shift
            ;;
            -o | --output)
                outf=$2
                shift
            ;;
            -c | --clean)
                find ~/.oskar -name '*.log' -type f -delete
                return 0
            ;;
            \?)
                
            ;;
        esac
        shift
    done

    if [ $gflag -eq 1 ]; then
        file="$prog.ini"

        cd ~/.oskar
    fi

    echo $file

    singularity exec --nv --bind $PWD --cleanenv --home $PWD ~/.oskar/OSKAR-2.8.3-Python3.sif $prog $file

    cd $prevd
    cp -rf ~/.oskar/output $outf
}

# Run beamformer
oskar -l -b -f oskar_beam.ini

# Define the presets
test_presets=("yuxiang1" "yuxiang2")

# Interferometer and image
for preset in ${test_presets[@]}; do
    # Move OSM folder to directory
    cp -r "../regrid/${preset}_osm" "test_${preset}_osm"

    # Get first OSM file name
    cd "test_${preset}_osm"
    files=(*)
    cd ..

    for file in ${files[@]}; do
        # Generate the INI files
        cp ../regrid/test_intif_inis/test_intif_gen.ini "test_intif_${preset}.ini"

        # Replace sky model location in INI file
        ofname="oskar_sky_model\/file=test_${preset}_osm\/${file}"
        sed -i "s/^preset.*/${ofname}/" "test_intif_${preset}.ini"

        # Replace frequency bin
        word=$(sed '5!d' "test_${preset}_osm/$file")
        ofname="start_frequency_hz=${word:45:9}"
        sed -i "s/^fset.*/${ofname}/" "test_intif_${preset}.ini"

        # Then, run the interferometry and imager simulations
        oskar -l -i -f "test_intif_${preset}.ini"
        oskar -l -I -f oskar_imager.ini

        # Copy output data to regrid folder
        cp -rf output/sim_image_I.fits "../regrid/test_output/${preset}_fits/${file}.fits"

        # Clear OSM folder and INI from directory
        rm "test_intif_$preset.ini"
    done

    rm -r "test_${preset}_osm"
done

cd ..