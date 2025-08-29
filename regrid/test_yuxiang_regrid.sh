#!/bin/bash
# Regrd yuxiang's data
cd ./regrid
python3 oskar-regrid.py

# Run a test of the OSKAR regridder
cd ./test_intifs

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

    if [ $gflag ]; then
        file="$prog.ini"

        cd ~/.oskar
    fi

    singularity exec --nv --bind $PWD --cleanenv --home $PWD ~/.oskar/OSKAR-2.8.3-Python3.sif $prog $file

    cd $prevd
    cp -rf ~/.oskar/output $outf
}

# First, run the beamformer
oskar -l -b -f test_beam.ini

# Create output directory
mkdir -p test_yuxiang_output

# Define the presets
test_presets=("yuxiang1" "yuxiang2")

# Interferometer and image
for preset in ${test_presets[@]}; do
    # Get first OSM file name
    cd "../$preset_osm"
    files=(*)
    file=${files[0]}
    cd "../test_intfs"

    # Then, run the interferometry simulation
    cp test_intf_gen.ini "test_intf_$preset.ini"
    ofname="oskar_sky_model\/file=..\/${preset}_osm\/${file}"
    sed -i "s/^preset.*/${ofname}/" "test_intf_$preset.ini"

    oskar -l -i -f "test_intf_$preset.ini"
    oskar -l -I -f test_image.ini
    
    cp -r output/sim.ms "test_output/sim_$preset.ms"
    cp output/sim_image_I.fits "test_output/sim_image_${preset}_I.fits"
    rm "test_intf_$preset.ini"
    rm -r output
done

cd ..
cd ..