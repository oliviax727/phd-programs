# Run a test of the OSKAR regridder
cd ./test_intifs

function oskar() {
    echo "Running $2 on $4"
}

# First, run the beamformer
oskar -l -b -f test_beam.ini

# Create output directory
mkdir -p test_output

# Define the presets
test_presets=("gaussian" "sinusoid" "flat" "point")

# Interferometer and image
for preset in ${test_presets[@]}; do
    # Then, run the interferometry simulation
    cp test_intf_gen.ini "test_intf_$preset.ini"
    ofname="oskar_sky_model\/file=test_${preset}_osm\/reformatted_no.1_177.500MHz.osm"
    sed -i "s/^preset.*/${ofname}/" "test_intf_$preset.ini"
    oskar -l -i -f "test_intf_$preset.ini"
    oskar -l -I -f test_image.ini
    cp -r output/sim.ms "test_output/sim_$preset.ms"
    cp output/sim_image.fits "test_output/sim_image_$preset.fits"
    rm "test_intf_$preset.ini"
    rm -r output
done

cd ..