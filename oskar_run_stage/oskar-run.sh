#!/bin/bash

source /software/projects/mwaeor/ohrw/install-scripts/oskar-install.sh

/software/projects/mwaeor/ohrw/.oskar/bin/oskar_sim_interferometer oskar_sim_interferometer.ini > test-intif.out 2>&1

/software/projects/mwaeor/ohrw/.oskar/bin/oskar_imager oskar_imager.ini > test-img.out 2>&1
