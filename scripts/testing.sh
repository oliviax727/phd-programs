#!/bin/bash

export OUTFILE="/scratch/mwaeor/ohrw/test.out"

(
	> $OUTFILE
	./scripts/chips-testing.py 2>&1 > $OUTFILE
	./scripts/cube-testing.py 2>&1 > $OUTFILE
)
