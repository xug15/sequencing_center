#!/bin/bash
SHOREmap backcross --chrsizes ~/genome4ann/chrSizes.txt --marker wt/1_converted_variant.txt --consen wt/1_converted_consen.txt --folder SHOREmap_analysisi.wt/ -plot-bc --marker-score 20 --marker-freq 0.2222 --min-coverage 20 --max-coverage 120 --marker-hit 1 --bg bg/1_converted_variant.txt --bg-score 1 --bg-freq 0.0 --bg-cov 1 -non-EMS  --cluster 1 -verbose
