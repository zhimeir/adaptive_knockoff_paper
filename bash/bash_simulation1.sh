#!/bin/bash
for SEED in {1..100}
do
    Rscript ../simulations/simulation1.R $SEED
done