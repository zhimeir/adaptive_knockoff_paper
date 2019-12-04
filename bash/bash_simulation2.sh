#!/bin/bash
for SEED in {1..100}
do
    Rscript ../simulations/simulation2.R $SEED
done