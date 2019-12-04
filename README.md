# adaptive_knockoff_paper

This is a repository for exactly reproducing the simulation results in the paper Adaptive Knockoffs. The scripts `simulaition1.R` and `simulation2.R` in the `R`folder correspond to the simulations in Section 5.2 and Section 5.3 respectively. These R scripts can be executed on a personal computer and each script generates one realization of the simulation.

In the paper, each experiment is repeated for 100 times. To reproduce that, one needs to use the bash files in the `bash` folder. The linux commands are  `bash bash_simulation1.sh` and `bash bash_simulation2.sh`.

We also provide a R package `adaptiveKnockoff`, which will be constantly updated. To use the lastst version of the package, visit <https://github.com/zhimeir/akn>.
