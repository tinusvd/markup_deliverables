# Markup languages and reproducible programming in statistics

Deliverable 1: Reproducible research compendium.

This research compendium allows for replication of Bayes Factors and Posterior Model Probabilities that were used in one of my assignments. While creating this compendium, I aimed for making the replication as easy as possible. First of all, the user can find the to-be-replicated table in the `output` folder. There, a file named `Table_2.pdf` holds the table. The output file also holds a `renv` `.lock` file which should be loaded into the R-environment for optimal replicability.

## Steps
1. Open the file named `gibbs_sampler.R`. Run the first line. 
2. Read the message given by the `restore` function of `renv` Select `1` in the R console to replicate the environment used to generate the results. Select `2` or `3` to continue with your current libraries. Note: selecting `2` or `3` might lead to failure of exact replication.
3. Run the remainder of the `gibbs_sampler.R` script.
4. Check the results.

If you have any questions regarding this project, please contact me at: m.l.g.vandam@uu.nl. 

