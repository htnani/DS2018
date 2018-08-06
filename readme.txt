Instructions to reproduce the experiments in Predicting Rice Phenotypes with Meta-Learning.
Please address all queries to: oghenejokpeme.orhobor *at* manchester *dot* ac *dot* uk

The three folders correspond to the three main experiments.
1. base_experiments -> base case (group and non-group)
2. framework_a -> proposed meta-learning framework A
3. framework_b -> proposed meta-learning framework B

Steps.
1. First download the genotype (imputed) and phenotype files from http://dx.doi.org/10.17632/86ygms76pb.1
2. Copy the genotype and phenotype data to their corresponding folder for the experiment one is interested in replicating.
3. Run the R scripts in the "code" folder for the experiment one is interested in replicating. *See Notes below*
The results will be in "/logs/results/" for the framework A and B experiments and in "/output/" for the base experiments.

Notes.
1. Ensure the libraries used in each file are installed.
2. Run the scripts for each experiment in the order in which they've been numbered. Especially for the meta-learning experiments.
3. The experiments were written to be run on a system which supports parallel compute. One can change the number of cores by searching for "registerDoMC" in each script before execution.