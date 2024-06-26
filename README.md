# Network Analysis 
This repository contains the code used for data analysis in the manuscript titled *"Longitudinal evolution of the transdiagnostic prodrome to severe mental disorders: a dynamic temporal network analysis informed by natural language processing and electronic health records"* by Arribas et al (under review, pre-print: https://www.medrxiv.org/content/10.1101/2024.03.08.24303965v1). The purpose of this project is to provide a **transparent** and **reproducible** method for analyzing the data presented in our publication, as well as making our code publically avaiable for researchers interested in using network analyses with their own databases.

Network analysis allows us to statistically model the relationships between nodes connected by edges (directed and undirected relationships) within a network and how these change over time. Weak, sparsely connected networks are more modifiable, while stronger, densely connected networks resist change and need greater disruption through higher intensity interventions to alter the state. Furthermore, edge estimates for temporal nodes may indicate directed causality between network nodes. Similarly, the centrality of a node, which represents the relative strength of connections in and out of the node, may signal the importance of that node to the system. 

Please refer to https://doi.org/10.1007/s11336-020-09697-3 for further information on this statistical approach. 

### Data Sharing
The data accessed by CRIS remain within an NHS firewall and governance is provided by a patient-led oversight committee. However, due to data sharing policy (see https://doi.org/10.1186/1471-244X-9-51 for futher details) we are unable to share the raw data. 

### Data Files
All the data is generated synthetically in the first script "0. Synth_data.R". The output file "dynamic_syn.csv" has exactly the same structure as our original dataset. 
Therefore, the analyses can be replicated with the synthetic dataset but the results will differ from our original results. 

Instead, we have provided actual data matrices in "features_covariances.RData" and "features_boot_temporal.RDS", as these are not individual-level data and can be made public. These files are loaded into "2. Network.R" and "3. Bootstrapping.R" respectively, so that the actual models from our manuscript can be visualised and the same figures can be replicated. 

## Getting Started
These instructions will guide you through setting up the project and running the analysis on your own machine.

### Installation
The scripts are all written in R. Please make sure you have installed the libraries which are loaded at the start of each script. 

### Running the Analysis
The code is divided into 6 parts. Each script will depend on data files outputted from previous scripts:

- "0. Synth_data.R" : code generates synthetic data and outputs it as "dynamic_syn.csv" file
- "1. Pre-processing.R" : requires the "dynamic_syn.csv" file. Pre-processes the data and outputs it as "residuals.csv" and "total_2yr.csv"
- "2. Network.R": builds the models using the "residuals.csv" file and visualises actual networks using "features_covariances.RData" 
- "3. Bootstrapping.R": model bootstrapping using the "total_2yr.csv" file. Actual model vs bootstrapped estimates can be visualised using "features_boot_temporal.RDS"
- "4. Permutation.R": requires "total_2yr.csv" and files provided in the /data/perm folder. Permuted results can be visualised as heatmaps.
- "5. Community.R": requires files provided in the /data/perm folder. Node-node covariance across communitues can be visualised as heatmaps.

## Contributing
Since this project is associated with a submitted manuscript, modifications might be limited. However, if you find any bugs or have suggestions, please open an issue in the repository.

## Manuscript Citation
Please cite our manuscript as follows if you use our code or analysis:

Arribas M, Barnby J, Patel R, et al. Longitudinal evolution of the transdiagnostic prodrome to severe mental disorders: a dynamic temporal network analysis informed by natural language processing and electronic health records. Published online March 9, 2024:2024.03.08.24303965. doi:10.1101/2024.03.08.24303965

## Authors
Maite Arribas - Initial Analysis and Development - m-arribas

## License
This project is licensed under the Apache License - see the LICENSE.md file for details.

## Acknowledgments
Special thanks to josephmbarnby for their repository [ParanoiaLongitudinalNetworkAnalysis]. This project's code has been largely adapted (with permission) from his work. Check out the original repository here: https://github.com/josephmbarnby/ParanoiaLongitudinalNetworkAnalysis.git. His contributions have been invaluable in the development of this project.

MA is supported by the UK Medical Research Council (MR/N013700/1) and King’s College London member of the MRC Doctoral Training Partnership in Biomedical Sciences. JMB has received funding from the Wellcome Trust (WT228268/Z/23/Z). 
