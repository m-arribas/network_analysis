<h1>Manuscript Data Analysis Project</h1>
This repository contains the code used for data analysis in the manuscript titled *"[Your Manuscript Title]"*. The purpose of this project is to provide a **transparent** and **reproducible** method for analyzing the data presented in our publication, as well as making code publically avaiable for researchers interested in using network analyses with their own databases.

Network analysis allows us to statistically model the relationships between nodes (prodromal features) connected by edges (directed and undirected relationships) within a network (prodrome) and how these change over time. Weak, sparsely connected networks are more modifiable, while stronger, densely connected networks resist change and need greater disruption through higher intensity interventions to alter the state (e.g. prevent disorder onset). Furthermore, edge estimates for temporal nodes may indicate directed causality between network features, which may aid understanding of the development of these disorders. Similarly, the centrality of a node, which represents the relative strength of connections in and out of the node, may signal the importance of that node to the disorder and as a potential intervention target.

Please refer to xxx for further information on this statistical approach. 

## Getting Started
These instructions will guide you through setting up the project and running the analysis on your own machine.

### Prerequisites

The code is divided into 4 parts. Each script will depend on data files outputted from previous scripts:

- "0. Synth_data.R" : code generates synthetic data and outputs the "dynamic_syn.csv" file
- "1. Pre-processing.R" : requires the "dynamic_syn.csv" file. Pre-processes the data and outputs it as "residuals.csv" and "total_2yr.csv"
- "2. Network.R": requires the "residuals.csv" and "features_covariances.RData" files
- "3. Bootstrapping.R": requires the "total_2yr.csv" and "features_boot_temporal.RDS" files

Before running the analysis, ensure you have the following software/tools installed:

### Installation
The scripts are all written in R. Please make sure you have installed the libraries which are loaded at the start of each script. 

### Running the Analysis
To replicate the analysis presented in the manuscript, follow these steps:

### Data Sharing
The data accessed by CRIS remain within an NHS firewall and governance is provided by a patient-led oversight committee. However, due to data sharing policy (see https://doi.org/10.1186/1471-244X-9-51 for futher details) we are unable to share the raw data. 

### Data Files
All the data is generated synthetically in the first script ""0. Synth_data.R". The output file "dynamic_syn.csv" has exactly the same structure as our original dataset. 
Therefore, the analyses can be replicated with the synthetic dataset but the results will differ from our original results. 

Instead, we have provided actual data matrices in "features_covariances.RData" and "features_boot_temporal.RDS", as these are not individual-level data and can be made public. These files are loaded into "2. Network.R" and "3. Bootstrapping.R" respectively, so that the actual models from our manuscript can be visualised and the same figures can be replicated. 

## Contributing
Since this project is associated with a published manuscript, modifications might be limited. However, if you find any bugs or have suggestions, please open an issue in the repository.

## Manuscript Citation
Please cite our manuscript as follows if you use our code or analysis:

## Authors
Maite Arribas - Initial Analysis and Development - m-arribas

License
This project is licensed under the MIT License - see the LICENSE.md file for details.

## Acknowledgments
Special thanks to josephmbarnby for their repository [ParanoiaLongitudinalNetworkAnalysis]. This project's code has been largely adapted (with permission) from their work. Check out the original repository here: https://github.com/josephmbarnby/ParanoiaLongitudinalNetworkAnalysis.git. Their contributions have been invaluable in the development of this project.

MA is supported by the UK Medical Research Council (MR/N013700/1) and King’s College London member of the MRC Doctoral Training Partnership in Biomedical Sciences. JMB has received funding from the Wellcome Trust (WT228268/Z/23/Z). 

