# novel proteins in PPI data
The code in this repo was created for the identification altProts in high throughput AP-MS data (BioPlex 2.0) and relevant network analysis. It accompanies our article: [Newfound coding potential of transcripts unveils missing members of human protein communities](https://www.biorxiv.org/content/10.1101/2020.12.02.406710v1).

## analytical pipeline
First the raw MS data was reanalyzed using the SearchGUI PeptideShaker tools as described in the manuscript. The code here starts with the hierarchical reports outputed by PeptideShaker.
1. Parse hierarchical reports<br/>
The notebook parse_psm_reports.ipynb aggregates all peptide-spectrum matches (PSMs) from the PeptideShaker hierarchical reports and produces the csv file for input to CompPASS.
3. Run CompPASS<br/>
PSM counts from previous step are ran through CompPASS with R.
5. Compute CompPASS Plus features<br/>
The ouput from CompPASS is used to compute 9 features representing each candidate interaction (compute_CompPASS_Plus_features.ipynb). Satistical filters are also applied, as described in [Huttlin et.al.](https://doi.org/10.1016/j.cell.2015.06.043).
7. run CompPASS Plus<br/>
The features are used to train naive bayes classifiers in cross-validation splited by batch of AP-MS experiments. (run_CompPASS_Plus.ipynb)
9. Assemble network<br/>
The scored interactions are filtered, assembled into a network and compared with the BioPlex networks. (assemble_network.ipynb)
11. Network topological analysis<br/>
Topological features or the resulting network are computed and visualized (degree distribution, average shortest paths, eigenvector centralites etc.). (full_network_features.ipynb)
13. Clustering and functional analysis<br/>
The network is partitioned with the markov clustering algorithm. Clusters are analyzed for enrchiment of Gene Ontology terms (clustering_GO.ipynb). Disease association are also computed for each cluster (disease_associations.ipynb).

## Final notes
These scripts and notebooks were run on a highperformance computing platform with 24 cores and 256G of memory. Some adjustments may be necessary when less resources are available.
