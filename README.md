# novel proteins in PPI data
The code in this repo was created for the identification altProts in high throughput AP-MS data (BioPlex 2.0) and relevant network analysis. It accompanies our article: [Newfound coding potential of transcripts unveils missing members of human protein communities](https://www.biorxiv.org/content/10.1101/2020.12.02.406710v1).

## analytical pipeline
First the raw MS data was reanalyzed using the SearchGUI PeptideShaker tools as described in the manuscript. The code here starts with the hierarchical reports outputed by PeptideShaker.
1. Parse hierarchical reports ()
2. Run CompPASS
3. Compute CompPASS Plus features
4. run CompPASS Plus
5. Assemble network
6. Network topological analysis
7. Clustering and functional analysis

## Final notes
These scripts and notebooks were run on a highperformance computing platform with 24 cores and 256G of memory. Some adjustments may be necessary when less resources are available.
