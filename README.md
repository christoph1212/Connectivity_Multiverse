# 🚧 Welcome 🚧
**This Repo is under construction**

This is the github repository for the study [Robustness of EEG Functional Brain Networks Associated With Fluid Intelligence: A Multiverse Analysis of Connectivity and Thresholding Methods](). The preregistration of this study can be found [here]().

To reproduce the results from the paper, simply clone the repository and request the data from the authors.

## Dependencies & Credits
This project uses functions and code from the following third-party toolboxes. Please refer to the respective repositories for license details. Code from the Orthogonal Minimum Spanning Tree repository (OMST; Dimitriadis et al., 2017) is licenced under GPL v3, which requires this project to adopt the same license.
- Brain Connectivity Toolbox (Rubinov & Sporns, 2010). see [Homepage](https://sites.google.com/site/bctnet/)
- EEGLab (Delorme & Makeig, 2004). see [Git Repo](https://github.com/sccn/eeglab)
- Efficiency Cost Optimization (ECO; De Vico Fallani et al., 2017). see [Git Repo](https://github.com/devuci/3n)
- Orthogonal Minimum Spanning Tree (OMST; Dimitriadis et al., 2017). see [Git Repo](https://github.com/stdimitr/multi-group-analysis-OMST-GDD)
- SmallWorldNess (Humphries & Gurney, 2008). see [Git Repo](https://github.com/mdhumphries/SmallWorldNess)
- Wavelet Enhanced ICA (wICA; Castellanos & Makarov, 2006). see [Git Repo](https://github.com/Masoud-Ghodrati/wICA)

## Analysis Pipeline
Make sure to stay in the directory of the script you run (e.g. Matlab script $\rightarrow$ Matlab directory).
1. Navigate to the `/Matlab` folder. Adapt the folder paths according to your system. You can change the analysis settings if you want to reproduce the code. Leave it unchanged to match with the study. This script automatically preprocesses and epochs the data, calculates connectivity measures, applied thresholding, and calculates graph-theoretic metrics. Depending on your machine and the number of cores, this will take quite a while.
2. Change to the `/R` folder. Run `XYZ.R` for statistical analysis

## References
Castellanos, N. P. & Makarov, V. A. (2006). Recovering EEG brain signals: Artifact suppression with wavelet enhanced independent component analysis. *Journal Of Neuroscience Methods, 158*(2), 300–312. https://doi.org/10.1016/j.jneumeth.2006.05.033

Delorme, A., & Makeig, S. (2004). EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. *Journal of Neuroscience Methods, 134*(1), 9–21. https://doi.org/10.1016/j.jneumeth.2003.10.009

De Vico Fallani, F., Latora, V. & Chavez, M. (2017). A Topological Criterion for Filtering Information in Complex Brain Networks. *PLoS Computational Biology, 13*(1), e1005305. https://doi.org/10.1371/journal.pcbi.1005305

Dimitriadis, S. I., Salis, C., Tarnanas, I., & Linden, D. E. (2017). Topological Filtering of Dynamic Functional Brain Networks Unfolds Informative Chronnectomics: A Novel Data-Driven Thresholding Scheme Based on Orthogonal Minimal Spanning Trees (OMSTs). *Frontiers in Neuroinformatics, 11*, 28. https://doi.org/10.3389/fninf.2017.00028

Humphries, M. D. & Gurney, K. (2008). Network ‘Small-World-Ness’: A Quantitative Method for Determining Canonical Network Equivalence. *PLoS ONE, 3*(4), e0002051. https://doi.org/10.1371/journal.pone.0002051

Rubinov, M. & Sporns, O. (2010). Complex network measures of brain connectivity: Uses and interpretations. *NeuroImage, 52*(3), 1059–1069. https://doi.org/10.1016/j.neuroimage.2009.10.003
