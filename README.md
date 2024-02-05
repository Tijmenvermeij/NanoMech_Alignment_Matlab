## NanoMechanical Alignment Framework
A nanomechanical alignment framework for alignment of SEM-DIC strain data and microstructure data

![SSLIP](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs11340-022-00884-0/MediaObjects/11340_2022_884_Fig3_HTML.png)

# Introduction to Alignment Framework

This Matlab framework allows alignment between different datasets, for example: EBSD data, BSE/SE data, DIC data, etc. The framework functions through selection of commong features in two datasets, followed by fitting of a distortion field to align the data.
More details can be found in [**this paper**](https://doi.org/10.1007/s11340-022-00884-0).

The **Alignment framework** is written in [**MATLAB**](https://mathworks.com/products/matlab.html) and uses several functionalities of the MATLAB-based crystallographic toolbox [**MTEX**](https://mtex-toolbox.github.io).

The **Alignment Methodology** and plotting functionalities are highlighted in two example scripts that showcase how the functions work and what their output comprises. The data needs to be downloaded and put in the data folders.
The Zn data can be found [here](https://www.dropbox.com/scl/fo/q9mi4eoujh25emnq31d41/h?rlkey=90a7turjxldqsrmtebs9y5rbf&dl=1).
The Ni based super alloy data can be found [here](https://www.dropbox.com/scl/fo/q9mi4eoujh25emnq31d41/h?rlkey=90a7turjxldqsrmtebs9y5rbf&dl=1).

Make sure to follow the **sections** in the example script step by step, because options need to be specified at various locations in the script.

Please report any bugs you encounter.

# Authors
**The alignment framework** has mainly been conceptualized and created by [**Tijmen Vermeij**](https://www.tue.nl/en/research/researchers/tijmen-vermeij/), under supervision of **Johan Hoefnagels**. Other contributors to the code: Jorn Verstijnen, Gert-Jan Slokker and Casper Mornout.

# How to cite this work
If you have applied the alignment framework to your research, please cite this open-access paper as your reference:
[**T. Vermeij, J.A.C. Verstijnen, T.J.J. Ramirez y Cantador, B. Blaysat, J. Neggers, J.P.M. Hoefnagels. A Nanomechanical Testing Framework Yielding Front&Rear-Sided, High-Resolution, Microstructure-Correlated SEM-DIC Strain Fields, Experimental Mechanics: 62, 1625-1646. (2022)**](https://doi.org/10.1007/s11340-022-00884-0).
