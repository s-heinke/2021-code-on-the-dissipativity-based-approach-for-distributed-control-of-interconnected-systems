This folder contains the files to reproduce the results presented in

S. Heinke and H. Werner, "On the Dissipativity Based Approach for 
Distributed Control of Interconnected Systems", 2021.

The data and figures can be reproduced using the file 'main.m'. 
The functions h2** design continuous time H2 output feedback controllers
and the functions **2MIMO convert the distributed models into lumped MIMO
models, where ** stands for
- AG: Arbitrary graph
- SI: Spatially invariant
- AD: Alpha-heterogeneous decomposable
- AS: Alpha-heterogeneous string interconnected

In order to run the code, the following Matlab packages are required:
- yalmip (https://yalmip.github.io/download/)

The code has been tested using:
- Matlab R2020a
- yalmip version 21-Nov-2017

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.