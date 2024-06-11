# Name

Geometry-based kinetic simulation of DNA-gold nanoparticle motor
https://github.com/TakanoriHarashima/DNA_NP_motor_simulation/blob/main/DNAmotor_simu_v5.02.py

# Manuscript

bioRxiv: https://www.biorxiv.org/content/10.1101/2024.05.23.595615v2

# Features

Details of the algorithm and the simulation parameters are described in the manuscript (https://www.biorxiv.org/content/10.1101/2024.05.23.595615v2). First, a two-dimensional pixel matrix was defined to model the RNAs on the glass surface. Pixel size was normalized by the DNA density on the AuNP. RNAs were randomly distributed with a ratio of the RNA density to the DNA density. Reaction at a single RNA site is assumed to contain three sequential elementary steps: DNA/RNA hybridization, RNase H binding to DNA/RNA duplex, and RNA hydrolysis. The reaction proceeds within the accessible area of DNA with the radius of 28.0 nm. Rate constants of the three elementary steps, konDNA/RNA, k E, and kcatE are the simulation parameters. Dwell time for each elementary step at each RNA site is determined by a random sampling from the exponential distribution: P(τ) = kexp(-kt), and the site with the shortest dwell time is changed to the next state. Simulation steps proceeded in units of reaction events. Constraint of the motor position was considered by introducing mobile region of each DNA/RNA hybrid (25.2 nm). The motor position was determined randomly within the region where all mobile regions overlap. Unless all mobile regions overlapped, the motor position was fixed. 

# Enviroment

python3.9

# Requirement

Package                       Version
------------------            ---------
pandas                        1.4.4
numpy                         1.21.5
scipy                         1.9.1
opencv-python                 4.7.0.68
matplotlib                    3.5.2
psutil                        5.9.0
glob2                         0.7

# Author

* Takanori Harashima
* Institute for Molecular Science
* harashimna@ims.ac.jp

# License
It is made available under a CC-BY-NC-ND 4.0 International license.

