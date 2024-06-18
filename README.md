# Name

Geometry-based kinetic simulation of DNA-gold nanoparticle motor
https://github.com/TakanoriHarashima/DNA_NP_motor_simulation/blob/main/DNAmotor_simu_v5.02.py

# Manuscript

bioRxiv: https://www.biorxiv.org/content/10.1101/2024.05.23.595615v2

# Features

Details of the algorithm and the simulation parameters are described in the manuscript (https://www.biorxiv.org/content/10.1101/2024.05.23.595615v2). First, a two-dimensional pixel matrix was defined to model the RNAs on the glass surface. Pixel size was normalized by the DNA density on the AuNP. RNAs were randomly distributed with a ratio of the RNA density to the DNA density. Reaction at a single RNA site is assumed to contain three sequential elementary steps: DNA/RNA hybridization, RNase H binding to DNA/RNA duplex, and RNA hydrolysis. The reaction proceeds within the accessible area of DNA with the radius of 28.0 nm. Rate constants of the three elementary steps, konDNA/RNA, k E, and kcatE are the simulation parameters. Dwell time for each elementary step at each RNA site is determined by a random sampling from the exponential distribution: P(τ) = kexp(-kt), and the site with the shortest dwell time is changed to the next state. Simulation steps proceeded in units of reaction events. Constraint of the motor position was considered by introducing mobile region of each DNA/RNA hybrid (25.2 nm). The motor position was determined randomly within the region where all mobile regions overlap. Unless all mobile regions overlapped, the motor position was fixed. 

# Enviroment Tested
  - Windows 10 Pro
  - Anaconda 3
  - python 3.9
  - Spyder IDE 5.4.1

# Requirement

| Package  | Version |
| ------------- | ------------- |
| pandas  | 1.4.4  |
| numpy  | 1.21.5  |
| scipy  | 1.9.1  |
| opencv-python  | 4.7.0.68  |
| matplotlib  | 3.5.2  |
| psutil  | 5.9.0  |
| glob2  | 0.7  |

# Installation Guide
Whole installation procedure takes approximately 1 hour on a normal desktop computer.
1. Install Anaconda:
  - https://www.anaconda.com/download

2. Open Spyder:
  - To run the bundled version of Spyder after installing it with Anaconda, the recommended method on Windows is to launch it via the Start menu shortcut. On other platforms, open Anaconda Navigator, scroll to Spyder under Home and click Launch.
  - https://docs.spyder-ide.org/current/installation.html

3. Install required python modules:
  - Open a console and install required python modules
  - pip install -r requirements.txt

5. Run the simulation:
  - Download DNAmotor_simu_v5.02.py from Github.
  - Spyder -> File -> Open -> Select DNAmotor_simu_v5.02.py
  - Set parameters.
  - Run File (push F5 key) and simulation starts.

# Demo
1. Setting
See the manuscript and parameters are described in Supplementary Table 3.  
Just for a demo, here is the recommended parameters  
```
# Basic parameters
globalfol = r'Directory\to\perform\DNAmotor\simulation'
date = '20240618'
RNA_size = 2000                        # Full range of RNA substrate
tmax_simu = 100000                     # Uplimit for simulation (sec) 
N_simu = 5                             # Nunber of trajectory generate
frame_per_event = 1000                 # Span frame to check the progress of simulation
foli=0                                 # ID of the condition of kinetic parameters
# Kinetic parameters
khyb_list = np.array([0.3])                  # DNA/RNA Hybridization rate [s-1]
konE_list = np.array([1.0]) *10**6           # RNase H binding rate [M-1 s-1]
kcat_list = np.array([4])                    # RNA hydrolysis rate [s-1] 
RNaseH_list = np.array([36])                 # RNase H condition [nM]
```
2. run the simulation
push F5 key to run the simulation
3. Expected output
The simulation make folders as follow.
001_khyb=0.30_kcatE=4.0_konE=1.0x106.
└─20240618_{'N_simu'= 5, 'tmax'= 100000, 'RNaseH'= 36, 'frame_per_event'= 1000}
    └─progress
        ├─000
        ├─001
        ├─002
        ├─003
        └─004

# Author

* Takanori Harashima
* Institute for Molecular Science, National Institutes of Natural Sciences, Okazaki, Aichi 444-8787, Japan
* harashima@ims.ac.jp

# License

Our simulation code used in this study is provided under the MIT License.

