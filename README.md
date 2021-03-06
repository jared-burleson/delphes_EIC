# Quick Note

delphes_EIC is a project for using the fast detector-response simulator DELPHES to model potential EIC studies. The main github for the project is located at https://github.com/miguelignacio/delphes_EIC and https://github.com/stephensekula/delphes_EIC/.

My work is specifically looking at subjet studies, particular with the exist methods inside of DELPHES and perhaps applying to a multivariate algorithm in the future. I will use this repository to post code updates/changes to both DELPHES files and my own additions for gathering data. In addition, I have added a folder where I will upload my jupyter notebook python scripts that cover the data analysis portion of my work.

# delphes_EIC

Install Delphes3 following:
https://github.com/delphes/delphes

The detector card contains an EIC detector based on the EIC detector handbook v1.2
http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf

So far it incorporates tracking, EMCAL and HCAL but lacks implementation of PID (it can be done though, following the LHCb card example)

Magnetic field: 1.5 T, Solenoid length: 2.0 m, Tracker radius: 80 cm. 

You can run Pythia8 within Delphes. The command file shown here is suitable for DIS at EIC. 

Run generation command:
`./DelphesPythia8 cards/delphes_card_EIC.tcl examples/Pythia8/DIS.cmnd out.root`

You can see examples of analysis code in the Delphes page above

Run visualization command:
 `root -l examples/EventDisplay.C'("cards/delphes_card_EIC.tcl","out.root")'`
 
The two examples shown here are for neutral-current and charged-current event 
for beam energies of 10 GeV electron on 100 GeV proton (63 GeV center-of-mass energy). 


## Setting up the code


1. Install LHAPDF
   * https://lhapdf.hepforge.org/
   * Using LHAPDF 6.2.3
   * Download the tarball and unpack it
   * BUGFIX in 6.2.3: there is a python coding error in bin/lhapdf. Edit this file and find and replace "add_add_mutually_exclusive_group" with "add_mutually_exclusive_group"
   * Configure it for local installation in your work area, e.g. ```./configure --prefix=/users/ssekula/scratch/EIC/```
   * Build it, ```make -j```,
   * Install it, ```make install```,
   * Make sure the environment is set properly to find the binaries, libraries, and python code (c.f. https://lhapdf.hepforge.org/install.html#lxplus for examples)
1. Install PYTHIA8,
   * http://home.thep.lu.se/~torbjorn/Pythia.html,
   * Download the tarball and unpack it. ,
   * Configure it for local installation in your work area, e.g. ```./configure --prefix=/users/ssekula/scratch/EIC/ --with-lhapdf6=/scratch/users/ssekula/EIC/```,
   * Build it, ```make -j```,
   * Install it, ```make install```,
   * Make sure the work area binary directory is in your PATH: ```PATH=/users/ssekula/scratch/EIC/bin:${PATH}```,
1. Install Delphes,
   * https://github.com/delphes/delphes,
   * Clone the project and make sure you are on the master branch,
   * Make sure ROOT is available in your path, e.g. ```lsetup \"root 6.18.04-x86_64-centos7-gcc8-opt\"```,
   * Compile with PYTHIA8: ```HAS_PYTHIA8=true PYTHIA8=/users/ssekula/scratch/EIC ./configure --prefix=/users/ssekula/scratch/EIC/```,
   * Build: ```make -j```,
   * Install: ```make install```,
1. Get the Delphes/EIC code for simulation and analysis of a detector baseline/configuration.,
   * https://github.com/miguelignacio/delphes_EIC,
   * Clone the repository locally,
   * Follow the instructions to run the example and generate a ROOT file.

## Bug Fixing

There are occasionally bugs or other easily fixable issues that we have come across in our work that if not addressed can cause issues when using DELPHES.

The first is the make the following adjustment to the ```BeamRemnants.cc``` file inside of ```/pythiaXXXX/src/```. Add the following loop after line 978:

```if (iLepScat > (event.size()-1)) { return false; }```

The issue is that occasionally pythia produces a messed-up beam remnant and this small addition will require pythia to try again rather than causing a run-time error in the simulation generation.

The second has to do with the Jet Trimming Algorithm in DELPHES. DELPHES has FastJet implemented and it calls FastJet and performs a trimming analysis in ```/Delphes-X.X.X/modules/FastJetFinder.cc```. However, based on my studies I found that the trimmed subjets only ever contain a single constituent regardless of radius or pT cuts. Upon investigation, I determined that the following line of code is extraneous and should be commented out or removed from the file in order to allow for the correct trimming algorithm to be employed. This is located around line 490 in the fComputeTrimming section:

```fastjet::PseudoJet trimmed_jet = join(trimmed_jet.constituents());```

Once this line is commented out, the trimmed_jet can be analyzed for pieces and have substructure that contains clustering of more than one individual constituent.

UPDATE: On 07/22/2020 DELPHES accepted the pull request addressing the issue with FastJetFinder.cc and the trimming section. If downloading DELPHES straight from the github after the date, the above issue should be fixed.

## Running Monte Carlo Production

The following command line will run these options:

* e-p beam energies of 10 and 275, respectively
* LHAPDF6:CT18NNLO PDF set
* Output will be written to CC_DIS_e10_p275_CT18NNLO/0/ (where 0 is the task ID, specified using SLURM_ARRAY_TASK_ID)
* 100,000 events generated in a single run of the command

```
SLURM_ARRAY_TASK_ID=0 ./run_study.py --template delphes_card_EIC.tcl --commands CC_DIS_template.cmnd -p '{"PARAM_NEVENTS": 100000, "PARAM_HADBEAM_ENERGY": 275, "PARAM_EBEAM_ENERGY": 10, "PARAM_PDFSET": "LHAPDF6:CT18NNLO"}' -n CC_DIS_e10_p275_CT18NNLO
```

## EIC Collider Variations

Beam energy recommended benchmarking points are (the order is hadron on lepton):

* 275 on 10 GeV (what we are currently using) 
* 275 on 18 GeV
* 100 on 10 GeV
* 100 on 5 GeV
* 41 on 5 GeV


## Running the "SimpleAnalysis" framework

SimpleAnalysis is a basic C++ framework that can operate on the ROOT files produced by Delphes. To compile:

```
cd SimpleAnalysis/
export DELPHES_PATH=<PATH TO DELPHES INSTALLATION>
make
```

This builds ```SimpleAnalysis.exe```. You can then run this on files produced above using the ```run_study.py``` script:

```
./SimpleAnalysis.exe --input_dir '../CC_DIS_e10_p275_CT18NNLO/0/out.root' --output_file test.root --module_sequence KaonPIDModule,ElectronPIDModule,MuonPIDModule,TaggingModule --nevents=100000
```

This runs on a single Delphes ROOT file and produces a new output file, test.root, containing the results of running four modules in sequence:

* KaonPIDModule: uses a very basic PID model (efficiency and K/pi separation) to "identify" kaons from tracks
* ElectronPIDModule: same as kaon PID, but for electrons
* MuonPIDModule: same as kaon PID, but for muons
* TaggingModule: Uses the lists of kaon, muon, and electron candidates provided by the above modules, as well as an implementation of signed-high-impact-parameter track finding, to tag jets.
* SubJetModule: work in progress, is capable of reporting important kinematics for analysis on the subjets for different algorithm/method types and the jets that they are constructed into.

The output file contains one entry per jet studied with a few basic jet variables. These can be processed using the scripts in ```SimpleAnalysis/scripts``` to make some plots.


