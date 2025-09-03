# Tree Producer for H(ZZ→4ℓ)+c Analysis

This repository provides tools to produce analysis trees for the **H(ZZ→4ℓ)+c** study starting from official **NanoAOD** samples.

## Setup

Follow the steps below to set up the environment in your **AFS area**.

### 1. Set up CMSSW

```bash
cmsrel CMSSW_13_3_3
cd CMSSW_13_3_3/src
cmsenv
```

### 2. Install required packages

```bash
# NanoAODTools and patches
git cms-addpkg PhysicsTools/NanoAODTools
git fetch https://github.com/namapane/cmssw.git NAT-dev2:namapane_NAT-dev2
git cherry-pick aa9ecbd04d6 98f8692142f

# NanoAOD and PR backport
git cms-addpkg PhysicsTools/NanoAOD
git cms-cherry-pick-pr 43536 CMSSW_13_0_X

# Additional fixes
git fetch https://github.com/namapane/cmssw.git NAT-dev:namapane_NAT-dev
git cherry-pick 3e73ca4c2f8

# NATModules
git clone https://github.com/cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
cd PhysicsTools/NATModules
git checkout -b from-37a092e 37a092e
```
### 3. Install this repository

```bash
cd CMSSW_13_3_3/src
git clone git@github.com:H-charm/NanoHc.git PhysicsTools/NanoHc
```
### 4. Select branch

Choose the branch corresponding to your use case and scram:
```bash
cd PhysicsTools/NanoHc
git checkout [branch]
cd CMSSW_13_3_3/src
scram b -j8
```

* **Run3** – Standard Run 3 analysis
* **Run2** – Standard Run 2 analysis
* **Run3\_CR** – Control region (Run 3)
* **Zpeak** – Validation of the Z peak
* **Zpeak\_CR** – Validation of the Z peak using the Z+l CR


Navigate to the corresponding branch for instructions on how to run and details regarding the producers.
