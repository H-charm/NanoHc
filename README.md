Tree producer for H(ZZ->4l)+c analysis starting from official NANOAOD  

Setup  
-----  
In your **AFS area**, install the required CMSSW:  
```  
cmsrel CMSSW_13_3_3
cd CMSSW_13_3_3/src
cmsenv
```
Then, install the required packages:
```
git cms-addpkg PhysicsTools/NanoAODTools
git fetch https://github.com/namapane/cmssw.git NAT-dev2:namapane_NAT-dev2
git cherry-pick aa9ecbd04d6 98f8692142f

git cms-addpkg PhysicsTools/NanoAOD
git cms-cherry-pick-pr 43536 CMSSW_13_0_X

git fetch https://github.com/namapane/cmssw.git NAT-dev:namapane_NAT-dev
git cherry-pick 3e73ca4c2f8

git clone https://github.com/cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
cd PhysicsTools/NATModules; git checkout -b from-37a092e 37a092e
```

Finally, install this repository:
```
git clone git@github.com:H-charm/NanoHc.git PhysicsTools/NanoHc
```
and checkout to branch:
```
git checkout [branch]
```
Where 
- `Run3`:
- `Run2`:
- `Run3_CR`:
- `Zpeak`: 
- `Zpeak_CR`:
