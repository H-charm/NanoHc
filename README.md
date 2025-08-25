Tree producer for H+c analysis starting from official NANOAOD  

Setup  
-----  
In your **AFS area**:  
```  
cmsrel CMSSW_13_3_3
cd CMSSW_13_3_3/src
cmsenv
```

```
git clone https://github.com/H-charm/NanoHc.git PhysicsTools/NanoHc

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