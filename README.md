# pOAnalysis

## Installation instructions

To install the package, execute the following in your work area.
If compilation fails for some reason, repeat the scram b...

```
#Setup CMSSW (done once, on el9)
SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_15_0_0_pre2
cd CMSSW_15_0_0_pre2/src
cmsenv

#Clone this package
cd $CMSSW_BASE/src
git clone https://github.com/michael-pitt/pOAnalysis.git
cd pOAnalysis/Analyzer/

scram b -j 8
```

## Running the analysis

Analysis are performed in two steps
   1. process data/MC files from MiniAOD format
   2. Ntuplize output files - at this stage trackID is computed together with some high-level variables.
   
### Processing the data

To process a single data file (stored in `/eos/cms/store/group/phys_heavyions/mpitt/pO_miniAOD`) run the following command:
```
cmsRun $CMSSW_BASE/src/pOAnalysis/Analyzer/python/ConfFile_cfg.py inputFiles=file:$FILENAME applyFilt=True
```
You can set `applyFilt=False` of no event pre-selection should be applied (usually used to debug). The _FILENAME_ argument can be any file from the EOS folder. For example:

```
FILENAME=/eos/cms/store/group/phys_heavyions/mpitt/pO_miniAOD/MinBias_pythia_Op/miniAOD_0487.root
```

### Running ntuplizer (IN DEVELOPMENT)

At this step, track identification is executed. The trackID is based on dEdX of the reconstructed track [[ref](https://indico.cern.ch/event/1154003/#4-ntuple-production-for-glueba)].

To execute this step, run the following command:
```
analysisWrapper --in output.root --out tracks.root --method TrackAnalysis
```
Where `output.root` is the output file from the previous step. 

## Producing MC samples

### MinimumBias

The Pythia fragment of MinimumBias events can be found in `/afs/cern.ch/user/m/mpitt/public/oxygen_run/fragments/pythia_Op_fragment.py`. 

To generate GEN samples, the following code can be executed using the `cmsDriver.py` command:

```
#copy the fragment to Configuration/GenProduction/python/MinBias_fragment.py
#produce SIM (50 events, using a random seed = 12345)
cmsDriver.py Configuration/GenProduction/python/MinBias_fragment.py --python_filename step1_cfg.py \
--eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM \
--fileout file:stepSIM.root --141X_mcRun3_2024_realistic_HI_v13 --beamspot DBrealistic \
--customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32(12345)" \
--step GEN,SIM --geometry DB:Extended --era Run3_2025_OXY --runUnscheduled --mc -n 50
``` 

### SIM -> AOD step

The following code can be executed using the `cmsDriver.py` command to generate AOD output from the `stepSIM.root` file:

```
cmsDriver.py --python_filename step2_cfg.py --fileout file:stepDR.root --eventcontent RAWSIM \
 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI-RAW \
 --conditions 141X_mcRun3_2024_realistic_HI_v13 --step DIGI,L1,DIGI2RAW,HLT:HIon \
 --customise_commands "process.HcalTPGCoderULUT.FG_HF_thresholds = [16, 19]\n process.RAWSIMoutput.outputCommands.extend(['keep *_mix_MergedTrackTruth_*', 'keep *Link*_simSiPixelDigis__*', 'keep *Link*_simSiStripDigis__*'])" \
 --geometry DB:Extended --filein file:stepSIM.root --era Run3_2025_OXY --mc -n -1
``` 

```
cmsDriver.py --python_filename step3_cfg.py --eventcontent AODSIM --fileout file:stepAOD.root \
 --datatier AODSIM --conditions 141X_mcRun3_2024_realistic_HI_v13 --step RAW2DIGI,L1Reco,RECO,RECOSIM \
 --geometry DB:Extended --filein file:stepDR.root --era Run3_2025_OXY --mc -n -1
```

### miniAOD step

The following code can be executed using the `cmsDriver.py` command to generate MINIAOD output from the `stepAOD.root` file:

```
cmsDriver.py  --python_filename step4_cfg.py --fileout file:miniAOD.root --eventcontent MINIAODSIM \
 --datatier MINIAODSIM --conditions 141X_mcRun3_2024_realistic_HI_v13 --step PAT \
 --geometry DB:Extended --filein file:stepAOD.root --era Run3_2025_OXY --mc -n -1
```

