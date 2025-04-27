import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('lumiJson', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'apply this lumi json file'
                 )
options.register('applyFilt', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'Apply filters'
                 )
options.parseArguments()

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Analysis",eras.Run2_2018_highBetaStar)

# import of standard configurations
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('ERROR')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )

#configure for data/MC based on lumi json input (if passed in command line)
from Configuration.AlCa.GlobalTag import GlobalTag
if options.lumiJson:
  print('Running on data file:')
  process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v11"
  print('Lumi sections will be selected with',options.lumiJson)
  from FWCore.PythonUtilities.LumiList import LumiList
  myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
  process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
  process.source.lumisToProcess.extend(myLumis)
else: process.GlobalTag.globaltag = "141X_mcRun3_2024_realistic_HI_v13"

process.TFileService = cms.Service("TFileService",
                    fileName = cms.string('output.root')
					)
					
process.analysis = cms.EDAnalyzer('Analyzer',
       tracks = cms.untracked.InputTag('tracks'),
       pfCands = cms.untracked.InputTag('packedPFCandidates'),
#       DeDxData = cms.untracked.InputTag('energyLossProducer','energyLossAllHits'),
#       vertices = cms.InputTag('offlinePrimaryVertices'),	 
#       vertices = cms.InputTag('offlinePrimaryVerticesWithBS'),	 
#       rpTrackTag = cms.InputTag('ctppsLocalTrackLiteProducer'),
#       genParticles = cms.InputTag('genParticles'),
       applyFilt = cms.bool(options.applyFilt),
)

process.anal = cms.Path(process.analysis)


###############################################################################
# Schedule
process.schedule = cms.Schedule(process.anal)

