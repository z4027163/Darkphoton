import FWCore.ParameterSet.Config as cms

# Set parameters externally 
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterMuons', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to contain at least two muons'
)

params.register(
    'reducedInfo', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to store just the reduced information'
)

params.register(
    'trigProcess', 
    'HLT', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData', 
    '101X_dataRun2_Prompt_v11', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagMC', 
    '101X_dataRun2_Prompt_v11', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'xsec', 
    0.001, 
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000000) )

# Input EDM files
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring([
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_700.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_701.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_702.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_703.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_704.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_705.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_706.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_707.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_708.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_709.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_70.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_710.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_711.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_712.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_713.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_714.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_715.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_716.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_717.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_718.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_719.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_71.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_720.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_721.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_722.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_723.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_724.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_725.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_726.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_727.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_728.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_729.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_72.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_730.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_731.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_732.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_733.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_734.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_735.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_736.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_737.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_738.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_739.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_73.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_740.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_741.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_742.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_743.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_744.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_745.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_746.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_747.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_748.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_749.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_74.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_750.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_751.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_752.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_753.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_754.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_755.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_756.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_757.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_758.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_759.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_75.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_760.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_761.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_762.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_763.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_764.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_765.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_766.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_767.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_768.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_769.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_76.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_770.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_771.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_772.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_773.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_774.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_775.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_776.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_777.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_778.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_779.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_77.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_780.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_781.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_782.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_783.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_784.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_785.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_786.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_787.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_788.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_789.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_78.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_790.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_791.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_792.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_793.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_794.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_795.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_796.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_797.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_798.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_799.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_79.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_7.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_800.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_801.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_802.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_803.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_804.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_805.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_806.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_807.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_808.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_809.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_80.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_810.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_811.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_812.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_813.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_814.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_815.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_816.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_817.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_818.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_819.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_81.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_820.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_821.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_822.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_823.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_824.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_825.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_826.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_827.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_828.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_829.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_82.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_830.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_831.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_832.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_833.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_834.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_835.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_836.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_837.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_838.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_839.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_83.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_840.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_841.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_842.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_843.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_844.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_845.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_846.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_847.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_848.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_849.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_84.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_850.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_851.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_852.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_853.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_854.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_855.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_856.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_857.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_858.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_859.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_85.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_860.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_861.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_862.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_863.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_864.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_865.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_866.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_867.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_868.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_869.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_86.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_870.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_871.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_872.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_873.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_874.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_875.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_876.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_877.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_878.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_879.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_87.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_880.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_881.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_882.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_883.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_884.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_885.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_886.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_887.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_888.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_889.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_88.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_890.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_891.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_892.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_893.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_894.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_895.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_896.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_897.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_898.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_899.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_89.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_8.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_900.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_901.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_902.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_903.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_904.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_905.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_906.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_907.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_908.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_909.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_90.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_910.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_911.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_912.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_913.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_914.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_915.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_916.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_918.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_919.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_91.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_920.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_921.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_922.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_923.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_924.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_925.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_926.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_927.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_928.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_929.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_92.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_930.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_931.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_932.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_933.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_934.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_935.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_936.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_937.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_938.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_939.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_93.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_940.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_941.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_942.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_943.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_944.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_945.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_946.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_947.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_948.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_949.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_94.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_950.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_951.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_952.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_953.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_954.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_955.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_956.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_957.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_958.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_959.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_95.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_960.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_961.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_962.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_963.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_964.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_965.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_966.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_967.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_968.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_969.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_96.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_970.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_971.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_972.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_973.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_974.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_975.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_976.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_977.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_978.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_979.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_97.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_980.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_981.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_982.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_983.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_984.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_985.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_986.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_987.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_988.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_989.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_98.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_990.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_991.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_992.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_993.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_994.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_995.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_996.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_997.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_998.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_999.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_99.root',
'file:/eos/user/w/wangz/aod_sample/lowDY/aod_9.root'
	])
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC : 
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("output.root")
)

# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheInfo = cms.InputTag("externalLHEProducer"),
    genInfo = cms.InputTag("generator"),
    useLHEWeights = cms.bool(params.useWeights)
)

from DarkPhotonAnalysisV2.DimuonAnalysis.TriggerPaths2018_cfi import getL1Conf

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingTreeMaker2017',
	isMC             = cms.bool(params.isMC),
	xsec             = cms.double(params.xsec),
	useLHEWeights    = cms.bool(params.useWeights),
	applyHLTFilter   = cms.bool(params.filterTrigger),
	require2Muons    = cms.bool(params.filterMuons),
	storeReducedInfo = cms.bool(params.reducedInfo),
    	triggerresults   = cms.InputTag("TriggerResults", "", params.trigProcess),
	doL1 = cms.bool(True),
        triggerConfiguration = cms.PSet(
    		hltResults            = cms.InputTag('TriggerResults','','HLT'),
    		l1tResults            = cms.InputTag(''),
    		daqPartitions         = cms.uint32(1),
    		l1tIgnoreMaskAndPrescale = cms.bool(False),
    		#l1tIgnoreMask         = cms.bool(False),
   		#l1techIgnorePrescales = cms.bool(False),
    		throw                 = cms.bool(False)
  	),
	ReadPrescalesFromFile = cms.bool( False ),
        AlgInputTag = cms.InputTag("gtStage2Digis"),
        l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1Seeds           = cms.vstring(getL1Conf()),
#	vertices          = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
	vertices          = cms.InputTag("hltScoutingMuonPackerCalo","displacedVtx"),
#	beamspot          = cms.InputTag("")
	muons             = cms.InputTag("hltScoutingMuonPackerCalo"),
	#pfcands          = cms.InputTag("hltScoutingPFPacker"),
	rho               = cms.InputTag("hltScoutingCaloPacker", "rho"),
    	pileupinfo        = cms.InputTag("addPileupInfo"),
    	geneventinfo      = cms.InputTag("generator"),
    	gens              = cms.InputTag("genParticles", "", "HLT")
)

# Analysis path
process.p = cms.Path(  process.gtStage2Digis + process.mmtree)
#if params.isMC : 
#    process.p = cms.Path(process.gentree + process.mmtree)
#else : 
#    process.p = cms.Path(  process.gtStage2Digis + process.mmtree)

