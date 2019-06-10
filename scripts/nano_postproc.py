#!/usr/bin/env python

#NanoAOD postprocessing for the Caltech Hmumu analysis
#Run as ./nano_postproc.py --help

import os
import argparse

#load the NanoAOD tools
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016All, jetmetUncertainties2017All, jetmetUncertainties2018All

def parse_args():
    parser = argparse.ArgumentParser(description='Do NanoAOD post-processing by computing PU weights etc on the given NanoAOD MC files')
    parser.add_argument('--data_period', type=str, choices=["2016", "2017", "2018"],
        help="The data-taking period that corresponds to the MC file being processed", default="2017"
    )
    parser.add_argument('--outdir', type=str, required=True,
        help="The output directory where to write the *_Friend.root files, must be writable (not /mnt/hadoop!)"
    )
    parser.add_argument('input_files', nargs=argparse.REMAINDER, help="List of all input files, must be readable using ROOT")
   
    args = parser.parse_args()

    #example input files just to test
    if len(args.input_files) == 0:
    	args.input_files = [
            "/mnt/hadoop/store/mc/RunIIFall17NanoAOD/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/F2ABC928-DAA5-E811-8A4A-842B2B689031.root",
            "/mnt/hadoop/store/mc/RunIIFall17NanoAOD/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/F80A3F6D-43A5-E811-BE12-7CD30AD08ED2.root",
            "/mnt/hadoop/store/mc/RunIIFall17NanoAOD/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FE35702D-DAA5-E811-9355-008CFA1113D8.root",
    	]

    return args

#generated using HmmAna/scripts/makeDataPileupHist.sh on lxplus
pileup_files_data = {
    "2016": "RunII_2016_data.root", 
    "2017": "RunII_2017_data.root", 
    "2018": "RunII_2018_data.root", 
}


jetmet_uncertainties = {
    "2016": jetmetUncertainties2016All, 
    "2017": jetmetUncertainties2017All, 
    "2018": jetmetUncertainties2018All
}

def main(args):

    #where to find our user-generated data PU files
    datadir = os.path.join(os.environ["CMSSW_BASE"], "src/MyAnalysis/HmmAna/data")

    #configure the puWeightProducer to extract the MC distribution directly from the current NanoAOD file ("auto")
    puAutoWeight = lambda : puWeightProducer("auto", os.path.join(datadir, "pileup", pileup_files_data[args.data_period]), "pu_mc", "pileup", verbose=True)
   
    #book all the different things we want to compute in the postprocessor
    modules = [puAutoWeight(), jetmet_uncertainties[args.data_period]()]
    
    #friend: if True, create a new tree that contains *only* the output weights, without copying all the inputs. Much faster.
    p=PostProcessor(args.outdir, args.input_files, friend=True, modules=modules)

    #actually run the code, creating one output file per each input file
    p.run()

if __name__ == "__main__":
    args = parse_args()
    main(args)
