To set up this code

cmsrel CMSSW_9_4_9

cd CMSSW_9_4_9/src

mkdir HmmAna

cd HmmAna

git clone git@github.com:irenedutta23/HmmAna master

make


To run the code

./analyzeHmm runList.txt out.root mc F 2016

OR 

./analyzeHmm runList.txt out.root data T 2016


============================================

CONDOR scripts (not yet adapted to work for three years)

===========================================

Before running condor, make sure to change the following items:
1. give appropriate addresses in the run_myprog.sh, proto_condor_submit and makecondorsubmit.py
2. When running over data, use the option "T" for IsData in makeCondorSubmit.py. While running over mc, use "F"
3. to initiate condor jobs, do python makeCondoSubmit.py
4. For every data[0] name in the makeCondorSubmit, make sure to have the appropriate runlist file in the condor directory and the latest version of analyzeHmm (executable)
5. Also create the following directories before submitting jobs

  a. condor/condor_output/condor_logs
  
  b. condor/condor_submit
 
6. Make sure to have the RoccoR2017.txt in condor/condor_output/condor_logs directory.
7. You can resubmit failed jobs by doing condor_submit condor_submit/submit_condor_*job* (appropriate job name)
