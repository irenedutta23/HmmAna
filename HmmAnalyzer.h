 
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 27 15:49:34 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/GluGluHToMuMu_M-125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/AA4BF847-3785-E811-B130-001E67DFFF5F.root
//////////////////////////////////////////////////////////

#ifndef HmmAnalyzer_h
#define HmmAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <map>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "TRandom.h"
#include "MainEvent.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
#include "RoccoR.h"
#include "LeptonEfficiencyCorrector.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif
// Header file for the classes stored in the TTree if any.

class HmmAnalyzer : public MainEvent {
 public :
   HmmAnalyzer(const TString &inputFileList="foo.txt", const char *outFileName="histo.root", TString dataset="data",const char *isData="F", TString year_num="2017");
   virtual ~HmmAnalyzer();
   void Analyze(bool isData, int option, string outputFileName, string label);
   
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  float getPileupWeight(int);
  float getPileupWeightUp(int);
  float getPileupWeightDown(int);
  void CorrectPtRoch( const RoccoR& _calib, const bool _doSys, const TLorentzVector _mu_vec,
                    float& _pt, float& _ptErr, float& _pt_sys_up, float& _pt_sys_down,
                    const int _charge, const int _trk_layers, const float _GEN_pt, const bool _isData );
  void     EventLoop(const char *, const char *);
  //declare any specific function required
  
  void clearTreeVectors();
  void BookTreeBranches();
  bool DataIs;
  TString year;
  std::string yearst;
  std::map<std::string,float> muon_pt_cut;
  std::map<std::string,float> btag_cut;

  TH1D *h_sumOfgw = new TH1D("h_sumOfgenWeight","h_sumOfgenWeight",1,0,1);
  TH1D *h_sumOfgpw = new TH1D("h_sumOfgenpuWeight","h_sumOfgenpuWeight",1,0,1);
  RoccoR _Roch_calib;
  
  TFile *pileupWeightFile;
  TH1F *pileupWeightHist, *pileupWeightSysUpHist, *pileupWeightSysDownHist;

  std::vector<std::string> muon_effSF_TRIG_files, muon_effSF_ID_files, muon_effSF_ISO_files;
  std::vector<std::string> histo_names_TRIG, histo_names_ID, histo_names_ID_stat, histo_names_ID_syst, histo_names_ISO, histo_names_ISO_stat, histo_names_ISO_syst;
  LeptonEfficiencyCorrector Mu_eff_SF_TRIG;
  LeptonEfficiencyCorrector Mu_eff_SF_ID;
  LeptonEfficiencyCorrector Mu_eff_SF_ID_stat;
  LeptonEfficiencyCorrector Mu_eff_SF_ID_syst;
  LeptonEfficiencyCorrector Mu_eff_SF_ISO;
  LeptonEfficiencyCorrector Mu_eff_SF_ISO_stat;
  LeptonEfficiencyCorrector Mu_eff_SF_ISO_syst;

  TFile *oFile;
  //TFile *ohistFile;
  TTree* tree;
  uint          t_run;
  uint          t_luminosityBlock;
  ulong       t_event;
  float       t_genWeight;
  float       t_puWeight;
  float       t_puWeightUp;
  float       t_puWeightDown;
  float       t_PrefireWeight;
  float       t_PrefireWeight_Up;
  float       t_PrefireWeight_Down;
  int         t_mu1;
  int         t_mu2;
  int         t_index_trigm_mu;
  std::vector<int>           *t_El_genPartIdx;
  std::vector<UChar_t>       *t_El_genPartFlav;
  std::vector <int>          *t_El_charge;
  std::vector<float>         *t_El_pt;
  std::vector<float>         *t_El_phi;
  std::vector<float>         *t_El_eta;   
  std::vector<float>         *t_El_mass;
  std::vector <int>          *t_El_cutBased;   
  std::vector <int>          *t_El_tightCharge;   
  std::vector<bool>          *t_El_cutBased_HEEP;   
  std::vector<bool>          *t_El_isPFcand;   
  std::vector<float>         *t_El_pfRelIso03_all;   
  std::vector<float>         *t_El_pfRelIso03_chg;      
  std::vector<float>         *t_El_miniPFRelIso_all;
  std::vector<float>         *t_El_miniPFRelIso_chg;
  std::vector<float>         *t_El_dxy;   
  std::vector<float>         *t_El_dxyErr;   
  std::vector<float>         *t_El_dz;   
  std::vector<float>         *t_El_dzErr; 
  std::vector<float>         *t_El_sip3d; 
  std::vector<float>         *t_Electron_mvaFall17Iso; 
  std::vector<bool>         *t_Electron_mvaFall17Iso_WP80;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17Iso_WP90;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17Iso_WPL;   //[nElectron]
  std::vector<float>        *t_Electron_mvaFall17noIso;
  std::vector<bool>         *t_Electron_mvaFall17noIso_WP80;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17noIso_WP90;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17noIso_WPL;   //[nElectron]
  
  std::vector<int>           *t_Mu_genPartIdx;
  std::vector<UChar_t>       *t_Mu_genPartFlav;
  std::vector<int>           *t_Mu_charge;  
  std::vector<float>         *t_Mu_EffSF_TRIG;
  std::vector<float>         *t_Mu_EffSFErr_TRIG;
  std::vector<float>         *t_Mu_EffSF_ID;
  std::vector<float>         *t_Mu_EffSF_ID_stat;
  std::vector<float>         *t_Mu_EffSF_ID_syst;
  std::vector<float>         *t_Mu_EffSF_ISO;
  std::vector<float>         *t_Mu_EffSF_ISO_stat;
  std::vector<float>         *t_Mu_EffSF_ISO_syst;
  std::vector<float>         *t_Mu_pt;   
  std::vector<float>         *t_Mu_ptErr;   
  std::vector<float>         *t_Mu_phi;   
  std::vector<float>         *t_Mu_eta;   
  std::vector<float>         *t_Mu_mass;  
  std::vector<float>         *t_Mu_dxy;   
  std::vector<float>         *t_Mu_dxyErr;   
  std::vector<float>         *t_Mu_dz;   
  std::vector<float>         *t_Mu_dzErr;
  std::vector<float>         *t_Mu_sip3d; 
  std::vector<float>         *t_Mu_pfRelIso03_all;   
  std::vector<float>         *t_Mu_pfRelIso03_chg;   
  std::vector<float>         *t_Mu_pfRelIso04_all;   
  std::vector<float>         *t_Mu_miniPFRelIso_all;
  std::vector<float>         *t_Mu_miniPFRelIso_chg;
  std::vector<int>           *t_Mu_tightCharge;   
  std::vector<bool>          *t_Mu_isPFcand;  
  std::vector<bool>          *t_Mu_istracker;
  std::vector<bool>          *t_Mu_isglobal; 
  std::vector<bool>          *t_Mu_mediumId;   
  std::vector<bool>          *t_Mu_softId;   
  std::vector<bool>          *t_Mu_tightId;    
  std::vector<int>           *t_Mu_nStations;   
  std::vector<int>           *t_Mu_nTrackerLayers;   

  float t_diMuon_pt;
  float t_diMuon_eta;
  float t_diMuon_phi;
  float t_diMuon_mass;
  
  float         t_MET_phi;
  float         t_MET_pt;
  float         t_MET_sumEt;

  std::vector<float>         *t_FatJet_area;  
  std::vector<float>         *t_FatJet_btagCMVA;  
  std::vector<float>         *t_FatJet_btagCSVV2;  
  std::vector<float>         *t_FatJet_btagDeepB;  
  std::vector<float>         *t_FatJet_eta;  
  std::vector<float>         *t_FatJet_mass;  
  std::vector<float>         *t_FatJet_msoftdrop;  
  std::vector<float>         *t_FatJet_n2b1;  
  std::vector<float>         *t_FatJet_n3b1;  
  std::vector<float>         *t_FatJet_phi;  
  std::vector<float>         *t_FatJet_pt;  
  std::vector<float>         *t_FatJet_tau1;  
  std::vector<float>         *t_FatJet_tau2;  
  std::vector<float>         *t_FatJet_tau3;  
  std::vector<float>         *t_FatJet_tau4;  
  std::vector<int>           *t_FatJet_jetId;  
  std::vector<int>           *t_FatJet_subJetIdx1;  
  std::vector<int>           *t_FatJet_subJetIdx2;  

  std::vector<float>         *t_SubJet_btagCMVA;   
  std::vector<float>         *t_SubJet_btagCSVV2;   
  std::vector<float>         *t_SubJet_btagDeepB;   
  std::vector<float>         *t_SubJet_eta;   
  std::vector<float>         *t_SubJet_mass;   
  std::vector<float>         *t_SubJet_n2b1;   
  std::vector<float>         *t_SubJet_n3b1;   
  std::vector<float>         *t_SubJet_phi;   
  std::vector<float>         *t_SubJet_pt;   
  std::vector<float>         *t_SubJet_tau1;   
  std::vector<float>         *t_SubJet_tau2;   
  std::vector<float>         *t_SubJet_tau3;   
  std::vector<float>         *t_SubJet_tau4;   
  int t_nJet;
  std::vector<float>         *t_Jet_area;   
  //std::vector<float>         *t_Jet_btagCMVA;   
  //std::vector<float>         *t_Jet_btagCSVV2;   
  std::vector<float>         *t_Jet_btagDeepB;   
  //std::vector<float>         *t_Jet_btagDeepC;   
  //std::vector<float>         *t_Jet_btagDeepFlavB;   
  std::vector<float>         *t_Jet_chEmEF;   
  std::vector<float>         *t_Jet_chHEF;   
  std::vector<float>         *t_Jet_eta;   
  std::vector<float>         *t_Jet_mass;   
  std::vector<float>         *t_Jet_neEmEF;   
  std::vector<float>         *t_Jet_neHEF;   
  std::vector<float>         *t_Jet_phi;   
  std::vector<float>         *t_Jet_pt;   
  std::vector<float>         *t_Jet_qgl;   
  std::vector<int>           *t_Jet_jetId;   
  std::vector<int>           *t_Jet_nConstituents;   
  std::vector<int>           *t_Jet_nElectrons;   
  std::vector<int>           *t_Jet_nMuons;   
  std::vector<int>           *t_Jet_puId;   

  int   t_SoftActivityJetNjets5;
  float t_diJet_pt;
  float t_diJet_eta;
  float t_diJet_phi;
  float t_diJet_mass;
  float t_diJet_mass_mo;

  int t_nbJet;
  std::vector<float>         *t_bJet_area;   
  std::vector<float>         *t_bJet_btagCMVA;   
  std::vector<float>         *t_bJet_btagCSVV2;   
  std::vector<float>         *t_bJet_btagDeepB;   
  std::vector<float>         *t_bJet_btagDeepC;   
  std::vector<float>         *t_bJet_btagDeepFlavB;   
  std::vector<float>         *t_bJet_chEmEF;   
  std::vector<float>         *t_bJet_chHEF;   
  std::vector<float>         *t_bJet_eta;   
  std::vector<float>         *t_bJet_mass;   
  std::vector<float>         *t_bJet_neEmEF;   
  std::vector<float>         *t_bJet_neHEF;   
  std::vector<float>         *t_bJet_phi;   
  std::vector<float>         *t_bJet_pt;   
  std::vector<float>         *t_bJet_qgl;   
  std::vector<int>           *t_bJet_jetId;   
  std::vector<int>           *t_bJet_nConstituents;   
  std::vector<int>           *t_bJet_nElectrons;   
  std::vector<int>           *t_bJet_nMuons;   
  std::vector<int>           *t_bJet_puId;   
  std::vector<double>         *t_bJet_SF;
  std::vector<double>         *t_bJet_SFup;
  std::vector<double>         *t_bJet_SFdown;

  float      t_PV_ndof;
  float      t_PV_x;
  float      t_PV_y;
  float      t_PV_z;
  int        t_PV_npvs;
  int        t_PV_npvsGood;

  UInt_t     t_nLHEPdfWeight;
  UInt_t     t_nLHEScaleWeight;
  UInt_t     t_nPSWeight;

  std::vector<float>         *t_LHEPdfWeight;
  std::vector<float>         *t_LHEScaleWeight;
  std::vector<float>         *t_PSWeight;

  std::vector<float>         *t_GenPart_eta;
  std::vector<float>         *t_GenPart_mass;
  std::vector<float>         *t_GenPart_phi;
  std::vector<float>         *t_GenPart_pt;
  std::vector<int>           *t_GenPart_genPartIdxMother;
  std::vector<int>           *t_GenPart_pdgId;
  std::vector<int>           *t_GenPart_status;

  std::vector<float>         *t_GenJet_eta;
  std::vector<float>         *t_GenJet_mass;
  std::vector<float>         *t_GenJet_phi;
  std::vector<float>         *t_GenJet_pt;
};

#endif

#ifdef HmmAnalyzer_cxx
HmmAnalyzer::HmmAnalyzer(const TString &inputFileList, const char *outFileName, TString dataset, const char *isData, TString year_num) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if(*isData!='T') DataIs = false;
  else DataIs = true;
  year = year_num;
  yearst = std::string(year.Data());
  std::string path_RochCor = "data/Rocco/RoccoR"+yearst+".txt";
  std::cout << "Rochester correction files: " << path_RochCor << std::endl;
  _Roch_calib.init(path_RochCor);

  h_sumOfgw->SetBinContent(1,0.0);
  h_sumOfgpw->SetBinContent(1,0.0);

  //muon pT selection
  muon_pt_cut["2016"] = 26.0;
  muon_pt_cut["2017"] = 29.0;
  muon_pt_cut["2018"] = 26.0;   
  //b-tag score selection
  btag_cut["2016"] = 0.6321; 
  btag_cut["2017"] = 0.4941;
  btag_cut["2018"] = 0.4184;

  //muon eff SFs
  muon_effSF_TRIG_files.clear();
  muon_effSF_ID_files.clear();
  muon_effSF_ISO_files.clear();
  histo_names_TRIG.clear();
  histo_names_ID.clear();
  histo_names_ISO.clear();
  if(year=="2016"){
  	std::string Mu_ID_file1 = "data/leptonSF/"+yearst+"/RunBCDEF_SF_ID.root";
 	std::string Mu_Iso_file1 = "data/leptonSF/"+yearst+"/RunBCDEF_SF_ISO.root";
  	std::string Mu_ID_file2 = "data/leptonSF/"+yearst+"/RunGH_SF_ID.root";
  	std::string Mu_Iso_file2 = "data/leptonSF/"+yearst+"/RunGH_SF_ISO.root";
  	std::string Mu_Trg_file1 = "data/leptonSF/"+yearst+"/EfficienciesAndSF_RunBtoF.root";
  	std::string Mu_Trg_file2 = "data/leptonSF/"+yearst+"/EfficienciesAndSF_RunGtoH.root";
  	std::cout << "Muon ID correction files: " << Mu_ID_file1 << " "<<Mu_ID_file2<<std::endl;
  	std::cout << "Muon Isolation correction files: " << Mu_Iso_file1 <<" "<<Mu_Iso_file2 << std::endl;
  	std::cout << "Muon triger SF files: " << Mu_Trg_file1 <<" "<<Mu_Trg_file2<< std::endl;

  	muon_effSF_TRIG_files.push_back(Mu_Trg_file1);
  	muon_effSF_TRIG_files.push_back(Mu_Trg_file2);
  	muon_effSF_ID_files.push_back(Mu_ID_file1);
  	muon_effSF_ID_files.push_back(Mu_ID_file2);
  	muon_effSF_ISO_files.push_back(Mu_Iso_file1);
  	muon_effSF_ISO_files.push_back(Mu_Iso_file2);
  	std::string Mu_Trg_name = "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio";
  	std::string Mu_ID_name = "NUM_MediumID_DEN_genTracks_eta_pt";
  	std::string Mu_Iso_name = "NUM_LooseRelIso_DEN_MediumID_eta_pt";
        histo_names_TRIG.push_back(Mu_Trg_name);
        histo_names_ID.push_back(Mu_ID_name);
        histo_names_ISO.push_back(Mu_Iso_name);
        histo_names_TRIG.push_back(Mu_Trg_name);
        histo_names_ID.push_back(Mu_ID_name);
        histo_names_ISO.push_back(Mu_Iso_name);
  }
  else if(year=="2017"){
        std::string Mu_Trg_file = "data/leptonSF/"+yearst+"/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root";
  	std::string Mu_ID_file = "data/leptonSF/"+yearst+"/RunBCDEF_SF_ID_syst.root";
  	std::string Mu_Iso_file = "data/leptonSF/"+yearst+"/RunBCDEF_SF_ISO_syst.root";
  	muon_effSF_TRIG_files.push_back(Mu_Trg_file);
  	muon_effSF_ID_files.push_back(Mu_ID_file);
  	muon_effSF_ISO_files.push_back(Mu_Iso_file);
  	std::string Mu_Trg_name = "IsoMu27_PtEtaBins/pt_abseta_ratio";
  	std::string Mu_ID_name = "NUM_MediumID_DEN_genTracks_pt_abseta";
  	std::string Mu_Iso_name = "NUM_LooseRelIso_DEN_MediumID_pt_abseta";
        histo_names_TRIG.push_back(Mu_Trg_name);
        histo_names_ID.push_back(Mu_ID_name);
        histo_names_ISO.push_back(Mu_Iso_name);
  }
  else{
        std::string Mu_Trg_file = "data/leptonSF/"+yearst+"/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root";
  	std::string Mu_ID_file = "data/leptonSF/"+yearst+"/RunABCD_SF_ID.root";
  	std::string Mu_Iso_file = "data/leptonSF/"+yearst+"/RunABCD_SF_ISO.root";
  	muon_effSF_TRIG_files.push_back(Mu_Trg_file);
  	muon_effSF_ID_files.push_back(Mu_ID_file);
  	muon_effSF_ISO_files.push_back(Mu_Iso_file);
        std::string Mu_Trg_name = "IsoMu24_PtEtaBins/pt_abseta_ratio";
  	std::string Mu_ID_name = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
  	std::string Mu_Iso_name = "NUM_LooseRelIso_DEN_MediumID_pt_abseta";
        histo_names_TRIG.push_back(Mu_Trg_name);
        histo_names_ID.push_back(Mu_ID_name);
        histo_names_ISO.push_back(Mu_Iso_name);
  }

  Mu_eff_SF_TRIG.init(muon_effSF_TRIG_files,histo_names_TRIG);
  Mu_eff_SF_ID.init(muon_effSF_ID_files,histo_names_ID);
  Mu_eff_SF_ISO.init(muon_effSF_ISO_files,histo_names_ISO);

  std::string Mu_ID_name_stat = "NUM_MediumID_DEN_genTracks_eta_pt_stat";
  std::string Mu_ID_name_syst = "NUM_MediumID_DEN_genTracks_eta_pt_syst";
  std::string Mu_Iso_name_stat = "NUM_LooseRelIso_DEN_MediumID_eta_pt_stat";
  std::string Mu_Iso_name_syst = "NUM_LooseRelIso_DEN_MediumID_eta_pt_syst";
  if(year=="2017"){
      Mu_ID_name_stat = "NUM_MediumID_DEN_genTracks_pt_abseta_stat";
      Mu_ID_name_syst = "NUM_MediumID_DEN_genTracks_pt_abseta_syst";
      Mu_Iso_name_stat = "NUM_LooseRelIso_DEN_MediumID_pt_abseta_stat";
      Mu_Iso_name_syst = "NUM_LooseRelIso_DEN_MediumID_pt_abseta_syst";
  }
  else if(year=="2018"){
      Mu_ID_name_stat = "NUM_MediumID_DEN_TrackerMuons_pt_abseta_stat";
      Mu_ID_name_syst = "NUM_MediumID_DEN_TrackerMuons_pt_abseta_syst";
      Mu_Iso_name_stat = "NUM_LooseRelIso_DEN_MediumID_pt_abseta_stat";
      Mu_Iso_name_syst = "NUM_LooseRelIso_DEN_MediumID_pt_abseta_syst";
  }
  histo_names_ID_stat.push_back(Mu_ID_name_stat);
  histo_names_ID_syst.push_back(Mu_ID_name_syst);
  histo_names_ISO_stat.push_back(Mu_Iso_name_stat);
  histo_names_ISO_syst.push_back(Mu_Iso_name_syst);
  if(year=="2016"){
     histo_names_ID_stat.push_back(Mu_ID_name_stat);
     histo_names_ID_syst.push_back(Mu_ID_name_syst);
     histo_names_ISO_stat.push_back(Mu_Iso_name_stat);
     histo_names_ISO_syst.push_back(Mu_Iso_name_syst);
  }
  Mu_eff_SF_ID_stat.init(muon_effSF_ID_files,histo_names_ID_stat);
  Mu_eff_SF_ID_syst.init(muon_effSF_ID_files,histo_names_ID_syst);
  Mu_eff_SF_ISO_stat.init(muon_effSF_ISO_files,histo_names_ISO_stat);
  Mu_eff_SF_ISO_syst.init(muon_effSF_ISO_files,histo_names_ISO_syst);

  TChain *tree = new TChain("Events");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  MainEvent::Init(tree);
  oFile = new TFile(outFileName, "recreate");
  TString histname(outFileName);
  //ohistFile = new TFile("hist_"+histname, "recreate");
  BookTreeBranches();
}

float HmmAnalyzer::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

float HmmAnalyzer::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

float HmmAnalyzer::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

void HmmAnalyzer::CorrectPtRoch( const RoccoR& _calib, const bool _doSys, const TLorentzVector _mu_vec,
                    float& _pt, float& _ptErr, float& _pt_sys_up, float& _pt_sys_down,
                    const int _charge, const int _trk_layers, const float _GEN_pt, const bool _isData ) {

  _pt = _mu_vec.Pt();
  _pt_sys_up = -999;
  _pt_sys_down = -999;
  float q_term     = 1.0;
  float q_term_sys = -99;

  float fRand_1 = gRandom->Rndm();
  //float fRand_2 = gRandom->Rndm();

  if (_isData){
              q_term = _calib.kScaleDT( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), 0, 0 );
  }
  else if (_GEN_pt > 0) { 
             q_term = _calib.kSpreadMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), _GEN_pt, 0, 0 );
             //q_term = _calib.kScaleFromGenMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),_trk_layers, _GEN_pt, fRand_1, 0, 0 );
  } else {      
            q_term = _calib.kSmearMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), _trk_layers,fRand_1, 0, 0 );
            //q_term = _calib.kScaleAndSmearMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), _trk_layers, fRand_1, fRand_2, 0, 0 );
  }
  if ( fabs(q_term - 1.0) > 0.4 ) {
    std::cout << "\n*** BIZZARELY HIGH QTERM ***" << std::endl;
    std::cout << "GEN pT = " << _GEN_pt << ", RECO pT = " << _mu_vec.Pt() << ", Q term = " << q_term
              << ", fRand_1 = " << fRand_1 << std::endl;
    std::cout << "Layers = " << _trk_layers << ", charge = " << _charge
              << ", eta = " << _mu_vec.Eta() << ", phi = " << _mu_vec.Phi() << std::endl;
  }


  int nUp   = 0;
  int nDown = 0;
  double sum_sq_up   = 0;
  double sum_sq_down = 0;

  /*
  for (int i = 0; i < 100; i++) {
    if (!_doSys) break;

    if (_isData)          q_term_sys = _calib.kScaleDT( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), 1, i );
    else if (_GEN_pt > 0) q_term_sys = _calib.kScaleFromGenMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
                                                               _trk_layers, _GEN_pt, fRand_1, 1, i );
    else                  q_term_sys = _calib.kScaleAndSmearMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
                                                                _trk_layers, fRand_1, fRand_2, 1, i );
    if ( q_term_sys >= q_term ) {
      nUp   += 1;
      sum_sq_up   += pow( q_term_sys - q_term, 2 );
    } else {
      nDown += 1;
      sum_sq_down += pow( q_term_sys - q_term, 2 );
    }
  }
  */

  _pt          = _mu_vec.Pt() * q_term;
  _ptErr       = _ptErr * q_term; // Account for shift in pT scale
  if (!_isData)  _ptErr *= std::fmax(q_term, 1./q_term); // Account for smearing in MC
  _pt_sys_up   = ( _doSys ? _pt * sqrt(sum_sq_up   / nUp)   : -999 );
  _pt_sys_down = ( _doSys ? _pt * sqrt(sum_sq_down / nDown) : -999 );

}

bool HmmAnalyzer::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;

  while ( getline (infile,buffer) )
  {
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
    if(!DataIs){
      std::string friend_buffer = "";

      std::string delimiter = "/";

      size_t pos = 0;
      std::string buffer_f = "/storage/user/nlu/Hmm/puJEC/"+yearst+"/17July_v1/"; //"/mnt/hadoop/store/user/nlu/Hmm/ntuple/2016/Nano14Dec2018/MC/puJEC/16June_v1/";
      std::string buffer_tmp = buffer;
      int dataname_pos = 7;
      if(buffer.find("storage")!=std::string::npos) dataname_pos = 8;
      int cout=0;
      while ((pos = buffer_tmp.find(delimiter)) != std::string::npos) {
         cout++;
         friend_buffer = buffer_tmp.substr(0, pos);
         std::cout <<"friend_buffer: "<<friend_buffer << std::endl;
         if(cout==dataname_pos){ //8 for storage
                buffer_f +=friend_buffer;
         }
         buffer_tmp.erase(0, pos + delimiter.length());
      }
      std::cout << friend_buffer << std::endl;
      std::cout << "buffer_tmp: "<<buffer_tmp << std::endl;

      buffer_f += buffer_tmp;
      std::cout << "Adding friend tree from " << buffer_f.c_str() << std::endl;
      chain->AddFriend("Friends",buffer_f.c_str());
    }
  }
  infile.close();

  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

 HmmAnalyzer::~HmmAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   h_sumOfgw->Write();
   h_sumOfgpw->Write();
   oFile->Write();
   oFile->Close();
}

Long64_t HmmAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HmmAnalyzer::clearTreeVectors(){
  t_run=0;
  t_luminosityBlock=0;
  t_event=0;
  t_genWeight = -999.;
  t_puWeight = -999.;
  t_puWeightUp = -999.;
  t_puWeightDown = -999.;
  t_PrefireWeight = 1.0;
  t_PrefireWeight_Up = 1.0;
  t_PrefireWeight_Down = 1.0;
  t_mu1=-999; 
  t_mu2=-999;
  t_index_trigm_mu=-999;
  t_nJet=0;
  t_nbJet=0;
  t_SoftActivityJetNjets5=-1000;
  t_El_genPartIdx->clear();
  t_El_genPartFlav->clear();
  t_El_charge->clear();
  t_El_pt->clear();
  t_El_phi->clear();
  t_El_eta->clear();   
  t_El_mass->clear();
  t_El_cutBased->clear();   
  t_El_tightCharge->clear();   
  t_El_cutBased_HEEP->clear();   
  t_El_isPFcand->clear();   
  t_El_pfRelIso03_all->clear();   
  t_El_pfRelIso03_chg->clear();      
  t_El_miniPFRelIso_all->clear();
  t_El_miniPFRelIso_chg->clear();
  t_El_dxy->clear();   
  t_El_dxyErr->clear();   
  t_El_dz->clear();   
  t_El_dzErr->clear();  
  t_El_sip3d->clear();
  t_Electron_mvaFall17Iso->clear(); 
  t_Electron_mvaFall17Iso_WP80->clear();
  t_Electron_mvaFall17Iso_WP90->clear();
  t_Electron_mvaFall17Iso_WPL->clear();
  t_Electron_mvaFall17noIso->clear();
  t_Electron_mvaFall17noIso_WP80->clear();
  t_Electron_mvaFall17noIso_WP90->clear();
  t_Electron_mvaFall17noIso_WPL->clear();

  t_Mu_genPartIdx->clear();
  t_Mu_genPartFlav->clear();  
  t_Mu_charge->clear();
  t_Mu_EffSF_TRIG->clear();
  t_Mu_EffSFErr_TRIG->clear();
  t_Mu_EffSF_ID->clear();
  t_Mu_EffSF_ID_stat->clear();
  t_Mu_EffSF_ID_syst->clear();
  t_Mu_EffSF_ISO->clear();
  t_Mu_EffSF_ISO_stat->clear();
  t_Mu_EffSF_ISO_syst->clear();
  t_Mu_pt->clear();   
  t_Mu_ptErr->clear();   
  t_Mu_phi->clear();   
  t_Mu_eta->clear();   
  t_Mu_mass->clear();  
  t_Mu_dxy->clear();   
  t_Mu_dxyErr->clear();   
  t_Mu_dz->clear();   
  t_Mu_dzErr->clear(); 
  t_Mu_sip3d->clear();
  t_Mu_pfRelIso03_all->clear();   
  t_Mu_pfRelIso03_chg->clear();   
  t_Mu_pfRelIso04_all->clear();   
  t_Mu_miniPFRelIso_all->clear();
  t_Mu_miniPFRelIso_chg->clear();
  t_Mu_tightCharge->clear();   
  t_Mu_isPFcand->clear();   
  t_Mu_isglobal->clear();
  t_Mu_istracker->clear();
  t_Mu_mediumId->clear();   
  t_Mu_softId->clear();   
  t_Mu_tightId->clear();    
  t_Mu_nStations->clear();   
  t_Mu_nTrackerLayers->clear();   

  t_diMuon_pt=-1000;
  t_diMuon_eta=-1000;
  t_diMuon_phi=-1000;
  t_diMuon_mass=-1000;

  t_MET_phi=-1000;
  t_MET_pt=-1000;
  t_MET_sumEt=-1000;

  t_FatJet_area->clear();  
  t_FatJet_btagCMVA->clear();  
  t_FatJet_btagCSVV2->clear();  
  t_FatJet_btagDeepB->clear();  
  t_FatJet_eta->clear();  
  t_FatJet_mass->clear();  
  t_FatJet_msoftdrop->clear();  
  t_FatJet_n2b1->clear();  
  t_FatJet_n3b1->clear();  
  t_FatJet_phi->clear();  
  t_FatJet_pt->clear();  
  t_FatJet_tau1->clear();  
  t_FatJet_tau2->clear();  
  t_FatJet_tau3->clear();  
  t_FatJet_tau4->clear();  
  t_FatJet_jetId->clear();  
  t_FatJet_subJetIdx1->clear();  
  t_FatJet_subJetIdx2->clear();  

  t_SubJet_btagCMVA->clear();   
  t_SubJet_btagCSVV2->clear();   
  t_SubJet_btagDeepB->clear();   
  t_SubJet_eta->clear();   
  t_SubJet_mass->clear();   
  t_SubJet_n2b1->clear();   
  t_SubJet_n3b1->clear();   
  t_SubJet_phi->clear();   
  t_SubJet_pt->clear();   
  t_SubJet_tau1->clear();   
  t_SubJet_tau2->clear();   
  t_SubJet_tau3->clear();   
  t_SubJet_tau4->clear();   
  
  t_Jet_area->clear();   
  //t_Jet_btagCMVA->clear();   
  //t_Jet_btagCSVV2->clear();   
  t_Jet_btagDeepB->clear();   
  //t_Jet_btagDeepC->clear();   
  //t_Jet_btagDeepFlavB->clear();   
  t_Jet_chEmEF->clear();   
  t_Jet_chHEF->clear();   
  t_Jet_eta->clear();   
  t_Jet_mass->clear();   
  t_Jet_neEmEF->clear();   
  t_Jet_neHEF->clear();   
  t_Jet_phi->clear();   
  t_Jet_pt->clear();   
  t_Jet_qgl->clear();   
  t_Jet_jetId->clear();   
  t_Jet_nConstituents->clear();   
  t_Jet_nElectrons->clear();   
  t_Jet_nMuons->clear();   
  t_Jet_puId->clear();   

  t_diJet_pt=-1000;
  t_diJet_eta=-1000;
  t_diJet_phi=-1000;
  t_diJet_mass=-1000;
  t_diJet_mass_mo=-1000;

  t_bJet_area->clear();   
  t_bJet_btagCMVA->clear();   
  t_bJet_btagCSVV2->clear();   
  t_bJet_btagDeepB->clear();   
  t_bJet_btagDeepC->clear();   
  t_bJet_btagDeepFlavB->clear();   
  t_bJet_chEmEF->clear();   
  t_bJet_chHEF->clear();   
  t_bJet_eta->clear();   
  t_bJet_mass->clear();   
  t_bJet_neEmEF->clear();   
  t_bJet_neHEF->clear();   
  t_bJet_phi->clear();   
  t_bJet_pt->clear();   
  t_bJet_qgl->clear();   
  t_bJet_jetId->clear();   
  t_bJet_nConstituents->clear();   
  t_bJet_nElectrons->clear();   
  t_bJet_nMuons->clear();   
  t_bJet_puId->clear();   
  t_bJet_SF->clear();
  t_bJet_SFup->clear();
  t_bJet_SFdown->clear();

  //SoftActivityJetHT5
  //SoftActivityJetNjets5
  
  t_PV_ndof=-1000;
  t_PV_x=-1000;
  t_PV_y=-1000;
  t_PV_z=-1000;
  t_PV_npvs=-1000;  
  t_PV_npvsGood=-1000;

  t_nLHEPdfWeight = 0;
  t_nLHEScaleWeight = 0;
  t_nPSWeight = 0;
  t_LHEPdfWeight->clear();
  t_LHEScaleWeight->clear();
  t_PSWeight->clear();

  t_GenPart_eta->clear();
  t_GenPart_mass->clear();
  t_GenPart_phi->clear();
  t_GenPart_pt->clear();
  t_GenPart_genPartIdxMother->clear();
  t_GenPart_pdgId->clear();
  t_GenPart_status->clear();

  t_GenJet_eta->clear();
  t_GenJet_mass->clear();
  t_GenJet_phi->clear();
  t_GenJet_pt->clear();
}

void HmmAnalyzer::BookTreeBranches(){
  tree = new TTree("tree","tree");
  tree->SetAutoSave(10000);

  tree->Branch("t_run", &t_run,"t_run/i");
  tree->Branch("t_luminosityBlock", &t_luminosityBlock,"t_luminosityBlock/i");
  tree->Branch("t_event", &t_event,"t_event/l");
  tree->Branch("t_genWeight", &t_genWeight,"t_genWeight/F");
  tree->Branch("t_puWeight", &t_puWeight,"t_puWeight/F");
  tree->Branch("t_puWeightUp", &t_puWeightUp,"t_puWeightUp/F");
  tree->Branch("t_puWeightDown", &t_puWeightDown,"t_puWeightDown/F");
  //tree->Branch("t_pileupWeight", &t_pileupWeight,"t_pileupWeight/F");
  //tree->Branch("t_pileupupWeight", &t_pileupupWeight,"t_pileupupWeight/F");
  //tree->Branch("t_pileupdnWeight", &t_pileupdnWeight,"t_pileupdnWeight/F");
  tree->Branch("t_PrefireWeight", &t_PrefireWeight, "t_PrefireWeight/F");
  tree->Branch("t_PrefireWeight_Up", &t_PrefireWeight_Up, "t_PrefireWeight_Up/F");
  tree->Branch("t_PrefireWeight_Down", &t_PrefireWeight_Down, "t_PrefireWeight_Down/F");
  tree->Branch("t_mu1", &t_mu1,"t_mu1/I");
  tree->Branch("t_mu2", &t_mu2,"t_mu2/I");
  tree->Branch("t_index_trigm_mu", &t_index_trigm_mu, "t_index_trigm_mu/I");

  t_El_genPartIdx= new std::vector<int>();
  t_El_genPartFlav= new std::vector<UChar_t>();
  t_El_charge= new std::vector<int>();
  t_El_pt= new std::vector<float>();
  t_El_phi= new std::vector<float>();
  t_El_eta= new std::vector<float>();   
  t_El_mass= new std::vector<float>();
  t_El_cutBased= new std::vector<int>();   
  t_El_tightCharge= new std::vector<int>();   
  t_El_cutBased_HEEP= new std::vector<bool>();   
  t_El_isPFcand= new std::vector<bool>();   
  t_El_pfRelIso03_all= new std::vector<float>();   
  t_El_pfRelIso03_chg= new std::vector<float>();      
  t_El_miniPFRelIso_all= new std::vector<float>();
  t_El_miniPFRelIso_chg= new std::vector<float>();
  t_El_dxy= new std::vector<float>();   
  t_El_dxyErr= new std::vector<float>();   
  t_El_dz= new std::vector<float>();   
  t_El_dzErr= new std::vector<float>();  
  t_El_sip3d= new std::vector<float>();
  t_Electron_mvaFall17Iso= new std::vector<float>(); 
  t_Electron_mvaFall17Iso_WP80= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17Iso_WP90= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17Iso_WPL= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso= new std::vector<float>();
  t_Electron_mvaFall17noIso_WP80= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso_WP90= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso_WPL= new std::vector<bool>();   //[nElectron]
 
  tree->Branch("t_El_genPartIdx",   "vector<int>",        &t_El_genPartIdx);
  tree->Branch("t_El_genPartFlav",   "vector<UChar_t>",       &t_El_genPartFlav); 
  tree->Branch("t_El_charge",        "vector<int>"         ,&t_El_charge);
  tree->Branch("t_El_pt",        "vector<float>"         ,&t_El_pt);
  tree->Branch("t_El_phi",        "vector<float>"         ,&t_El_phi);
  tree->Branch("t_El_eta",        "vector<float>"         ,&t_El_eta);   
  tree->Branch("t_El_mass",        "vector<float>"         ,&t_El_mass);
  tree->Branch("t_El_cutBased",        "vector<int>"         ,&t_El_cutBased);   
  tree->Branch("t_El_tightCharge",        "vector<int>"         ,&t_El_tightCharge);   
  tree->Branch("t_El_cutBased_HEEP",        "vector<bool>"         ,&t_El_cutBased_HEEP);   
  tree->Branch("t_El_isPFcand",        "vector<bool>"         ,&t_El_isPFcand);   
  tree->Branch("t_El_pfRelIso03_all",        "vector<float>"         ,&t_El_pfRelIso03_all);   
  tree->Branch("t_El_pfRelIso03_chg",        "vector<float>"         ,&t_El_pfRelIso03_chg);      
  tree->Branch("t_El_miniPFRelIso03_all",        "vector<float>"         ,&t_El_miniPFRelIso_all);
  tree->Branch("t_El_miniPFRelIso03_chg",        "vector<float>"         ,&t_El_miniPFRelIso_chg);
  tree->Branch("t_El_dxy",        "vector<float>"         ,&t_El_dxy);   
  tree->Branch("t_El_dxyErr",        "vector<float>"         ,&t_El_dxyErr);   
  tree->Branch("t_El_dz",        "vector<float>"         ,&t_El_dz);   
  tree->Branch("t_El_dzErr",        "vector<float>"         ,&t_El_dzErr);   
  tree->Branch("t_El_sip3d",        "vector<float>"         ,&t_El_sip3d);
  tree->Branch("t_Electron_mvaFall17Iso",        "vector<float>",         &t_Electron_mvaFall17Iso); 
  tree->Branch("t_Electron_mvaFall17Iso_WP80",        "vector<bool>",         &t_Electron_mvaFall17Iso_WP80);
  tree->Branch("t_Electron_mvaFall17Iso_WP90",        "vector<bool>",         &t_Electron_mvaFall17Iso_WP90);
  tree->Branch("t_Electron_mvaFall17Iso_WPL",        "vector<bool>",         &t_Electron_mvaFall17Iso_WPL);
  tree->Branch("t_Electron_mvaFall17noIso",        "vector<float>",         &t_Electron_mvaFall17noIso);
  tree->Branch("t_Electron_mvaFall17noIso_WP80",        "vector<bool>",         &t_Electron_mvaFall17noIso_WP80);
  tree->Branch("t_Electron_mvaFall17noIso_WP90",        "vector<bool>",         &t_Electron_mvaFall17noIso_WP90);
  tree->Branch("t_Electron_mvaFall17noIso_WPL",        "vector<bool>",         &t_Electron_mvaFall17noIso_WPL);

  t_LHEPdfWeight = new std::vector<float>();
  t_LHEScaleWeight = new std::vector<float>();
  t_PSWeight = new std::vector<float>(); 
  tree->Branch("t_nLHEPdfWeight", &t_nLHEPdfWeight, "t_nLHEPdfWeight/i");
  tree->Branch("t_LHEPdfWeight","vector<float>", &t_LHEPdfWeight);
  tree->Branch("t_nLHEScaleWeight", &t_nLHEScaleWeight, "t_nLHEScaleWeight/i");
  tree->Branch("t_LHEScaleWeight", "vector<float>", &t_LHEScaleWeight);
  tree->Branch("t_nPSWeight", &t_nPSWeight, "t_nPSWeight/i");
  tree->Branch("t_PSWeight", "vector<float>", &t_PSWeight);

 
  t_Mu_genPartIdx= new std::vector<int>();
  t_Mu_genPartFlav= new std::vector<UChar_t>(); 
  t_Mu_charge= new std::vector<int>();   
  t_Mu_EffSF_TRIG= new std::vector<float>();
  t_Mu_EffSFErr_TRIG= new std::vector<float>();  
  t_Mu_EffSF_ID= new std::vector<float>();
  t_Mu_EffSF_ID_stat= new std::vector<float>();
  t_Mu_EffSF_ID_syst= new std::vector<float>();
  t_Mu_EffSF_ISO= new std::vector<float>();
  t_Mu_EffSF_ISO_stat= new std::vector<float>();
  t_Mu_EffSF_ISO_syst= new std::vector<float>();
  t_Mu_pt= new std::vector<float>();   
  t_Mu_ptErr= new std::vector<float>();   
  t_Mu_phi= new std::vector<float>();   
  t_Mu_eta= new std::vector<float>();   
  t_Mu_mass= new std::vector<float>();  
  t_Mu_dxy= new std::vector<float>();   
  t_Mu_dxyErr= new std::vector<float>();   
  t_Mu_dz= new std::vector<float>();   
  t_Mu_dzErr= new std::vector<float>(); 
  t_Mu_sip3d= new std::vector<float>();
  t_Mu_pfRelIso03_all= new std::vector<float>();   
  t_Mu_pfRelIso03_chg= new std::vector<float>();   
  t_Mu_pfRelIso04_all= new std::vector<float>();   
  t_Mu_miniPFRelIso_all= new std::vector<float>();
  t_Mu_miniPFRelIso_chg= new std::vector<float>();
  t_Mu_tightCharge= new std::vector<int>();   
  t_Mu_isPFcand= new std::vector<bool>();  
  t_Mu_isglobal= new std::vector<bool>();
  t_Mu_istracker= new std::vector<bool>(); 
  t_Mu_mediumId= new std::vector<bool>();   
  t_Mu_softId= new std::vector<bool>();   
  t_Mu_tightId= new std::vector<bool>();    
  t_Mu_nStations= new std::vector<int>();   
  t_Mu_nTrackerLayers= new std::vector<int>();   

  tree->Branch("t_Mu_genPartIdx"    , "vector<int>"     ,&t_Mu_genPartIdx);
  tree->Branch("t_Mu_genPartFlav",  "vector<UChar_t>",   &t_Mu_genPartFlav);
  tree->Branch("t_Mu_charge"    , "vector<int>"         ,&t_Mu_charge );
  tree->Branch("t_Mu_EffSF_TRIG",    "vector<float>"   ,&t_Mu_EffSF_TRIG);
  tree->Branch("t_Mu_EffSFErr_TRIG",    "vector<float>"   ,&t_Mu_EffSFErr_TRIG);
  tree->Branch("t_Mu_EffSF_ID",    "vector<float>"   ,&t_Mu_EffSF_ID);
  tree->Branch("t_Mu_EffSF_ID_stat",    "vector<float>"   ,&t_Mu_EffSF_ID_stat);
  tree->Branch("t_Mu_EffSF_ID_syst",    "vector<float>"   ,&t_Mu_EffSF_ID_syst);
  tree->Branch("t_Mu_EffSF_ISO",    "vector<float>"   ,&t_Mu_EffSF_ISO);
  tree->Branch("t_Mu_EffSF_ISO_stat",    "vector<float>"   ,&t_Mu_EffSF_ISO_stat);
  tree->Branch("t_Mu_EffSF_ISO_syst",    "vector<float>"   ,&t_Mu_EffSF_ISO_syst);
  tree->Branch("t_Mu_pt"    , "vector<float>"         ,&t_Mu_pt );   
  tree->Branch("t_Mu_ptErr"    , "vector<float>"         ,&t_Mu_ptErr );   
  tree->Branch("t_Mu_phi"    , "vector<float>"         ,&t_Mu_phi );   
  tree->Branch("t_Mu_eta"    , "vector<float>"         ,&t_Mu_eta );   
  tree->Branch("t_Mu_mass"    , "vector<float>"         ,&t_Mu_mass );  
  tree->Branch("t_Mu_dxy"    , "vector<float>"         ,&t_Mu_dxy );   
  tree->Branch("t_Mu_dxyErr"    , "vector<float>"         ,&t_Mu_dxyErr );   
  tree->Branch("t_Mu_dz"    , "vector<float>"         ,&t_Mu_dz );   
  tree->Branch("t_Mu_dzErr"    , "vector<float>"         ,&t_Mu_dzErr );   
  tree->Branch("t_Mu_sip3d"   , "vector<float>"         ,&t_Mu_sip3d );
  tree->Branch("t_Mu_pfRelIso03_all"    , "vector<float>"         ,&t_Mu_pfRelIso03_all );   
  tree->Branch("t_Mu_pfRelIso03_chg"    , "vector<float>"         ,&t_Mu_pfRelIso03_chg );   
  tree->Branch("t_Mu_pfRelIso04_all"    , "vector<float>"         ,&t_Mu_pfRelIso04_all );   
  tree->Branch("t_Mu_miniPFRelIso_all"    , "vector<float>"         ,&t_Mu_miniPFRelIso_all );
  tree->Branch("t_Mu_miniPFRelIso_chg"    , "vector<float>"         ,&t_Mu_miniPFRelIso_chg );
  tree->Branch("t_Mu_tightCharge"    , "vector<int>"         ,&t_Mu_tightCharge );   
  tree->Branch("t_Mu_isPFcand"    , "vector<bool>"         ,&t_Mu_isPFcand );   
  tree->Branch("t_Mu_isglobal"    , "vector<bool>"         ,&t_Mu_isglobal );
  tree->Branch("t_Mu_istracker"    , "vector<bool>"         ,&t_Mu_istracker );
  tree->Branch("t_Mu_mediumId"    , "vector<bool>"         ,&t_Mu_mediumId );   
  tree->Branch("t_Mu_softId"    , "vector<bool>"         ,&t_Mu_softId );   
  tree->Branch("t_Mu_tightId"    , "vector<bool>"         ,&t_Mu_tightId );    
  tree->Branch("t_Mu_nStations"    , "vector<int>"         ,&t_Mu_nStations );   
  tree->Branch("t_Mu_nTrackerLayers"    , "vector<int>"         ,&t_Mu_nTrackerLayers );

  tree->Branch("t_diMuon_pt",   &t_diMuon_pt,"t_diMuon_pt/F");  
  tree->Branch("t_diMuon_eta",   &t_diMuon_eta,"t_diMuon_eta/F");
  tree->Branch("t_diMuon_phi",  &t_diMuon_phi,"t_diMuon_phi/F");
  tree->Branch("t_diMuon_mass",   &t_diMuon_mass,"t_diMuon_mass/F");

  tree->Branch("t_MET_phi",  &t_MET_phi,"t_MET_phi/F");
  tree->Branch("t_MET_pt",   &t_MET_pt,"t_MET_pt/F");
  tree->Branch("t_MET_sumEt",&t_MET_sumEt,"t_MET_sumEt/F");
  
  t_FatJet_area= new std::vector<float>();  
  t_FatJet_btagCMVA= new std::vector<float>();  
  t_FatJet_btagCSVV2= new std::vector<float>();  
  t_FatJet_btagDeepB= new std::vector<float>();  
  t_FatJet_eta= new std::vector<float>();  
  t_FatJet_mass= new std::vector<float>();  
  t_FatJet_msoftdrop= new std::vector<float>();  
  t_FatJet_n2b1= new std::vector<float>();  
  t_FatJet_n3b1= new std::vector<float>();  
  t_FatJet_phi= new std::vector<float>();  
  t_FatJet_pt= new std::vector<float>();  
  t_FatJet_tau1= new std::vector<float>();  
  t_FatJet_tau2= new std::vector<float>();  
  t_FatJet_tau3= new std::vector<float>();  
  t_FatJet_tau4= new std::vector<float>();  
  t_FatJet_jetId= new std::vector<int>();  
  t_FatJet_subJetIdx1= new std::vector<int>();  
  t_FatJet_subJetIdx2= new std::vector<int>();  

  tree->Branch("t_FatJet_area",         "vector<float>", &t_FatJet_area);  
  tree->Branch("t_FatJet_btagCMVA",         "vector<float>", &t_FatJet_btagCMVA);  
  tree->Branch("t_FatJet_btagCSVV2",         "vector<float>", &t_FatJet_btagCSVV2);  
  tree->Branch("t_FatJet_btagDeepB",         "vector<float>", &t_FatJet_btagDeepB);  
  tree->Branch("t_FatJet_eta",         "vector<float>", &t_FatJet_eta);  
  tree->Branch("t_FatJet_mass",         "vector<float>", &t_FatJet_mass);  
  tree->Branch("t_FatJet_msoftdrop",         "vector<float>", &t_FatJet_msoftdrop);  
  tree->Branch("t_FatJet_n2b1",         "vector<float>", &t_FatJet_n2b1);  
  tree->Branch("t_FatJet_n3b1",         "vector<float>", &t_FatJet_n3b1);  
  tree->Branch("t_FatJet_phi",         "vector<float>", &t_FatJet_phi);  
  tree->Branch("t_FatJet_pt",         "vector<float>", &t_FatJet_pt);  
  tree->Branch("t_FatJet_tau1",         "vector<float>", &t_FatJet_tau1);  
  tree->Branch("t_FatJet_tau2",         "vector<float>", &t_FatJet_tau2);  
  tree->Branch("t_FatJet_tau3",         "vector<float>", &t_FatJet_tau3);  
  tree->Branch("t_FatJet_tau4",         "vector<float>", &t_FatJet_tau4);  
  tree->Branch("t_FatJet_jetId",         "vector<int>", &t_FatJet_jetId);  
  tree->Branch("t_FatJet_subJetIdx1",         "vector<int>", &t_FatJet_subJetIdx1);  
  tree->Branch("t_FatJet_subJetIdx2",         "vector<int>", &t_FatJet_subJetIdx2);  
  
  t_SubJet_btagCMVA= new std::vector<float>();   
  t_SubJet_btagCSVV2= new std::vector<float>();   
  t_SubJet_btagDeepB= new std::vector<float>();   
  t_SubJet_eta= new std::vector<float>();   
  t_SubJet_mass= new std::vector<float>();   
  t_SubJet_n2b1= new std::vector<float>();   
  t_SubJet_n3b1= new std::vector<float>();   
  t_SubJet_phi= new std::vector<float>();   
  t_SubJet_pt= new std::vector<float>();   
  t_SubJet_tau1= new std::vector<float>();   
  t_SubJet_tau2= new std::vector<float>();   
  t_SubJet_tau3= new std::vector<float>();   
  t_SubJet_tau4= new std::vector<float>();   
  
  tree->Branch("t_SubJet_btagCMVA",        "vector<float>", &t_SubJet_btagCMVA);   
  tree->Branch("t_SubJet_btagCSVV2",        "vector<float>", &t_SubJet_btagCSVV2);   
  tree->Branch("t_SubJet_btagDeepB",        "vector<float>", &t_SubJet_btagDeepB);   
  tree->Branch("t_SubJet_eta",        "vector<float>", &t_SubJet_eta);   
  tree->Branch("t_SubJet_mass",        "vector<float>", &t_SubJet_mass);   
  tree->Branch("t_SubJet_n2b1",        "vector<float>", &t_SubJet_n2b1);   
  tree->Branch("t_SubJet_n3b1",        "vector<float>", &t_SubJet_n3b1);   
  tree->Branch("t_SubJet_phi",        "vector<float>", &t_SubJet_phi);   
  tree->Branch("t_SubJet_pt",        "vector<float>", &t_SubJet_pt);   
  tree->Branch("t_SubJet_tau1",        "vector<float>", &t_SubJet_tau1);   
  tree->Branch("t_SubJet_tau2",        "vector<float>", &t_SubJet_tau2);   
  tree->Branch("t_SubJet_tau3",        "vector<float>", &t_SubJet_tau3);   
  tree->Branch("t_SubJet_tau4",        "vector<float>", &t_SubJet_tau4);   
   

  t_Jet_area= new std::vector<float>();   
  //t_Jet_btagCMVA= new std::vector<float>();   
  //t_Jet_btagCSVV2= new std::vector<float>();   
  t_Jet_btagDeepB= new std::vector<float>();   
  //t_Jet_btagDeepC= new std::vector<float>();   
  //t_Jet_btagDeepFlavB= new std::vector<float>();   
  t_Jet_chEmEF= new std::vector<float>();   
  t_Jet_chHEF= new std::vector<float>();   
  t_Jet_eta= new std::vector<float>();   
  t_Jet_mass= new std::vector<float>();   
  t_Jet_neEmEF= new std::vector<float>();   
  t_Jet_neHEF= new std::vector<float>();   
  t_Jet_phi= new std::vector<float>();   
  t_Jet_pt= new std::vector<float>();   
  t_Jet_qgl= new std::vector<float>();   
  t_Jet_jetId= new std::vector<int>();   
  t_Jet_nConstituents= new std::vector<int>();   
  t_Jet_nElectrons= new std::vector<int>();   
  t_Jet_nMuons= new std::vector<int>();   
  t_Jet_puId= new std::vector<int>();   

  tree->Branch("t_SoftActivityJetNjets5"    , &t_SoftActivityJetNjets5, "t_SoftActivityJetNjets5/I");
  tree->Branch("t_nJet",  &t_nJet,"t_nJet/I");
  tree->Branch("t_Jet_area"    , "vector<float>"         ,&t_Jet_area);   
  //tree->Branch("t_Jet_btagCMVA"    , "vector<float>"         ,&t_Jet_btagCMVA);   
  //tree->Branch("t_Jet_btagCSVV2"    , "vector<float>"         ,&t_Jet_btagCSVV2);   
  tree->Branch("t_Jet_btagDeepB"    , "vector<float>"         ,&t_Jet_btagDeepB);   
  //tree->Branch("t_Jet_btagDeepC"    , "vector<float>"         ,&t_Jet_btagDeepC);   
  //tree->Branch("t_Jet_btagDeepFlavB"    , "vector<float>"         ,&t_Jet_btagDeepFlavB);   
  tree->Branch("t_Jet_chEmEF"    , "vector<float>"         ,&t_Jet_chEmEF);   
  tree->Branch("t_Jet_chHEF"    , "vector<float>"         ,&t_Jet_chHEF);   
  tree->Branch("t_Jet_eta"    , "vector<float>"         ,&t_Jet_eta);   
  tree->Branch("t_Jet_mass"    , "vector<float>"         ,&t_Jet_mass);   
  tree->Branch("t_Jet_neEmEF"    , "vector<float>"         ,&t_Jet_neEmEF);   
  tree->Branch("t_Jet_neHEF"    , "vector<float>"         ,&t_Jet_neHEF);   
  tree->Branch("t_Jet_phi"    , "vector<float>"         ,&t_Jet_phi);   
  tree->Branch("t_Jet_pt"    , "vector<float>"         ,&t_Jet_pt);   
  tree->Branch("t_Jet_qgl"    , "vector<float>"         ,&t_Jet_qgl);   
  tree->Branch("t_Jet_jetId"    , "vector<int>"         ,&t_Jet_jetId);   
  tree->Branch("t_Jet_nConstituents"    , "vector<int>"         ,&t_Jet_nConstituents);   
  tree->Branch("t_Jet_nElectrons"    , "vector<int>"         ,&t_Jet_nElectrons);   
  tree->Branch("t_Jet_nMuons"    , "vector<int>"         ,&t_Jet_nMuons);   
  tree->Branch("t_Jet_puId"    , "vector<int>"         ,&t_Jet_puId);   

  tree->Branch("t_diJet_pt",   &t_diJet_pt,"t_diJet_pt/F");  
  tree->Branch("t_diJet_eta",   &t_diJet_eta,"t_diJet_eta/F");
  tree->Branch("t_diJet_phi",  &t_diJet_phi,"t_diJet_phi/F");
  tree->Branch("t_diJet_mass",   &t_diJet_mass,"t_diJet_mass/F");
  tree->Branch("t_diJet_mass_mo",   &t_diJet_mass_mo,"t_diJet_mass_mo/F");

  tree->Branch("t_nbJet",  &t_nbJet,"t_nbJet/I");
  t_bJet_area= new std::vector<float>();   
  t_bJet_btagCMVA= new std::vector<float>();   
  t_bJet_btagCSVV2= new std::vector<float>();   
  t_bJet_btagDeepB= new std::vector<float>();   
  t_bJet_btagDeepC= new std::vector<float>();   
  t_bJet_btagDeepFlavB= new std::vector<float>();   
  t_bJet_chEmEF= new std::vector<float>();   
  t_bJet_chHEF= new std::vector<float>();   
  t_bJet_eta= new std::vector<float>();   
  t_bJet_mass= new std::vector<float>();   
  t_bJet_neEmEF= new std::vector<float>();   
  t_bJet_neHEF= new std::vector<float>();   
  t_bJet_phi= new std::vector<float>();   
  t_bJet_pt= new std::vector<float>();   
  t_bJet_qgl= new std::vector<float>();   
  t_bJet_jetId= new std::vector<int>();   
  t_bJet_nConstituents= new std::vector<int>();   
  t_bJet_nElectrons= new std::vector<int>();   
  t_bJet_nMuons= new std::vector<int>();   
  t_bJet_puId= new std::vector<int>();   
  t_bJet_SF= new std::vector<double>();
  t_bJet_SFup= new std::vector<double>();
  t_bJet_SFdown= new std::vector<double>();
  
  tree->Branch("t_bJet_area"    , "vector<float>"         ,&t_bJet_area);   
  tree->Branch("t_bJet_btagCMVA"    , "vector<float>"         ,&t_bJet_btagCMVA);   
  tree->Branch("t_bJet_btagCSVV2"    , "vector<float>"         ,&t_bJet_btagCSVV2);   
  tree->Branch("t_bJet_btagDeepB"    , "vector<float>"         ,&t_bJet_btagDeepB);   
  tree->Branch("t_bJet_btagDeepC"    , "vector<float>"         ,&t_bJet_btagDeepC);   
  tree->Branch("t_bJet_btagDeepFlavB"    , "vector<float>"         ,&t_bJet_btagDeepFlavB);   
  tree->Branch("t_bJet_chEmEF"    , "vector<float>"         ,&t_bJet_chEmEF);   
  tree->Branch("t_bJet_chHEF"    , "vector<float>"         ,&t_bJet_chHEF);   
  tree->Branch("t_bJet_eta"    , "vector<float>"         ,&t_bJet_eta);   
  tree->Branch("t_bJet_mass"    , "vector<float>"         ,&t_bJet_mass);   
  tree->Branch("t_bJet_neEmEF"    , "vector<float>"         ,&t_bJet_neEmEF);   
  tree->Branch("t_bJet_neHEF"    , "vector<float>"         ,&t_bJet_neHEF);   
  tree->Branch("t_bJet_phi"    , "vector<float>"         ,&t_bJet_phi);   
  tree->Branch("t_bJet_pt"    , "vector<float>"         ,&t_bJet_pt);   
  tree->Branch("t_bJet_qgl"    , "vector<float>"         ,&t_bJet_qgl);   
  tree->Branch("t_bJet_jetId"    , "vector<int>"         ,&t_bJet_jetId);   
  tree->Branch("t_bJet_nConstituents"    , "vector<int>"         ,&t_bJet_nConstituents);   
  tree->Branch("t_bJet_nElectrons"    , "vector<int>"         ,&t_bJet_nElectrons);   
  tree->Branch("t_bJet_nMuons"    , "vector<int>"         ,&t_bJet_nMuons);   
  tree->Branch("t_bJet_puId"    , "vector<int>"         ,&t_bJet_puId);   
  tree->Branch("t_bJet_SF"    , "vector<double>"         ,&t_bJet_SF);
  tree->Branch("t_bJet_SFup"    , "vector<double>"         ,&t_bJet_SFup);
  tree->Branch("t_bJet_SFdown"    , "vector<double>"     ,&t_bJet_SFdown);


  tree->Branch("t_PV_ndof", &t_PV_ndof, "t_PV_ndof/F");
  tree->Branch("t_PV_x", &t_PV_x,"t_PV_x/F");
  tree->Branch("t_PV_y", &t_PV_y,"t_PV_y/F");
  tree->Branch("t_PV_z", &t_PV_z,"t_PV_z/F");
  tree->Branch("t_PV_npvs",&t_PV_npvs,"t_PV_npvs/I");  
  tree->Branch("t_PV_npvsGood",&t_PV_npvsGood,"t_PV_npvsGood/I"); 

  
  t_GenPart_eta= new std::vector<float>();
  t_GenPart_mass= new std::vector<float>();
  t_GenPart_phi= new std::vector<float>();
  t_GenPart_pt= new std::vector<float>();
  t_GenPart_genPartIdxMother= new std::vector<int>();
  t_GenPart_pdgId= new std::vector<int>();
  t_GenPart_status= new std::vector<int>();

  t_GenJet_eta= new std::vector<float>();
  t_GenJet_mass= new std::vector<float>();
  t_GenJet_phi= new std::vector<float>();
  t_GenJet_pt= new std::vector<float>();

  tree->Branch("t_GenPart_eta",      "vector<float>", &t_GenPart_eta);
  tree->Branch("t_GenPart_mass",      "vector<float>", &t_GenPart_mass);
  tree->Branch("t_GenPart_phi",      "vector<float>", &t_GenPart_phi);
  tree->Branch("t_GenPart_pt",      "vector<float>", &t_GenPart_pt);
  tree->Branch("t_GenPart_genPartIdxMother",      "vector<int>", &t_GenPart_genPartIdxMother);
  tree->Branch("t_GenPart_pdgId",      "vector<int>", &t_GenPart_pdgId);
  tree->Branch("t_GenPart_status",      "vector<int>", &t_GenPart_status);
  
  tree->Branch("t_GenJet_eta",      "vector<float>", &t_GenJet_eta);
  tree->Branch("t_GenJet_mass",      "vector<float>", &t_GenJet_mass);
  tree->Branch("t_GenJet_phi",      "vector<float>", &t_GenJet_phi);
  tree->Branch("t_GenJet_pt",      "vector<float>", &t_GenJet_pt);

}
#endif // #ifdef HmmAnalyzer_cxx
