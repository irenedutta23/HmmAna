/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 10 04:03:58 2018 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: ../condor/condor_output/condor_logs/ElData16_job0.root
//////////////////////////////////////////////////////////

#ifndef HiggsMuMu_h
#define HiggsMuMu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "NtupleVariables.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
#include<map>
class HiggsMuMu:public NtupleVariables  {
public :

  HiggsMuMu(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *isData="F");
  virtual ~HiggsMuMu();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *);
  void     Categorization(const char *, const char *, float , float,double, TString procname);
  void     GenInfo(const char *, const char *);
  void     BookHistogram(const char *);
  void     clearTreeVectors();
 
  TFile *oFile;
  TTree *cattree;

  UInt_t    run;
  UInt_t    lumi;
  ULong64_t event;
  float genWeight = -9999.;
  float genWeight_Up = -9999.;
  float genWeight_Down = -9999.;
  int cat_index = -9999.;
  float Higgs_mass = -9999.;
  float Higgs_pt = -9999.;
  float Higgs_eta = -9999.;
  float extralep_pfRelIso03 = -9999.;
  float extralep_pt = -9999.;
  float extralep_eta = -9999.;
  float dRlepHiggs = -9999.;
  float dEtamm = -9999.;
  int softJet5=9999;
  float dRmm = -9999.;
  float dPhimm = -9999.;  
  float MET_phi = -9999.;
  float MET_pt = -9999.;
  float dEta_jj = -9999.;
  float M_jj= -9999.;
  float pt_jj=-9999;
  float eta_jj=-9999;
  float phi_jj= -9999;
  float M_mmjj= -9999.;
  float pt_mmjj=-9999;
  float eta_mmjj=-9999;
  float phi_mmjj= -9999;
  float Zep = -9999.;
  float dRmin_mj = -9999.;
  float dRmax_mj = -9999.;
  float dRmin_mmj = -9999.;
  float dRmax_mmj= -9999.;
  float leadingJet_pt = -9999.;
  float subleadingJet_pt = -9999.;
  float leadingJet_eta = -9999.;
  float leadingJet_qgl = -9999.;
  float subleadingJet_qgl = -9999.;
  float subleadingJet_eta = -9999.;
  float dPhijj=-9999.;
  float cthetaCS = -9999.;
  float BDT_incl = -9999.;
  std::vector<int>   *l1_index;
  std::vector<int>   *l2_index;

  //BDT vars
   float ll_mass =0.;
   float MqqLog=0.;
   float mumujj_pt=0.;
   float DeltaEtaQQ=0.;
   float softActivityEWK_njets5=0.;
   float ll_zstar=0.;
   float ll_pt=0.;
   float theta2=0.;
   float impulsoZ=0.;
   float maxAbsEta=0.;
   float qgl_2qAtanh=0.;

  //define histograms
  TH1D *catyield;
  TH1D *h_diMuon_mass_cat;
  TH1D *h_diMuon_mass_ggH;
  TH1D *h_diMuon_mass_VBF;
  TH1D *h_diMuon_mass_Up_VBF;
  TH1D *h_diMuon_mass_Down_VBF;
  TH1D *h_diMuon_mass_Z;
  TH1D *h_diMuon_mass_Up_Z;
  TH1D *h_diMuon_mass_Down_Z;
  TH1D *h_diMuon_mass_VHHad;
  TH1D *h_diMuon_mass_ZHll;
  TH1D *h_diMuon_mass_WHlv;
  TH1D *h_diMuon_mass_ttHLep;
  TH1D *h_diMuon_mass_ttHHad;  
  TH1D *h_diMuon_mass_ttHLoose;
  TH1D *h_dRlepH;
  TH1D *h_extralep1_pt;
  TH1D *h_extralep1_eta;
  TH1D *h_extralep_Electron_mvaFall17Iso;

  TH1D *h_diMuon_mass_110To140_VBF;
  TH1D *h_diMuon_mass_110To140_Up_VBF;
  TH1D *h_diMuon_mass_110To140_Down_VBF;
  TH1D *h_diMuon_mass_110To140_WHlv;
  TH1D *h_diMuon_mass_110To140_ttHLep;
  TH1D *h_diMuon_mass_110To140_ttHHad;  
  TH1D *h_diMuon_mass_110To140_ttHLoose;

  TH1D *h_mu1mu2dR_VBF;
  TH1D *h_mu1mu2dPhi_VBF;
  TH1D *h_mu1mu2dEta_VBF;
  TH1D *h_diMuon_pt_VBF;
  TH1D *h_diMuon_eta_VBF;
  TH1D *h_diMuon_phi_VBF;
  TH1D *h_dijet_pt_VBF;
  TH1D *h_dijet_eta_VBF;
  TH1D *h_dijet_phi_VBF;
  TH1D *h_Mjj_VBF;
  TH1D *h_dijet_dEta_VBF;
  TH1D *h_dijet_dPhijj_VBF;
  TH1D *h_jet1_pt_VBF;
  TH1D *h_jet2_pt_VBF;
  TH1D *h_jet1_eta_VBF;
  TH1D *h_jet2_eta_VBF;
  TH1D *h_cthetaCS_VBF;
  TH1D *h_dRmin_mj_VBF;
  TH1D *h_dRmax_mj_VBF;
  TH1D *h_dRmin_mmj_VBF;
  TH1D *h_dRmax_mmj_VBF;
  TH1D *h_Zep_VBF;
  TH1D *h_M_mmjj_VBF;
  TH1D *h_pt_mmjj_VBF;
  TH1D *h_eta_mmjj_VBF;
  TH1D *h_phi_mmjj_VBF;
  TH1D *h_leadingJet_qgl_VBF;
  TH1D *h_subleadingJet_qgl_VBF;
  TH1D *h_softJet5_VBF;
  TH2D *h_2D_ptjj_dRmm;
  TH2D *h_2D_CS_dEtamm;
  TH2D *h_2D_Hpt_ptjj;
  TH2D *h_2D_Hpt_dRmm;
  TH2D *h_2D_j1eta_dEtajj;
  TH2D *h_2D_j2eta_dEtajj;
  TH2D *h_2D_etammjj_etajj;
  TH2D *h_2D_j1pt_ptjj;

  TH1D *h_mu1mu2dR_Up_VBF;
  TH1D *h_mu1mu2dPhi_Up_VBF;
  TH1D *h_mu1mu2dEta_Up_VBF;
  TH1D *h_diMuon_pt_Up_VBF;
  TH1D *h_diMuon_eta_Up_VBF;
  TH1D *h_diMuon_phi_Up_VBF;
  TH1D *h_dijet_pt_Up_VBF;
  TH1D *h_dijet_eta_Up_VBF;
  TH1D *h_dijet_phi_Up_VBF;
  TH1D *h_Mjj_Up_VBF;
  TH1D *h_dijet_dEta_Up_VBF;
  TH1D *h_dijet_dPhijj_Up_VBF;
  TH1D *h_jet1_pt_Up_VBF;
  TH1D *h_jet2_pt_Up_VBF;
  TH1D *h_jet1_eta_Up_VBF;
  TH1D *h_jet2_eta_Up_VBF;
  TH1D *h_cthetaCS_Up_VBF;
  TH1D *h_dRmin_mj_Up_VBF;
  TH1D *h_dRmax_mj_Up_VBF;
  TH1D *h_dRmin_mmj_Up_VBF;
  TH1D *h_dRmax_mmj_Up_VBF;
  TH1D *h_Zep_Up_VBF;
  TH1D *h_M_mmjj_Up_VBF;
  TH1D *h_pt_mmjj_Up_VBF;
  TH1D *h_eta_mmjj_Up_VBF;
  TH1D *h_phi_mmjj_Up_VBF;
  TH1D *h_leadingJet_qgl_Up_VBF;
  TH1D *h_subleadingJet_qgl_Up_VBF;
  TH1D *h_softJet5_Up_VBF;


  TH1D *h_mu1mu2dR_Down_VBF;
  TH1D *h_mu1mu2dPhi_Down_VBF;
  TH1D *h_mu1mu2dEta_Down_VBF;
  TH1D *h_diMuon_pt_Down_VBF;
  TH1D *h_diMuon_eta_Down_VBF;
  TH1D *h_diMuon_phi_Down_VBF;
  TH1D *h_dijet_pt_Down_VBF;
  TH1D *h_dijet_eta_Down_VBF;
  TH1D *h_dijet_phi_Down_VBF;
  TH1D *h_Mjj_Down_VBF;
  TH1D *h_dijet_dEta_Down_VBF;
  TH1D *h_dijet_dPhijj_Down_VBF;
  TH1D *h_jet1_pt_Down_VBF;
  TH1D *h_jet2_pt_Down_VBF;
  TH1D *h_jet1_eta_Down_VBF;
  TH1D *h_jet2_eta_Down_VBF;
  TH1D *h_cthetaCS_Down_VBF;
  TH1D *h_dRmin_mj_Down_VBF;
  TH1D *h_dRmax_mj_Down_VBF;
  TH1D *h_dRmin_mmj_Down_VBF;
  TH1D *h_dRmax_mmj_Down_VBF;
  TH1D *h_Zep_Down_VBF;
  TH1D *h_M_mmjj_Down_VBF;
  TH1D *h_pt_mmjj_Down_VBF;
  TH1D *h_eta_mmjj_Down_VBF;
  TH1D *h_phi_mmjj_Down_VBF;
  TH1D *h_leadingJet_qgl_Down_VBF;
  TH1D *h_subleadingJet_qgl_Down_VBF;
  TH1D *h_softJet5_Down_VBF;

  TH1D *h_mu1mu2dR_Z;
  TH1D *h_mu1mu2dPhi_Z;
  TH1D *h_mu1mu2dEta_Z;
  TH1D *h_diMuon_pt_Z;
  TH1D *h_diMuon_eta_Z;
  TH1D *h_diMuon_phi_Z;
  TH1D *h_dijet_pt_Z;
  TH1D *h_dijet_eta_Z;
  TH1D *h_dijet_phi_Z;
  TH1D *h_Mjj_Z;
  TH1D *h_dijet_dEta_Z;
  TH1D *h_dijet_dPhijj_Z;
  TH1D *h_jet1_pt_Z;
  TH1D *h_jet2_pt_Z;
  TH1D *h_jet1_eta_Z;
  TH1D *h_jet2_eta_Z;
  TH1D *h_cthetaCS_Z;
  TH1D *h_dRmin_mj_Z;
  TH1D *h_dRmax_mj_Z;
  TH1D *h_dRmin_mmj_Z;
  TH1D *h_dRmax_mmj_Z;
  TH1D *h_Zep_Z;
  TH1D *h_M_mmjj_Z;
  TH1D *h_pt_mmjj_Z;
  TH1D *h_eta_mmjj_Z;
  TH1D *h_phi_mmjj_Z;
  TH1D *h_leadingJet_qgl_Z;
  TH1D *h_subleadingJet_qgl_Z;
  TH1D *h_softJet5_Z;


  TH1D *h_mu1mu2dR_Up_Z;
  TH1D *h_mu1mu2dPhi_Up_Z;
  TH1D *h_mu1mu2dEta_Up_Z;
  TH1D *h_diMuon_pt_Up_Z;
  TH1D *h_diMuon_eta_Up_Z;
  TH1D *h_diMuon_phi_Up_Z;
  TH1D *h_dijet_pt_Up_Z;
  TH1D *h_dijet_eta_Up_Z;
  TH1D *h_dijet_phi_Up_Z;
  TH1D *h_Mjj_Up_Z;
  TH1D *h_dijet_dEta_Up_Z;
  TH1D *h_dijet_dPhijj_Up_Z;
  TH1D *h_jet1_pt_Up_Z;
  TH1D *h_jet2_pt_Up_Z;
  TH1D *h_jet1_eta_Up_Z;
  TH1D *h_jet2_eta_Up_Z;
  TH1D *h_cthetaCS_Up_Z;
  TH1D *h_dRmin_mj_Up_Z;
  TH1D *h_dRmax_mj_Up_Z;
  TH1D *h_dRmin_mmj_Up_Z;
  TH1D *h_dRmax_mmj_Up_Z;
  TH1D *h_Zep_Up_Z;
  TH1D *h_M_mmjj_Up_Z;
  TH1D *h_pt_mmjj_Up_Z;
  TH1D *h_eta_mmjj_Up_Z;
  TH1D *h_phi_mmjj_Up_Z;
  TH1D *h_leadingJet_qgl_Up_Z;
  TH1D *h_subleadingJet_qgl_Up_Z;
  TH1D *h_softJet5_Up_Z;

  TH1D *h_mu1mu2dR_Down_Z;
  TH1D *h_mu1mu2dPhi_Down_Z;
  TH1D *h_mu1mu2dEta_Down_Z;
  TH1D *h_diMuon_pt_Down_Z;
  TH1D *h_diMuon_eta_Down_Z;
  TH1D *h_diMuon_phi_Down_Z;
  TH1D *h_dijet_pt_Down_Z;
  TH1D *h_dijet_eta_Down_Z;
  TH1D *h_dijet_phi_Down_Z;
  TH1D *h_Mjj_Down_Z;
  TH1D *h_dijet_dEta_Down_Z;
  TH1D *h_dijet_dPhijj_Down_Z;
  TH1D *h_jet1_pt_Down_Z;
  TH1D *h_jet2_pt_Down_Z;
  TH1D *h_jet1_eta_Down_Z;
  TH1D *h_jet2_eta_Down_Z;
  TH1D *h_cthetaCS_Down_Z;
  TH1D *h_dRmin_mj_Down_Z;
  TH1D *h_dRmax_mj_Down_Z;
  TH1D *h_dRmin_mmj_Down_Z;
  TH1D *h_dRmax_mmj_Down_Z;
  TH1D *h_Zep_Down_Z;
  TH1D *h_M_mmjj_Down_Z;
  TH1D *h_pt_mmjj_Down_Z;
  TH1D *h_eta_mmjj_Down_Z;
  TH1D *h_phi_mmjj_Down_Z;
  TH1D *h_leadingJet_qgl_Down_Z;
  TH1D *h_subleadingJet_qgl_Down_Z;
  TH1D *h_softJet5_Down_Z;

  TH1D *h_leading_bJet_pt_ttHLoose;
  TH1D *h_leading_bJet_eta_ttHLoose;
  TH1D *h_leading_bJet_phi_ttHLoose;
  TH1D *h_mu1mu2dR_ttHLoose;
  TH1D *h_mu1mu2dPhi_ttHLoose;
  TH1D *h_mu1mu2dEta_ttHLoose;
  TH1D *h_diMuon_pt_ttHLoose;
  TH1D *h_diMuon_eta_ttHLoose;
  TH1D *h_diMuon_phi_ttHLoose;

  TH1D *h_leading_bJet_pt_ttHHad;
  TH1D *h_leading_bJet_eta_ttHHad;
  TH1D *h_leading_bJet_phi_ttHHad;
  TH1D *h_mu1mu2dR_ttHHad;
  TH1D *h_mu1mu2dPhi_ttHHad;
  TH1D *h_mu1mu2dEta_ttHHad;
  TH1D *h_diMuon_pt_ttHHad;
  TH1D *h_diMuon_eta_ttHHad;
  TH1D *h_diMuon_phi_ttHHad;
  
  TH1D *h_leading_bJet_pt_ttHLep;
  TH1D *h_leading_bJet_eta_ttHLep;
  TH1D *h_leading_bJet_phi_ttHLep;
  TH1D *h_mu1mu2dR_ttHLep;
  TH1D *h_mu1mu2dPhi_ttHLep;
  TH1D *h_mu1mu2dEta_ttHLep;
  TH1D *h_diMuon_pt_ttHLep;
  TH1D *h_diMuon_eta_ttHLep;
  TH1D *h_diMuon_phi_ttHLep;

  TH1D *h_MET_pt_WHTolv;
  TH1D *h_mu1mu2dR_WHTolv;
  TH1D *h_mu1mu2dPhi_WHTolv;
  TH1D *h_mu1mu2dEta_WHTolv;
  TH1D *h_diMuon_pt_WHTolv;
  TH1D *h_diMuon_eta_WHTolv;
  TH1D *h_diMuon_phi_WHTolv;

  TH1D *h_gen_diMuon_m;
  TH1D *h_gen_extralep;
  TH1D *h_gen_dRlepH;
  TH1D *h_gen_extralep1_pt;
  TH1D *h_gen_extralep1_eta;

  TH1D *h_mu1pt;
  TH1D *h_mu2pt;
  TH1D *h_mu1eta;
  TH1D *h_mu2eta;
  TH1D *h_mu1phi;
  TH1D *h_mu2phi;
  TH1D *h_mu1mu2dR;
  TH1D *h_mu1mu2dPhi;
  TH1D *h_diMuon_pt;
  TH1D *h_diMuon_eta;
  TH1D *h_diMuon_phi;
  TH1D *h_diMuon_mass;
  TH1D *h_diMuon_mass_SR;
  TH1D *h_diMuon_mass_110To120;
  TH1D *h_diMuon_mass_130To150;
  TH1D *h_diMuon_mass_110To150;
  TH1D *h_j1pt;
  TH1D *h_j1phi;
  TH1D *h_j1eta;

  TH1D *h_j2pt;
  TH1D *h_j2phi;
  TH1D *h_j2eta;

  TH1D *h_Njet;
  TH1D *h_Nbjet;
  TH1D *h_j1j2dR;
  TH1D *h_j1j2dPhi;

  TH1D *h_dijet_pt;
  TH1D *h_dijet_eta;
  TH1D *h_dijet_phi;
  TH1D *h_Mjj;
  
  TH1D *h_MET_pt;
  TH1D *h_METphi;
  TH1D *h_MET_sumEt;
};

#endif

#ifdef HiggsMuMu_cxx

void HiggsMuMu::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;0

 

  oFile = new TFile(outFileName, "recreate");
  //oFile->mkdir("Cutflow");
  //oFile->cd("Cutflow");
  cattree =  new TTree("cattree","cattree");
  l1_index = new std::vector<int>();
  l2_index= new std::vector<int>(); 
  cattree->Branch("run", &run,"run/i");
  cattree->Branch("lumi", &lumi, "lumi/i");
  cattree->Branch("event", &event,"event/l");
  cattree->Branch("cat_index", &cat_index, "cat_index/i");
  cattree->Branch("genWeight", &genWeight,"genWeight/F");
  cattree->Branch("genWeight_Up", &genWeight_Up,"genWeight_Up/F");
  cattree->Branch("genWeight_Down", &genWeight_Down,"genWeight_Down/F");
  cattree->Branch("Higgs_mass", &Higgs_mass,"Higgs_mass/F");
  cattree->Branch("Higgs_pt", &Higgs_pt,"Higgs_pt/F");
  cattree->Branch("Higgs_eta", &Higgs_eta,"Higgs_eta/F");
  cattree->Branch("extralep_pfRelIso03", &extralep_pfRelIso03, "extralep_pfRelIso03/F");
  cattree->Branch("extralep_pt", &extralep_pt,"extralep_pt/F");
  cattree->Branch("extralep_eta", &extralep_eta,"extralep_eta/F");
  cattree->Branch("dRlepHiggs", &dRlepHiggs,"dRlepHiggs/F");
  cattree->Branch("softJet5", &softJet5,"softJet5/I");
  cattree->Branch("dRmm", &dRmm,"dRmm/F");
  cattree->Branch("dEtamm", &dEtamm,"dEtamm/F");
  cattree->Branch("dPhimm", &dPhimm,"dPhimm/F");  
  cattree->Branch("l1_index", "vector<int>", &l1_index);
  cattree->Branch("l2_index", "vector<int>", &l2_index);
  cattree->Branch("MET_phi", &MET_phi, "MET_phi/F");
  cattree->Branch("MET_pt", &MET_pt, "MET_pt/F");
  cattree->Branch("dEta_jj", &dEta_jj,"dEta_jj/F");
  cattree->Branch("M_jj", &M_jj,"M_jj/F");
  cattree->Branch("pt_jj", &pt_jj,"pt_jj/F");
  cattree->Branch("eta_jj", &eta_jj,"eta_jj/F");
  cattree->Branch("phi_jj", &phi_jj,"phi_jj/F");
  cattree->Branch("M_mmjj", &M_mmjj,"M_mmjj/F");
  cattree->Branch("pt_mmjj", &pt_mmjj,"pt_mmjj/F");
  cattree->Branch("eta_mmjj", &eta_mmjj,"eta_mmjj/F");
  cattree->Branch("phi_mmjj", &phi_mmjj,"phi_mmjj/F");
  cattree->Branch("Zep", &Zep,"Zep/F");
  cattree->Branch("dRmin_mj", &dRmin_mj,"dRmin_mj/F");
  cattree->Branch("dRmax_mj", &dRmax_mj,"dRmax_mj/F");
  cattree->Branch("dRmin_mmj", &dRmin_mmj,"dRmin_mmj/F");
  cattree->Branch("dRmax_mmj", &dRmax_mmj,"dRmax_mmj/F");
  cattree->Branch("leadingJet_pt",&leadingJet_pt,"leadingJet_pt/F");
  cattree->Branch("subleadingJet_pt",&subleadingJet_pt,"subleadingJet_pt/F");
  cattree->Branch("leadingJet_eta",&leadingJet_eta,"leadingJet_eta/F");
  cattree->Branch("subleadingJet_eta",&subleadingJet_eta,"subleadingJet_eta/F");
  cattree->Branch("leadingJet_qgl",&leadingJet_qgl,"leadingJet_qgl/F");
  cattree->Branch("subleadingJet_qgl",&subleadingJet_qgl,"subleadingJet_qgl/F");
  cattree->Branch("dPhijj", &dPhijj,"dPhijj/F");
  cattree->Branch("cthetaCS",&cthetaCS,"cthetaCS/F");
  cattree->Branch("BDT_incl", &BDT_incl,"BDT_incl/F");
  cattree->Branch("ll_mass", &ll_mass,"ll_mass/F");
  cattree->Branch("MqqLog", &MqqLog,"MqqLog/F");
  cattree->Branch("mumujj_pt", &mumujj_pt,"mumujj_pt/F");
  cattree->Branch("DeltaEtaQQ", &DeltaEtaQQ,"DeltaEtaQQ/F"); 
  cattree->Branch("softActivityEWK_njets5", &softActivityEWK_njets5,"softActivityEWK_njets5/F"); 
  cattree->Branch("ll_zstar", &ll_zstar,"ll_zstar/F");
  cattree->Branch("ll_pt", &ll_pt,"ll_pt/F");
  cattree->Branch("theta2", &theta2,"theta2/F");
  cattree->Branch("impulsoZ", &impulsoZ,"impulsoZ/F");
  cattree->Branch("maxAbsEta", &maxAbsEta,"maxAbsEta/F");
  cattree->Branch("qgl_2qAtanh", &qgl_2qAtanh,"qgl_2qAtanh/F");
  h_diMuon_mass_ggH = new TH1D("h_diMuon_mass_ggH","diMuon_mass_ggH",10,120,130);
  h_diMuon_mass_VHHad = new TH1D("h_diMuon_mass_VHHad","diMuon_mass_VHHad",10,120,130);
  h_diMuon_mass_ZHll = new TH1D("h_diMuon_mass_ZHll","diMuon_mass_ZHll",10,120,130);

  //==========================================================================//
  oFile->mkdir("Categorization_Hists");
  oFile->cd("Categorization_Hists");
  catyield = new TH1D("h_category_yield","h_category_yield",10,0,10);
  h_diMuon_mass_cat=new TH1D("diMuon_mass_cat","Mass for dimuon system",30,110.,140.);
  //==========================================================================//
  oFile->mkdir("Gen_Hists");
  oFile->cd("Gen_Hists");
  h_gen_diMuon_m = new TH1D("h_gen_higgs_mass","Higgs mass",10,110,150);
  h_gen_extralep = new TH1D("h_gen_extralep","number of extra lepton",5,0,5);
  h_gen_dRlepH = new TH1D("h_gen_dRlepH","gen_deltaR(extra lepton and Higgs)",10,0,10);
  h_gen_extralep1_pt = new TH1D("h_gen_extralep1_pt","gen_extralep1_pt",10,0,100);
  h_gen_extralep1_eta = new TH1D("h_gen_extralep1_eta","gen_extralep1_eta",10,-5,5);
  //========================================================================================//
  oFile->mkdir("WHTolv");
  oFile->cd("WHTolv");
  h_dRlepH = new TH1D("h_dRlepH","deltaR(extra lepton and Higgs)",10,0,10);
  
  h_diMuon_mass_WHlv = new TH1D("h_diMuon_mass_WHlv","diMuon_mass_WHlv",10,120,130);
  h_diMuon_mass_110To140_WHlv = new TH1D("h_diMuon_mass_110To140_WHlv","diMuon_mass_WHlv",15,110.,140.);
  h_extralep1_pt = new TH1D("h_extralep1_pt","extralep1_pt",10,0,100);
  h_extralep1_eta = new TH1D("h_extralep1_eta","extralep1_eta",10,-5,5);
  h_extralep_Electron_mvaFall17Iso =new TH1D("Electron_mvaFall17Iso","Electron_mvaFall17Iso",100,0.,1.);
  h_MET_pt_WHTolv=new TH1D("MET_pt_WHlv","MET pT",25,0.,200.);
  h_mu1mu2dR_WHTolv=new TH1D("mu1mu2dR_WHlv","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_WHTolv=new TH1D("mu1mu2dPhi_WHlv","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_WHTolv=new TH1D("mu1mu2dEta_WHlv","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_WHTolv=new TH1D("diMuon_pt_WHlv","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_WHTolv=new TH1D("diMuon_phi_WHlv","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_WHTolv=new TH1D("diMuon_eta_WHlv","#eta for dimuon system",40,-4.,4.);
  
  //========================================================================================//
  oFile->mkdir("VBF");
  oFile->cd("VBF");
  
  h_diMuon_mass_VBF = new TH1D("h_diMuon_mass_VBF","diMuon_mass_VBF",10,120,130);
  
  h_diMuon_mass_110To140_VBF = new TH1D("h_diMuon_mass_110To150_VBF","diMuon_mass_VBF",20,110.,150.);
  h_mu1mu2dR_VBF=new TH1D("mu1mu2dR_VBF","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_VBF=new TH1D("mu1mu2dPhi_VBF","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_VBF=new TH1D("mu1mu2dEta_VBF","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_VBF=new TH1D("diMuon_pt_VBF","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_VBF=new TH1D("diMuon_phi_VBF","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_VBF=new TH1D("diMuon_eta_VBF","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_VBF=new TH1D("dijet_pt_VBF","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_VBF=new TH1D("dijet_phi_VBF","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_VBF=new TH1D("dijet_eta_VBF","#eta for dijet system",40,-4.,4.);
  h_Mjj_VBF=new TH1D("diJet_mass_VBF","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_VBF=new TH1D("dijet_dEta_VBF","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_VBF=new TH1D("dijet_dPhijj_VBF","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_VBF=new TH1D("jet1_pt_VBF","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_VBF=new TH1D("jet2_pt_VBF","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_VBF=new TH1D("jet1_eta_VBF","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_VBF=new TH1D("jet2_eta_VBF","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_VBF=new TH1D("cthetaCS_VBF","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_VBF=new TH1D("dRmin_mj_VBF","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_VBF=new TH1D("dRmax_mj_VBF","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_VBF=new TH1D("dRmin_mmj_VBF","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_VBF=new TH1D("dRmax_mmj_VBF","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_VBF=new TH1D("Zep_VBF","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_VBF = new TH1D("di_muon_jet_mass_VBF","diMuon_mass_VBF",100,500.,2000.);
  h_pt_mmjj_VBF=new TH1D("di_muon_jet_pt_VBF","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_VBF=new TH1D("di_muon_jet_eta_VBF","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_VBF=new TH1D("di_muon_jet_phi_VBF","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_VBF=new TH1D("jet1_qgl_VBF","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_VBF=new TH1D("jet2_qgl_VBF","QGL for subleading jet",20,0.,1.);
  h_softJet5_VBF=new TH1D("softJet5_VBF","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);
  h_2D_ptjj_dRmm=new TH2D("2D_ptjj_dRmm_VBF","P_{T} for dijet vs #Delta R between muons",100,0.0,1000.,15,0.,4.);
  h_2D_CS_dEtamm=new TH2D("2D_CS_dEtamm_VBF","cos(#theta^{*}) vs #Delta #eta for two leading muons",20,0.,1.,40,-4.,4.);
  h_2D_Hpt_ptjj=new TH2D("2D_diMuonpt_ptjj_VBF","P_{T} for dimuon system vs P_{T} for dijet",100,0.0,1000.,100,0.0,1000.);
  h_2D_Hpt_dRmm=new TH2D("2D_diMuonpt_dRmm_VBF","P_{T} for dimuon system vs #Delta R between muons",100,0.0,1000.,15,0.,4.);
  h_2D_j1eta_dEtajj=new TH2D("2D_j1eta_dEtajj_VBF","#eta for leading jet vs #Delta #eta for two leading jets",40,-4.,4.,40,-4.,4.);
  h_2D_j2eta_dEtajj=new TH2D("2D_j2eta_dEtajj_VBF","#eta for sub-leading jet vs #Delta #eta for two leading jets",40,-4.,4.,40,-4.,4.);
  h_2D_etammjj_etajj=new TH2D("2D_etammjj_etajj_VBF","#eta for dimuon+dijet system vs #Delta #eta for two leading jets",40,-4.,4.,40,-4.,4.);
  h_2D_j1pt_ptjj=new TH2D("2D_j1pt_ptjj_VBF","P_{T} for leading jet  vs P_{T} for dijet",100,0.0,1000., 100,0.0,1000.);


  h_diMuon_mass_Up_VBF = new TH1D("h_diMuon_mass_Up_VBF","diMuon_mass_Up_VBF",10,120,130);
  
  h_diMuon_mass_110To140_Up_VBF = new TH1D("h_diMuon_mass_110To150_Up_VBF","diMuon_mass_Up_VBF",20,110.,150.);
  h_mu1mu2dR_Up_VBF=new TH1D("mu1mu2dR_Up_VBF","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_Up_VBF=new TH1D("mu1mu2dPhi_Up_VBF","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_Up_VBF=new TH1D("mu1mu2dEta_Up_VBF","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_Up_VBF=new TH1D("diMuon_pt_Up_VBF","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_Up_VBF=new TH1D("diMuon_phi_Up_VBF","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_Up_VBF=new TH1D("diMuon_eta_Up_VBF","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_Up_VBF=new TH1D("dijet_pt_Up_VBF","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_Up_VBF=new TH1D("dijet_phi_Up_VBF","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_Up_VBF=new TH1D("dijet_eta_Up_VBF","#eta for dijet system",40,-4.,4.);
  h_Mjj_Up_VBF=new TH1D("diJet_mass_Up_VBF","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_Up_VBF=new TH1D("dijet_dEta_Up_VBF","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_Up_VBF=new TH1D("dijet_dPhijj_Up_VBF","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_Up_VBF=new TH1D("jet1_pt_Up_VBF","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_Up_VBF=new TH1D("jet2_pt_Up_VBF","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_Up_VBF=new TH1D("jet1_eta_Up_VBF","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_Up_VBF=new TH1D("jet2_eta_Up_VBF","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_Up_VBF=new TH1D("cthetaCS_Up_VBF","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_Up_VBF=new TH1D("dRmin_mj_Up_VBF","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_Up_VBF=new TH1D("dRmax_mj_Up_VBF","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_Up_VBF=new TH1D("dRmin_mmj_Up_VBF","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_Up_VBF=new TH1D("dRmax_mmj_Up_VBF","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_Up_VBF=new TH1D("Zep_Up_VBF","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_Up_VBF = new TH1D("di_muon_jet_mass_Up_VBF","diMuon_mass_Up_VBF",100,500.,2000.);
  h_pt_mmjj_Up_VBF=new TH1D("di_muon_jet_pt_Up_VBF","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_Up_VBF=new TH1D("di_muon_jet_eta_Up_VBF","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_Up_VBF=new TH1D("di_muon_jet_phi_Up_VBF","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_Up_VBF=new TH1D("jet1_qgl_Up_VBF","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_Up_VBF=new TH1D("jet2_qgl_Up_VBF","QGL for subleading jet",20,0.,1.);
  h_softJet5_Up_VBF=new TH1D("softJet5_Up_VBF","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);

  h_diMuon_mass_Down_VBF = new TH1D("h_diMuon_mass_Down_VBF","diMuon_mass_Down_VBF",10,120,130);
  
  h_diMuon_mass_110To140_Down_VBF = new TH1D("h_diMuon_mass_110To150_Down_VBF","diMuon_mass_Down_VBF",20,110.,150.);
  h_mu1mu2dR_Down_VBF=new TH1D("mu1mu2dR_Down_VBF","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_Down_VBF=new TH1D("mu1mu2dPhi_Down_VBF","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_Down_VBF=new TH1D("mu1mu2dEta_Down_VBF","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_Down_VBF=new TH1D("diMuon_pt_Down_VBF","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_Down_VBF=new TH1D("diMuon_phi_Down_VBF","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_Down_VBF=new TH1D("diMuon_eta_Down_VBF","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_Down_VBF=new TH1D("dijet_pt_Down_VBF","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_Down_VBF=new TH1D("dijet_phi_Down_VBF","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_Down_VBF=new TH1D("dijet_eta_Down_VBF","#eta for dijet system",40,-4.,4.);
  h_Mjj_Down_VBF=new TH1D("diJet_mass_Down_VBF","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_Down_VBF=new TH1D("dijet_dEta_Down_VBF","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_Down_VBF=new TH1D("dijet_dPhijj_Down_VBF","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_Down_VBF=new TH1D("jet1_pt_Down_VBF","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_Down_VBF=new TH1D("jet2_pt_Down_VBF","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_Down_VBF=new TH1D("jet1_eta_Down_VBF","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_Down_VBF=new TH1D("jet2_eta_Down_VBF","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_Down_VBF=new TH1D("cthetaCS_Down_VBF","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_Down_VBF=new TH1D("dRmin_mj_Down_VBF","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_Down_VBF=new TH1D("dRmax_mj_Down_VBF","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_Down_VBF=new TH1D("dRmin_mmj_Down_VBF","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_Down_VBF=new TH1D("dRmax_mmj_Down_VBF","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_Down_VBF=new TH1D("Zep_Down_VBF","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_Down_VBF = new TH1D("di_muon_jet_mass_Down_VBF","diMuon_mass_Down_VBF",100,500.,2000.);
  h_pt_mmjj_Down_VBF=new TH1D("di_muon_jet_pt_Down_VBF","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_Down_VBF=new TH1D("di_muon_jet_eta_Down_VBF","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_Down_VBF=new TH1D("di_muon_jet_phi_Down_VBF","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_Down_VBF=new TH1D("jet1_qgl_Down_VBF","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_Down_VBF=new TH1D("jet2_qgl_Down_VBF","QGL for subleading jet",20,0.,1.);
  h_softJet5_Down_VBF=new TH1D("softJet5_Down_VBF","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);


  //========================================================================================//
  oFile->mkdir("Z");
  oFile->cd("Z");
  h_diMuon_mass_Z = new TH1D("h_diMuon_mass_Z","diMuon_mass_Z",20,76,106);
  h_mu1mu2dR_Z=new TH1D("mu1mu2dR_Z","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_Z=new TH1D("mu1mu2dPhi_Z","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_Z=new TH1D("mu1mu2dEta_Z","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_Z=new TH1D("diMuon_pt_Z","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_Z=new TH1D("diMuon_phi_Z","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_Z=new TH1D("diMuon_eta_Z","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_Z=new TH1D("dijet_pt_Z","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_Z=new TH1D("dijet_phi_Z","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_Z=new TH1D("dijet_eta_Z","#eta for dijet system",40,-4.,4.);
  h_Mjj_Z=new TH1D("diJet_mass_Z","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_Z=new TH1D("dijet_dEta_Z","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_Z=new TH1D("dijet_dPhijj_Z","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_Z=new TH1D("jet1_pt_Z","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_Z=new TH1D("jet2_pt_Z","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_Z=new TH1D("jet1_eta_Z","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_Z=new TH1D("jet2_eta_Z","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_Z=new TH1D("cthetaCS_Z","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_Z=new TH1D("dRmin_mj_Z","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_Z=new TH1D("dRmax_mj_Z","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_Z=new TH1D("dRmin_mmj_Z","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_Z=new TH1D("dRmax_mmj_Z","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_Z=new TH1D("Zep_Z","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_Z = new TH1D("di_muon_jet_mass_Z","diMuon_mass_Z",100,500.,2000.);
  h_pt_mmjj_Z=new TH1D("di_muon_jet_pt_Z","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_Z=new TH1D("di_muon_jet_eta_Z","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_Z=new TH1D("di_muon_jet_phi_Z","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_Z=new TH1D("jet1_qgl_Z","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_Z=new TH1D("jet2_qgl_Z","QGL for subleading jet",20,0.,1.);
  h_softJet5_Z=new TH1D("softJet5_Z","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);





  h_diMuon_mass_Up_Z = new TH1D("h_diMuon_mass_Up_Z","diMuon_mass_Up_Z",20,76,106);
  h_mu1mu2dR_Up_Z=new TH1D("mu1mu2dR_Up_Z","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_Up_Z=new TH1D("mu1mu2dPhi_Up_Z","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_Up_Z=new TH1D("mu1mu2dEta_Up_Z","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_Up_Z=new TH1D("diMuon_pt_Up_Z","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_Up_Z=new TH1D("diMuon_phi_Up_Z","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_Up_Z=new TH1D("diMuon_eta_Up_Z","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_Up_Z=new TH1D("dijet_pt_Up_Z","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_Up_Z=new TH1D("dijet_phi_Up_Z","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_Up_Z=new TH1D("dijet_eta_Up_Z","#eta for dijet system",40,-4.,4.);
  h_Mjj_Up_Z=new TH1D("diJet_mass_Up_Z","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_Up_Z=new TH1D("dijet_dEta_Up_Z","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_Up_Z=new TH1D("dijet_dPhijj_Up_Z","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_Up_Z=new TH1D("jet1_pt_Up_Z","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_Up_Z=new TH1D("jet2_pt_Up_Z","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_Up_Z=new TH1D("jet1_eta_Up_Z","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_Up_Z=new TH1D("jet2_eta_Up_Z","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_Up_Z=new TH1D("cthetaCS_Up_Z","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_Up_Z=new TH1D("dRmin_mj_Up_Z","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_Up_Z=new TH1D("dRmax_mj_Up_Z","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_Up_Z=new TH1D("dRmin_mmj_Up_Z","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_Up_Z=new TH1D("dRmax_mmj_Up_Z","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_Up_Z=new TH1D("Zep_Up_Z","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_Up_Z = new TH1D("di_muon_jet_mass_Up_Z","diMuon_mass_Up_Z",100,500.,2000.);
  h_pt_mmjj_Up_Z=new TH1D("di_muon_jet_pt_Up_Z","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_Up_Z=new TH1D("di_muon_jet_eta_Up_Z","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_Up_Z=new TH1D("di_muon_jet_phi_Up_Z","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_Up_Z=new TH1D("jet1_qgl_Up_Z","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_Up_Z=new TH1D("jet2_qgl_Up_Z","QGL for subleading jet",20,0.,1.);
  h_softJet5_Up_Z=new TH1D("softJet5_Up_Z","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);
  

  h_diMuon_mass_Down_Z = new TH1D("h_diMuon_mass_Down_Z","diMuon_mass_Down_Z",20,76,106);
  h_mu1mu2dR_Down_Z=new TH1D("mu1mu2dR_Down_Z","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_Down_Z=new TH1D("mu1mu2dPhi_Down_Z","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_Down_Z=new TH1D("mu1mu2dEta_Down_Z","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_Down_Z=new TH1D("diMuon_pt_Down_Z","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_Down_Z=new TH1D("diMuon_phi_Down_Z","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_Down_Z=new TH1D("diMuon_eta_Down_Z","#eta for dimuon system",40,-4.,4.);
  h_dijet_pt_Down_Z=new TH1D("dijet_pt_Down_Z","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_phi_Down_Z=new TH1D("dijet_phi_Down_Z","#phi for dijet system",32,-3.2,3.2);
  h_dijet_eta_Down_Z=new TH1D("dijet_eta_Down_Z","#eta for dijet system",40,-4.,4.);
  h_Mjj_Down_Z=new TH1D("diJet_mass_Down_Z","Mass for dijet system",100,60.,1000.);
  h_dijet_dEta_Down_Z=new TH1D("dijet_dEta_Down_Z","#Delta #eta for two leading jets",80,-8.,8.);
  h_dijet_dPhijj_Down_Z=new TH1D("dijet_dPhijj_Down_Z","#Delta #phi for dijet system",32,-3.2,3.2);
  h_jet1_pt_Down_Z=new TH1D("jet1_pt_Down_Z","P_{T} for leading jet",100,0.0,1000.);
  h_jet2_pt_Down_Z=new TH1D("jet2_pt_Down_Z","P_{T} for sub-leading jet",100,0.0,1000.);
  h_jet1_eta_Down_Z=new TH1D("jet1_eta_Down_Z","#eta for leading jet",40,-4.,4.);
  h_jet2_eta_Down_Z=new TH1D("jet2_eta_Down_Z","#eta for sub-leading jet",40,-4.,4.);
  h_cthetaCS_Down_Z=new TH1D("cthetaCS_Down_Z","cos(#theta^{*})",20,0.,1.);
  h_dRmin_mj_Down_Z=new TH1D("dRmin_mj_Down_Z","Min. #Delta R between muon and jet",15,0.,4.);
  h_dRmax_mj_Down_Z=new TH1D("dRmax_mj_Down_Z","Max. #Delta R between muon and jet",15,0.,4.);
  h_dRmin_mmj_Down_Z=new TH1D("dRmin_mmj_Down_Z","Min. #Delta R between dimuon pair and jet",15,0.,4.);
  h_dRmax_mmj_Down_Z=new TH1D("dRmax_mmj_Down_Z","Max. #Delta R between dimuon pair and jet",15,0.,4.);
  h_Zep_Down_Z=new TH1D("Zep_Down_Z","Zeppenfeld variable",50,-5.,5.);
  h_M_mmjj_Down_Z = new TH1D("di_muon_jet_mass_Down_Z","diMuon_mass_Down_Z",100,500.,2000.);
  h_pt_mmjj_Down_Z=new TH1D("di_muon_jet_pt_Down_Z","P_{T} for dimuon+dijet system",100,0.0,1000.);
  h_eta_mmjj_Down_Z=new TH1D("di_muon_jet_eta_Down_Z","#eta for dimuon+dijet system",40,-4.,4.);
  h_phi_mmjj_Down_Z=new TH1D("di_muon_jet_phi_Down_Z","#phi for dimuon+dijet system",32,-3.2,3.2);
  h_leadingJet_qgl_Down_Z=new TH1D("jet1_qgl_Down_Z","QGL for leading jet",20,0.,1.);
  h_subleadingJet_qgl_Down_Z=new TH1D("jet2_qgl_Down_Z","QGL for subleading jet",20,0.,1.);
  h_softJet5_Down_Z=new TH1D("softJet5_Down_Z","Number of soft EWK jets with pt > 5 GeV",20,0.,20.);
  
  //=============================================================================================//
  oFile->mkdir("ttHHad");
  oFile->cd("ttHHad");
  
  h_diMuon_mass_ttHHad = new TH1D("h_diMuon_mass_ttHHad","diMuon_mass_ttHHad",10,120,130);
  h_diMuon_mass_110To140_ttHHad = new TH1D("h_diMuon_mass_110To140_ttHHad","diMuon_mass_ttHHad",15,110.,140.);
  h_mu1mu2dR_ttHHad=new TH1D("mu1mu2dR_ttHHad","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_ttHHad=new TH1D("mu1mu2dPhi_ttHHad","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_ttHHad=new TH1D("mu1mu2dEta_ttHHad","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_ttHHad=new TH1D("diMuon_pt_ttHHad","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_ttHHad=new TH1D("diMuon_phi_ttHHad","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_ttHHad=new TH1D("diMuon_eta_ttHHad","#eta for dimuon system",40,-4.,4.);

  h_leading_bJet_pt_ttHHad=new TH1D("leading_bJet_pt_ttHHad","P_{T} for leading b-jet system",100,0.0,1000.);
  h_leading_bJet_phi_ttHHad=new TH1D("leading_bJet_phi_ttHHad","#phi for leading b-jet system",32,-3.2,3.2);
  h_leading_bJet_eta_ttHHad=new TH1D("leading_bJet_eta_ttHHad","#eta for leading b-jet system",40,-4.,4.);

  //========================================================================================//
  oFile->mkdir("ttHLep");
  oFile->cd("ttHLep");
  
  h_diMuon_mass_ttHLep = new TH1D("h_diMuon_mass_ttHLep","diMuon_mass_ttHLep",10,120,130);
  h_diMuon_mass_110To140_ttHLep = new TH1D("h_diMuon_mass_110To140_TTHLep","diMuon_mass_ttHLep",15,110.,140.);
  h_mu1mu2dR_ttHLep=new TH1D("mu1mu2dR_ttHLep","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_ttHLep=new TH1D("mu1mu2dPhi_ttHLep","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_ttHLep=new TH1D("mu1mu2dEta_ttHLep","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_ttHLep=new TH1D("diMuon_pt_ttHLep","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_ttHLep=new TH1D("diMuon_phi_ttHLep","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_ttHLep=new TH1D("diMuon_eta_ttHLep","#eta for dimuon system",40,-4.,4.);

  h_leading_bJet_pt_ttHLep=new TH1D("leading_bJet_pt_ttHLep","P_{T} for leading b-jet system",100,0.0,1000.);
  h_leading_bJet_phi_ttHLep=new TH1D("leading_bJet_phi_ttHLep","#phi for leading b-jet system",32,-3.2,3.2);
  h_leading_bJet_eta_ttHLep=new TH1D("leading_bJet_eta_ttHLep","#eta for leading b-jet system",40,-4.,4.);

  //========================================================================================//
  oFile->mkdir("ttHLoose");
  oFile->cd("ttHLoose");
  
  h_diMuon_mass_ttHLoose = new TH1D("h_diMuon_mass_ttHLoose","diMuon_mass_ttHLoose",10,120,130);
  h_diMuon_mass_110To140_ttHLoose = new TH1D("h_diMuon_mass_110To140_ttHLoose","diMuon_mass_ttHLoose",15,110.,140.);
  h_mu1mu2dR_ttHLoose=new TH1D("mu1mu2dR_ttHLoose","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi_ttHLoose=new TH1D("mu1mu2dPhi_ttHLoose","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_mu1mu2dEta_ttHLoose=new TH1D("mu1mu2dEta_ttHLoose","#Delta #eta for two leading muons",80,-8.,8.);
  h_diMuon_pt_ttHLoose=new TH1D("diMuon_pt_ttHLoose","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi_ttHLoose=new TH1D("diMuon_phi_ttHLoose","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta_ttHLoose=new TH1D("diMuon_eta_ttHLoose","#eta for dimuon system",40,-4.,4.);

  h_leading_bJet_pt_ttHLoose=new TH1D("leading_bJet_pt_ttHLoose","P_{T} for leading b-jet system",100,0.0,1000.);
  h_leading_bJet_phi_ttHLoose=new TH1D("leading_bJet_phi_ttHLoose","#phi for leading b-jet system",32,-3.2,3.2);
  h_leading_bJet_eta_ttHLoose=new TH1D("leading_bJet_eta_ttHLoose","#eta for leading b-jet system",40,-4.,4.);

  //==================================================================================//
  oFile->mkdir("Inclusive");
  oFile->cd("Inclusive");
  h_mu1pt=new TH1D("mu1pt","P_{T} for leading muon",100,0.0,1000.);
  h_mu1phi=new TH1D("mu1phi","#phi for leading muon",32,-3.2,3.2);
  h_mu1eta=new TH1D("mu1eta","#eta for leading muon",40,-4.,4.);

  h_mu2pt=new TH1D("mu2pt","P_{T} for sub-leading muon",100,0.0,1000.);
  h_mu2phi=new TH1D("mu2phi","#phi for sub-leading muon",32,-3.2,3.2);
  h_mu2eta=new TH1D("mu2eta","#eta for sub-leading muon",40,-4.,4.);
  h_mu1mu2dR=new TH1D("mu1mu2dR","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi=new TH1D("mu1mu2dPhi","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_diMuon_pt=new TH1D("diMuon_pt","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_phi=new TH1D("diMuon_phi","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_eta=new TH1D("diMuon_eta","#eta for dimuon system",40,-4.,4.);
  h_diMuon_mass=new TH1D("diMuon_mass","Mass for dimuon system",70,60.,200.);
  h_diMuon_mass_SR=new TH1D("diMuon_mass_SR","Mass for dimuon system",70,60.,200.);
  h_diMuon_mass_110To120=new TH1D("diMuon_mass_110To120","Mass for dimuon system",5,110.,120.);
  h_diMuon_mass_130To150=new TH1D("diMuon_mass_130To150","Mass for dimuon system",10,130.,150.);
  h_diMuon_mass_110To150=new TH1D("diMuon_mass_110To150","Mass for dimuon system",20,110.,150.);

  h_j1pt=new TH1D("j1pt","P_{T} for leading jet",100,0.0,1000.);
  h_j1phi=new TH1D("j1phi","#phi for leading jet",32,-3.2,3.2);
  h_j1eta=new TH1D("j1eta","#eta for leading jet",40,-4.,4.);

  h_j2pt=new TH1D("j2pt","P_{T} for sub-leading jet",100,0.0,1000.);
  h_j2phi=new TH1D("j2phi","#phi for sub-leading jet",32,-3.2,3.2);
  h_j2eta=new TH1D("j2eta","#eta for sub-leading jet",40,-4.,4.);

  h_Njet=new TH1D("Njet","Number of jets",11,0.0,11.);
  h_Nbjet=new TH1D("Nbjet","Number of b-jets",11,0.0,11.);
  h_j1j2dR=new TH1D("j1j2dR","#Delta R between two leading jets",15,0.,4.);
  h_j1j2dPhi=new TH1D("j1j2dPhi","#Delta #phi for two leading jets",32,-3.2,3.2);

  h_dijet_pt=new TH1D("dijet_pt","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_eta=new TH1D("dijet_phi","#phi for dijet system",32,-3.2,3.2);
  h_dijet_phi=new TH1D("dijet_eta","#eta for dijet system",40,-4.,4.);
  h_Mjj=new TH1D("diJet_mass","Mass for dijet system",100,60.,1000.);
  
  h_MET_pt=new TH1D("MET_pt","MET P_{T}",100,0.0,1000.);
  h_METphi=new TH1D("MET_phi","MET #phi",32,-3.2,3.2);
  h_MET_sumEt=new TH1D("MET_sumEt","MET Sum E_{T}",100,0.,1000.);
}

void HiggsMuMu::clearTreeVectors(){
   run = -9999;
   lumi = -9999;
   event = -9999;
   genWeight = -9999.;
   genWeight_Up = -9999.;
   genWeight_Down = -9999.;
   cat_index = -9999.;
   Higgs_mass = -9999.;
   Higgs_pt = -9999.;
   Higgs_eta = -9999.;
   extralep_pfRelIso03 = -9999.;
   extralep_pt = -9999.;
   extralep_eta = -9999.;
   dRlepHiggs = -9999.;
   dEtamm = -9999.;
   dRmm = -9999.;
   dPhimm = -9999.;
   MET_phi = -9999.;
   MET_pt = -9999.;
   dEta_jj = -9999.;
   M_jj= -9999.;
   pt_jj=-9999;
   eta_jj=-9999;
   phi_jj= -9999;
   M_mmjj= -9999.;
   pt_mmjj=-9999;
   eta_mmjj=-9999;
   phi_mmjj= -9999;
   Zep = -9999.;
   dRmin_mj = -9999.;
   dRmax_mj = -9999.;
   dRmin_mmj = -9999.;
   dRmax_mmj= -9999.;
   leadingJet_pt = -9999.;
   subleadingJet_pt = -9999.;
   leadingJet_eta = -9999.;
   leadingJet_qgl = -9999.;
   subleadingJet_qgl = -9999.;
   subleadingJet_eta = -9999.;
   dPhijj=-9999.;
   cthetaCS = -9999.;
   BDT_incl = -9999.;
   l1_index->clear();
   l2_index->clear();
   //BDT vars
   ll_mass =0.;
   MqqLog=0.;
   mumujj_pt=0.;
   DeltaEtaQQ=0.;
   softActivityEWK_njets5=0.;
   ll_zstar=0.;
   ll_pt=0.;
   theta2=0.;
   impulsoZ=0.;
   maxAbsEta=0.;
   qgl_2qAtanh=0.;
}

HiggsMuMu::HiggsMuMu(const TString &inputFileList, const char *outFileName, const char* dataset, const char *isData)
{
TChain *tree = new TChain("tree");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  NtupleVariables::Init(tree);
  BookHistogram(outFileName);
  }
Bool_t HiggsMuMu::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

HiggsMuMu::~HiggsMuMu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   oFile->Write();
   oFile->Close();
   
}

Long64_t HiggsMuMu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      // Notify();
   }
   return centry;
}
#endif // #ifdef HiggsMuMu_cxx
