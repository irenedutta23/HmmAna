//////////////////////////////////////////////////////////
// 
// Original Author : Irene Dutta
//                   Caltech
// Date Created    : Mon 27 Aug, 2018
//////////////////////////////////////////////////////////
#ifndef NtupleVariables_h
#define NtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TSelector.h>
#include "TH1F.h"
#include "TH2F.h"

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif

using namespace std;
class NtupleVariables : public TSelector {
public :

   NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   void    Init(TTree *tree);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   double  DeltaPhi(double, double);
   double  DeltaR(double eta1, double phi1, double eta2, double phi2);
   double  HT (std::vector<TLorentzVector>);
   std::pair<double,double> helicityAngles( TLorentzVector v1, TLorentzVector v2 );
   std::pair<double,double> CSAngles( TLorentzVector v1, TLorentzVector v2, int charge);
   TLorentzVector MHT(std::vector<TLorentzVector>);
   
   // void LeptHist(TH1D *h[20],const char *leptype, int Zcut);
   //ClassDef(NtupleVariables,0);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
   UInt_t          t_run;
   UInt_t          t_luminosityBlock;
   ULong64_t       t_event;
   Float_t         t_genWeight;
   Float_t         t_puWeight;
   Float_t         t_puWeightUp;
   Float_t         t_puWeightDown;
   Float_t         t_PrefireWeight;
   Float_t         t_PrefireWeight_Up;
   Float_t         t_PrefireWeight_Down;
   UInt_t          t_mu1;
   UInt_t          t_mu2;
   UInt_t          t_index_trigm_mu;
   vector<int>     *t_El_charge;
   vector<float>   *t_El_pt;
   vector<float>   *t_El_phi;
   vector<float>   *t_El_eta;
   vector<float>   *t_El_mass;
   vector<int>     *t_El_cutBased;
   vector<int>     *t_El_tightCharge;
   vector<bool>    *t_El_cutBased_HEEP;
   vector<bool>    *t_El_isPFcand;
   vector<float>   *t_El_pfRelIso03_all;
   vector<float>   *t_El_pfRelIso03_chg;
   vector<float>   *t_El_miniPFRelIso03_all;
   vector<float>   *t_El_dxy;
   vector<float>   *t_El_dxyErr;
   vector<float>   *t_El_dz;
   vector<float>   *t_El_dzErr;
   vector<float>   *t_El_sip3d;
   vector<float>   *t_Electron_mvaFall17Iso;
   vector<bool>    *t_Electron_mvaFall17Iso_WP80;
   vector<bool>    *t_Electron_mvaFall17Iso_WP90;
   vector<bool>    *t_Electron_mvaFall17Iso_WPL;
   vector<bool>    *t_Electron_mvaFall17noIso_WP80;
   vector<bool>    *t_Electron_mvaFall17noIso_WP90;
   vector<bool>    *t_Electron_mvaFall17noIso_WPL;
   vector<int>     *t_Mu_charge;
   vector<float>   *t_Mu_EffSF_TRIG;
   vector<float>   *t_Mu_EffSFErr_TRIG;
   vector<float>   *t_Mu_EffSF_ID;
   //vector<float>   *t_Mu_EffSFErr_ID;
   vector<float>   *t_Mu_EffSF_ISO;
   //vector<float>   *t_Mu_EffSFErr_ISO; 
   vector<float>   *t_Mu_pt;
   vector<float>   *t_Mu_ptErr;
   vector<float>   *t_Mu_phi;
   vector<float>   *t_Mu_eta;
   vector<float>   *t_Mu_mass;
   vector<float>   *t_Mu_dxy;
   vector<float>   *t_Mu_dxyErr;
   vector<float>   *t_Mu_dz;
   vector<float>   *t_Mu_dzErr;
   vector<float>   *t_Mu_sip3d;
   vector<float>   *t_Mu_pfRelIso03_all;
   vector<float>   *t_Mu_pfRelIso03_chg;
   vector<float>   *t_Mu_pfRelIso04_all;
   vector<float>   *t_Mu_miniPFRelIso_all;
   vector<int>     *t_Mu_tightCharge;
   vector<bool>    *t_Mu_isPFcand;
   vector<bool>    *t_Mu_isglobal;
   vector<bool>    *t_Mu_istracker;
   vector<bool>    *t_Mu_mediumId;
   vector<bool>    *t_Mu_softId;
   vector<bool>    *t_Mu_tightId;
   vector<int>     *t_Mu_nStations;
   vector<int>     *t_Mu_nTrackerLayers;
   Float_t         t_diMuon_pt;
   Float_t         t_diMuon_eta;
   Float_t         t_diMuon_phi;
   Float_t         t_diMuon_mass;
   Float_t         t_diJet_mass_mo;
   Float_t         t_MET_phi;
   Float_t         t_MET_pt;
   Float_t         t_MET_sumEt;
   vector<float>   *t_FatJet_area;
   vector<float>   *t_FatJet_btagCMVA;
   vector<float>   *t_FatJet_btagCSVV2;
   vector<float>   *t_FatJet_btagDeepB;
   vector<float>   *t_FatJet_eta;
   vector<float>   *t_FatJet_mass;
   vector<float>   *t_FatJet_msoftdrop;
   vector<float>   *t_FatJet_n2b1;
   vector<float>   *t_FatJet_n3b1;
   vector<float>   *t_FatJet_phi;
   vector<float>   *t_FatJet_pt;
   vector<float>   *t_FatJet_tau1;
   vector<float>   *t_FatJet_tau2;
   vector<float>   *t_FatJet_tau3;
   vector<float>   *t_FatJet_tau4;
   vector<int>     *t_FatJet_jetId;
   vector<int>     *t_FatJet_subJetIdx1;
   vector<int>     *t_FatJet_subJetIdx2;
   vector<float>   *t_SubJet_btagCMVA;
   vector<float>   *t_SubJet_btagCSVV2;
   vector<float>   *t_SubJet_btagDeepB;
   vector<float>   *t_SubJet_eta;
   vector<float>   *t_SubJet_mass;
   vector<float>   *t_SubJet_n2b1;
   vector<float>   *t_SubJet_n3b1;
   vector<float>   *t_SubJet_phi;
   vector<float>   *t_SubJet_pt;
   vector<float>   *t_SubJet_tau1;
   vector<float>   *t_SubJet_tau2;
   vector<float>   *t_SubJet_tau3;
   vector<float>   *t_SubJet_tau4;
   Int_t           t_nJet;
   vector<float>   *t_Jet_area;
   vector<float>   *t_Jet_btagCMVA;
   vector<float>   *t_Jet_btagCSVV2;
   vector<float>   *t_Jet_btagDeepB;
   vector<float>   *t_Jet_btagDeepC;
   vector<float>   *t_Jet_btagDeepFlavB;
   vector<float>   *t_Jet_chEmEF;
   vector<float>   *t_Jet_chHEF;
   vector<float>   *t_Jet_eta;
   vector<float>   *t_Jet_mass;
   vector<float>   *t_Jet_neEmEF;
   vector<float>   *t_Jet_neHEF;
   vector<float>   *t_Jet_phi;
   vector<float>   *t_Jet_pt;
   vector<float>   *t_Jet_qgl;
   vector<int>     *t_Jet_jetId;
   vector<int>     *t_Jet_nConstituents;
   vector<int>     *t_Jet_nElectrons;
   vector<int>     *t_Jet_nMuons;
   vector<int>     *t_Jet_puId;
   Float_t         t_diJet_pt;
   Float_t         t_diJet_eta;
   Float_t         t_diJet_phi;
   Float_t         t_diJet_mass;
   Int_t           t_nbJet;
   vector<float>   *t_bJet_area;
   vector<float>   *t_bJet_btagCMVA;
   vector<float>   *t_bJet_btagCSVV2;
   vector<float>   *t_bJet_btagDeepB;
   vector<float>   *t_bJet_btagDeepC;
   vector<float>   *t_bJet_btagDeepFlavB;
   vector<float>   *t_bJet_chEmEF;
   vector<float>   *t_bJet_chHEF;
   vector<float>   *t_bJet_eta;
   vector<float>   *t_bJet_mass;
   vector<float>   *t_bJet_neEmEF;
   vector<float>   *t_bJet_neHEF;
   vector<float>   *t_bJet_phi;
   vector<float>   *t_bJet_pt;
   vector<float>   *t_bJet_qgl;
   vector<int>     *t_bJet_jetId;
   vector<int>     *t_bJet_nConstituents;
   vector<int>     *t_bJet_nElectrons;
   vector<int>     *t_bJet_nMuons;
   vector<int>     *t_bJet_puId;
   vector<double>   *t_bJet_SF;
   vector<double>   *t_bJet_SFup;
   vector<double>   *t_bJet_SFdown;

   Float_t         t_PV_ndof;
   Float_t         t_PV_x;
   Float_t         t_PV_y;
   Float_t         t_PV_z;
   Int_t           t_PV_npvs;
   Int_t           t_PV_npvsGood;
   vector<float>   *t_GenPart_eta;
   vector<float>   *t_GenPart_mass;
   vector<float>   *t_GenPart_phi;
   vector<float>   *t_GenPart_pt;
   vector<int>     *t_GenPart_genPartIdxMother;
   vector<int>     *t_GenPart_pdgId;
   vector<int>     *t_GenPart_status;
   Int_t           t_SoftActivityJetNjets5;
   vector<float>   *t_GenJet_eta;
   vector<float>   *t_GenJet_mass;
   vector<float>   *t_GenJet_phi;
   vector<float>   *t_GenJet_pt;
    
   // List of branches
   TBranch        *b_t_run;   //!
   TBranch        *b_t_luminosityBlock;   //!
   TBranch        *b_t_event;   //!
   TBranch        *b_t_genWeight;   //!
   TBranch        *b_t_puWeight;
   TBranch        *b_t_puWeightUp;
   TBranch        *b_t_puWeightDown;
   TBranch        *b_t_PrefireWeight;
   TBranch        *b_t_PrefireWeight_Up;
   TBranch        *b_t_PrefireWeight_Down;
   TBranch        *b_t_mu1;   //!
   TBranch        *b_t_mu2;   //!
   TBranch        *b_t_index_trigm_mu;
   TBranch        *b_t_El_charge;   //!
   TBranch        *b_t_El_pt;   //!
   TBranch        *b_t_El_phi;   //!
   TBranch        *b_t_El_eta;   //!
   TBranch        *b_t_El_mass;   //!
   TBranch        *b_t_El_cutBased;   //!
   TBranch        *b_t_El_tightCharge;   //!
   TBranch        *b_t_El_cutBased_HEEP;   //!
   TBranch        *b_t_El_isPFcand;   //!
   TBranch        *b_t_El_pfRelIso03_all;   //!
   TBranch        *b_t_El_pfRelIso03_chg;   //!
   TBranch        *b_t_El_miniPFRelIso03_all; //!
   TBranch        *b_t_El_dxy;   //!
   TBranch        *b_t_El_dxyErr;   //!
   TBranch        *b_t_El_dz;   //!
   TBranch        *b_t_El_dzErr;   //!
   TBranch        *b_t_El_sip3d;   //!
   TBranch        *b_t_Electron_mvaFall17Iso;
   TBranch        *b_t_Electron_mvaFall17Iso_WP80;   //!
   TBranch        *b_t_Electron_mvaFall17Iso_WP90;   //!
   TBranch        *b_t_Electron_mvaFall17Iso_WPL;   //!
   TBranch        *b_t_Electron_mvaFall17noIso_WP80;   //!
   TBranch        *b_t_Electron_mvaFall17noIso_WP90;   //!
   TBranch        *b_t_Electron_mvaFall17noIso_WPL;   //!
   TBranch        *b_t_Mu_charge;   //!
   TBranch        *b_t_Mu_EffSF_TRIG;
   TBranch        *b_t_Mu_EffSFErr_TRIG;
   TBranch        *b_t_Mu_EffSF_ID;
   //TBranch        *b_t_Mu_EffSFErr_ID;
   TBranch        *b_t_Mu_EffSF_ISO;
   //TBranch        *b_t_Mu_EffSFErr_ISO;
   TBranch        *b_t_Mu_pt;   //!
   TBranch        *b_t_Mu_ptErr;   //!
   TBranch        *b_t_Mu_phi;   //!
   TBranch        *b_t_Mu_eta;   //!
   TBranch        *b_t_Mu_mass;   //!
   TBranch        *b_t_Mu_dxy;   //!
   TBranch        *b_t_Mu_dxyErr;   //!
   TBranch        *b_t_Mu_dz;   //!
   TBranch        *b_t_Mu_dzErr;   //!
   TBranch        *b_t_Mu_sip3d;   //!
   TBranch        *b_t_Mu_pfRelIso03_all;   //!
   TBranch        *b_t_Mu_pfRelIso03_chg;   //!
   TBranch        *b_t_Mu_pfRelIso04_all;   //!
   TBranch        *b_t_Mu_miniPFRelIso_all; //!
   TBranch        *b_t_Mu_tightCharge;   //!
   TBranch        *b_t_Mu_isPFcand;   //!
   TBranch        *b_t_Mu_isglobal;   //!
   TBranch        *b_t_Mu_istracker;   //!
   TBranch        *b_t_Mu_mediumId;   //!
   TBranch        *b_t_Mu_softId;   //!
   TBranch        *b_t_Mu_tightId;   //!
   TBranch        *b_t_Mu_nStations;   //!
   TBranch        *b_t_Mu_nTrackerLayers;   //!
   TBranch        *b_t_diMuon_pt;   //!
   TBranch        *b_t_diMuon_eta;   //!
   TBranch        *b_t_diMuon_phi;   //!
   TBranch        *b_t_diMuon_mass;   //!
   TBranch        *b_t_diJet_mass_mo;   //!
   TBranch        *b_t_MET_phi;   //!
   TBranch        *b_t_MET_pt;   //!
   TBranch        *b_t_MET_sumEt;   //!
   TBranch        *b_t_FatJet_area;   //!
   TBranch        *b_t_FatJet_btagCMVA;   //!
   TBranch        *b_t_FatJet_btagCSVV2;   //!
   TBranch        *b_t_FatJet_btagDeepB;   //!
   TBranch        *b_t_FatJet_eta;   //!
   TBranch        *b_t_FatJet_mass;   //!
   TBranch        *b_t_FatJet_msoftdrop;   //!
   TBranch        *b_t_FatJet_n2b1;   //!
   TBranch        *b_t_FatJet_n3b1;   //!
   TBranch        *b_t_FatJet_phi;   //!
   TBranch        *b_t_FatJet_pt;   //!
   TBranch        *b_t_FatJet_tau1;   //!
   TBranch        *b_t_FatJet_tau2;   //!
   TBranch        *b_t_FatJet_tau3;   //!
   TBranch        *b_t_FatJet_tau4;   //!
   TBranch        *b_t_FatJet_jetId;   //!
   TBranch        *b_t_FatJet_subJetIdx1;   //!
   TBranch        *b_t_FatJet_subJetIdx2;   //!
   TBranch        *b_t_SubJet_btagCMVA;   //!
   TBranch        *b_t_SubJet_btagCSVV2;   //!
   TBranch        *b_t_SubJet_btagDeepB;   //!
   TBranch        *b_t_SubJet_eta;   //!
   TBranch        *b_t_SubJet_mass;   //!
   TBranch        *b_t_SubJet_n2b1;   //!
   TBranch        *b_t_SubJet_n3b1;   //!
   TBranch        *b_t_SubJet_phi;   //!
   TBranch        *b_t_SubJet_pt;   //!
   TBranch        *b_t_SubJet_tau1;   //!
   TBranch        *b_t_SubJet_tau2;   //!
   TBranch        *b_t_SubJet_tau3;   //!
   TBranch        *b_t_SubJet_tau4;   //!
   TBranch        *b_t_nJet;   //!
   TBranch        *b_t_Jet_area;   //!
   TBranch        *b_t_Jet_btagCMVA;   //!
   TBranch        *b_t_Jet_btagCSVV2;   //!
   TBranch        *b_t_Jet_btagDeepB;   //!
   TBranch        *b_t_Jet_btagDeepC;   //!
   TBranch        *b_t_Jet_btagDeepFlavB;   //!
   TBranch        *b_t_Jet_chEmEF;   //!
   TBranch        *b_t_Jet_chHEF;   //!
   TBranch        *b_t_Jet_eta;   //!
   TBranch        *b_t_Jet_mass;   //!
   TBranch        *b_t_Jet_neEmEF;   //!
   TBranch        *b_t_Jet_neHEF;   //!
   TBranch        *b_t_Jet_phi;   //!
   TBranch        *b_t_Jet_pt;   //!
   TBranch        *b_t_Jet_qgl;   //!
   TBranch        *b_t_Jet_jetId;   //!
   TBranch        *b_t_Jet_nConstituents;   //!
   TBranch        *b_t_Jet_nElectrons;   //!
   TBranch        *b_t_Jet_nMuons;   //!
   TBranch        *b_t_Jet_puId;   //!
   TBranch        *b_t_diJet_pt;   //!
   TBranch        *b_t_diJet_eta;   //!
   TBranch        *b_t_diJet_phi;   //!
   TBranch        *b_t_diJet_mass;   //!
   TBranch        *b_t_nbJet;   //!
   TBranch        *b_t_bJet_area;   //!
   TBranch        *b_t_bJet_btagCMVA;   //!
   TBranch        *b_t_bJet_btagCSVV2;   //!
   TBranch        *b_t_bJet_btagDeepB;   //!
   TBranch        *b_t_bJet_btagDeepC;   //!
   TBranch        *b_t_bJet_btagDeepFlavB;   //!
   TBranch        *b_t_bJet_chEmEF;   //!
   TBranch        *b_t_bJet_chHEF;   //!
   TBranch        *b_t_bJet_eta;   //!
   TBranch        *b_t_bJet_mass;   //!
   TBranch        *b_t_bJet_neEmEF;   //!
   TBranch        *b_t_bJet_neHEF;   //!
   TBranch        *b_t_bJet_phi;   //!
   TBranch        *b_t_bJet_pt;   //!
   TBranch        *b_t_bJet_qgl;   //!
   TBranch        *b_t_bJet_jetId;   //!
   TBranch        *b_t_bJet_nConstituents;   //!
   TBranch        *b_t_bJet_nElectrons;   //!
   TBranch        *b_t_bJet_nMuons;   //!
   TBranch        *b_t_bJet_puId;   //!
   TBranch        *b_t_bJet_SF;
   TBranch        *b_t_bJet_SFup;
   TBranch        *b_t_bJet_SFdown;
   TBranch        *b_t_PV_ndof;   //!
   TBranch        *b_t_PV_x;   //!
   TBranch        *b_t_PV_y;   //!
   TBranch        *b_t_PV_z;   //!
   TBranch        *b_t_PV_npvs;   //!
   TBranch        *b_t_PV_npvsGood;   //!
   TBranch        *b_t_GenPart_eta;   //!
   TBranch        *b_t_GenPart_mass;   //!
   TBranch        *b_t_GenPart_phi;   //!
   TBranch        *b_t_GenPart_pt;   //!
   TBranch        *b_t_GenPart_genPartIdxMother;   //!
   TBranch        *b_t_GenPart_pdgId;   //!
   TBranch        *b_t_GenPart_status;   //!
   TBranch        *b_t_SoftActivityJetNjets5;
   TBranch        *b_t_GenJet_eta;   //!
   TBranch        *b_t_GenJet_mass;   //!
   TBranch        *b_t_GenJet_phi;   //!
   TBranch        *b_t_GenJet_pt;   //!

};
#endif

#ifdef NtupleVariables_cxx
void NtupleVariables::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  // Set object pointer
  

 // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // Set object pointer
   t_El_charge = 0;
   t_El_pt = 0;
   t_El_phi = 0;
   t_El_eta = 0;
   t_El_mass = 0;
   t_El_cutBased = 0;
   t_El_tightCharge = 0;
   t_El_cutBased_HEEP = 0;
   t_El_isPFcand = 0;
   t_El_pfRelIso03_all = 0;
   t_El_pfRelIso03_chg = 0;
   t_El_miniPFRelIso03_all = 0;
   t_El_dxy = 0;
   t_El_dxyErr = 0;
   t_El_dz = 0;
   t_El_dzErr = 0;
   t_El_sip3d = 0;
   t_Electron_mvaFall17Iso=0;
   t_Electron_mvaFall17Iso_WP80 = 0;
   t_Electron_mvaFall17Iso_WP90 = 0;
   t_Electron_mvaFall17Iso_WPL = 0;
   t_Electron_mvaFall17noIso_WP80 = 0;
   t_Electron_mvaFall17noIso_WP90 = 0;
   t_Electron_mvaFall17noIso_WPL = 0;
   t_Mu_charge = 0;
   t_Mu_EffSF_TRIG=0;
   t_Mu_EffSFErr_TRIG=0;
   t_Mu_EffSF_ID=0;
   //t_Mu_EffSFErr_ID=0;
   t_Mu_EffSF_ISO=0;
   //t_Mu_EffSFErr_ISO=0;
   t_Mu_pt = 0;
   t_Mu_ptErr = 0;
   t_Mu_phi = 0;
   t_Mu_eta = 0;
   t_Mu_mass = 0;
   t_Mu_dxy = 0;
   t_Mu_dxyErr = 0;
   t_Mu_dz = 0;
   t_Mu_dzErr = 0;
   t_Mu_sip3d = 0;
   t_Mu_pfRelIso03_all = 0;
   t_Mu_pfRelIso03_chg = 0;
   t_Mu_pfRelIso04_all = 0;
   t_Mu_miniPFRelIso_all = 0;
   t_Mu_tightCharge = 0;
   t_Mu_isPFcand = 0;
   t_Mu_isglobal = 0;
   t_Mu_istracker = 0;
   t_Mu_mediumId = 0;
   t_Mu_softId = 0;
   t_Mu_tightId = 0;
   t_Mu_nStations = 0;
   t_Mu_nTrackerLayers = 0;
   t_FatJet_area = 0;
   t_FatJet_btagCMVA = 0;
   t_FatJet_btagCSVV2 = 0;
   t_FatJet_btagDeepB = 0;
   t_FatJet_eta = 0;
   t_FatJet_mass = 0;
   t_FatJet_msoftdrop = 0;
   t_FatJet_n2b1 = 0;
   t_FatJet_n3b1 = 0;
   t_FatJet_phi = 0;
   t_FatJet_pt = 0;
   t_FatJet_tau1 = 0;
   t_FatJet_tau2 = 0;
   t_FatJet_tau3 = 0;
   t_FatJet_tau4 = 0;
   t_FatJet_jetId = 0;
   t_FatJet_subJetIdx1 = 0;
   t_FatJet_subJetIdx2 = 0;
   t_SubJet_btagCMVA = 0;
   t_SubJet_btagCSVV2 = 0;
   t_SubJet_btagDeepB = 0;
   t_SubJet_eta = 0;
   t_SubJet_mass = 0;
   t_SubJet_n2b1 = 0;
   t_SubJet_n3b1 = 0;
   t_SubJet_phi = 0;
   t_SubJet_pt = 0;
   t_SubJet_tau1 = 0;
   t_SubJet_tau2 = 0;
   t_SubJet_tau3 = 0;
   t_SubJet_tau4 = 0;
   t_Jet_area = 0;
   t_Jet_btagCMVA = 0;
   t_Jet_btagCSVV2 = 0;
   t_Jet_btagDeepB = 0;
   t_Jet_btagDeepC = 0;
   t_Jet_btagDeepFlavB = 0;
   t_Jet_chEmEF = 0;
   t_Jet_chHEF = 0;
   t_Jet_eta = 0;
   t_Jet_mass = 0;
   t_Jet_neEmEF = 0;
   t_Jet_neHEF = 0;
   t_Jet_phi = 0;
   t_Jet_pt = 0;
   t_Jet_qgl = 0;
   t_Jet_jetId = 0;
   t_Jet_nConstituents = 0;
   t_Jet_nElectrons = 0;
   t_Jet_nMuons = 0;
   t_Jet_puId = 0;
   t_bJet_area = 0;
   t_bJet_btagCMVA = 0;
   t_bJet_btagCSVV2 = 0;
   t_bJet_btagDeepB = 0;
   t_bJet_btagDeepC = 0;
   t_bJet_btagDeepFlavB = 0;
   t_bJet_chEmEF = 0;
   t_bJet_chHEF = 0;
   t_bJet_eta = 0;
   t_bJet_mass = 0;
   t_bJet_neEmEF = 0;
   t_bJet_neHEF = 0;
   t_bJet_phi = 0;
   t_bJet_pt = 0;
   t_bJet_qgl = 0;
   t_bJet_jetId = 0;
   t_bJet_nConstituents = 0;
   t_bJet_nElectrons = 0;
   t_bJet_nMuons = 0;
   t_bJet_puId = 0;
   t_bJet_SF = 0;
   t_bJet_SFup = 0;
   t_bJet_SFdown = 0;
   t_GenPart_eta = 0;
   t_GenPart_mass = 0;
   t_GenPart_phi = 0;
   t_GenPart_pt = 0;
   t_GenPart_genPartIdxMother = 0;
   t_GenPart_pdgId = 0;
   t_GenPart_status = 0;
   t_GenJet_eta = 0;
   t_GenJet_mass = 0;
   t_GenJet_phi = 0;
   t_GenJet_pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("t_run", &t_run, &b_t_run);
   fChain->SetBranchAddress("t_luminosityBlock", &t_luminosityBlock, &b_t_luminosityBlock);
   fChain->SetBranchAddress("t_event", &t_event, &b_t_event);
   fChain->SetBranchAddress("t_genWeight", &t_genWeight, &b_t_genWeight);
   fChain->SetBranchAddress("t_puWeight", &t_puWeight, &b_t_puWeight);
   fChain->SetBranchAddress("t_puWeightUp", &t_puWeightUp, &b_t_puWeightUp);
   fChain->SetBranchAddress("t_puWeightDown", &t_puWeightDown, &b_t_puWeightDown);
   fChain->SetBranchAddress("t_PrefireWeight", &t_PrefireWeight, &b_t_PrefireWeight);
   fChain->SetBranchAddress("t_PrefireWeight_Up", &t_PrefireWeight_Up, &b_t_PrefireWeight_Up);
   fChain->SetBranchAddress("t_PrefireWeight_Down", &t_PrefireWeight_Down, &b_t_PrefireWeight_Down);
   fChain->SetBranchAddress("t_mu1", &t_mu1, &b_t_mu1);
   fChain->SetBranchAddress("t_mu2", &t_mu2, &b_t_mu2);
   fChain->SetBranchAddress("t_index_trigm_mu", &t_index_trigm_mu, &b_t_index_trigm_mu);
   fChain->SetBranchAddress("t_El_charge", &t_El_charge, &b_t_El_charge);
   fChain->SetBranchAddress("t_El_pt", &t_El_pt, &b_t_El_pt);
   fChain->SetBranchAddress("t_El_phi", &t_El_phi, &b_t_El_phi);
   fChain->SetBranchAddress("t_El_eta", &t_El_eta, &b_t_El_eta);
   fChain->SetBranchAddress("t_El_mass", &t_El_mass, &b_t_El_mass);
   fChain->SetBranchAddress("t_El_cutBased", &t_El_cutBased, &b_t_El_cutBased);
   fChain->SetBranchAddress("t_El_tightCharge", &t_El_tightCharge, &b_t_El_tightCharge);
   fChain->SetBranchAddress("t_El_cutBased_HEEP", &t_El_cutBased_HEEP, &b_t_El_cutBased_HEEP);
   fChain->SetBranchAddress("t_El_isPFcand", &t_El_isPFcand, &b_t_El_isPFcand);
   fChain->SetBranchAddress("t_El_pfRelIso03_all", &t_El_pfRelIso03_all, &b_t_El_pfRelIso03_all);
   fChain->SetBranchAddress("t_El_pfRelIso03_chg", &t_El_pfRelIso03_chg, &b_t_El_pfRelIso03_chg);
   fChain->SetBranchAddress("t_El_miniPFRelIso03_all", &t_El_miniPFRelIso03_all, &b_t_El_miniPFRelIso03_all);
   fChain->SetBranchAddress("t_El_dxy", &t_El_dxy, &b_t_El_dxy);
   fChain->SetBranchAddress("t_El_dxyErr", &t_El_dxyErr, &b_t_El_dxyErr);
   fChain->SetBranchAddress("t_El_dz", &t_El_dz, &b_t_El_dz);
   fChain->SetBranchAddress("t_El_dzErr", &t_El_dzErr, &b_t_El_dzErr);
   fChain->SetBranchAddress("t_El_sip3d", &t_El_sip3d, &b_t_El_sip3d);
   fChain->SetBranchAddress("t_Electron_mvaFall17Iso", &t_Electron_mvaFall17Iso, &b_t_Electron_mvaFall17Iso);
   fChain->SetBranchAddress("t_Electron_mvaFall17Iso_WP80", &t_Electron_mvaFall17Iso_WP80, &b_t_Electron_mvaFall17Iso_WP80);
   fChain->SetBranchAddress("t_Electron_mvaFall17Iso_WP90", &t_Electron_mvaFall17Iso_WP90, &b_t_Electron_mvaFall17Iso_WP90);
   fChain->SetBranchAddress("t_Electron_mvaFall17Iso_WPL", &t_Electron_mvaFall17Iso_WPL, &b_t_Electron_mvaFall17Iso_WPL);
   fChain->SetBranchAddress("t_Electron_mvaFall17noIso_WP80", &t_Electron_mvaFall17noIso_WP80, &b_t_Electron_mvaFall17noIso_WP80);
   fChain->SetBranchAddress("t_Electron_mvaFall17noIso_WP90", &t_Electron_mvaFall17noIso_WP90, &b_t_Electron_mvaFall17noIso_WP90);
   fChain->SetBranchAddress("t_Electron_mvaFall17noIso_WPL", &t_Electron_mvaFall17noIso_WPL, &b_t_Electron_mvaFall17noIso_WPL);
   fChain->SetBranchAddress("t_Mu_charge", &t_Mu_charge, &b_t_Mu_charge);
   fChain->SetBranchAddress("t_Mu_EffSF_TRIG", &t_Mu_EffSF_TRIG, &b_t_Mu_EffSF_TRIG);
   fChain->SetBranchAddress("t_Mu_EffSFErr_TRIG", &t_Mu_EffSFErr_TRIG, &b_t_Mu_EffSFErr_TRIG);
   fChain->SetBranchAddress("t_Mu_EffSF_ID", &t_Mu_EffSF_ID, &b_t_Mu_EffSF_ID);
   //fChain->SetBranchAddress("t_Mu_EffSFErr_ID", &t_Mu_EffSFErr_ID, &b_t_Mu_EffSFErr_ID);
   fChain->SetBranchAddress("t_Mu_EffSF_ISO", &t_Mu_EffSF_ISO, &b_t_Mu_EffSF_ISO);
   //fChain->SetBranchAddress("t_Mu_EffSFErr_ISO", &t_Mu_EffSFErr_ISO, &b_t_Mu_EffSFErr_ISO);
   fChain->SetBranchAddress("t_Mu_pt", &t_Mu_pt, &b_t_Mu_pt);
   fChain->SetBranchAddress("t_Mu_ptErr", &t_Mu_ptErr, &b_t_Mu_ptErr);
   fChain->SetBranchAddress("t_Mu_phi", &t_Mu_phi, &b_t_Mu_phi);
   fChain->SetBranchAddress("t_Mu_eta", &t_Mu_eta, &b_t_Mu_eta);
   fChain->SetBranchAddress("t_Mu_mass", &t_Mu_mass, &b_t_Mu_mass);
   fChain->SetBranchAddress("t_Mu_dxy", &t_Mu_dxy, &b_t_Mu_dxy);
   fChain->SetBranchAddress("t_Mu_dxyErr", &t_Mu_dxyErr, &b_t_Mu_dxyErr);
   fChain->SetBranchAddress("t_Mu_dz", &t_Mu_dz, &b_t_Mu_dz);
   fChain->SetBranchAddress("t_Mu_dzErr", &t_Mu_dzErr, &b_t_Mu_dzErr);
   fChain->SetBranchAddress("t_Mu_sip3d", &t_Mu_sip3d, &b_t_Mu_sip3d);
   fChain->SetBranchAddress("t_Mu_pfRelIso03_all", &t_Mu_pfRelIso03_all, &b_t_Mu_pfRelIso03_all);
   fChain->SetBranchAddress("t_Mu_pfRelIso03_chg", &t_Mu_pfRelIso03_chg, &b_t_Mu_pfRelIso03_chg);
   fChain->SetBranchAddress("t_Mu_pfRelIso04_all", &t_Mu_pfRelIso04_all, &b_t_Mu_pfRelIso04_all);
   fChain->SetBranchAddress("t_Mu_miniPFRelIso_all", &t_Mu_miniPFRelIso_all, &b_t_Mu_miniPFRelIso_all);
   fChain->SetBranchAddress("t_Mu_tightCharge", &t_Mu_tightCharge, &b_t_Mu_tightCharge);
   fChain->SetBranchAddress("t_Mu_isPFcand", &t_Mu_isPFcand, &b_t_Mu_isPFcand);
   fChain->SetBranchAddress("t_Mu_isglobal", &t_Mu_isglobal, &b_t_Mu_isglobal);
   fChain->SetBranchAddress("t_Mu_istracker", &t_Mu_istracker, &b_t_Mu_istracker);
   fChain->SetBranchAddress("t_Mu_mediumId", &t_Mu_mediumId, &b_t_Mu_mediumId);
   fChain->SetBranchAddress("t_Mu_softId", &t_Mu_softId, &b_t_Mu_softId);
   fChain->SetBranchAddress("t_Mu_tightId", &t_Mu_tightId, &b_t_Mu_tightId);
   fChain->SetBranchAddress("t_Mu_nStations", &t_Mu_nStations, &b_t_Mu_nStations);
   fChain->SetBranchAddress("t_Mu_nTrackerLayers", &t_Mu_nTrackerLayers, &b_t_Mu_nTrackerLayers);
   fChain->SetBranchAddress("t_diMuon_pt", &t_diMuon_pt, &b_t_diMuon_pt);
   fChain->SetBranchAddress("t_diMuon_eta", &t_diMuon_eta, &b_t_diMuon_eta);
   fChain->SetBranchAddress("t_diMuon_phi", &t_diMuon_phi, &b_t_diMuon_phi);
   fChain->SetBranchAddress("t_diMuon_mass", &t_diMuon_mass, &b_t_diMuon_mass);
   fChain->SetBranchAddress("t_MET_phi", &t_MET_phi, &b_t_MET_phi);
   fChain->SetBranchAddress("t_MET_pt", &t_MET_pt, &b_t_MET_pt);
   fChain->SetBranchAddress("t_MET_sumEt", &t_MET_sumEt, &b_t_MET_sumEt);
   fChain->SetBranchAddress("t_FatJet_area", &t_FatJet_area, &b_t_FatJet_area);
   fChain->SetBranchAddress("t_FatJet_btagCMVA", &t_FatJet_btagCMVA, &b_t_FatJet_btagCMVA);
   fChain->SetBranchAddress("t_FatJet_btagCSVV2", &t_FatJet_btagCSVV2, &b_t_FatJet_btagCSVV2);
   fChain->SetBranchAddress("t_FatJet_btagDeepB", &t_FatJet_btagDeepB, &b_t_FatJet_btagDeepB);
   fChain->SetBranchAddress("t_FatJet_eta", &t_FatJet_eta, &b_t_FatJet_eta);
   fChain->SetBranchAddress("t_FatJet_mass", &t_FatJet_mass, &b_t_FatJet_mass);
   fChain->SetBranchAddress("t_FatJet_msoftdrop", &t_FatJet_msoftdrop, &b_t_FatJet_msoftdrop);
   fChain->SetBranchAddress("t_FatJet_n2b1", &t_FatJet_n2b1, &b_t_FatJet_n2b1);
   fChain->SetBranchAddress("t_FatJet_n3b1", &t_FatJet_n3b1, &b_t_FatJet_n3b1);
   fChain->SetBranchAddress("t_FatJet_phi", &t_FatJet_phi, &b_t_FatJet_phi);
   fChain->SetBranchAddress("t_FatJet_pt", &t_FatJet_pt, &b_t_FatJet_pt);
   fChain->SetBranchAddress("t_FatJet_tau1", &t_FatJet_tau1, &b_t_FatJet_tau1);
   fChain->SetBranchAddress("t_FatJet_tau2", &t_FatJet_tau2, &b_t_FatJet_tau2);
   fChain->SetBranchAddress("t_FatJet_tau3", &t_FatJet_tau3, &b_t_FatJet_tau3);
   fChain->SetBranchAddress("t_FatJet_tau4", &t_FatJet_tau4, &b_t_FatJet_tau4);
   fChain->SetBranchAddress("t_FatJet_jetId", &t_FatJet_jetId, &b_t_FatJet_jetId);
   fChain->SetBranchAddress("t_FatJet_subJetIdx1", &t_FatJet_subJetIdx1, &b_t_FatJet_subJetIdx1);
   fChain->SetBranchAddress("t_FatJet_subJetIdx2", &t_FatJet_subJetIdx2, &b_t_FatJet_subJetIdx2);
   fChain->SetBranchAddress("t_SubJet_btagCMVA", &t_SubJet_btagCMVA, &b_t_SubJet_btagCMVA);
   fChain->SetBranchAddress("t_SubJet_btagCSVV2", &t_SubJet_btagCSVV2, &b_t_SubJet_btagCSVV2);
   fChain->SetBranchAddress("t_SubJet_btagDeepB", &t_SubJet_btagDeepB, &b_t_SubJet_btagDeepB);
   fChain->SetBranchAddress("t_SubJet_eta", &t_SubJet_eta, &b_t_SubJet_eta);
   fChain->SetBranchAddress("t_SubJet_mass", &t_SubJet_mass, &b_t_SubJet_mass);
   fChain->SetBranchAddress("t_SubJet_n2b1", &t_SubJet_n2b1, &b_t_SubJet_n2b1);
   fChain->SetBranchAddress("t_SubJet_n3b1", &t_SubJet_n3b1, &b_t_SubJet_n3b1);
   fChain->SetBranchAddress("t_SubJet_phi", &t_SubJet_phi, &b_t_SubJet_phi);
   fChain->SetBranchAddress("t_SubJet_pt", &t_SubJet_pt, &b_t_SubJet_pt);
   fChain->SetBranchAddress("t_SubJet_tau1", &t_SubJet_tau1, &b_t_SubJet_tau1);
   fChain->SetBranchAddress("t_SubJet_tau2", &t_SubJet_tau2, &b_t_SubJet_tau2);
   fChain->SetBranchAddress("t_SubJet_tau3", &t_SubJet_tau3, &b_t_SubJet_tau3);
   fChain->SetBranchAddress("t_SubJet_tau4", &t_SubJet_tau4, &b_t_SubJet_tau4);
   fChain->SetBranchAddress("t_nJet", &t_nJet, &b_t_nJet);
   fChain->SetBranchAddress("t_Jet_area", &t_Jet_area, &b_t_Jet_area);
   fChain->SetBranchAddress("t_Jet_btagCMVA", &t_Jet_btagCMVA, &b_t_Jet_btagCMVA);
   fChain->SetBranchAddress("t_Jet_btagCSVV2", &t_Jet_btagCSVV2, &b_t_Jet_btagCSVV2);
   fChain->SetBranchAddress("t_Jet_btagDeepB", &t_Jet_btagDeepB, &b_t_Jet_btagDeepB);
   fChain->SetBranchAddress("t_Jet_btagDeepC", &t_Jet_btagDeepC, &b_t_Jet_btagDeepC);
   fChain->SetBranchAddress("t_Jet_btagDeepFlavB", &t_Jet_btagDeepFlavB, &b_t_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("t_Jet_chEmEF", &t_Jet_chEmEF, &b_t_Jet_chEmEF);
   fChain->SetBranchAddress("t_Jet_chHEF", &t_Jet_chHEF, &b_t_Jet_chHEF);
   fChain->SetBranchAddress("t_Jet_eta", &t_Jet_eta, &b_t_Jet_eta);
   fChain->SetBranchAddress("t_Jet_mass", &t_Jet_mass, &b_t_Jet_mass);
   fChain->SetBranchAddress("t_Jet_neEmEF", &t_Jet_neEmEF, &b_t_Jet_neEmEF);
   fChain->SetBranchAddress("t_Jet_neHEF", &t_Jet_neHEF, &b_t_Jet_neHEF);
   fChain->SetBranchAddress("t_Jet_phi", &t_Jet_phi, &b_t_Jet_phi);
   fChain->SetBranchAddress("t_Jet_pt", &t_Jet_pt, &b_t_Jet_pt);
   fChain->SetBranchAddress("t_Jet_qgl", &t_Jet_qgl, &b_t_Jet_qgl);
   fChain->SetBranchAddress("t_Jet_jetId", &t_Jet_jetId, &b_t_Jet_jetId);
   fChain->SetBranchAddress("t_Jet_nConstituents", &t_Jet_nConstituents, &b_t_Jet_nConstituents);
   fChain->SetBranchAddress("t_Jet_nElectrons", &t_Jet_nElectrons, &b_t_Jet_nElectrons);
   fChain->SetBranchAddress("t_Jet_nMuons", &t_Jet_nMuons, &b_t_Jet_nMuons);
   fChain->SetBranchAddress("t_Jet_puId", &t_Jet_puId, &b_t_Jet_puId);
   fChain->SetBranchAddress("t_diJet_pt", &t_diJet_pt, &b_t_diJet_pt);
   fChain->SetBranchAddress("t_diJet_eta", &t_diJet_eta, &b_t_diJet_eta);
   fChain->SetBranchAddress("t_diJet_phi", &t_diJet_phi, &b_t_diJet_phi);
   fChain->SetBranchAddress("t_diJet_mass", &t_diJet_mass, &b_t_diJet_mass);
   fChain->SetBranchAddress("t_diJet_mass_mo", &t_diJet_mass_mo, &b_t_diJet_mass_mo);
   fChain->SetBranchAddress("t_nbJet", &t_nbJet, &b_t_nbJet);
   fChain->SetBranchAddress("t_bJet_area", &t_bJet_area, &b_t_bJet_area);
   fChain->SetBranchAddress("t_bJet_btagCMVA", &t_bJet_btagCMVA, &b_t_bJet_btagCMVA);
   fChain->SetBranchAddress("t_bJet_btagCSVV2", &t_bJet_btagCSVV2, &b_t_bJet_btagCSVV2);
   fChain->SetBranchAddress("t_bJet_btagDeepB", &t_bJet_btagDeepB, &b_t_bJet_btagDeepB);
   fChain->SetBranchAddress("t_bJet_btagDeepC", &t_bJet_btagDeepC, &b_t_bJet_btagDeepC);
   fChain->SetBranchAddress("t_bJet_btagDeepFlavB", &t_bJet_btagDeepFlavB, &b_t_bJet_btagDeepFlavB);
   fChain->SetBranchAddress("t_bJet_chEmEF", &t_bJet_chEmEF, &b_t_bJet_chEmEF);
   fChain->SetBranchAddress("t_bJet_chHEF", &t_bJet_chHEF, &b_t_bJet_chHEF);
   fChain->SetBranchAddress("t_bJet_eta", &t_bJet_eta, &b_t_bJet_eta);
   fChain->SetBranchAddress("t_bJet_mass", &t_bJet_mass, &b_t_bJet_mass);
   fChain->SetBranchAddress("t_bJet_neEmEF", &t_bJet_neEmEF, &b_t_bJet_neEmEF);
   fChain->SetBranchAddress("t_bJet_neHEF", &t_bJet_neHEF, &b_t_bJet_neHEF);
   fChain->SetBranchAddress("t_bJet_phi", &t_bJet_phi, &b_t_bJet_phi);
   fChain->SetBranchAddress("t_bJet_pt", &t_bJet_pt, &b_t_bJet_pt);
   fChain->SetBranchAddress("t_bJet_qgl", &t_bJet_qgl, &b_t_bJet_qgl);
   fChain->SetBranchAddress("t_bJet_jetId", &t_bJet_jetId, &b_t_bJet_jetId);
   fChain->SetBranchAddress("t_bJet_nConstituents", &t_bJet_nConstituents, &b_t_bJet_nConstituents);
   fChain->SetBranchAddress("t_bJet_nElectrons", &t_bJet_nElectrons, &b_t_bJet_nElectrons);
   fChain->SetBranchAddress("t_bJet_nMuons", &t_bJet_nMuons, &b_t_bJet_nMuons);
   fChain->SetBranchAddress("t_bJet_puId", &t_bJet_puId, &b_t_bJet_puId);
   fChain->SetBranchAddress("t_bJet_SF", &t_bJet_SF, &b_t_bJet_SF);
   fChain->SetBranchAddress("t_bJet_SFup", &t_bJet_SFup, &b_t_bJet_SFup);
   fChain->SetBranchAddress("t_bJet_SFdown", &t_bJet_SFdown, &b_t_bJet_SFdown);
   fChain->SetBranchAddress("t_PV_ndof", &t_PV_ndof, &b_t_PV_ndof);
   fChain->SetBranchAddress("t_PV_x", &t_PV_x, &b_t_PV_x);
   fChain->SetBranchAddress("t_PV_y", &t_PV_y, &b_t_PV_y);
   fChain->SetBranchAddress("t_PV_z", &t_PV_z, &b_t_PV_z);
   fChain->SetBranchAddress("t_PV_npvs", &t_PV_npvs, &b_t_PV_npvs);
   fChain->SetBranchAddress("t_PV_npvsGood", &t_PV_npvsGood, &b_t_PV_npvsGood);
   fChain->SetBranchAddress("t_GenPart_eta", &t_GenPart_eta, &b_t_GenPart_eta);
   fChain->SetBranchAddress("t_GenPart_mass", &t_GenPart_mass, &b_t_GenPart_mass);
   fChain->SetBranchAddress("t_GenPart_phi", &t_GenPart_phi, &b_t_GenPart_phi);
   fChain->SetBranchAddress("t_GenPart_pt", &t_GenPart_pt, &b_t_GenPart_pt);
   fChain->SetBranchAddress("t_GenPart_genPartIdxMother", &t_GenPart_genPartIdxMother, &b_t_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("t_GenPart_pdgId", &t_GenPart_pdgId, &b_t_GenPart_pdgId);
   fChain->SetBranchAddress("t_GenPart_status", &t_GenPart_status, &b_t_GenPart_status);
   fChain->SetBranchAddress("t_SoftActivityJetNjets5", &t_SoftActivityJetNjets5, &b_t_SoftActivityJetNjets5);
   fChain->SetBranchAddress("t_GenJet_eta", &t_GenJet_eta, &b_t_GenJet_eta);
   fChain->SetBranchAddress("t_GenJet_mass", &t_GenJet_mass, &b_t_GenJet_mass);
   fChain->SetBranchAddress("t_GenJet_phi", &t_GenJet_phi, &b_t_GenJet_phi);
   fChain->SetBranchAddress("t_GenJet_pt", &t_GenJet_pt, &b_t_GenJet_pt);
   Notify();
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NtupleVariables_cxx
