#define HmmAnalyzer_cxx
#include "HmmAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

#ifdef MAKECINT
#pragma link C++ class vector<float>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<int>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<bool>+;
#endif

int main(int argc, char* argv[])
{

  if(argc < 4) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  HmmAnalyzer Hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << endl;
  Hmm.EventLoop(data, isData);

  return 0;
}

void HmmAnalyzer::EventLoop(const char *data,const char *isData)
{ 
  if (fChain == 0) return;
  //clearTreeVectors();
  //cout<<"cleared tree vectors\n";
  //BookTreeBranches();
  //cout<<"booked tree branches\n";
  float muon_mass = 0.1056583745;
  bool datafile = true;
  if(*isData=='F') datafile = false;
 
  //btag SF
  BTagCalibration calib("deepcsv","data/btagSF/DeepCSV_94XSF_V3_B_F.csv");
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,  // operating point
			       "central",             // central sys type
			       {"up", "down"});      // other sys types

  reader.load(calib,                // calibration instance
	      BTagEntry::FLAV_B,    // btag flavour
	      "comb")       ;        // measurement type
  
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 6;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%5000==0) cout <<"entry: "<<jentry<<endl;

      //if(event!=709267) continue;

      clearTreeVectors();

      //sum of genWeight
      float value_h_sumOfgw = h_sumOfgw->GetBinContent(1);
      if(*isData=='F')   value_h_sumOfgw = value_h_sumOfgw + genWeight;
      else value_h_sumOfgw = value_h_sumOfgw + 1.0;
      h_sumOfgw->SetBinContent(1,value_h_sumOfgw);

      //sum of genWeight and pileupweight
      float value_h_sumOfgpw = h_sumOfgpw->GetBinContent(1);
      if(*isData=='F')   value_h_sumOfgpw = value_h_sumOfgpw + genWeight*puWeight;
      else value_h_sumOfgpw = value_h_sumOfgpw + 1.0;
      h_sumOfgpw->SetBinContent(1,value_h_sumOfgpw);

      bool trig_decision = false;
      if( HLT_IsoMu24==1 || HLT_IsoTkMu24==1 ) trig_decision =true;

      int index_mu1(-999), index_mu2(-999); 
      bool run_muChecks =false; 
      if(nMuon>=2 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && Flag_BadChargedCandidateFilter && trig_decision && PV_npvsGood>0) run_muChecks =true;

      vector<float> mu_pt_Roch_corr, mu_ptErr_Roch_corr, mu_Iso_Roch_corr;
      mu_pt_Roch_corr.clear(), mu_ptErr_Roch_corr.clear(), mu_Iso_Roch_corr.clear();

      float pt_Roch, ptErr_Roch, pt_Roch_sys_up, pt_Roch_sys_down;
      pt_Roch = 0, ptErr_Roch = 0, pt_Roch_sys_up = 0, pt_Roch_sys_down =0;
      bool Event_sel= false;
      bool trig_match = false;  
      if(run_muChecks){
        for(int i=0;i<nMuon;i++){
            TLorentzVector mu_raw;
            mu_raw.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],muon_mass); //Muon_mass[i]);
            pt_Roch = 0, ptErr_Roch = 0, pt_Roch_sys_up = 0, pt_Roch_sys_down =0;
            //cout <<"pt "<<Muon_pt[i]<<" Err "<<Muon_ptErr[i] <<endl;

            float gen_pt = -999.;
            if(*isData=='F'){
              for(int j=0;j<nGenPart;j++){
                if(Muon_charge[i]==-1 && GenPart_pdgId[j]==13 && DeltaR(Muon_eta[i], Muon_phi[i], GenPart_eta[j], GenPart_phi[j]) <0.1){ gen_pt = GenPart_pt[j]; break;}
                else if(Muon_charge[i]==1 && GenPart_pdgId[j]==-13 && DeltaR(Muon_eta[i], Muon_phi[i], GenPart_eta[j], GenPart_phi[j]) <0.1){ gen_pt = GenPart_pt[j]; break;}
              }
            }
            CorrectPtRoch( _Roch_calib, false, mu_raw,
                   pt_Roch, ptErr_Roch, pt_Roch_sys_up, pt_Roch_sys_down,
                    Muon_charge[i], Muon_nTrackerLayers[i], gen_pt, datafile );
            //cout <<"pt_Roch "<<pt_Roch<<endl;
            mu_pt_Roch_corr.push_back(pt_Roch);
            mu_ptErr_Roch_corr.push_back(ptErr_Roch); 
        }
        
        for(int i=0;i<nMuon;i++){
          if( Muon_isGlobal[i] && Muon_isTracker[i] && mu_pt_Roch_corr[i]>26. && Muon_mediumId[i] && fabs(Muon_eta[i])<2.4 && Muon_pfRelIso04_all[i]<0.25){
              for(int j=i+1;j<nMuon;j++){
                 if( Muon_isGlobal[j] && Muon_isTracker[j] && Muon_charge[i]*Muon_charge[j]== -1 && mu_pt_Roch_corr[j]>20. && Muon_mediumId[j] && fabs(Muon_eta[j])<2.4 && Muon_pfRelIso04_all[j]<0.25){
                    Event_sel =true; 
                    index_mu1 = i; 
                    index_mu2 = j;
                    break;
                 }
              }
              if(Event_sel =true) break;
          }
       }
      
       for(int i=0; i<nTrigObj; i++){
         float dR_TrigObj = 999.;
         if(TrigObj_id[i]==13){
              dR_TrigObj = DeltaR(Muon_eta[index_mu1], Muon_phi[index_mu1], TrigObj_eta[i], TrigObj_phi[i]);
              if(dR_TrigObj<0.1 && mu_pt_Roch_corr[index_mu1]>26.){ 
                 trig_match = true;
                 t_index_trigm_mu = 1;
                 break;
              }
              else{
                if(mu_pt_Roch_corr[index_mu2]>26.){
                  dR_TrigObj = DeltaR(Muon_eta[index_mu2], Muon_phi[index_mu2], TrigObj_eta[i], TrigObj_phi[i]);
                  if(dR_TrigObj<0.1){ 
                    trig_match = true;
                    t_index_trigm_mu = 2;
                    break;
                  }
                }
             } 
          }
       }//end of triger match, end of loop over trigger objects
      
      }
      if(Event_sel && trig_match){
	t_run =run;
	t_luminosityBlock=luminosityBlock;
	t_event=event;
	//cout<<jentry<<" : "<<t_event<<"-------------------\n";
        int t_index_mu1 = -999;
        int t_index_mu2 = -999;
        int t_index = 0;
        /*
        if(*isData=='F'){
               t_pileupWeight = getPileupWeight(Pileup_nPU);
               t_pileupupWeight = getPileupWeightUp(Pileup_nPU);
               t_pileupdnWeight = getPileupWeightDown(Pileup_nPU);
        }
        */
	for(int i=0;i<nMuon;i++){
          if(fabs(Muon_eta[i])<2.4 && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.25){
            if(i==index_mu1) t_index_mu1 = t_index;
            if(i==index_mu2) t_index_mu2 = t_index;
            if(*isData=='F'){
               t_Mu_EffSF_TRIG->push_back(Mu_eff_SF_TRIG.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSFErr_TRIG->push_back(Mu_eff_SF_TRIG.getSFErr(13,Muon_pt[i],Muon_eta[i]));
               t_Mu_EffSF_ID->push_back(Mu_eff_SF_ID.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSF_ID_stat->push_back(Mu_eff_SF_ID_stat.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSF_ID_syst->push_back(Mu_eff_SF_ID_syst.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSF_ISO->push_back(Mu_eff_SF_ISO.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSF_ISO_stat->push_back(Mu_eff_SF_ISO_stat.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
               t_Mu_EffSF_ISO_syst->push_back(Mu_eff_SF_ISO_syst.getSFAve(11,Muon_pt[i],Muon_eta[i],0.5548));
            }
	    t_Mu_charge->push_back(Muon_charge[i]);   
	    t_Mu_pt->push_back(mu_pt_Roch_corr[i]);   
	    t_Mu_ptErr->push_back(mu_ptErr_Roch_corr[i]);   
	    t_Mu_phi->push_back(Muon_phi[i]);   
	    t_Mu_eta->push_back(Muon_eta[i]);   
	    t_Mu_mass->push_back(muon_mass); //Muon_mass[i]);  
	    t_Mu_dxy->push_back(Muon_dxy[i]);   
	    t_Mu_dxyErr->push_back(Muon_dxyErr[i]);   
	    t_Mu_dz->push_back(Muon_dz[i]);   
	    t_Mu_dzErr->push_back(Muon_dzErr[i]); 
	    t_Mu_pfRelIso03_all->push_back(Muon_pfRelIso03_all[i]);   
	    t_Mu_pfRelIso03_chg->push_back(Muon_pfRelIso03_chg[i]);   
	    t_Mu_pfRelIso04_all->push_back(Muon_pfRelIso04_all[i]);   
            t_Mu_miniPFRelIso_all->push_back(Muon_miniPFRelIso_all[i]);
            t_Mu_miniPFRelIso_chg->push_back(Muon_miniPFRelIso_chg[i]);
	    t_Mu_tightCharge->push_back(Muon_tightCharge[i]);   
	    t_Mu_isPFcand->push_back(Muon_isPFcand[i]);
            //t_Mu_isglobal->push_back(Muon_isglobal[i]);
            //t_Mu_istracker->push_back(Muon_istracker[i]);   
	    t_Mu_mediumId->push_back(Muon_mediumId[i]);   
	    t_Mu_softId->push_back(Muon_softId[i]);   
	    t_Mu_tightId->push_back(Muon_tightId[i]);    
	    t_Mu_nStations->push_back(Muon_nStations[i]);   
	    t_Mu_nTrackerLayers->push_back(Muon_nTrackerLayers[i]);   
            t_index++;
	  }
	}
        t_mu1 = t_index_mu1;
        t_mu2 = t_index_mu2;
        //if(t_index_trigm_mu ==1) t_index_trigm_mu = t_mu1;
        //else t_index_trigm_mu = t_mu2;

	TLorentzVector dimu, mu1,mu2;
	mu1.SetPtEtaPhiM((mu_pt_Roch_corr)[index_mu1],(Muon_eta)[index_mu1],(Muon_phi)[index_mu1],muon_mass);
	mu2.SetPtEtaPhiM((mu_pt_Roch_corr)[index_mu2],(Muon_eta)[index_mu2],(Muon_phi)[index_mu2],muon_mass);
	dimu=mu1+mu2;
	t_diMuon_pt = dimu.Pt();
	t_diMuon_eta= dimu.Eta();
	t_diMuon_phi= dimu.Phi();
	t_diMuon_mass= dimu.M();
        t_SoftActivityJetNjets5 = SoftActivityJetNjets5;

	for (int j =0;j<nJet;j++){
	  bool matched_mu=false;
	  for(int i=0;i<t_Mu_pt->size();i++){
	    double dR= DeltaR((*t_Mu_eta)[i],(*t_Mu_phi)[i],Jet_eta[j],Jet_phi[j]);
	    if(dR<0.4){
	      matched_mu=true;
	      break;
	    }
	  }
	  if(!matched_mu && Jet_pt[j]>25. && fabs(Jet_eta[j])<4.7 && Jet_jetId[j]>=2 && Jet_puId[j]>=1){
	    t_nJet++;
	    t_Jet_area->push_back(Jet_area[j]);
	    //t_Jet_btagCMVA->push_back(Jet_btagCMVA[j]);   
	    //t_Jet_btagCSVV2->push_back(Jet_btagCSVV2[j]);   
	    t_Jet_btagDeepB->push_back(Jet_btagDeepB[j]);   
	    t_Jet_btagDeepC->push_back(Jet_btagDeepC[j]);   
	    t_Jet_btagDeepFlavB->push_back(Jet_btagDeepFlavB[j]);   
	    t_Jet_chEmEF->push_back(Jet_chEmEF[j]);   
	    t_Jet_chHEF->push_back(Jet_chHEF[j]);   
	    t_Jet_eta->push_back(Jet_eta[j]);   
	    t_Jet_mass->push_back(Jet_mass[j]);   
	    t_Jet_neEmEF->push_back(Jet_neEmEF[j]);   
	    t_Jet_neHEF->push_back(Jet_neHEF[j]);   
	    t_Jet_phi->push_back(Jet_phi[j]);   
	    t_Jet_pt->push_back(Jet_pt[j]);   
	    t_Jet_qgl->push_back(Jet_qgl[j]);   
	    t_Jet_jetId->push_back(Jet_jetId[j]);   
	    t_Jet_nConstituents->push_back(Jet_nConstituents[j]);   
	    t_Jet_nElectrons->push_back(Jet_nElectrons[j]);   
	    t_Jet_nMuons->push_back(Jet_nMuons[j]);   
	    t_Jet_puId->push_back(Jet_puId[j]);   
            if(Jet_btagDeepB[j]>0.6321){ //medium WP
              if(*isData!='T'){
		double jet_scalefactor    = reader.eval_auto_bounds(
								    "central", 
								    BTagEntry::FLAV_B, 
								    fabs(Jet_eta[j]), // absolute value of eta
								    Jet_pt[j]
								    ); 
		double jet_scalefactor_up = reader.eval_auto_bounds(
								    "up", BTagEntry::FLAV_B, fabs(Jet_eta[j]), Jet_pt[j]);
		double jet_scalefactor_do = reader.eval_auto_bounds(
								    "down", BTagEntry::FLAV_B, fabs(Jet_eta[j]), Jet_pt[j]); 
		//cout<<jet_scalefactor<<" "<<jet_scalefactor_up<<" "<<jet_scalefactor_do<<endl;
		t_bJet_SF->push_back(jet_scalefactor);
		t_bJet_SFup->push_back(jet_scalefactor_up);
		t_bJet_SFdown->push_back(jet_scalefactor_do);
	      }
	      t_nbJet++;
	      t_bJet_area->push_back(Jet_area[j]);
	      //t_bJet_btagCMVA->push_back(Jet_btagCMVA[j]);   
	      //t_bJet_btagCSVV2->push_back(Jet_btagCSVV2[j]);   
	      t_bJet_btagDeepB->push_back(Jet_btagDeepB[j]);   
	      t_bJet_btagDeepC->push_back(Jet_btagDeepC[j]);   
	      t_bJet_btagDeepFlavB->push_back(Jet_btagDeepFlavB[j]);   
	      t_bJet_chEmEF->push_back(Jet_chEmEF[j]);   
	      t_bJet_chHEF->push_back(Jet_chHEF[j]);   
	      t_bJet_eta->push_back(Jet_eta[j]);   
	      t_bJet_mass->push_back(Jet_mass[j]);   
	      t_bJet_neEmEF->push_back(Jet_neEmEF[j]);   
	      t_bJet_neHEF->push_back(Jet_neHEF[j]);   
	      t_bJet_phi->push_back(Jet_phi[j]);   
	      t_bJet_pt->push_back(Jet_pt[j]);   
	      t_bJet_qgl->push_back(Jet_qgl[j]);   
	      t_bJet_jetId->push_back(Jet_jetId[j]);   
	      t_bJet_nConstituents->push_back(Jet_nConstituents[j]);   
	      t_bJet_nElectrons->push_back(Jet_nElectrons[j]);   
	      t_bJet_nMuons->push_back(Jet_nMuons[j]);   
	      t_bJet_puId->push_back(Jet_puId[j]);   

	  
	    }
 
	  }
	}
	if(t_Jet_pt->size()>=2){
          for(int k=0;k<t_Jet_pt->size();k++){
            for(int m=k+1; m<t_Jet_pt->size();m++){
              TLorentzVector j1,j2, jj;
              if(k==0 && m==1){
          	  j1.SetPtEtaPhiM((*t_Jet_pt)[0], (*t_Jet_eta)[0],(*t_Jet_phi)[0],(*t_Jet_mass)[0]);
	          j2.SetPtEtaPhiM((*t_Jet_pt)[1], (*t_Jet_eta)[1],(*t_Jet_phi)[1],(*t_Jet_mass)[1]);
	          jj=j1+j2;	
	          t_diJet_pt = jj.Pt();
	          t_diJet_eta=jj.Eta();
	          t_diJet_phi=jj.Phi();
	          t_diJet_mass=jj.M();
                  t_diJet_mass_mo=jj.M();
              }
              else{
                  j1.SetPtEtaPhiM((*t_Jet_pt)[k], (*t_Jet_eta)[k],(*t_Jet_phi)[k],(*t_Jet_mass)[k]);
                  j2.SetPtEtaPhiM((*t_Jet_pt)[m], (*t_Jet_eta)[m],(*t_Jet_phi)[m],(*t_Jet_mass)[m]);
                  jj=j1+j2;
                  if(t_diJet_mass_mo<jj.M()) t_diJet_mass_mo=jj.M();

             }
            }
          }
	}

	for(int i=0;i<nElectron;i++){
          //t_El_genPartIdx->push_back(Electron_genPartIdx[i]);
          //t_El_genPartFlav->push_back(Electron_genPartFlav[i]);
	  t_El_charge->push_back(Electron_charge[i]);
	  t_El_pt->push_back(Electron_pt[i]);
	  t_El_phi->push_back(Electron_phi[i]);
	  t_El_eta->push_back(Electron_eta[i]);   
	  t_El_mass->push_back(Electron_mass[i]);
	  t_El_cutBased->push_back(Electron_cutBased[i]);   
	  t_El_tightCharge->push_back(Electron_tightCharge[i]);   
	  t_El_cutBased_HEEP->push_back(Electron_cutBased_HEEP[i]);   
	  t_El_isPFcand->push_back(Electron_isPFcand[i]);   
	  t_El_pfRelIso03_all->push_back(Electron_pfRelIso03_all[i]);   
	  t_El_pfRelIso03_chg->push_back(Electron_pfRelIso03_chg[i]);      
          t_El_miniPFRelIso_all->push_back(Electron_miniPFRelIso_all[i]);
          t_El_miniPFRelIso_chg->push_back(Electron_miniPFRelIso_chg[i]);
	  t_El_dxy->push_back(Electron_dxy[i]);   
	  t_El_dxyErr->push_back(Electron_dxyErr[i]);   
	  t_El_dz->push_back(Electron_dz[i]);   
	  t_El_dzErr->push_back(Electron_dzErr[i]);   
          t_Electron_mvaFall17Iso->push_back(Electron_mvaFall17V2Iso[i]);
	  t_Electron_mvaFall17Iso_WP80->push_back(Electron_mvaFall17V2Iso_WP80[i]);
	  t_Electron_mvaFall17Iso_WP90->push_back(Electron_mvaFall17V2Iso_WP90[i]);
	  t_Electron_mvaFall17Iso_WPL->push_back(Electron_mvaFall17V2Iso_WPL[i]);
          t_Electron_mvaFall17noIso->push_back(Electron_mvaFall17V2noIso[i]);
	  t_Electron_mvaFall17noIso_WP80->push_back(Electron_mvaFall17V2noIso_WP80[i]);
	  t_Electron_mvaFall17noIso_WP90->push_back(Electron_mvaFall17V2noIso_WP90[i]);
	  t_Electron_mvaFall17noIso_WPL->push_back(Electron_mvaFall17V2noIso_WPL[i]);
	}

	t_MET_pt = MET_pt;
	t_MET_phi = MET_phi;
	t_MET_sumEt  = MET_sumEt;
 
	for(int i=0;i<nFatJet;i++){
	  t_FatJet_area->push_back(FatJet_area[i]);  
	  t_FatJet_btagCMVA->push_back(FatJet_btagCMVA[i]);  
	  t_FatJet_btagCSVV2->push_back(FatJet_btagCSVV2[i]);  
	  t_FatJet_btagDeepB->push_back(FatJet_btagDeepB[i]);  
	  t_FatJet_eta->push_back(FatJet_eta[i]);  
	  t_FatJet_mass->push_back(FatJet_mass[i]);  
	  t_FatJet_msoftdrop->push_back(FatJet_msoftdrop[i]);  
	  t_FatJet_n2b1->push_back(FatJet_n2b1[i]);  
	  t_FatJet_n3b1->push_back(FatJet_n3b1[i]);  
	  t_FatJet_phi->push_back(FatJet_phi[i]);  
	  t_FatJet_pt->push_back(FatJet_pt[i]);  
	  t_FatJet_tau1->push_back(FatJet_tau1[i]);  
	  t_FatJet_tau2->push_back(FatJet_tau2[i]);  
	  t_FatJet_tau3->push_back(FatJet_tau3[i]);  
	  t_FatJet_tau4->push_back(FatJet_tau4[i]);  
	  t_FatJet_jetId->push_back(FatJet_tau4[i]);  
	  t_FatJet_subJetIdx1->push_back(FatJet_subJetIdx1[i]);  
	  t_FatJet_subJetIdx2->push_back(FatJet_subJetIdx2[i]);  
	} 
	for(int i=0;i<nSubJet;i++){
	  t_SubJet_btagCMVA->push_back(SubJet_btagCMVA[i]);   
	  t_SubJet_btagCSVV2->push_back(SubJet_btagCSVV2[i]);   
	  t_SubJet_btagDeepB->push_back(SubJet_btagDeepB[i]);   
	  t_SubJet_eta->push_back(SubJet_eta[i]);   
	  t_SubJet_mass->push_back(SubJet_mass[i]);   
	  t_SubJet_n2b1->push_back(SubJet_n2b1[i]);   
	  t_SubJet_n3b1->push_back(SubJet_n3b1[i]);   
	  t_SubJet_phi->push_back(SubJet_phi[i]);   
	  t_SubJet_pt->push_back(SubJet_pt[i]);   
	  t_SubJet_tau1->push_back(SubJet_tau1[i]);   
	  t_SubJet_tau2->push_back(SubJet_tau2[i]);   
	  t_SubJet_tau3->push_back(SubJet_tau3[i]);   
	  t_SubJet_tau4->push_back(SubJet_tau4[i]);   
	}

	t_PV_ndof = PV_ndof;
	t_PV_x = PV_x;
	t_PV_y = PV_y;
	t_PV_z = PV_z;
	t_PV_npvs = PV_npvs;
	t_PV_npvs = PV_npvsGood;

	if(*isData=='F'){
          t_genWeight = genWeight;
          t_puWeight = puWeight;
          t_puWeightUp = puWeightUp;
          t_puWeightDown = puWeightDown;
	  for(int i=0;i<nGenPart;i++){
	  
	    if((abs(GenPart_pdgId[i])>=11 && abs(GenPart_pdgId[i])<=16) || (abs(GenPart_pdgId[i])>=23 && abs(GenPart_pdgId[i])<=25) || (abs(GenPart_genPartIdxMother[i])>=23 && abs(GenPart_genPartIdxMother[i])<=25) || (abs(GenPart_genPartIdxMother[i])>=11 && abs(GenPart_genPartIdxMother[i])<=16) ){
	      t_GenPart_eta->push_back(GenPart_eta[i]);
	      t_GenPart_mass->push_back(GenPart_mass[i]);
	      t_GenPart_phi->push_back(GenPart_phi[i]);
	      t_GenPart_pt->push_back(GenPart_pt[i]);
	      t_GenPart_genPartIdxMother->push_back(GenPart_genPartIdxMother[i]);
	      t_GenPart_pdgId->push_back(GenPart_pdgId[i]);
	      t_GenPart_status->push_back(GenPart_status[i]);
	    }
	  }
          for(int i=0;i<nGenJet;i++){
              t_GenJet_eta->push_back(GenJet_eta[i]);
              t_GenJet_mass->push_back(GenJet_mass[i]);
              t_GenJet_phi->push_back(GenJet_phi[i]);
              t_GenJet_pt->push_back(GenJet_pt[i]);
          }
          t_nLHEPdfWeight = nLHEPdfWeight;
          for(int i=0; i<nLHEPdfWeight; i++){
              t_LHEPdfWeight->push_back(LHEPdfWeight[i]);
          }
         
          t_nLHEScaleWeight = nLHEScaleWeight;
          for(int i=0; i<nLHEScaleWeight; i++){
             t_LHEScaleWeight->push_back(LHEScaleWeight[i]);
          }
          t_nPSWeight = nPSWeight;
          for(int i=0; i<nPSWeight; i++){
             t_PSWeight->push_back(PSWeight[i]);
          }
	}
      
      //cout <<run<<" "<<luminosityBlock<<" "<<event<<" "<<mu_pt_Roch_corr[t_index_mu1]<<" "<<Muon_eta[t_index_mu1]<<mu_pt_Roch_corr[t_index_mu2]<<" "<<Muon_eta[t_index_mu2]<<" "<<t_diMuon_mass<<" "<<t_nJet<<" "<<t_diJet_mass_mo<<" "<<t_nbJet<<endl;

      tree->Fill();
      }
   }
}
