#define HiggsMuMu_cxx
#include "HiggsMuMu.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <map>
#include <vector>
#include <cstring>
#include<string>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
int main(int argc, char* argv[])
{

  if(argc < 3) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  /*
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
    TString scales        = argv[5];
  HiggsMuMu hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << scales<<endl;
  float scalef = scales.Atof();
  */
    
    const char *inputFileList = argv[1];
    const char *outFileName   = argv[2];
    const char *data          = argv[3];
    const char *isData        = argv[4];
    HiggsMuMu hmm(inputFileList, outFileName, data, isData);
    cout << "dataset " << data << " " << endl;
    
    map<TString,double> proc_scale;
    double lumi = 41.529*1000.;
    double lumi_18 = 59.74*1000.;
    double lumi_16 = 35.9*1000.;
   proc_scale["VBFHToMuMu_2017"]=0.000823*lumi/5913454.85352;
   proc_scale["VBFHToMuMu_2018"]=0.000823*lumi_18/3858502.5;
   proc_scale["VBFHToMuMu_2016"]=0.000823*lumi_16/8342013.06055;
   proc_scale["ZH_2018"]=0.000192*lumi_18/379921.21875;
   proc_scale["WplusH_2018"]=0.000183*lumi_18/260275.25;
   proc_scale["WminusH_2018"]=0.000116*lumi_18/162347.503906;
   proc_scale["ZH_2016"]=0.000192*lumi_16/413010.311768;
   proc_scale["WplusH_2016"]=0.000183*lumi_16/360018.046387;//105908.515625;
   proc_scale["WminusH_2016"]=0.000116*lumi_16/158522.949219;//66702.0078125;
   proc_scale["ttH_2018"]=0.000110*lumi_18/250980.332031;
   proc_scale["ZH_2017"]=0.000192*lumi/234457.452209;
   proc_scale["WplusH_2017"]=0.000183*lumi/260182.667969;
   proc_scale["WminusH_2017"]=0.000116*lumi/162194.871094;
   proc_scale["ttH_2017"]=0.000110*lumi/155117.609375;
   proc_scale["ttH_2016"]=0.000110*lumi_16/258743.374512;//258860.732483;
   proc_scale["DYJetsToLL"]=6225.42*lumi/(3258462550016.0+492179082112.0);
   proc_scale["DYJetsToLL_VBFfilter_2018"]=2.02*lumi_18/3173402777.5;//(3355073909.5);
   proc_scale["DYJetsToLL_VBFfilter_2017"]=2.02*lumi/(2259008890.6);
   proc_scale["DYJetsToLL_VBFfilter_2016"]=2.02*lumi_16/2570758278.5;
   proc_scale["DYJetsToLL_M105To160_incl_2018"]=46.9479*lumi_18/292139596.5;//(3355073909.5);
   proc_scale["DYJetsToLL_M105To160_incl_2017"]=46.9479*lumi/(7139863856.25);
   proc_scale["DYJetsToLL_M105To160_incl_2016"]=46.9479*lumi_16/7223964604.28;//(2699711552.0);
   proc_scale["ttTosemileptonic"]=6.871e+02*lumi/11784986264.000000;
   proc_scale["ttsl_2018"]=6.871e+02*lumi_18/29967356136.0;
   proc_scale["ttTo2l2v_2018"]=87.31*lumi_18/4600510308.0;
   proc_scale["ttsl_2017"]=6.871e+02*lumi/12941416132.0;
   proc_scale["ttTo2l2v_2017"]=87.31*lumi/5593093066.0;
   //proc_scale["ttsl_2016"]=6.871e+02*lumi_16/12941416132.0;
   proc_scale["ttTo2l2v_2016"]=87.31*lumi_16/79083400.0;//64203460.0;
   proc_scale["ttsl_2016"]=6.871e+02*lumi_16/150930664.0;
   proc_scale["ttJets_DiLept_2018"]=85.656*lumi_18/28671097.0469;// 287030939.5;
   proc_scale["ttJets_DiLept_2017"]=85.656*lumi/28349068.0;
   proc_scale["ttJets_DiLept_2016"]=85.656*lumi_16/24622652.0;
   proc_scale["EWK_2016"]=1.608*lumi_16/1576200.0;
   proc_scale["EWK_2017"]=1.608*lumi/3676958.2182;//3678099.0;
   proc_scale["EWK_2018"]=1.608*lumi_18/1978654.78125;//1978878.0;
   proc_scale["ggH_2018"]=0.010571*lumi_18/216906552.0;
   proc_scale["ggH_2017"]=0.010571*lumi/381714429.875;
   proc_scale["ggH_2016"]=0.010571*lumi_16/996618.0;//58565928.4755;
   proc_scale["ttTo2l2v"]=85.656*lumi/(623402174.0+4782395097.687500+199762.000000);
   proc_scale["TTTo2L"]=85.656*lumi/(623402174.0+4782395097.687500);
   proc_scale["WZTo1L1Nu2Q"]=1.161e+01*lumi/352741934.218750;
   proc_scale["WZTo3LNu_2016"]=4.42965*lumi_16/97457043.1172;//95696911.5444;
   proc_scale["WZTo2L2Q_2016"]=5.595*lumi_16/76043637.625;//211724236.125;
   proc_scale["WZTo3LNu_2017"]=4.42965*lumi/4571207.61914;
   proc_scale["WZTo2L2Q_2017"]=5.595*lumi/266683248.0;
   proc_scale["WZTo3LNu_2018"]=4.42965*lumi_18/187398698.125;//9218653.66309;
   proc_scale["WZTo2L2Q_2018"]=5.595*lumi_18/253329022.092;
   proc_scale["ZZTo4L_2016"]=1.256*lumi_16/82437095.0;//82474616.0;//5728400.0;
   proc_scale["ZZTo2L2Q_2016"]=3.22*lumi_16/496436.0;
   proc_scale["ZZTo2L2Nu_2016"]=0.564*lumi_16/8931750.0;//7530975.0;
   proc_scale["ZZTo4L_2017"]=1.256*lumi/159446233.461;
   proc_scale["ZZTo2L2Q_2017"]=3.22*lumi/157683326.25;
   proc_scale["ZZTo2L2Nu_2017"]=0.564*lumi/5303017.48047;
   proc_scale["ZZTo4L_2018"]=1.256*lumi_18/8919844.51562;
   proc_scale["ZZTo2L2Q_2018"]=3.22*lumi_18/157419820.886;
   proc_scale["ZZTo2L2Nu_2018"]=0.564*lumi_18/5097461.55273;
   
   proc_scale["WWTo2L2Nu_2016"]=12.46*lumi_16/1999000.0;//1578014.0;
   proc_scale["WWToLNuQQ_2016"]=4.599e+01*lumi_16/1999200.0;
   proc_scale["WWTo2L2Nu_2017"]=12.46*lumi/43938136.25;
   proc_scale["WWToLNuQQ_2017"]=4.599e+01*lumi/1285749509.5;
   proc_scale["WWTo2L2Nu_2018"]=12.46*lumi_18/85338022.0947;
   proc_scale["WWToLNuQQ_2018"]=4.599e+01*lumi_18/643636627.5;
   proc_scale["WWW_4F_2016"]=0.2086*lumi_16/50012.1679688;//50040.8632812;
   proc_scale["WZZ_2016"]=0.05565*lumi_16/13735.4013062;//11341.4560547;
   proc_scale["ZZZ_2016"]=0.01398*lumi_16/3499.77664185;//3500.18438721;
   proc_scale["WWW_4F_2017"]=0.2086*lumi/49913.048584;
   proc_scale["WZZ_2017"]=0.05565*lumi/13943.3115234;
   proc_scale["ZZZ_2017"]=0.01398*lumi/3676.50537109;
   proc_scale["WWZ_4F_2017"]=0.1651*lumi/41844.3398438;
   proc_scale["WWW_4F_2018"]=0.2086*lumi_18/51529.8018799;
   proc_scale["WZZ_2018"]=0.05565*lumi_18/14236.171875;
   proc_scale["ZZZ_2018"]=0.01398*lumi_18/7366.66020203;//3683.33010101;
   proc_scale["WWZ_4F_2018"]=0.1651*lumi/41737.7509766;
   proc_scale["DYJetsToLL_2017"]=5765.4*lumi/3.74020946509e+12;
   proc_scale["DYJetsToLL_2018"]=5765.4*lumi_18/2.83840120718e+12;//100113543.375;//17846205568.0;
   proc_scale["DYJetsToLL_2016"]=5765.4*lumi_16/1.86045481216e+12;//1856607046400.000000000;//1.8938682816e+12;
   proc_scale["DYJetsToLL_ext"]=5765.4*lumi/3258462550016.0;
   proc_scale["DYJetsToLL_small"]=5765.4*lumi/492179082112.0;
   proc_scale["TTTo2L_small"]=85.656*lumi/623402174.0;
   proc_scale["TTTo2L2Nu_TuneCP5_PSweights"]=85.656*lumi/4782395097.687500;
   proc_scale["TTZToLLNuNu_2017"]=0.2529*lumi/1837346.4375;
   proc_scale["TTWJetsToLNu_2017"]=0.2001*lumi/1680699.34375;
   proc_scale["TTZToLLNuNu_2016"]=0.2529*lumi_16/1473794.00171;
   proc_scale["TTWJetsToLNu_2016"]=0.2001*lumi_16/1786523.84839;
   proc_scale["TTZToLLNuNu_2018"]=0.2529*lumi_18/3231284.61328;
   proc_scale["TTWJetsToLNu_2018"]=0.2001*lumi_18/1689122.54346;
   proc_scale["DY0J"]=5409.0*lumi/525496638080.000000;
   proc_scale["DY1J"]=937.4*lumi/587697890042.468750;
   proc_scale["DY2J"]=291.0*lumi/137566456384.000000;
   
    TString procname   = argv[3];
    if(*isData=='F') cout <<"process name: "<<procname<<" scale: "<<proc_scale[procname]<<endl;
   
  hmm.Categorization(data, isData, 110, 150,proc_scale[procname],procname);
    
  /*
  if(*isData=='F'){
        cout <<"process name: "<<data<<" scale: "<<scalef<<endl;//<<proc_scale[procname]<<endl;
        hmm.Categorization(data, isData, 95, 160, scalef);
  }
  else hmm.Categorization(data, isData, 95, 160, 1.0);
  
  if(*isData=='T') hmm.EventLoop(data, isData, 1.0);
  else hmm.EventLoop(data, isData, scalef);
  */
  //hmm.GenInfo(data, isData, proc_scale[procname]);
  return 0;
}

void HiggsMuMu::GenInfo(const char *data,const char *isData, double scale)
{  
   

   float muon_mass = 0.1056583745;
   float el_mass = 0.000511;

   float mlo = 110., mhi = 140.;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%5000 == 0) cout <<"jentry: "<<jentry<<endl;
      if(t_mu1>1000 || t_mu2>1000) continue;
      clearTreeVectors();

      //cout <<"event number: "<<jentry<<endl;
      int index_mu1 = -999, index_mu2= -999;
      vector<int> el;
      vector<int> mu;
      vector<int> extralep;
      el.clear();
      mu.clear();
      extralep.clear();
      TLorentzVector mu1, mu2, reco_mu1,reco_mu2;
      
      double evt_wt = 1., tmp_wt = 1.;
      if(*isData=='T'){
         evt_wt=1.;
         tmp_wt=1.;
      }
      else{
         evt_wt=t_genWeight*scale;
         tmp_wt=t_genWeight;
      }

      reco_mu1.SetPtEtaPhiM((*t_Mu_pt)[t_mu1],(*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],muon_mass);
      reco_mu2.SetPtEtaPhiM((*t_Mu_pt)[t_mu2],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],muon_mass);
      reco_dimu_mass = (reco_mu1 + reco_mu2).M();
      reco_mu1_E = reco_mu1.Energy();
      reco_mu2_E = reco_mu2.Energy();
      cos12 = TMath::Cos(reco_mu1.Angle(reco_mu2.Vect()));
      double lepSF = 1.0;
      if((reco_dimu_mass<mlo) || (reco_dimu_mass>mhi)) continue;
      if((*isData=='T') && reco_dimu_mass<130. && reco_dimu_mass>120. ) continue;
      if(*isData=='F'){
            //if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_TRIG)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
            //else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_TRIG)[t_mu2]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
            //2016,2018
            if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
            else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];

            evt_wt*=lepSF;
            //evt_wt*=t_PrefireWeight;
            tmp_wt*=lepSF;
            //tmp_wt*=t_PrefireWeight;
      }
 
      run  =  t_run;
      event = t_event;
      lumi = t_luminosityBlock;
      evt_weight = evt_wt;
      genWeight = t_genWeight;
      reco_mu1_pt = (*t_Mu_pt)[t_mu1];
      reco_mu2_pt = (*t_Mu_pt)[t_mu2];
      reco_mu1_eta = (*t_Mu_eta)[t_mu1];
      reco_mu2_eta = (*t_Mu_eta)[t_mu2];
      reco_mu1_phi = (*t_Mu_phi)[t_mu1];
      reco_mu2_phi = (*t_Mu_phi)[t_mu2];
      //Gen part 
      int n_GenPart = t_GenPart_pdgId->size();
      for( int i=0; i<n_GenPart; i++){
       //cout <<"pdgID: "<<i<<" "<<(*t_GenPart_pdgId)[i]<<" t_GenPart_status "<<(*t_GenPart_status)[i]<<endl;
       if(abs((*t_GenPart_pdgId)[i]) == 13){
        //cout <<"loop first lepton "<<endl;
        //cout <<"i "<<i<<" t_GenPart_status "<<(*t_GenPart_status)[i]<<" genPartIdxMother "<<(*t_GenPart_genPartIdxMother)[i]<<" pt "<<(*t_GenPart_pt)[i]<<" eta "<<(*t_GenPart_eta)[i] <<endl;
        for( int j=i+1; j<n_GenPart; j++){
           if(abs((*t_GenPart_pdgId)[j]) == 13){
             //cout <<"j "<<j<<" t_GenPart_status "<<(*t_GenPart_status)[j]<<" genPartIdxMother "<<(*t_GenPart_genPartIdxMother)[j]<<" pt "<<(*t_GenPart_pt)[j]<<" eta "<<(*t_GenPart_eta)[j] <<endl;
         if((*t_GenPart_status)[i] ==23 && (*t_GenPart_status)[j] ==23){
             mu1.SetPtEtaPhiM((*t_GenPart_pt)[i],(*t_GenPart_eta)[i],(*t_GenPart_phi)[i],muon_mass);
             mu2.SetPtEtaPhiM((*t_GenPart_pt)[j],(*t_GenPart_eta)[j],(*t_GenPart_phi)[j],muon_mass);
             float gen_diMuon_mass = (mu1 + mu2).M();
             //cout <<"Higgs mass "<<diMuon_mass<<endl;
                  if(fabs(gen_dimu_mass-125.0) > fabs(gen_diMuon_mass-125.0)){
                      gen_dimu_mass = gen_diMuon_mass;
                      index_mu1 = i; index_mu2 = j;
                      //cout <<"lepton indexes "<<i<<" "<<j<<endl;
             } 
           }
        }
       }
      }
     }
     if(index_mu1 == -999 && index_mu2 == -999){

        for( int i=0; i<n_GenPart; i++){
           if(abs((*t_GenPart_pdgId)[i]) == 13){ 
            for( int j=i+1; j<n_GenPart; j++){
                if(abs((*t_GenPart_pdgId)[j]) == 13){
                      if((*t_GenPart_status)[i] ==1 && (*t_GenPart_status)[j] ==1){
                           mu1.SetPtEtaPhiM((*t_GenPart_pt)[i],(*t_GenPart_eta)[i],(*t_GenPart_phi)[i],muon_mass);
                           mu2.SetPtEtaPhiM((*t_GenPart_pt)[j],(*t_GenPart_eta)[j],(*t_GenPart_phi)[j],muon_mass);
                           float gen_diMuon_mass = (mu1 + mu2).M();
                           if(fabs(gen_dimu_mass-125.0) > fabs(gen_diMuon_mass-125.0)){
                             gen_dimu_mass = gen_diMuon_mass;
                             index_mu1 = i; index_mu2 = j;
                           }              
                      }
                }
            }
        }
       }
      }
      if(index_mu1 != -999 && index_mu2 != -999){
        float dR1 = DeltaR(t_GenPart_eta->at(index_mu1),t_GenPart_phi->at(index_mu1),reco_mu1_eta,reco_mu1_phi);
        float dR2 = DeltaR(t_GenPart_eta->at(index_mu1),t_GenPart_phi->at(index_mu1),reco_mu2_eta,reco_mu2_phi);  
        if(dR1<dR2){
           gen_mu1_pt = (*t_GenPart_pt)[index_mu1];
           gen_mu2_pt = (*t_GenPart_pt)[index_mu2];
           gen_mu1_eta = (*t_GenPart_eta)[index_mu1];
           gen_mu2_eta = (*t_GenPart_eta)[index_mu2];
           gen_mu1_phi = (*t_GenPart_phi)[index_mu1];
           gen_mu2_phi = (*t_GenPart_phi)[index_mu2];
           mu1.SetPtEtaPhiM((*t_GenPart_pt)[index_mu1],(*t_GenPart_eta)[index_mu1],(*t_GenPart_phi)[index_mu1],muon_mass);
           mu2.SetPtEtaPhiM((*t_GenPart_pt)[index_mu2],(*t_GenPart_eta)[index_mu2],(*t_GenPart_phi)[index_mu2],muon_mass);
           gen_mu1_E = mu1.Energy();
           gen_mu2_E = mu2.Energy();
           mu1_E_err = (reco_mu1_E - gen_mu1_E)/gen_mu1_E;
           mu2_E_err = (reco_mu2_E - gen_mu2_E)/gen_mu2_E;
        }
        else{
           gen_mu1_pt = (*t_GenPart_pt)[index_mu2];
           gen_mu2_pt = (*t_GenPart_pt)[index_mu1];
           gen_mu1_eta = (*t_GenPart_eta)[index_mu2];
           gen_mu2_eta = (*t_GenPart_eta)[index_mu1];
           gen_mu1_phi = (*t_GenPart_phi)[index_mu2];
           gen_mu2_phi = (*t_GenPart_phi)[index_mu1];
           mu1.SetPtEtaPhiM((*t_GenPart_pt)[index_mu2],(*t_GenPart_eta)[index_mu2],(*t_GenPart_phi)[index_mu2],muon_mass);
           mu2.SetPtEtaPhiM((*t_GenPart_pt)[index_mu1],(*t_GenPart_eta)[index_mu1],(*t_GenPart_phi)[index_mu1],muon_mass);
           gen_mu1_E = mu2.Energy();
           gen_mu2_E = mu1.Energy();
           mu1_E_err = (reco_mu1_E - gen_mu2_E)/gen_mu2_E;
           mu2_E_err = (reco_mu2_E - gen_mu1_E)/gen_mu1_E;
        }
      }
      /*
      h_gen_mu1mu2dR->Fill(dR,evt_wt);
      h_gen_mu1mu2dPhi->Fill(dPhi,evt_wt);
      h_gen_diMuon_pt->Fill(diMuon_pt,evt_wt);
      h_gen_diMuon_eta->Fill(diMuon_eta,evt_wt);
      */
      /*
      for( int i=0; i<n_GenPart; i++){
        if(i!=index_mu1 && i!=index_mu2 && (*t_GenPart_pt)[i] > 10. &&( abs((*t_GenPart_pdgId)[i])== 13 || abs((*t_GenPart_pdgId)[i]) == 11) && (*t_GenPart_status)[i] ==1 ){
        bool overlap=false;
        for(int j=0; j<extralep.size();j++){
           float dRlepH = DeltaR(t_GenPart_eta->at(i),t_GenPart_phi->at(i),t_GenPart_eta->at(extralep.at(j)), t_GenPart_phi->at(extralep.at(j)));
           if(dRlepH<0.1){
              overlap = true;
              break;
           }
           if(overlap) break;
           else extralep.push_back(i);
        }

       
        h_gen_extralep1_pt->Fill(t_GenPart_pt->at(i),evt_wt);
        h_gen_extralep1_eta->Fill(t_GenPart_eta->at(i),evt_wt);
        float dRlepH = DeltaR(t_GenPart_eta->at(i),t_GenPart_phi->at(i),(mu1+mu2).Eta(), (mu1+mu2).Phi());
        h_gen_dRlepH->Fill(dRlepH,evt_wt);
  
      }
     
     }
     */
     //h_gen_extralep->Fill(extralep.size(),evt_wt);
     h_gen_diMuon_m->Fill(gen_dimu_mass);
     cattree->Fill();
     //cout <<"gen Higgs mass: "<<gen_dimu_mass<<endl;
  }
}

void HiggsMuMu::Categorization(const char *data,const char *isData, float mlo, float mhi, double scale, TString procname)
{
   if (fChain == 0) return;
   cout <<"mass window: "<<mlo<<" - "<<mhi<<endl;
   float muon_mass = 0.1056583745;
   float el_mass = 0.000511;
    
   //Pisa VBF BDT
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    
    reader->AddVariable( "ll_mass", &ll_mass );
    reader->AddVariable( "MqqLog", &MqqLog );
    reader->AddVariable( "mumujj_pt", &mumujj_pt );
    reader->AddVariable( "DeltaEtaQQ", &DeltaEtaQQ );
    reader->AddVariable( "softActivityEWK_njets5",&softActivityEWK_njets5);
    reader->AddVariable( "ll_zstar", &ll_zstar );
    reader->AddVariable( "ll_pt", &ll_pt );
    reader->AddVariable( "theta2", &theta2 );
    reader->AddVariable( "impulsoZ", &impulsoZ );
    reader->AddVariable( "maxAbsEta", &maxAbsEta );
    reader->AddVariable( "qgl_2qAtanh", &qgl_2qAtanh );
    TString methodName = "Classification_BDTG.variables__ll_mass__MqqLog__mumujj_pt__DeltaEtaQQ__softActivityEWK_njets5__ll_zstar__ll_pt__theta2__impulsoZ__maxAbsEta__qgl_2qAtanh";
    TString weightfile = "Classification_BDTG.variables__ll_mass__MqqLog__mumujj_pt__DeltaEtaQQ__softActivityEWK_njets5__ll_zstar__ll_pt__theta2__impulsoZ__maxAbsEta__qgl_2qAtanh.xml";
    reader->BookMVA( methodName, weightfile );
    
    
   TMVA::Reader *reader2016 = new TMVA::Reader( "!Color:!Silent" );

   reader2016->AddVariable( "dimu_pt", &dimu_pt );
   reader2016->AddVariable( "dimu_eta", &dimu_eta );
      reader2016->AddVariable( "dimu_abs_dEta", &dimu_abs_dEta );
      reader2016->AddVariable( "dimu_abs_dPhi", &dimu_abs_dPhi );
      reader2016->AddVariable( "jet1_eta", &jet1_eta);
      reader2016->AddVariable( "jet2_eta", &jet2_eta);
      reader2016->AddVariable( "dijet1_mass", &dijet1_mass);
      reader2016->AddVariable( "dijet2_mass", &dijet2_mass);
      reader2016->AddVariable( "dijet1_abs_dEta", &dijet1_abs_dEta);
      reader2016->AddVariable( "dijet2_abs_dEta", &dijet2_abs_dEta);
      reader2016->AddVariable( "nJetsCent", &nJetsCent);
      reader2016->AddVariable( "nJetsFwd", &nJetsFwd);
      reader2016->AddVariable( "nBMed", &nBMed);
      reader2016->AddVariable( "MET", &MET);
      reader2016->AddSpectator( "samp_ID", &samp_ID);
      reader2016->AddSpectator( "samp_wgt", &samp_wgt);
      reader2016->AddSpectator( "res_wgt", &res_wgt);
      reader2016->AddSpectator( "LHE_HT", &LHE_HT);
      reader2016->AddSpectator( "dimu_mass_Roch", &dimu_mass_Roch);
      reader2016->AddSpectator( "BASE_cat", &BASE_cat);
      cout <<"p6" <<endl;
      reader->AddSpectator( "BASE_cat", &BASE_cat);
      TString methodName2016 = "f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1";
      TString weightfile2016 = "f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml";
      reader2016->BookMVA( methodName2016, weightfile2016 );
    
   cout <<"p7" <<endl;
   TMVA::Reader *reader2j = new TMVA::Reader( "!Color:!Silent" );
   reader2j->AddVariable( "hmmpt", &dimu_pt );
   reader2j->AddVariable( "hmmrap", &hmmrap );
      reader2j->AddVariable( "hmmthetacs", &dtheta );
      reader2j->AddVariable( "hmmphics", &dphi );
      reader2j->AddVariable( "j1pt", &jet1_pt);
      reader2j->AddVariable( "j1eta", &jet1_eta);
      reader2j->AddVariable( "j2pt", &jet2_pt);
      reader2j->AddVariable( "detajj", &dijet1_abs_dEta);
      reader2j->AddVariable( "dphijj", &dphijj);
      reader2j->AddVariable( "mjj", &dijet_mass);
      reader2j->AddVariable( "met", &MET);
      reader2j->AddVariable( "zepen", &zepen);
      reader2j->AddVariable( "njets", &njets);
      reader2j->AddVariable( "drmj", &drmj); 
      reader2j->AddVariable( "m1ptOverMass", &m1ptOverMass);
      reader2j->AddVariable( "m2ptOverMass", &m2ptOverMass);
      reader2j->AddVariable( "m1eta", &m1eta);
      reader2j->AddVariable( "m2eta", &m2eta);
      reader2j->AddSpectator( "hmerr", &hmerr);
      reader2j->AddSpectator( "weight", &weight);
      reader2j->AddSpectator( "hmass", &hmass);
      reader2j->AddSpectator( "nbjets", &nBMed);
      reader2j->AddSpectator( "bdtucsd_inclusive", &bdtucsd_inclusive);
      reader2j->AddSpectator( "bdtucsd_01jet", &bdtucsd_01jet);
      reader2j->AddSpectator( "bdtucsd_2jet", &bdtucsd_2jet);
      TString methodName2j = "TMVAClassification_BDTG.weights.2jet_bveto";
      TString weightfile2j = "TMVAClassification_BDTG.weights.2jet_bveto.xml";
      reader2j->BookMVA( methodName2j, weightfile2j );

   TMVA::Reader *reader01j = new TMVA::Reader( "!Color:!Silent" );
   reader01j->AddVariable( "hmmpt", &dimu_pt );
   reader01j->AddVariable( "hmmrap", &hmmrap );
      reader01j->AddVariable( "hmmthetacs", &dtheta );
      reader01j->AddVariable( "hmmphics", &dphi );
      reader01j->AddVariable( "j1pt", &jet1_pt);
      reader01j->AddVariable( "j1eta", &jet1_eta);
      reader01j->AddVariable( "met", &MET);
      reader01j->AddVariable( "nbjets", &nBMed);
      reader01j->AddVariable( "drmj", &drmj);
      reader01j->AddVariable( "njets", &njets);
      reader01j->AddVariable( "m1ptOverMass", &m1ptOverMass);
      reader01j->AddVariable( "m2ptOverMass", &m2ptOverMass);
      reader01j->AddVariable( "m1eta", &m1eta);
      reader01j->AddVariable( "m2eta", &m2eta);
      reader01j->AddSpectator( "hmerr", &hmerr);
      reader01j->AddSpectator( "weight", &weight);
      reader01j->AddSpectator( "hmass", &hmass);
      reader01j->AddSpectator( "bdtucsd_inclusive", &bdtucsd_inclusive);
      reader01j->AddSpectator( "bdtucsd_01jet", &bdtucsd_01jet);
      reader01j->AddSpectator( "bdtucsd_2jet", &bdtucsd_2jet);
      TString methodName01j = "TMVAClassification_BDTG.weights.01jet";
      TString weightfile01j = "TMVAClassification_BDTG.weights.01jet.xml";
      reader01j->BookMVA( methodName01j, weightfile01j );


   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 2;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      clearTreeVectors();
      if(jentry%50000 == 0) cout <<"jentry: "<<jentry<<endl;
      if(t_mu1>1000 || t_mu2>1000) continue;
       vector <int > genJet_idx;
       genJet_idx.clear();
       //cout<<procname<<endl;
       TString s1 = "DYJetsToLL_M105To160_incl_2017" ;
       TString s2 = "DYJetsToLL_M105To160_incl_2016" ;
       TString s3 = "DYJetsToLL_M105To160_incl_2018" ;
       TString s4 = "DYJetsToLL_2017" ;
       TString s5 = "DYJetsToLL_2018" ;
       TString s6 = "DYJetsToLL_2016" ;
       TString s7 = "DYJetsToLL_VBFfilter_2017" ;
       TString s8 = "DYJetsToLL_VBFfilter_2018" ;
       TString s9 = "DYJetsToLL_VBFfilter_2016" ;
       if (procname == s1 || procname == s2 || procname == s3){
           //cout<<"Entering!!!!\n";
           if(t_GenJet_pt->size() >= 2){
               for(int i=0; i<t_GenJet_pt->size(); i++){
                   bool foundGenLep=false;
                   for(int k=0; k<t_GenPart_pt->size(); k++){
                       
                       if(abs(t_GenPart_pdgId->at(k))>=11 && abs(t_GenPart_pdgId->at(k))<=16){
                           if(DeltaR((*t_GenJet_eta)[i],(*t_GenJet_phi)[i],(*t_GenPart_eta)[k],(*t_GenPart_phi)[k])<0.3){foundGenLep=true;break;}
                       }
                       
                   }
                   if(!foundGenLep)genJet_idx.push_back(i);
               }
           }
           
       }
      int index_mu1 = -999, index_mu2= -999;
      double evt_wt = 1., evt_wt_Up=1,evt_wt_Down=1, tmp_wt = 1.;
      if(*isData=='F'){
         evt_wt=t_genWeight*scale;
         tmp_wt=t_genWeight;
      }
      vector<int> el;
      vector<int> mu;
      el.clear();
      mu.clear();

      if(t_Mu_pt->size() > 2){
         for(int i=0; i<t_Mu_pt->size(); i++){
             if(i!=t_mu1 && i!=t_mu2){
                 if(t_Mu_pt->at(i)>20. && t_Mu_miniPFRelIso_all->at(i)<0.4) mu.push_back(i);
             }
         }  
      }

      if(t_El_pt->size() > 0){
         for(int i=0; i<t_El_pt->size(); i++){
            if(t_El_miniPFRelIso03_all->at(i)<0.4 &&  t_El_pt->at(i)>20. && fabs(t_El_eta->at(i))<2.5 && t_El_cutBased->at(i)==2){
                el.push_back(i);
            }
         }  
      }

      TLorentzVector dimu, reco_mu1, reco_mu2, mu1, mu2;
      reco_mu1.SetPtEtaPhiM((*t_Mu_pt)[t_mu1],(*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],muon_mass);
      reco_mu2.SetPtEtaPhiM((*t_Mu_pt)[t_mu2],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],muon_mass);
      dimu = reco_mu1 + reco_mu2;
      float diMuon_mass = dimu.M();
      //SR: 120-130, sideband: 100-120, 130-150
      if((diMuon_mass>mlo && diMuon_mass<mhi) || (diMuon_mass>76. && diMuon_mass<106.)){
          if((*isData=='T') && diMuon_mass<130. && diMuon_mass>120. ) continue;
          float diMuon_pt = dimu.Pt();
          float diMuon_eta = dimu.Eta();
          float diMuon_phi = dimu.Phi();
          float dR= DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2]);
          float dEta = (*t_Mu_eta)[t_mu1] - (*t_Mu_eta)[t_mu2];
          float dPhi= DeltaPhi((*t_Mu_phi)[t_mu2],(*t_Mu_phi)[t_mu1]);
          m1eta = (*t_Mu_eta)[t_mu1];
          m2eta = (*t_Mu_eta)[t_mu2];
          m1ptOverMass = (*t_Mu_pt)[t_mu1]/diMuon_mass;
          m2ptOverMass = (*t_Mu_pt)[t_mu2]/diMuon_mass;
          std::pair<double,double>  angles = CSAngles(reco_mu1,reco_mu2,t_Mu_charge->at(t_mu1));
          //SR: 120-130, sideband: 100-120, 130-150
       if(*isData=='F'){
         double lepSF = 1.0;
         if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_TRIG)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
         else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_TRIG)[t_mu2]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
         evt_wt*=lepSF;
         evt_wt*=t_puWeight;
         evt_wt*=t_PrefireWeight;
      
           if(procname == s2  || procname == s3  || procname == s5 || procname == s8 || procname == s1  ||  procname == s4   || procname == s7  ){
               if(diMuon_pt<10.)evt_wt*=0.910385;
               else if(diMuon_pt<20.)evt_wt*=1.13543;
               else if(diMuon_pt<30.)evt_wt*=1.10441;
               else if(diMuon_pt<40.)evt_wt*=1.01315;
               else if(diMuon_pt<50.)evt_wt*=0.982598;
               else if(diMuon_pt<60.)evt_wt*=0.980697;
               else if(diMuon_pt<70.)evt_wt*= 0.972673;
               else if(diMuon_pt<80.)evt_wt*=0.972325;
               else if(diMuon_pt<100.)evt_wt*=0.966127;
               else if(diMuon_pt<150.)evt_wt*=0.953262;
               else if(diMuon_pt<200.)evt_wt*= 0.933403;
               else if(diMuon_pt<1000.)evt_wt*=0.904518;
           }
           if(procname == s6 || procname == s9){
               if(diMuon_pt<10.)evt_wt*=1.05817;
               else if(diMuon_pt<20.)evt_wt*=0.994488;
               else if(diMuon_pt<30.)evt_wt*=0.930056;
               else if(diMuon_pt<40.)evt_wt*=0.925206;
               else if(diMuon_pt<50.)evt_wt*=0.946403;
               else if(diMuon_pt<60.)evt_wt*=0.962136;
               else if(diMuon_pt<70.)evt_wt*= 0.965316;
               else if(diMuon_pt<80.)evt_wt*=0.978209;
               else if(diMuon_pt<100.)evt_wt*=0.988761;
               else if(diMuon_pt<150.)evt_wt*=0.982497;
               else if(diMuon_pt<200.)evt_wt*= 0.971749;
               else if(diMuon_pt<1000.)evt_wt*=0.914429;
           }
           
           evt_wt_Up=evt_wt*t_puWeightUp*t_PrefireWeight_Up;
           evt_wt_Down=evt_wt*t_puWeightDown*t_PrefireWeight_Down;
           
            //from Gen info
            int n_GenPart = t_GenPart_pdgId->size();
            for( int i=0; i<n_GenPart; i++){
             if(abs((*t_GenPart_pdgId)[i]) == 13){
              for( int j=i+1; j<n_GenPart; j++){
                 if(abs((*t_GenPart_pdgId)[j]) == 13){
                    if((*t_GenPart_status)[i] ==23 && (*t_GenPart_status)[j] ==23){
                        mu1.SetPtEtaPhiM((*t_GenPart_pt)[i],(*t_GenPart_eta)[i],(*t_GenPart_phi)[i],muon_mass);
                        mu2.SetPtEtaPhiM((*t_GenPart_pt)[j],(*t_GenPart_eta)[j],(*t_GenPart_phi)[j],muon_mass);
                        float gen_diMuon_mass = (mu1 + mu2).M();
                        if(fabs(gen_dimu_mass-125.0) > fabs(gen_diMuon_mass-125.0)){
                               gen_dimu_mass = gen_diMuon_mass;
                               index_mu1 = i; index_mu2 = j;
                        }
                     }
                 }
             }
            }
          }
          if(index_mu1 == -999 && index_mu2 == -999){

           for( int i=0; i<n_GenPart; i++){
            if(abs((*t_GenPart_pdgId)[i]) == 13){
              for( int j=i+1; j<n_GenPart; j++){
                if(abs((*t_GenPart_pdgId)[j]) == 13){
                      if((*t_GenPart_status)[i] ==1 && (*t_GenPart_status)[j] ==1){
                           mu1.SetPtEtaPhiM((*t_GenPart_pt)[i],(*t_GenPart_eta)[i],(*t_GenPart_phi)[i],muon_mass);
                           mu2.SetPtEtaPhiM((*t_GenPart_pt)[j],(*t_GenPart_eta)[j],(*t_GenPart_phi)[j],muon_mass);
                           float gen_diMuon_mass = (mu1 + mu2).M();
                           if(fabs(gen_dimu_mass-125.0) > fabs(gen_diMuon_mass-125.0)){
                             gen_dimu_mass = gen_diMuon_mass;
                             index_mu1 = i; index_mu2 = j;
                           }
                      }
                }
              }
            }
          }
         }
      if(index_mu1 != -999 && index_mu2 != -999){
        float dR1 = DeltaR(t_GenPart_eta->at(index_mu1),t_GenPart_phi->at(index_mu1),reco_mu1_eta,reco_mu1_phi);
        float dR2 = DeltaR(t_GenPart_eta->at(index_mu1),t_GenPart_phi->at(index_mu1),reco_mu2_eta,reco_mu2_phi);
        if(dR1<dR2){
           gen_mu1_pt = (*t_GenPart_pt)[index_mu1];
           gen_mu2_pt = (*t_GenPart_pt)[index_mu2];
           gen_mu1_eta = (*t_GenPart_eta)[index_mu1];
           gen_mu2_eta = (*t_GenPart_eta)[index_mu2];
           gen_mu1_phi = (*t_GenPart_phi)[index_mu1];
           gen_mu2_phi = (*t_GenPart_phi)[index_mu2];
           mu1.SetPtEtaPhiM((*t_GenPart_pt)[index_mu1],(*t_GenPart_eta)[index_mu1],(*t_GenPart_phi)[index_mu1],muon_mass);
           mu2.SetPtEtaPhiM((*t_GenPart_pt)[index_mu2],(*t_GenPart_eta)[index_mu2],(*t_GenPart_phi)[index_mu2],muon_mass);
           gen_mu1_E = mu1.Energy();
           gen_mu2_E = mu2.Energy();
           mu1_E_err = (reco_mu1_E - gen_mu1_E)/gen_mu1_E;
           mu2_E_err = (reco_mu2_E - gen_mu2_E)/gen_mu2_E;
        }
        else{
           gen_mu1_pt = (*t_GenPart_pt)[index_mu2];
           gen_mu2_pt = (*t_GenPart_pt)[index_mu1];
           gen_mu1_eta = (*t_GenPart_eta)[index_mu2];
           gen_mu2_eta = (*t_GenPart_eta)[index_mu1];
           gen_mu1_phi = (*t_GenPart_phi)[index_mu2];
           gen_mu2_phi = (*t_GenPart_phi)[index_mu1];
           mu1.SetPtEtaPhiM((*t_GenPart_pt)[index_mu2],(*t_GenPart_eta)[index_mu2],(*t_GenPart_phi)[index_mu2],muon_mass);
           mu2.SetPtEtaPhiM((*t_GenPart_pt)[index_mu1],(*t_GenPart_eta)[index_mu1],(*t_GenPart_phi)[index_mu1],muon_mass);
           gen_mu1_E = mu2.Energy();
           gen_mu2_E = mu1.Energy();
           mu1_E_err = (reco_mu1_E - gen_mu2_E)/gen_mu2_E;
           mu2_E_err = (reco_mu2_E - gen_mu1_E)/gen_mu1_E;
        }
        reco_gen_dimu_mass = diMuon_mass - gen_dimu_mass;
        r_reco_gen_dimu_mass = (diMuon_mass - gen_dimu_mass)/diMuon_mass;
        fabs_reco_gen_dimu_mass = fabs(diMuon_mass - gen_dimu_mass);
        r_fabs_reco_gen_dimu_mass = fabs(diMuon_mass - gen_dimu_mass)/diMuon_mass;
         
      }

          //from Gen info
          reco_mu1_pt = (*t_Mu_pt)[t_mu1];
          reco_mu2_pt = (*t_Mu_pt)[t_mu2];
          reco_mu1_eta = (*t_Mu_eta)[t_mu1];
          reco_mu2_eta = (*t_Mu_eta)[t_mu2];
          reco_mu1_phi = (*t_Mu_phi)[t_mu1];
          reco_mu2_phi = (*t_Mu_phi)[t_mu2];

          max_reco_mu_eta = reco_mu1_eta;
          min_reco_mu_eta = reco_mu2_eta;
          if(fabs(reco_mu1_eta) < fabs(reco_mu2_eta)){
             max_reco_mu_eta = reco_mu2_eta;
             min_reco_mu_eta = reco_mu1_eta;
          }
       }
          
          //inclusive BDT
          ///vars
          dimu_pt = diMuon_pt;
          dimu_eta = diMuon_eta;
          hmmrap = dimu.Rapidity();
          dimu_abs_dEta = fabs(dEta);
          dimu_abs_dPhi = fabs(dPhi);
          njets = (float)t_nJet;
          
          if(t_nJet>0){
            drmj = 1000.;
            for(int i=0; i<t_nJet;i++){
                  float tmp = DeltaR((*t_Jet_eta)[i],(*t_Jet_phi)[i],(*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1] );
                  if(tmp < drmj) drmj = tmp;
                  tmp = DeltaR((*t_Jet_eta)[i],(*t_Jet_phi)[i],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2] );
                  if(tmp < drmj) drmj = tmp;
            }
              if(drmj == 1000){ drmj = 0; cout <<"drmj error"<< endl;}
                
            jet1_pt = (*t_Jet_pt)[0];
            jet1_eta = (*t_Jet_eta)[0];
            jet1_qgl = (*t_Jet_qgl)[0];
            TLorentzVector jet1;
            jet1.SetPtEtaPhiM(jet1_pt, jet1_eta, (*t_Jet_phi)[0], (*t_Jet_mass)[0]);
            dRJ1Higgs = DeltaR(jet1_eta,(*t_Jet_phi)[0],diMuon_pt, diMuon_phi);
            m_mmj = (jet1 +dimu).M();
            pt_mmj = (jet1 +dimu).Pt();
                
            if(t_nJet>1){
              jet2_pt = (*t_Jet_pt)[1];
              jet2_eta = (*t_Jet_eta)[1];
              jet2_qgl = (*t_Jet_qgl)[1];
              dRJ2Higgs = DeltaR(jet2_eta,(*t_Jet_phi)[1],diMuon_pt, diMuon_phi);
           }

      if(t_nJet>=2){
         dijet_mass = t_diJet_mass;
         for(int k=0; k<t_nJet; k++){
           float j1_pt = (*t_Jet_pt)[k];
           float j1_eta = (*t_Jet_eta)[k];
           float j1_phi = (*t_Jet_phi)[k];
           float j1_mass = (*t_Jet_mass)[k];
           for(int j=k+1; j<t_nJet; j++){
             float j2_pt = (*t_Jet_pt)[j];
             float j2_eta = (*t_Jet_eta)[j];
             float j2_phi = (*t_Jet_phi)[j];
             float j2_mass = (*t_Jet_mass)[j];
             TLorentzVector j1,j2, jj;
             j1.SetPtEtaPhiM(j1_pt, j1_eta, j1_phi, j1_mass);
             j2.SetPtEtaPhiM(j2_pt, j2_eta, j2_phi, j2_mass);
             jj=j1+j2;
             if(dijet1_mass<jj.M()){
                   zepen = dimu.Rapidity() - 0.5*(reco_mu1.Rapidity()+reco_mu2.Rapidity());
                   dphijj = DeltaPhi(j1_phi,j2_phi); 
                   dijet1_mass = jj.M();
                   dijet1_abs_dEta = fabs(j1_eta-j2_eta);
                   m_mmjj = (jj +dimu).M();
                   pT_mmjj = (jj +dimu).Pt();
             }
             else if(dijet2_mass<jj.M()){
                   dijet2_mass = jj.M();
                   dijet2_abs_dEta = fabs(j1_eta-j2_eta);
             }
           }
         }
        }
      }//close t_nJet>0
              
      for(int k=0; k<t_nJet; k++){
          float jet_eta_tmp = (*t_Jet_eta)[k];
          if(fabs(jet_eta_tmp)>2.4) nJetsFwd = nJetsFwd + 1.0;
          else nJetsCent = nJetsCent + 1.0;
      }

      nBMed = (float)t_nbJet;
      MET = t_MET_pt;
      MET_phi = t_MET_phi;

          double binv_inc = catyield->GetBinContent(20);
          binv_inc = binv_inc + evt_wt;
          catyield->SetBinContent(20,binv_inc);
          h_diMuon_mass_cat->Fill(diMuon_mass,evt_wt);
        
          run  =  t_run;
          event = t_event;
          lumi = t_luminosityBlock;
          genWeight = evt_wt;
          genWeight_Up = evt_wt_Up;
          genWeight_Down = evt_wt_Down;
          evt_weight = evt_wt;
          Higgs_mass = diMuon_mass;
          hmass = diMuon_mass;
          dtheta = angles.first;
          dphi = angles.second;
          Higgs_pt = diMuon_pt;
          Higgs_eta = diMuon_eta;
          dRmm = dR;
          dEtamm = dEta;
          dPhimm = dPhi;

          BDT_incl = reader2016->EvaluateMVA( methodName2016 );
          BDT_2j = reader2j->EvaluateMVA( methodName2j );
          BDT_01j = reader01j->EvaluateMVA( methodName01j );
          h_BDT->Fill(BDT_incl,evt_wt);
 
          //ttH
          if(t_nbJet>0){
	    //if(*isData=='F')evt_wt*=(*t_bJet_SF)[0];
              if(el.size()> 0  || mu.size()>0){
                  cat_index = 1;
                  double binv = catyield->GetBinContent(1);
                  binv = binv + evt_wt;
                  catyield->SetBinContent(1,binv);
                  if(*isData=='F'){
                     h_diMuon_mass_ttHLep->Fill(diMuon_mass,evt_wt);
                  }
		  else{
                     if(diMuon_mass<120. || diMuon_mass>130.){
                         h_diMuon_mass_ttHLep->Fill(diMuon_mass,evt_wt);
                     }
                  }
                  //h_diMuon_mass_110To140_ttHLep->Fill(diMuon_mass,evt_wt);
		  //else if(*isData=='F')h_diMuon_mass_110To140_ttHLep->Fill(diMuon_mass,evt_wt);
		  if(diMuon_mass<120. || diMuon_mass>130.){
                 
		    h_mu1mu2dR_ttHLep->Fill(dR,evt_wt);
		    h_mu1mu2dPhi_ttHLep->Fill(dPhi,evt_wt);
		    h_mu1mu2dEta_ttHLep->Fill(dEta,evt_wt);
		    h_diMuon_pt_ttHLep->Fill(diMuon_pt,evt_wt);
		    h_diMuon_eta_ttHLep->Fill(diMuon_eta,evt_wt);
		    h_diMuon_phi_ttHLep->Fill(diMuon_phi,evt_wt);
		    //h_Nbjet_ttHLep->Fill(t_nbJet);
		    h_leading_bJet_pt_ttHLep->Fill((*t_bJet_pt)[0],evt_wt);
		    h_leading_bJet_eta_ttHLep->Fill((*t_bJet_eta)[0],evt_wt);
		    h_leading_bJet_phi_ttHLep->Fill((*t_bJet_phi)[0],evt_wt);
		  }
	      }
	  }
          else if(t_nbJet>0 && t_nJet>1){
                  cat_index = 2;
                  double binv = catyield->GetBinContent(2);
                  binv = binv + evt_wt;
                  catyield->SetBinContent(2,binv);
                  if(*isData=='F'){
                     h_diMuon_mass_ttHHad->Fill(diMuon_mass,evt_wt);
                  }
                  else{
                     if(diMuon_mass<120. || diMuon_mass>130.){
                         h_diMuon_mass_ttHHad->Fill(diMuon_mass,evt_wt);
                     }
                  }

		  //if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_ttHHad->Fill(diMuon_mass,evt_wt);
		  //else if(*isData=='F')h_diMuon_mass_110To140_ttHHad->Fill(diMuon_mass,evt_wt);

		  if(diMuon_mass<120. || diMuon_mass>130.){
		    h_mu1mu2dR_ttHHad->Fill(dR,evt_wt);
		    h_mu1mu2dPhi_ttHHad->Fill(dPhi,evt_wt);
		    h_mu1mu2dEta_ttHHad->Fill(dEta,evt_wt);
		    h_diMuon_pt_ttHHad->Fill(diMuon_pt,evt_wt);
		    h_diMuon_eta_ttHHad->Fill(diMuon_eta,evt_wt);
		    h_diMuon_phi_ttHHad->Fill(diMuon_phi,evt_wt);
		    h_leading_bJet_pt_ttHHad->Fill((*t_bJet_pt)[0],evt_wt);
		    h_leading_bJet_eta_ttHHad->Fill((*t_bJet_eta)[0],evt_wt);
		    h_leading_bJet_phi_ttHHad->Fill((*t_bJet_phi)[0],evt_wt);
		  }
        }
             /*
              else{
                  cat_index = 3;
                  double binv = catyield->GetBinContent(3);
                  binv = binv + evt_wt;
                  catyield->SetBinContent(3,binv);
                  if(*isData=='F'){
                     h_diMuon_mass_ttHLoose->Fill(diMuon_mass,evt_wt);
                  }
                  else{
                     if(diMuon_mass<120. || diMuon_mass>130.){
                         h_diMuon_mass_ttHLoose->Fill(diMuon_mass,evt_wt);
                     }
                  }
                

		  if(diMuon_mass<120. || diMuon_mass>130.){
		    h_mu1mu2dR_ttHLoose->Fill(dR,evt_wt);
		    h_mu1mu2dPhi_ttHLoose->Fill(dPhi,evt_wt);
		    h_mu1mu2dEta_ttHLoose->Fill(dEta,evt_wt);
		    h_diMuon_pt_ttHLoose->Fill(diMuon_pt,evt_wt);
		    h_diMuon_eta_ttHLoose->Fill(diMuon_eta,evt_wt);
		    h_diMuon_phi_ttHLoose->Fill(diMuon_phi,evt_wt);
		    h_leading_bJet_pt_ttHLoose->Fill((*t_bJet_pt)[0],evt_wt);
		    h_leading_bJet_eta_ttHLoose->Fill((*t_bJet_eta)[0],evt_wt);
		    h_leading_bJet_phi_ttHLoose->Fill((*t_bJet_phi)[0],evt_wt);
		  }
              }
        */
          //ZH mm
          else if(mu.size()>0){
              cat_index = 3;
              double binv = catyield->GetBinContent(3);
              binv = binv + evt_wt;
              catyield->SetBinContent(3,binv);
          }
          //ZH ee
          else if(el.size()>0){
              cat_index = 4;
              double binv = catyield->GetBinContent(4);
              binv = binv + evt_wt;
              catyield->SetBinContent(4,binv);
          }
          //VBF
          else if(t_diJet_mass>400. && ((*t_Jet_qgl)[1]!=-1 && (*t_Jet_qgl)[0]!=-1) && (*t_Jet_pt)[0]>30.){// && (*t_Jet_pt)[1]>30.){
              if(genJet_idx.size()>=2){
                  TLorentzVector di_gJ, gJ1,gJ2;
                  gJ1.SetPtEtaPhiM((*t_GenJet_pt)[genJet_idx[0]],(*t_GenJet_eta)[genJet_idx[0]],(*t_GenJet_phi)[genJet_idx[0]],(*t_GenJet_mass)[genJet_idx[0]]);
                  gJ2.SetPtEtaPhiM((*t_GenJet_pt)[genJet_idx[1]],(*t_GenJet_eta)[genJet_idx[1]],(*t_GenJet_phi)[genJet_idx[1]],(*t_GenJet_mass)[genJet_idx[1]]);
                  di_gJ=gJ1+gJ2;
                  //cout<<"GenJet dimass : "<<di_gJ.M()<<endl;
                  if(di_gJ.M()>350. && (procname == s1 || procname == s2 || procname == s3)){continue;}
                  //cout<<"GenJet dimass : "<<di_gJ.M()<<endl;
              }
              
              cat_index = 5;
              double binv = catyield->GetBinContent(5);
              binv = binv + evt_wt;
              catyield->SetBinContent(5,binv);
              ll_mass=diMuon_mass; //BDT var
              MqqLog=log(t_diJet_mass); //BDT_var
              if(diMuon_mass>108.){
                  if(*isData=='F'){
                      h_diMuon_mass_VBF->Fill(diMuon_mass,evt_wt);
                      h_diMuon_mass_Up_VBF->Fill(diMuon_mass,evt_wt_Up);
                      h_diMuon_mass_Down_VBF->Fill(diMuon_mass,evt_wt_Down);
                  }
                  if(diMuon_mass<120. || diMuon_mass>130.){
                      h_diMuon_mass_110To140_VBF->Fill(diMuon_mass,evt_wt);
                      h_diMuon_mass_110To140_Up_VBF->Fill(diMuon_mass,evt_wt_Up);
                      h_diMuon_mass_110To140_Down_VBF->Fill(diMuon_mass,evt_wt_Down);
                  }
                  else if(*isData=='F'){
                      h_diMuon_mass_110To140_VBF->Fill(diMuon_mass,evt_wt);
                      h_diMuon_mass_110To140_Up_VBF->Fill(diMuon_mass,evt_wt_Up);
                      h_diMuon_mass_110To140_Down_VBF->Fill(diMuon_mass,evt_wt_Down);
                  }
                  if(diMuon_mass<120. || diMuon_mass>130.){
                      h_diMuon_pt_VBF->Fill(diMuon_pt,evt_wt);
                      h_diMuon_eta_VBF->Fill(diMuon_eta,evt_wt);
                      h_diMuon_phi_VBF->Fill(diMuon_phi,evt_wt);
                      h_mu1mu2dR_VBF->Fill(dR,evt_wt);
                      h_mu1mu2dPhi_VBF->Fill(dPhi,evt_wt);
                      h_mu1mu2dEta_VBF->Fill(dEta,evt_wt);
                      h_softJet5_VBF->Fill(t_SoftActivityJetNjets5,evt_wt);
                      h_dijet_pt_VBF->Fill(t_diJet_pt,evt_wt);
                      h_dijet_eta_VBF->Fill(t_diJet_eta,evt_wt);
                      h_dijet_phi_VBF->Fill(t_diJet_phi,evt_wt);
                      h_Mjj_VBF->Fill(t_diJet_mass,evt_wt);
                      h_dijet_dEta_VBF->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt);
                      
                      h_diMuon_pt_Up_VBF->Fill(diMuon_pt,evt_wt_Up);
                      h_diMuon_eta_Up_VBF->Fill(diMuon_eta,evt_wt_Up);
                      h_diMuon_phi_Up_VBF->Fill(diMuon_phi,evt_wt_Up);
                      h_mu1mu2dR_Up_VBF->Fill(dR,evt_wt_Up);
                      h_mu1mu2dPhi_Up_VBF->Fill(dPhi,evt_wt_Up);
                      h_mu1mu2dEta_Up_VBF->Fill(dEta,evt_wt_Up);
                      h_softJet5_Up_VBF->Fill(t_SoftActivityJetNjets5,evt_wt_Up);
                      h_dijet_pt_Up_VBF->Fill(t_diJet_pt,evt_wt_Up);
                      h_dijet_eta_Up_VBF->Fill(t_diJet_eta,evt_wt_Up);
                      h_dijet_phi_Up_VBF->Fill(t_diJet_phi,evt_wt_Up);
                      h_Mjj_Up_VBF->Fill(t_diJet_mass,evt_wt_Up);
                      h_dijet_dEta_Up_VBF->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt_Up);
                      
                      h_diMuon_pt_Down_VBF->Fill(diMuon_pt,evt_wt_Down);
                      h_diMuon_eta_Down_VBF->Fill(diMuon_eta,evt_wt_Down);
                      h_diMuon_phi_Down_VBF->Fill(diMuon_phi,evt_wt_Down);
                      h_mu1mu2dR_Down_VBF->Fill(dR,evt_wt_Down);
                      h_mu1mu2dPhi_Down_VBF->Fill(dPhi,evt_wt_Down);
                      h_mu1mu2dEta_Down_VBF->Fill(dEta,evt_wt_Down);
                      h_softJet5_Down_VBF->Fill(t_SoftActivityJetNjets5,evt_wt_Down);
                      h_dijet_pt_Down_VBF->Fill(t_diJet_pt,evt_wt_Down);
                      h_dijet_eta_Down_VBF->Fill(t_diJet_eta,evt_wt_Down);
                      h_dijet_phi_Down_VBF->Fill(t_diJet_phi,evt_wt_Down);
                      h_Mjj_Down_VBF->Fill(t_diJet_mass,evt_wt_Down);
                      h_dijet_dEta_Down_VBF->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt_Down);
                  }
              }
              else{
                  h_diMuon_mass_Z->Fill(diMuon_mass,evt_wt);
                  h_diMuon_pt_Z->Fill(diMuon_pt,evt_wt);
                  h_diMuon_eta_Z->Fill(diMuon_eta,evt_wt);
                  h_diMuon_phi_Z->Fill(diMuon_phi,evt_wt);
                  h_mu1mu2dR_Z->Fill(dR,evt_wt);
                  h_mu1mu2dPhi_Z->Fill(dPhi,evt_wt);
                  h_mu1mu2dEta_Z->Fill(dEta,evt_wt);
                  h_softJet5_Z->Fill(t_SoftActivityJetNjets5,evt_wt);
                  h_dijet_pt_Z->Fill(t_diJet_pt,evt_wt);
                  h_dijet_eta_Z->Fill(t_diJet_eta,evt_wt);
                  h_dijet_phi_Z->Fill(t_diJet_phi,evt_wt);
                  h_Mjj_Z->Fill(t_diJet_mass,evt_wt);
                  h_dijet_dEta_Z->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt);
                  
                  
                  h_diMuon_mass_Up_Z->Fill(diMuon_mass,evt_wt_Up);
                  h_diMuon_pt_Up_Z->Fill(diMuon_pt,evt_wt_Up);
                  h_diMuon_eta_Up_Z->Fill(diMuon_eta,evt_wt_Up);
                  h_diMuon_phi_Up_Z->Fill(diMuon_phi,evt_wt_Up);
                  h_mu1mu2dR_Up_Z->Fill(dR,evt_wt_Up);
                  h_mu1mu2dPhi_Up_Z->Fill(dPhi,evt_wt_Up);
                  h_mu1mu2dEta_Up_Z->Fill(dEta,evt_wt_Up);
                  h_softJet5_Up_Z->Fill(t_SoftActivityJetNjets5,evt_wt_Up);
                  h_dijet_pt_Up_Z->Fill(t_diJet_pt,evt_wt_Up);
                  h_dijet_eta_Up_Z->Fill(t_diJet_eta,evt_wt_Up);
                  h_dijet_phi_Up_Z->Fill(t_diJet_phi,evt_wt_Up);
                  h_Mjj_Up_Z->Fill(t_diJet_mass,evt_wt_Up);
                  h_dijet_dEta_Up_Z->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt_Up);
                  
                  
                  
                  h_diMuon_mass_Down_Z->Fill(diMuon_mass,evt_wt_Down);
                  h_diMuon_pt_Down_Z->Fill(diMuon_pt,evt_wt_Down);
                  h_diMuon_eta_Down_Z->Fill(diMuon_eta,evt_wt_Down);
                  h_diMuon_phi_Down_Z->Fill(diMuon_phi,evt_wt_Down);
                  h_mu1mu2dR_Down_Z->Fill(dR,evt_wt_Down);
                  h_mu1mu2dPhi_Down_Z->Fill(dPhi,evt_wt_Down);
                  h_mu1mu2dEta_Down_Z->Fill(dEta,evt_wt_Down);
                  h_softJet5_Down_Z->Fill(t_SoftActivityJetNjets5,evt_wt_Down);
                  h_dijet_pt_Down_Z->Fill(t_diJet_pt,evt_wt_Down);
                  h_dijet_eta_Down_Z->Fill(t_diJet_eta,evt_wt_Down);
                  h_dijet_phi_Down_Z->Fill(t_diJet_phi,evt_wt_Down);
                  h_Mjj_Down_Z->Fill(t_diJet_mass,evt_wt_Down);
                  h_dijet_dEta_Down_Z->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt_Down);
              }
              softJet5=t_SoftActivityJetNjets5;
              dRmm = dR;
              dEtamm = dEta;
              dPhimm = dPhi;
              M_jj=t_diJet_mass;
              pt_jj=t_diJet_pt;
              eta_jj=t_diJet_eta;
              phi_jj=t_diJet_phi;
              TLorentzVector dimu,dijet,mmjj;
              dimu.SetPtEtaPhiM(diMuon_pt,diMuon_eta,diMuon_phi,diMuon_mass);
              dijet.SetPtEtaPhiM(t_diJet_pt,t_diJet_eta,t_diJet_phi,t_diJet_mass);
              mmjj=dimu+dijet;
              M_mmjj=mmjj.M();
              pt_mmjj=mmjj.Pt();
              eta_mmjj=mmjj.Eta();
              phi_mmjj=mmjj.Phi();
              dEta_jj=(*t_Jet_eta)[0]-(*t_Jet_eta)[1];
              Zep=(diMuon_eta-0.5*((*t_Jet_eta)[0]+(*t_Jet_eta)[1]));///fabs((*t_Jet_eta)[0]-(*t_Jet_eta)[1]));
              //====================BDT vars===============================
              mumujj_pt=log(pt_mmjj);
              DeltaEtaQQ=fabs((*t_Jet_eta)[0]-(*t_Jet_eta)[1]);
              softActivityEWK_njets5=t_SoftActivityJetNjets5;
              ll_zstar=Zep/((*t_Jet_eta)[0]-(*t_Jet_eta)[1]);
              ll_pt=dimu.Pt();
              
              TVector3 dimuon, j2;
              dimuon.SetPtEtaPhi(dimu.Pt(),dimu.Eta(),dimu.Phi());
              //cout<<dimuon.Pt()<<", "<<dimuon.Eta()<<", "<<dimuon.Phi()<<endl;
              j2.SetPtEtaPhi((*t_Jet_pt)[1],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
              //cout<<j2.Pt()<<", "<<j2.Eta()<<", "<<j2.Phi()<<endl;
              //cout<<j2.Dot(dimuon)<<","<<j2.Mag()<<", "<<dimuon.Mag()<<","<<theta2<<endl;
              theta2=j2.Dot(dimuon)/(j2.Mag()*dimuon.Mag());
              
              impulsoZ=log(fabs(mmjj.Pz()));
              if(fabs((*t_Jet_eta)[0])>fabs((*t_Jet_eta)[1]))
                  maxAbsEta=fabs((*t_Jet_eta)[0]);
              else maxAbsEta=fabs((*t_Jet_eta)[1]);
              //cout<<1.999995*((*t_Jet_qgl)[1]-0.5)<<", "<<(*t_Jet_qgl)[1]<<endl;
              qgl_2qAtanh=atanh(1.999995*((*t_Jet_qgl)[1]-0.5));
              
              BDT_Pisa = reader->EvaluateMVA( methodName );
              //==========================================================
              
              double dr[4];
              dr[0]=DeltaR((*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
              dr[1] = DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
              dr[2]=DeltaR((*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
              dr[3]=DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
              
              for (int c = 0 ; c < 3; c++)
              {
                  for (int d = 0 ; d < 3 - c; d++)
                  {
                      if (dr[d] > dr[d+1]) /* For decreasing order use < */
                      {
                          double swap = dr[d];
                          dr[d]   = dr[d+1];
                          dr[d+1] = swap;
                      }
                  }
              }
              dRmin_mj=dr[0];
              dRmax_mj=dr[3];
              
              dr[0]=DeltaR(diMuon_eta,diMuon_phi,(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
              dr[1]=DeltaR(diMuon_eta,diMuon_phi,(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
              
              if(dr[0]>dr[1]){
                  dRmin_mmj=dr[1];
                  dRmax_mmj=dr[0];
              }
              else{
                  dRmin_mmj=dr[0];
                  dRmax_mmj=dr[1];
              }
              
              dPhijj=DeltaPhi((*t_Jet_phi)[1],(*t_Jet_phi)[0]);
              leadingJet_pt=(*t_Jet_pt)[0];
              subleadingJet_pt=(*t_Jet_pt)[1];
              leadingJet_eta = (*t_Jet_eta)[0];
              subleadingJet_eta = (*t_Jet_eta)[1];
              leadingJet_qgl=(*t_Jet_qgl)[0];
              subleadingJet_qgl=(*t_Jet_qgl)[1];
              cthetaCS=2*(mu2.E()*mu1.Pz()-mu1.E()*mu2.Pz())/(diMuon_mass*sqrt(pow(diMuon_mass,2)+pow(diMuon_pt,2)));
              //cout<<event<<" "<<cthetaCS<<" "<<Zep<<endl;
              cattree->Fill();
              if((diMuon_mass<120. || diMuon_mass>130.) && diMuon_mass >108.){
                  h_dijet_dPhijj_VBF->Fill(dPhijj,evt_wt);
                  h_jet1_pt_VBF->Fill(leadingJet_pt,evt_wt);
                  h_jet2_pt_VBF->Fill(subleadingJet_pt,evt_wt);
                  h_jet1_eta_VBF->Fill(leadingJet_eta,evt_wt);
                  h_jet2_eta_VBF->Fill(subleadingJet_eta,evt_wt);
                  h_cthetaCS_VBF->Fill(cthetaCS,evt_wt);
                  h_dRmin_mj_VBF->Fill(dRmin_mj,evt_wt);
                  h_dRmax_mj_VBF->Fill(dRmax_mj,evt_wt);
                  h_dRmin_mmj_VBF->Fill(dRmin_mmj,evt_wt);
                  h_dRmax_mmj_VBF->Fill(dRmax_mmj,evt_wt);
                  h_Zep_VBF->Fill(Zep,evt_wt);
                  h_M_mmjj_VBF->Fill(M_mmjj,evt_wt);
                  h_pt_mmjj_VBF->Fill(pt_mmjj,evt_wt);
                  h_eta_mmjj_VBF->Fill(eta_mmjj,evt_wt);
                  h_phi_mmjj_VBF->Fill(phi_mmjj,evt_wt);
                  h_leadingJet_qgl_VBF->Fill(leadingJet_qgl,evt_wt);
                  h_subleadingJet_qgl_VBF->Fill(subleadingJet_qgl,evt_wt);
                  h_2D_ptjj_dRmm->Fill(pt_jj,dRmm,evt_wt);
                  h_2D_CS_dEtamm->Fill(M_mmjj,M_jj,evt_wt);
                  h_2D_Hpt_ptjj->Fill(diMuon_pt,pt_jj,evt_wt);
                  h_2D_Hpt_dRmm->Fill(diMuon_pt,dRmm,evt_wt);
                  h_2D_j1eta_dEtajj->Fill(leadingJet_eta,dEta_jj,evt_wt);
                  h_2D_j2eta_dEtajj->Fill(subleadingJet_eta,dEta_jj,evt_wt);
                  h_2D_etammjj_etajj->Fill(eta_mmjj,eta_jj,evt_wt);
                  h_2D_j1pt_ptjj->Fill(leadingJet_pt,pt_jj,evt_wt);
                  
                  h_dijet_dPhijj_Up_VBF->Fill(dPhijj,evt_wt_Up);
                  h_jet1_pt_Up_VBF->Fill(leadingJet_pt,evt_wt_Up);
                  h_jet2_pt_Up_VBF->Fill(subleadingJet_pt,evt_wt_Up);
                  h_jet1_eta_Up_VBF->Fill(leadingJet_eta,evt_wt_Up);
                  h_jet2_eta_Up_VBF->Fill(subleadingJet_eta,evt_wt_Up);
                  h_cthetaCS_Up_VBF->Fill(cthetaCS,evt_wt_Up);
                  h_dRmin_mj_Up_VBF->Fill(dRmin_mj,evt_wt_Up);
                  h_dRmax_mj_Up_VBF->Fill(dRmax_mj,evt_wt_Up);
                  h_dRmin_mmj_Up_VBF->Fill(dRmin_mmj,evt_wt_Up);
                  h_dRmax_mmj_Up_VBF->Fill(dRmax_mmj,evt_wt_Up);
                  h_Zep_Up_VBF->Fill(Zep,evt_wt_Up);
                  h_M_mmjj_Up_VBF->Fill(M_mmjj,evt_wt_Up);
                  h_pt_mmjj_Up_VBF->Fill(pt_mmjj,evt_wt_Up);
                  h_eta_mmjj_Up_VBF->Fill(eta_mmjj,evt_wt_Up);
                  h_phi_mmjj_Up_VBF->Fill(phi_mmjj,evt_wt_Up);
                  h_leadingJet_qgl_Up_VBF->Fill(leadingJet_qgl,evt_wt_Up);
                  h_subleadingJet_qgl_Up_VBF->Fill(subleadingJet_qgl,evt_wt_Up);
                  
                  
                  h_dijet_dPhijj_Down_VBF->Fill(dPhijj,evt_wt_Down);
                  h_jet1_pt_Down_VBF->Fill(leadingJet_pt,evt_wt_Down);
                  h_jet2_pt_Down_VBF->Fill(subleadingJet_pt,evt_wt_Down);
                  h_jet1_eta_Down_VBF->Fill(leadingJet_eta,evt_wt_Down);
                  h_jet2_eta_Down_VBF->Fill(subleadingJet_eta,evt_wt_Down);
                  h_cthetaCS_Down_VBF->Fill(cthetaCS,evt_wt_Down);
                  h_dRmin_mj_Down_VBF->Fill(dRmin_mj,evt_wt_Down);
                  h_dRmax_mj_Down_VBF->Fill(dRmax_mj,evt_wt_Down);
                  h_dRmin_mmj_Down_VBF->Fill(dRmin_mmj,evt_wt_Down);
                  h_dRmax_mmj_Down_VBF->Fill(dRmax_mmj,evt_wt_Down);
                  h_Zep_Down_VBF->Fill(Zep,evt_wt_Down);
                  h_M_mmjj_Down_VBF->Fill(M_mmjj,evt_wt_Down);
                  h_pt_mmjj_Down_VBF->Fill(pt_mmjj,evt_wt_Down);
                  h_eta_mmjj_Down_VBF->Fill(eta_mmjj,evt_wt_Down);
                  h_phi_mmjj_Down_VBF->Fill(phi_mmjj,evt_wt_Down);
                  h_leadingJet_qgl_Down_VBF->Fill(leadingJet_qgl,evt_wt_Down);
                  h_subleadingJet_qgl_Down_VBF->Fill(subleadingJet_qgl,evt_wt_Down);
              }
              else{
                  h_dijet_dPhijj_Z->Fill(dPhijj,evt_wt);
                  h_jet1_pt_Z->Fill(leadingJet_pt,evt_wt);
                  h_jet2_pt_Z->Fill(subleadingJet_pt,evt_wt);
                  h_jet1_eta_Z->Fill(leadingJet_eta,evt_wt);
                  h_jet2_eta_Z->Fill(subleadingJet_eta,evt_wt);
                  h_cthetaCS_Z->Fill(cthetaCS,evt_wt);
                  h_dRmin_mj_Z->Fill(dRmin_mj,evt_wt);
                  h_dRmax_mj_Z->Fill(dRmax_mj,evt_wt);
                  h_dRmin_mmj_Z->Fill(dRmin_mmj,evt_wt);
                  h_dRmax_mmj_Z->Fill(dRmax_mmj,evt_wt);
                  h_Zep_Z->Fill(Zep,evt_wt);
                  h_M_mmjj_Z->Fill(M_mmjj,evt_wt);
                  h_pt_mmjj_Z->Fill(pt_mmjj,evt_wt);
                  h_eta_mmjj_Z->Fill(eta_mmjj,evt_wt);
                  h_phi_mmjj_Z->Fill(phi_mmjj,evt_wt);
                  h_leadingJet_qgl_Z->Fill(leadingJet_qgl,evt_wt);
                  h_subleadingJet_qgl_Z->Fill(subleadingJet_qgl,evt_wt);
                  
                  
                  h_dijet_dPhijj_Up_Z->Fill(dPhijj,evt_wt_Up);
                  h_jet1_pt_Up_Z->Fill(leadingJet_pt,evt_wt_Up);
                  h_jet2_pt_Up_Z->Fill(subleadingJet_pt,evt_wt_Up);
                  h_jet1_eta_Up_Z->Fill(leadingJet_eta,evt_wt_Up);
                  h_jet2_eta_Up_Z->Fill(subleadingJet_eta,evt_wt_Up);
                  h_cthetaCS_Up_Z->Fill(cthetaCS,evt_wt_Up);
                  h_dRmin_mj_Up_Z->Fill(dRmin_mj,evt_wt_Up);
                  h_dRmax_mj_Up_Z->Fill(dRmax_mj,evt_wt_Up);
                  h_dRmin_mmj_Up_Z->Fill(dRmin_mmj,evt_wt_Up);
                  h_dRmax_mmj_Up_Z->Fill(dRmax_mmj,evt_wt_Up);
                  h_Zep_Up_Z->Fill(Zep,evt_wt_Up);
                  h_M_mmjj_Up_Z->Fill(M_mmjj,evt_wt_Up);
                  h_pt_mmjj_Up_Z->Fill(pt_mmjj,evt_wt_Up);
                  h_eta_mmjj_Up_Z->Fill(eta_mmjj,evt_wt_Up);
                  h_phi_mmjj_Up_Z->Fill(phi_mmjj,evt_wt_Up);
                  h_leadingJet_qgl_Up_Z->Fill(leadingJet_qgl,evt_wt_Up);
                  h_subleadingJet_qgl_Up_Z->Fill(subleadingJet_qgl,evt_wt_Up);
                  
                  h_dijet_dPhijj_Down_Z->Fill(dPhijj,evt_wt_Down);
                  h_jet1_pt_Down_Z->Fill(leadingJet_pt,evt_wt_Down);
                  h_jet2_pt_Down_Z->Fill(subleadingJet_pt,evt_wt_Down);
                  h_jet1_eta_Down_Z->Fill(leadingJet_eta,evt_wt_Down);
                  h_jet2_eta_Down_Z->Fill(subleadingJet_eta,evt_wt_Down);
                  h_cthetaCS_Down_Z->Fill(cthetaCS,evt_wt_Down);
                  h_dRmin_mj_Down_Z->Fill(dRmin_mj,evt_wt_Down);
                  h_dRmax_mj_Down_Z->Fill(dRmax_mj,evt_wt_Down);
                  h_dRmin_mmj_Down_Z->Fill(dRmin_mmj,evt_wt_Down);
                  h_dRmax_mmj_Down_Z->Fill(dRmax_mmj,evt_wt_Down);
                  h_Zep_Down_Z->Fill(Zep,evt_wt_Down);
                  h_M_mmjj_Down_Z->Fill(M_mmjj,evt_wt_Down);
                  h_pt_mmjj_Down_Z->Fill(pt_mmjj,evt_wt_Down);
                  h_eta_mmjj_Down_Z->Fill(eta_mmjj,evt_wt_Down);
                  h_phi_mmjj_Down_Z->Fill(phi_mmjj,evt_wt_Down);
                  h_leadingJet_qgl_Down_Z->Fill(leadingJet_qgl,evt_wt_Down);
                  h_subleadingJet_qgl_Down_Z->Fill(subleadingJet_qgl,evt_wt_Down);
              }
          }
          
          //VH, had
          else if(t_diJet_mass>70. && t_diJet_mass<100){
              cat_index = 6;
              double binv = catyield->GetBinContent(6);
              binv = binv +evt_wt;
              catyield->SetBinContent(6,binv);
              h_diMuon_mass_VHHad->Fill(diMuon_mass,evt_wt);
          }
          else if(t_nJet>1){
              cat_index = 7;
              double binv = catyield->GetBinContent(7);
              binv = binv + evt_wt;
              catyield->SetBinContent(7,binv);
          }
          else{
              cat_index = 8;
              double binv = catyield->GetBinContent(8);
              binv = binv + evt_wt;
              catyield->SetBinContent(8,binv);
              h_diMuon_mass_ggH->Fill(diMuon_mass,evt_wt);
         }
         cattree->Fill();
      }//mass cut
   }//event loop
}

void HiggsMuMu::EventLoop(const char *data,const char *isData, double scale)
{  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%5000 ==0) cout<<jentry<<endl;
      TLorentzVector m1,m2,m12;
      bool olp=false;
      int mu1idx=t_mu1;
      int mu2idx=t_mu2; 
      if((t_mu1<100) && (t_mu2<100))  olp = true;

      if(olp){
         TLorentzVector m1,m2,m12;
         bool olp=false;
         int mu1idx=t_mu1;
         int mu2idx=t_mu2;
         double evt_wt=1.;
         evt_wt=t_genWeight*scale;
         double lepSF = 1.0;
         if(*isData=='F'){
         if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_TRIG)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
         else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_TRIG)[t_mu2]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
          
         evt_wt*=lepSF;
         evt_wt*=t_puWeight;
         evt_wt*=t_PrefireWeight;
             
         float w_zpt = 1.0;
             //CP5
             if(m12.Pt()<10) w_zpt = 0.910385;
             else if(m12.Pt()>10. && m12.Pt()<20.) w_zpt = 1.13543;
             else if(m12.Pt()>20. && m12.Pt()<30.) w_zpt = 1.10441;
             else if(m12.Pt()>30. && m12.Pt()<40.) w_zpt = 1.01315;
             else if(m12.Pt()>40. && m12.Pt()<50.) w_zpt = 0.982598;
             else if(m12.Pt()>50. && m12.Pt()<60.) w_zpt = 0.980697;
             else if(m12.Pt()>60. && m12.Pt()<70.) w_zpt = 0.972673;
             else if(m12.Pt()>70. && m12.Pt()<80.) w_zpt = 0.972325;
             else if(m12.Pt()>80. && m12.Pt()<100.) w_zpt = 0.966127;
             else if(m12.Pt()>100. && m12.Pt()<150.) w_zpt = 0.953262;
             else if(m12.Pt()>150. && m12.Pt()<200.) w_zpt = 0.933403;
             else w_zpt = 0.904518;
             
             evt_wt*=w_zpt;
             
        }
	//if(t_diMuon_mass<120. || t_diMuon_mass>130.)h_diMuon_mass->Fill(t_diMuon_mass,evt_wt);
        m1.SetPtEtaPhiM((*t_Mu_pt)[t_mu1],(*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Mu_mass)[t_mu1]);
        m2.SetPtEtaPhiM((*t_Mu_pt)[t_mu2],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],(*t_Mu_mass)[t_mu2]);
        m12=m1+m2;
        
       
        /*
        if(m12.Pt()<10) w_zpt = 1.05817;
        else if(m12.Pt()>10. && m12.Pt()<20.) w_zpt = 0.994488;
        else if(m12.Pt()>20. && m12.Pt()<30.) w_zpt = 0.930056;
        else if(m12.Pt()>30. && m12.Pt()<40.) w_zpt = 0.925206;
        else if(m12.Pt()>40. && m12.Pt()<50.) w_zpt = 0.946403;
        else if(m12.Pt()>50. && m12.Pt()<60.) w_zpt = 0.962136;
        else if(m12.Pt()>60. && m12.Pt()<70.) w_zpt = 0.965316;
        else if(m12.Pt()>70. && m12.Pt()<80.) w_zpt = 0.978209;
        else if(m12.Pt()>80. && m12.Pt()<100.) w_zpt = 0.988761;
        else if(m12.Pt()>100. && m12.Pt()<150.) w_zpt = 0.982497;
        else if(m12.Pt()>150. && m12.Pt()<200.) w_zpt = 0.971749;
        else w_zpt = 0.914429; 
        */
         
	if(*isData=='T'){
            if((m12.M()<120.) || (m12.M()>130.)) h_diMuon_mass->Fill(m12.M(),evt_wt);
        }
        else h_diMuon_mass->Fill(m12.M(),evt_wt);


        if(m12.M()<106. && m12.M()>76.){
            h_Zpt->Fill(m12.Pt(),evt_wt);
            h_Zm->Fill(m12.M(),evt_wt);
        }

	if(*isData=='F') if(m12.M()>120. && m12.M()<130.)h_diMuon_mass_SR->Fill(m12.M(),evt_wt);
	if((m12.M()>110. && m12.M()<120.) || (m12.M()>130. && m12.M()<150.)){ //make plot in mass sidebind 110-120, 130-150
          h_sidept->Fill(m12.Pt(),evt_wt);
	  h_mu1pt->Fill((*t_Mu_pt)[mu1idx],evt_wt);
	  h_mu2pt->Fill((*t_Mu_pt)[mu2idx],evt_wt);
	  h_mu1eta->Fill((*t_Mu_eta)[mu1idx],evt_wt);
	  h_mu2eta->Fill((*t_Mu_eta)[mu2idx],evt_wt);
	  h_mu1phi->Fill((*t_Mu_phi)[mu1idx],evt_wt);
	  h_mu2phi->Fill((*t_Mu_phi)[mu2idx],evt_wt);
	  double dR= DeltaR((*t_Mu_eta)[mu1idx],(*t_Mu_phi)[mu1idx],(*t_Mu_eta)[mu2idx],(*t_Mu_phi)[mu2idx]);
	  h_mu1mu2dR->Fill(dR,evt_wt);
	  double dPhi= DeltaPhi((*t_Mu_phi)[mu2idx],(*t_Mu_phi)[mu1idx]);
	  h_mu1mu2dPhi->Fill(dPhi,evt_wt);
	  // if (Cut(ientry) < 0) continue;
	
	  //cout<<"done till dPhi\n";
		
	  h_diMuon_pt->Fill(m12.Pt(),evt_wt);
	  h_diMuon_eta->Fill(m12.Eta(),evt_wt);
	  h_diMuon_phi->Fill(m12.Phi(),evt_wt);
	  h_diMuon_mass_110To150->Fill(m12.M(),evt_wt);
	  if(m12.M()>110.&& m12.M()<120.)h_diMuon_mass_110To120->Fill(m12.M(),evt_wt);
	  if(m12.M()>130. && m12.M()<150.)h_diMuon_mass_130To150->Fill(m12.M(),evt_wt);
	  //cout<<"done all dimu\n";	
	  h_MET_pt->Fill(t_MET_pt,evt_wt);
	  h_METphi->Fill(t_MET_phi,evt_wt);
	  h_MET_sumEt->Fill(t_MET_sumEt,evt_wt);

	  h_Njet->Fill(t_nJet,evt_wt);
	  h_Nbjet->Fill(t_nbJet,evt_wt);
	  //cout<<t_nJet<<endl;
	  if(t_nJet>0){
	    //cout<<(*t_Jet_pt)[0]<<endl;
	    h_j1pt->Fill((*t_Jet_pt)[0],evt_wt);
	    //cout<<"pt\n";
	    h_j1eta->Fill((*t_Jet_eta)[0],evt_wt);
	    h_j1phi->Fill((*t_Jet_phi)[0],evt_wt);
	  }
	  //cout<<"done with jet 1\n";
	  if(t_nJet>1){
	    h_j2pt->Fill((*t_Jet_pt)[1],evt_wt);
	    h_j2eta->Fill((*t_Jet_eta)[1],evt_wt);
	    h_j2phi->Fill((*t_Jet_phi)[1],evt_wt);

	    dR= DeltaR((*t_Jet_eta)[0],(*t_Jet_phi)[0],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
	    h_j1j2dR->Fill(dR,evt_wt);
	    dPhi= DeltaPhi((*t_Jet_phi)[0],(*t_Jet_phi)[1]);
	    h_j1j2dPhi->Fill(dPhi,evt_wt);

	    h_dijet_pt->Fill(t_diJet_pt,evt_wt);
	    h_dijet_eta->Fill(t_diJet_eta,evt_wt);
	    h_dijet_phi->Fill(t_diJet_phi,evt_wt);
	    h_Mjj->Fill(t_diJet_mass,evt_wt);
	  }
	}
      }
      //cout<<"ok done\n";
   }
}
