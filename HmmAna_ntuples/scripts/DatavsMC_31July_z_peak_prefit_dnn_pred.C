#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"
using std::cout;
using std::endl;
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) ;
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame );
void histDraw(TPad *p,TH1F *h1a,TGraphAsymmErrors *h1b, THStack *hs,TH1F *vFrame );
void DatavsMC_31July_z_peak_prefit_dnn_pred()
{

  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_13TeV  = "35.9  fb^{-1}";
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=11; //11 = default: left-aligned
  //int iPos=33; // right-aligned
   //int iPos=22; // center-aligned
  //int iPos=0;//out of frame (in exceptional cases)
  //First declare the file names

    // references for T, B, L, R
  /* float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
  */
  TString f[77];
  f[0]="fitDiagnosticstest_dijet_inv_mass.root";
  f[1]="fitDiagnosticstest_h_peak_sb.root";
  //Declare other constants, strings that you might need here.

  //Now let us define the plot name to overlay
  TString plotname[33];
  plotname[0]= "data";
  plotname[1] = "dy_0j";
  plotname[2] = "dy_1j";
  plotname[3] = "dy_2j";
  plotname[4] ="ewk_lljj_mll50_mjj120";
  plotname[5] = "ggh";
  plotname[6]= "st_tw_antitop";
  plotname[7] = "st_tw_top";
  plotname[8] = "ttjets_dl";
  plotname[9] = "ttjets_sl";
  plotname[10] = "vbf";
  plotname[11] = "ww_2l2nu";
  plotname[12] = "wz_2l2q";
  plotname[13] = "wz_3lnu";
  plotname[14] = "zz";
  plotname[15] = "zzz";
  plotname[16] = "wzz";
  plotname[17] = "www";
  plotname[18] = "tth";
  plotname[19] = "wmh";
  plotname[20] = "wph";
  plotname[21] = "zh";
  gStyle->SetOptStat(0);
  //Also give fancy name for the axis titles
  TString xtitle= "P_{T} (GeV)";
  TString xtitle1= "#eta";
  TString xtitle2= "#phi";
  TString xtitle3= "H_{T} (GeV)";
  TString ytitle = "Number of Events"; // Or "Events"
  
  //Now let us open the files
  TFile *file[77];
  for(int i=0;i<2;i++) file[i] = new TFile(f[i]);
  cout<<"2 files\n";
  //change directory if required.
  TDirectory *directory_md0[77];
  TDirectory *Main_Dir[77];
  char dir_name[100];
  for(int i=0;i<2;i++) Main_Dir[i] = (TDirectory*) file[i]->Get("shapes_prefit");
 
  for(int j=0;j<2;j++){
    directory_md0[j] = (TDirectory*) Main_Dir[j]->Get("ch1");
  
  
  }
  TString chname[33];
  
 
  
  chname[0]= "z_peak_dnn_pred_prefit";
  chname[1]= "z_peak_pnn_pred_prefit";
  // cout<<directory[0]<<endl;

  Double_t xbins[35]={0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300,400,500,700,1000.};
  
  TString name_unc;
  double sum=0.0;
  char hist_sum[20];
  TLegend *lg[32];
  char name[100];
  THStack *hs[32];
  TCanvas *c[32];
  TH1F *vFrame[32];
  TPad *pad1[32];
  TPad *pad2[32];
  TPad *pad3[32];
  TH1F *hd[32];
  TLine *l[32];
  //Now open the respective histograms from the files
  TH1F *h1[32],*h2[32],*h4[32];
  TH1F *h5[32],*h6[32],*h7[32],*h8[32];
  TH1F *h9[32],*h10[32],*h11[32],*h12[32];
  TH1F *h13[32],*h14[32],*h15[32],*h16[32];
  TH1F *h17[32],*h18[32],*h19[32],*h20[32];
  TH1F *h21[32],*h22[32],*h23[32],*h24[32];
  TH1F *h25[32],*h26[32];
  TH1F *h[32], *h_err[32], *h_temp,*h_total[32];
  TH2F *h_covar[32];
  TGraphAsymmErrors *h3[32];
   TAxis *xaxis_temp;
   int nbins_temp;
   float xlow;
   float xup;

  double bin_err, stat_err,jec_err,jer_err;
  int nbins;
  TString bins_labels[]={"0-0.001","0.001-0.009","0.009-0.01","0.01-0.011","0.011-0.012","0.012-0.013","0.013-0.014","0.014-0.015","0.015-0.016","0.016-0.017","0.017-0.018","0.018-0.019","0.019-0.02","0.02-0.021","0.021-0.022","0.02-0.023","0.023-0.024","0.024-0.025","0.025-0.026","0.026-0.027","0.027-0.028","0.028-0.029","0.029-0.03","0.03-0.031","0.031-0.032","0.032-0.033","0.033-0.034","0.034-0.035","0.035-0.036","0.036-0.037","0.037-0.038","0.038-0.039","0.039-1."};
  double dnn_bins[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
  for(int j=0;j<2;j++){
    //dy 0j
    h1[j] = (TH1F*)directory_md0[j]->Get(plotname[1]);
    nbins = h1[j]->GetXaxis()->GetNbins();
    
    h1[j]->SetBins(nbins, dnn_bins);
    //h1[j]->GetXaxis()->SetRangeUser(0.001,1.);
    sprintf(name,"h%i",j);
    h[j] = (TH1F*)h1[j]->Clone(name); //For MC // doing it here so that I can adjust the color of h[j]
    decorate(h1[j],"",ytitle,"",kYellow,2,kYellow,20,1);

    //dy 1j
    h25[j] = (TH1F*)directory_md0[j]->Get(plotname[2]);
    xaxis_temp=h25[j]->GetXaxis();
    nbins_temp=h25[j]->GetNbinsX();
    xlow=xaxis_temp->GetBinLowEdge(1);
    xup=xaxis_temp->GetBinUpEdge(nbins_temp);
    decorate(h25[j],"",ytitle,"",kYellow,2,kYellow,21,1);
    h25[j]->SetBins(nbins, dnn_bins);
    // h25[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //dy 2j
    h23[j] = (TH1F*)directory_md0[j]->Get(plotname[3]);
    decorate(h23[j],"",ytitle,"",kYellow,2,kYellow,21,1);
    h23[j]->SetBins(nbins, dnn_bins);
    //h23[j]->GetXaxis()->SetRangeUser(0.001,1.);
    cout<<"dy\n";
    
    //Data
    //h3[j] = (TH1F*)directory_md0[j]->Get(plotname[0]);
    //decorate(h3[j],"",ytitle,"",kBlack,2,kBlack,8,0);
    h3[j] = (TGraphAsymmErrors*)directory_md0[j]->Get(plotname[0]);
    h3[j]->SetMarkerColor(kBlack);
    h3[j]->SetLineColor(kBlack); 
    h3[j]->SetLineWidth(2);
    h3[j]->SetMarkerStyle(8);
    h3[j]->SetMarkerSize(1.3);
    
      //decorate(h3[j],"",ytitle,"",kBlack,2,kBlack,8,0);
    //VBF
    h2[j] = (TH1F*)directory_md0[j]->Get(plotname[10]);
    decorate(h2[j],"",ytitle,"",kViolet+1,2,kViolet+1,21,0);
    cout<<"vbf\n";
    h2[j]->SetBins(nbins, dnn_bins);
    //h2[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //ggh
    h6[j] = (TH1F*)directory_md0[j]->Get(plotname[5]);
    decorate(h6[j],"",ytitle,"",kCyan+4,2,kCyan+4,21,0);
    h6[j]->SetBins(nbins, dnn_bins);
    // h6[j]->GetXaxis()->SetRangeUser(0.001,1.);
    cout<<"ggh\n";
    //ewk
    h24[j] = (TH1F*)directory_md0[j]->Get(plotname[4]);
    decorate(h24[j],"",ytitle,"",kCyan,2,kCyan,21,1);
    h24[j]->SetBins(nbins, dnn_bins);
    //h24[j]->GetXaxis()->SetRangeUser(0.001,1.);
    cout<<"ewk\n";
    //ttdl
    h5[j] = (TH1F*)directory_md0[j]->Get(plotname[8]);
    decorate(h5[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);
    h5[j]->SetBins(nbins, dnn_bins);
    //h5[j]->GetXaxis()->SetRangeUser(0.001,1.);
     
    //ttsl
    h4[j] = (TH1F*)directory_md0[j]->Get(plotname[9]);
    decorate(h4[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);
    h4[j]->SetBins(nbins, dnn_bins);
    // h4[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //ww to 2l2v
    h15[j] = (TH1F*)directory_md0[j]->Get(plotname[11]);
    decorate(h15[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    h15[j]->SetBins(nbins, dnn_bins);
    // h15[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //wz to 2l2q
    h14[j] = (TH1F*)directory_md0[j]->Get(plotname[12]);
    decorate(h14[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    h14[j]->SetBins(nbins, dnn_bins);
    //h14[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //wz to 3lnu
    h13[j] = (TH1F*)directory_md0[j]->Get(plotname[13]);
    decorate(h13[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    h13[j]->SetBins(nbins, dnn_bins);
    //h13[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //zz
    h11[j] = (TH1F*)directory_md0[j]->Get(plotname[14]);
    decorate(h11[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    h11[j]->SetBins(nbins, dnn_bins);
    //h11[j]->GetXaxis()->SetRangeUser(0.001,1.);
   
    //st_Tw_antittop
    h18[j] = (TH1F*)directory_md0[j]->Get(plotname[6]);
    decorate(h18[j],"",ytitle,"",kPink-2,2,kPink-2,21,1);
    h18[j]->SetBins(nbins, dnn_bins);
    //h18[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //st_Tw_top
    h12[j] = (TH1F*)directory_md0[j]->Get(plotname[7]);
    decorate(h12[j],"",ytitle,"",kPink-2,2,kPink-2,21,1);
    h12[j]->SetBins(nbins, dnn_bins);
    //h12[j]->GetXaxis()->SetRangeUser(0.001,1.);
    cout<<"top\n";
    //www
    h17[j] = (TH1F*)directory_md0[j]->Get(plotname[15]);
    decorate(h17[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
    h17[j]->SetBins(nbins, dnn_bins);
    //h17[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //wzz
    h21[j] = (TH1F*)directory_md0[j]->Get(plotname[16]);
    decorate(h21[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
    h21[j]->SetBins(nbins, dnn_bins);
    //h21[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //zzz
    h22[j] = (TH1F*)directory_md0[j]->Get(plotname[17]);
    decorate(h22[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
    h22[j]->SetBins(nbins, dnn_bins);
    // h22[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //wph
    h7[j] = (TH1F*)directory_md0[j]->Get(plotname[20]);
    decorate(h7[j],"",ytitle,"",kBlue,2,kBlue,21,0);
    h7[j]->SetBins(nbins, dnn_bins);
    h7[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //wmh
    h8[j] = (TH1F*)directory_md0[j]->Get(plotname[19]);
    decorate(h8[j],"",ytitle,"",kPink-7,2,kPink-7,21,0);
    h8[j]->SetBins(nbins, dnn_bins);
    //h8[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //zh
    h9[j] = (TH1F*)directory_md0[j]->Get(plotname[21]);
    decorate(h9[j],"",ytitle,"",kRed-2,2,kRed-2,21,0);
    h9[j]->SetBins(nbins, dnn_bins);
    //h9[j]->GetXaxis()->SetRangeUser(0.001,1.);
    //tth
    h10[j] = (TH1F*)directory_md0[j]->Get(plotname[18]);
    decorate(h10[j],"",ytitle,"",kRed,2,kRed,21,0);
    h10[j]->SetBins(nbins, dnn_bins);
    //h10[j]->GetXaxis()->SetRangeUser(0.001,1.);

    h_covar[j]=(TH2F*)directory_md0[j]->Get("total_covar");
 
    h_total[j]=(TH1F*)directory_md0[j]->Get("total");
  }

  cout<<"got all hists\n";

  TLegend *lg1 = new TLegend(0.35,0.75,0.85,0.9);
  lg1->SetTextFont(132);
  lg1->SetBorderSize(0);
  lg1->SetNColumns(3);
  gStyle->SetLegendTextSize(0.04);
  lg1->AddEntry(h3[0],"Data","p");
  
  lg1->AddEntry(h1[0],"DY+jets","f");
  lg1->AddEntry(h24[0],"EWK LL+jets","f");
  lg1->AddEntry(h5[0],"t#bar{t}","f");
  lg1->AddEntry(h11[0],"VV ","f");
   lg1->AddEntry(h12[0],"Single top","f");
  lg1->AddEntry(h17[0],"VVV","f");
  
  //lg1->AddEntry(h19[0],"ttX+Jets","f");
  

  lg1->AddEntry(h2[0],"VBF","l");
  lg1->AddEntry(h6[0],"ggH","l");
  lg1->AddEntry(h7[0],"W^{+}H","l");
  lg1->AddEntry(h8[0],"W^{-}H","l");
  lg1->AddEntry(h10[0],"ttH","l");
  lg1->AddEntry(h9[0],"ZH","l");
  for(int j=0;j<2;j++){
    sprintf(dir_name,"ch%i",j+1);
    sprintf(name,"hs%i",j);
    hs[j] = new THStack(name,dir_name); // Stack of all backgrounds
   
    
   
   
    //hs[j]->Add(h19[j]);//ttw
    //hs[j]->Add(h20[j]); //ttz
    hs[j]->Add(h17[j]); //www
    hs[j]->Add(h21[j]); //wzz
    hs[j]->Add(h22[j]); //zzz
    //hs[j]->Add(h26[j]);//wwz
    hs[j]->Add(h24[j]);//EWK
    //hs[j]->Add(h16[j]); //wwtolnuqq
    hs[j]->Add(h12[j]); //st_tw_top
    hs[j]->Add(h18[j]);//st_tw_antitop
    hs[j]->Add(h11[j]); // zz ##### zz to 4l
   
    hs[j]->Add(h13[j]);//wz to 3l nu
    
    hs[j]->Add(h14[j]); //wz to 2l2q
   
    hs[j]->Add(h15[j]); //ww to 2l2v
   
    hs[j]->Add(h4[j]); //ttsl
   
    hs[j]->Add(h5[j]); //tt dl ##### ttbar
    
    hs[j]->Add(h23[j]);//dy 2j
    hs[j]->Add(h25[j]);// dy 1j #####DY incl 105To160
    hs[j]->Add(h1[j]);
    

    sprintf(name,"c%i",j);
    c[j] = new TCanvas(name,name,800,800);
    /* c[j]->SetFillColor(0);
    c[j]->SetBorderMode(0);
    c[j]->SetFrameFillStyle(0);
    c[j]->SetFrameBorderMode(0);
    c[j]->SetLeftMargin( L/W );
    c[j]->SetRightMargin( R/W );
    c[j]->SetTopMargin( T/H );
    c[j]->SetBottomMargin( B/H );
    c[j]->SetTickx(0);
    c[j]->SetTicky(0);*/
    sprintf(name,"pad1%i",j);
    pad1[j]= new TPad(name, name, 0, 0.3, 1.0, 1.0);
    pad1[j]->SetBottomMargin(0.05); // Upper and lower plot are joined
    //pad1[j]->SetGridx();         // Vertical grid
    pad1[j]->Draw();             // Draw the Upper pad: pad1
    pad1[j]->cd();               // pad1 becomes the current pad
    pad1[j]->SetLogy();
    
    vFrame[j] = pad1[j]->DrawFrame(0, 0.3, 1.0, 700);
    histDraw(pad1[j],h1[j],h3[j],hs[j],vFrame[j]);
     TH1F *h_tempdata = new TH1F(name,name,nbins_temp,xlow,xup); // the histogram (you should set the number of bins, the title etc)
    auto nPoints = h3[j]->GetN(); // number of points in your TGraph
    for(int i=0; i < nPoints; ++i) {
      double x,y;
      h3[j]->GetPoint(i, x, y);
      //cout<<"Filling point:"<<i<<endl;
      h_tempdata->Fill(x,y); // ?
      double graph_err=h3[j]->GetErrorY(i);
      h_tempdata->SetBinError(i,graph_err);
    }
    h_tempdata->SetBins(nbins, dnn_bins);
    //h_tempdata->GetXaxis()->SetRangeUser(0.001,1.);
    h_tempdata->Draw("ePsames");  //data 
    h2[j]->Draw("hist same");     //signal sample
    h9[j]->Draw("hist same");
    h6[j]->Draw("hist same");
    h7[j]->Draw("hist same");
    h8[j]->Draw("hist same");
    h10[j]->Draw("hist same");
   
    
    
   
    
    
    h[j]->Sumw2();
    h[j]->Add(h12[j]);
    h[j]->Add(h5[j]);
    h[j]->Add(h13[j]);
    h[j]->Add(h17[j]);
    h[j]->Add(h14[j]);
    h[j]->Add(h15[j]);
    h[j]->Add(h18[j]);
    h[j]->Add(h11[j]);
    //h[j]->Add(h16[j]);
    h[j]->Add(h4[j]);
    h[j]->Add(h23[j]);
    //h[j]->Add(h20[j]);
    //h[j]->Add(h19[j]);
    h[j]->Add(h21[j]);
    h[j]->Add(h22[j]);
    h[j]->Add(h24[j]);
    h[j]->Add(h25[j]);
    //h[j]->Add(h26[j]);

    
    
    sprintf(name,"h_err%i",j);
    h_err[j] = (TH1F*)h[j]->Clone(name);
    Int_t n=h_err[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++){
      
	h_err[j]->SetBinContent(i,1.0);
    }
    h[j]->SetBins(nbins,dnn_bins);
    //h[j]->SetFillColorAlpha(40,0.35);
    //h[j]->SetMarkerStyle(20);
    h[j]->SetFillColor(12);
    //h[j]->SetLineWidth(20);
    h[j]->SetFillStyle(3018);
    h[j]->Draw("e2 same");
   
    vFrame[j]->GetYaxis()->SetTitle(ytitle);
    vFrame[j]->GetYaxis()->SetTitleSize(0.06);
    vFrame[j]->GetYaxis()->SetTitleOffset(0.8);
    vFrame[j]->GetYaxis()->SetLabelSize(0.05);
    vFrame[j]->GetXaxis()->SetLabelOffset(999);
    vFrame[j]->GetXaxis()->SetLabelSize(0);
    lg1->Draw();
    CMS_lumi( pad1[j], iPeriod, iPos );
    gPad->RedrawAxis();

    
    c[j]->cd();          // Go back to the main canvas before defining pad2

    sprintf(name,"pad2%i",j);
    pad2[j] = new TPad(name, name , 0, 0.0, 1.0, 0.3);
    pad2[j]->SetTopMargin(0.03);
    pad2[j]->SetBottomMargin(0.25);
    pad2[j]->SetGridx(); // vertical grid
    pad2[j]->SetGridy(); 
    pad2[j]->Draw();
    pad2[j]->cd();

   

    sprintf(name,"htemp%i",j);
   
   
    /*
    n=h_tempdata->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      if(h_tempdata->GetBinContent(i)!=0){
	bin_err =h_tempdata->GetBinError(i)/h_tempdata->GetBinContent(i);
	h_tempdata->SetBinError(i,bin_err);
      }
      }*/
    sprintf(name,"hd%i",j);
    hd[j] = (TH1F*)h_tempdata->Clone(name); // For data *************
    
    hd[j]->SetLineColor(kBlack);
    hd[j]->Sumw2();
    hd[j]->SetStats(0);// No statistics on lower plot
    hd[j]->SetMinimum(0.5);  // Define Y ..
    hd[j]->SetMaximum(1.5); // .. range
    hd[j]->Divide(h[j]); //Data/MC
    
    n=h_tempdata->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      if(h[j]->GetBinContent(i)!=0){
	bin_err =h_tempdata->GetBinError(i)/h[j]->GetBinContent(i);
	hd[j]->SetBinError(i,bin_err);
      }
      }
    hd[j]->SetMarkerStyle(21);
    hd[j]->Draw("ep");
    // hd[j]->GetXaxis()->SetRangeUser(0.001,1.);
    n=h[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      if(h[j]->GetBinContent(i)!=0){
	bin_err =h[j]->GetBinError(i)/h[j]->GetBinContent(i);
	h_err[j]->SetBinError(i,bin_err);
      }
    }
    h_err[j]->SetFillColor(40);
    //h[j]->SetLineWidth(20);
    h_err[j]->SetFillStyle(3001);
    h_err[j]->Draw("e2 same");
    hd[j]->Draw("ep same");

    
    hd[j]->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    
    hd[j]->GetYaxis()->SetTitle("Data/MC");
    hd[j]->GetYaxis()->SetNdivisions(310);
    hd[j]->GetYaxis()->SetTitleSize(30);
    hd[j]->GetYaxis()->SetTitleFont(43);
    hd[j]->GetYaxis()->SetTitleOffset(1.);
    hd[j]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[j]->GetYaxis()->SetLabelSize(20);

    // X axis ratio plot settings
    hd[j]->GetXaxis()->SetTitleSize(30);
    hd[j]->GetXaxis()->SetTitleFont(43);
    hd[j]->GetXaxis()->SetTitleOffset(2.8);
    hd[j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[j]->GetXaxis()->SetLabelSize(25);
   // X axis ratio plot settings
    hd[j]->GetXaxis()->SetTitleSize(30);
    hd[j]->GetXaxis()->SetTitleFont(43);
    hd[j]->GetXaxis()->SetTitleOffset(2.8);
    hd[j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[j]->GetXaxis()->SetLabelSize(25);

    h_err[j]->GetYaxis()->SetTitle("Data/MC");
    h_err[j]->GetYaxis()->SetNdivisions(310);
    h_err[j]->GetYaxis()->SetTitleSize(30);
    h_err[j]->GetYaxis()->SetTitleFont(43);
    h_err[j]->GetYaxis()->SetTitleOffset(1.);
    h_err[j]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_err[j]->GetYaxis()->SetLabelSize(20);

    // X axis ratio plot settings
    h_err[j]->GetXaxis()->SetTitleSize(20);
    h_err[j]->GetXaxis()->SetTitleFont(43);
    h_err[j]->GetXaxis()->SetTitleOffset(5.);
    h_err[j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_err[j]->GetXaxis()->SetLabelSize(25);
    h_err[j]->GetXaxis()->SetLabelOffset(0.02);
    h_err[j]->SetStats(0);// No statistics on lower plot
    h_err[j]->SetMinimum(0.5);  // Define Y ..
    h_err[j]->SetMaximum(1.5); // .. range
   sprintf(name,"l%i",j);
   
  }
  vFrame[0]->SetTitle("dijet mass ");
  hd[0]->GetXaxis()->SetTitle("Mjj (GeV)");
 
 
  vFrame[1]->SetTitle("pt balance ");
  //hd[1]->GetXaxis()->SetTitle(" P_{T} balance");
  hd[1]->GetXaxis()->SetTitle("DNN score");
  //vFrame[1]->GetXaxis()->SetRangeUser(0,10);
  for(int i=0;i<=32;i++){
    h_err[1]->GetXaxis()->SetBinLabel(i,bins_labels[i]);
    hd[1]->GetXaxis()->SetBinLabel(i,bins_labels[i]);
  }
  h_err[1]->GetXaxis()->LabelsOption("v");
  hd[1]->GetXaxis()->LabelsOption("v");
  //pad1[1]->SetLogx();
  //pad2[1]->SetLogx();
  //vFrame[1]->GetXaxis()->SetRangeUser(0.001,1.);
  //hd[1]->GetXaxis()->SetRangeUser(0.001,1.);
  //h_err[1]->GetXaxis()->SetRangeUser(0.001,1.);
  c[1]->Update();
  //hd[1]->GetXaxis()->SetLimits(0.,5.);
  /*
  vFrame[2]->SetTitle(" Leading Jet P_{T} ");
  hd[2]->GetXaxis()->SetTitle("j_{1} P_{T} (GeV) ");
  
  vFrame[3]->SetTitle("Dijet mass");
  hd[3]->GetXaxis()->SetTitle("M(jj) (GeV)");
 
  vFrame[4]->SetTitle("Number of jets ");
  hd[4]->GetXaxis()->SetTitle("N_{jet} ");
  
  vFrame[5]->SetTitle("Leading Jet #eta ");
  hd[5]->GetXaxis()->SetTitle("j_{1}  #eta");

  vFrame[6]->SetTitle("Sub-leading Jet #eta");
  hd[6]->GetXaxis()->SetTitle("j_{2}  #eta ");
  
  vFrame[7]->SetTitle("Sub-leading Jet P_{T}  ");
  hd[7]->GetXaxis()->SetTitle("j_{2} P_{T} (GeV) " );
  /*
  vFrame[8]->SetTitle("dilepton #phi ");
  hd[8]->GetXaxis()->SetTitle("#mu#mu #phi " );
  
  vFrame[9]->SetTitle("dilepton mass");
  hd[9]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");

  vFrame[10]->SetTitle("dilepton mass");
  hd[10]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");

  vFrame[11]->SetTitle("Dijet mass");
  vFrame[11]->GetXaxis()->SetRangeUser(400.,1000.);
  hd[11]->GetXaxis()->SetTitle("M(jj) (GeV)");
  hd[11]->GetXaxis()->SetRangeUser(400.,1000.);
  vFrame[12]->SetTitle("#Delta #eta between j_{1} and j_{2}");
  hd[12]->GetXaxis()->SetTitle("#Delta #eta between j_{1} and j_{2}");

  vFrame[13]->SetTitle("#Delta #phi between j_{1} and j_{2}");
  hd[13]->GetXaxis()->SetTitle("#Delta #phi between j_{1} and j_{2}");

  vFrame[14]->SetTitle(" Leading Jet P_{T} ");
  hd[14]->GetXaxis()->SetTitle("j_{1} P_{T} (GeV) ");
 
  vFrame[16]->SetTitle("Leading Jet #eta ");
  hd[16]->GetXaxis()->SetTitle("j_{1}  #eta");
  
  vFrame[15]->SetTitle(" Sub-leading Jet P_{T} ");
  hd[15]->GetXaxis()->SetTitle("j_{2} P_{T} (GeV) ");
 
  vFrame[17]->SetTitle("Sub-leading Jet #eta ");
  hd[17]->GetXaxis()->SetTitle("j_{2}  #eta");

  vFrame[18]->SetTitle("cos(#theta^{*})");
  hd[18]->GetXaxis()->SetTitle("cos(#theta^{*})");

  vFrame[19]->SetTitle("Min #Delta R between muon and jet ");
  hd[19]->GetXaxis()->SetTitle("Min #Delta R (#mu j)");

  vFrame[20]->SetTitle("Max #Delta R between muon and jet ");
  hd[20]->GetXaxis()->SetTitle("Max #Delta R (#mu j)");

  vFrame[21]->SetTitle("Min #Delta R between dimuon and jet ");
  hd[21]->GetXaxis()->SetTitle("Min #Delta R (#mu#mu j) ");

  vFrame[22]->SetTitle("Max #Delta R between dimuon and jet ");
  hd[22]->GetXaxis()->SetTitle("Max #Delta R (#mu#mu j)");

  vFrame[23]->SetTitle("Zeppenfeld variable");
  hd[23]->GetXaxis()->SetTitle("Zeppenfeld variable");
  
  vFrame[24]->SetTitle("Dijet+dimuon mass");
  hd[24]->GetXaxis()->SetTitle("M(#mu#mu+jj) (GeV)");

  vFrame[25]->SetTitle(" Dimuon+Dijet P_{T} ");
  hd[25]->GetXaxis()->SetTitle("P_{T} (#mu#mu+jj) (GeV)");
  
  vFrame[26]->SetTitle("#Delta #eta of dimuon+dijet");
  hd[26]->GetXaxis()->SetTitle("#Delta #eta (#mu#mu+jj) ");

  vFrame[27]->SetTitle("#Delta #phi of dimuon+dijet");
  hd[27]->GetXaxis()->SetTitle("#Delta #phi (#mu#mu+jj)");

  vFrame[28]->SetTitle(" Leading Jet QGL ");
  hd[28]->GetXaxis()->SetTitle("jet 1 qgl");

  vFrame[29]->SetTitle(" sub-leading Jet QGL ");
  hd[29]->GetXaxis()->SetTitle("jet 2 qgl");
  vFrame[30]->SetTitle("# of EWK jets with P_{T} > 5 GeV ");
  hd[30]->GetXaxis()->SetTitle("# of soft EWK jets with P_{T} > 5 GeV");
  */
  for(int j=0;j<2;j++){
    
    pad2[j]->cd();
    l[j]= new TLine(1,1.0,32,1.0);
    l[j]->SetLineColor(kBlack); l[j]->SetLineWidth(2); l[j]->SetLineStyle(1);
    l[j]->Draw("same");
   
   
    //c[j]->cd();
    //c[j]->Update();
    //c[j]->RedrawAxis();
    //c[j]->GetFrame()->Draw();
    /* name_unc=plotname[j]+".pdf";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".png";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".C";
    c[j]->SaveAs(name_unc);
    */

  }
  for(int j=0;j<2;j++){
    //c[j]->SaveAs(name_unc);
    name_unc=chname[j]+".png";
    c[j]->SaveAs(name_unc);
    //name_unc=plotname[j]+".pdf";
    //c[j]->SaveAs(name_unc);
    //name_unc=plotname[j]+".C";
    //c[j]->SaveAs(name_unc);
  }
  

}


//Define the function decorate for histograms
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);
  
  h->SetMarkerSize(1.3);
  h->SetTitle(title);
  h->SetMinimum(0.00001);
}
//Define the function decorate for legends
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(100);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //g->SetHeader(legendheader);
}
/*
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame )
{
  
  double maxCont = -1.0;
  maxCont=h1a->GetBinContent(h1a->GetMaximumBin());
  if(maxCont<h1b->GetBinContent(h1b->GetMaximumBin()))maxCont=h1b->GetBinContent(h1b->GetMaximumBin());
  maxCont*=100.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.01, maxRange, maxCont);
  hs->Draw("same hist");
  }*/
void histDraw(TPad *p,TH1F *h1a,TGraphAsymmErrors *h1b, THStack *hs,TH1F *vFrame )
{
  
  double maxCont = -1.0;
  maxCont=h1a->GetBinContent(h1a->GetMaximumBin());
  if(maxCont<h1b->GetMaximum())maxCont=h1b->GetMaximum();
  maxCont*=1000.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.01, maxRange, maxCont);
  hs->Draw("same hist");
}
// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
/*
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);

  return nevents;
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);

  return nevents;
}
*/
void h12ascii (TH1* h)
{
   Int_t n = h->GetNbinsX();
   
   for (Int_t i=1; i<=n; i++) {
      printf("%g %g\n",
             h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
             h->GetBinContent(i));
   }
}
