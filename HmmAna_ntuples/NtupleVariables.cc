//////////////////////////////////////////////////////////
// 
// Original Author : Irene Dutta
//                   Caltech
// Date Created    : Mon 27 Aug, 2018
//////////////////////////////////////////////////////////
#define NtupleVariables_cxx

#include "NtupleVariables.h"
#include <TH2.h>
#include <TStyle.h>
#include<iostream>

double NtupleVariables::DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double NtupleVariables::DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

double NtupleVariables::HT (std::vector<TLorentzVector> vjets) {
  double ht = 0.0;
  for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
    //if( vjets[ijet].Pt()>50.0 && std::abs(vjets[ijet].Eta())<2.5 ) 
      ht += vjets[ijet].Pt();
  }
  return ht;
}
TLorentzVector NtupleVariables::MHT(std::vector<TLorentzVector> vjets) {;
  TLorentzVector mht(0.0, 0.0, 0.0, 0.0);
  for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
    if( vjets[ijet].Pt()>30.0 && std::abs(vjets[ijet].Eta())<5.0 ) 
      mht -= vjets[ijet];
  }

  return mht;
}

std::pair<double,double> NtupleVariables::helicityAngles( TLorentzVector v1, TLorentzVector v2) {

      TLorentzVector H = v1 + v2;
      v1.Boost( -( H.BoostVector() ) ); // go to Higgs RFR
      v2.Boost( -( H.BoostVector() ) );

      TVector3 reco_M1 = v1.Vect().Unit();
      TVector3 reco_M2 = v2.Vect().Unit();

      //Costh*
      //double cthstr=reco_M1.Z();
      double phi = TMath::ATan2(reco_M1.Y(), reco_M1.X());
      double theta = TMath::ACos(reco_M1.Z());
      return std::pair<double,double>(theta,phi);
}

std::pair<double,double> NtupleVariables::CSAngles( TLorentzVector v1, TLorentzVector v2, int charge) {
//https://github.com/alisw/AliPhysics/blob/master/PWGDQ/dielectron/core/AliDielectronPair.cxx
    float fBeamEnergy = 6500.0;
    float proMass = 0.938272;
    TLorentzVector pro1(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+proMass*proMass));
    TLorentzVector pro2(0.,0., fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+proMass*proMass));
    
    TLorentzVector H = v1 + v2;
    v1.Boost( -( H.BoostVector() ) ); // go to Higgs RFR
    v2.Boost( -( H.BoostVector() ) );
    pro1.Boost( -( H.BoostVector() ) );
    pro2.Boost( -( H.BoostVector() ) );
    
    TVector3 yAxis = ((pro1.Vect()).Cross(pro2.Vect())).Unit();
    TVector3 zAxisCS = ((pro1.Vect()).Unit()-(pro2.Vect()).Unit()).Unit();
    TVector3 xAxisCS = (yAxis.Cross(zAxisCS)).Unit();
    
    TVector3 reco_M1 = v1.Vect().Unit();
    TVector3 reco_M2 = v2.Vect().Unit();
    double thetaCS = zAxisCS.Dot((v1.Vect()).Unit());
    double phiCS   = TMath::ATan2((v1.Vect()).Dot(yAxis), (v1.Vect()).Dot(xAxisCS));
    if(charge<0){
        thetaCS = zAxisCS.Dot((v2.Vect()).Unit());
        phiCS   = TMath::ATan2((v2.Vect()).Dot(yAxis), (v2.Vect()).Dot(xAxisCS));
    }
    return std::pair<double,double>(thetaCS,phiCS);
}
