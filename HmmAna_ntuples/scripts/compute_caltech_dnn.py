import ROOT as rt
import numpy as np
import math
from keras.models import Model
from keras.models import load_model
import os, sys
try:
    import setGPU
except:
    os.environ['KERAS_BACKEND'] = 'tensorflow'

def DeltaPhi(phi1,phi2):
    result = phi1 - phi2
    while(result > M_PI):  
        result -= 2 * M_PI
    while(result <= -M_PI):
        result += 2 * M_PI

    return result

def deltar(obj1, obj2):
    deta = obj1.Eta() - obj2.Eta()
    dphi = DeltaPhi(obj1.Phi(),obj2.Phi())

    dr = math.sqrt(deta**2 + dphi**2)
    
    return deta, dphi, dr 

def dnn_variables(leading_muon, subleading_muon, leading_jet, subleading_jet, nsoft, j1_qgl, j2_qgl): #leading_muon, subleading_muon, leading_jet, subleading_jet are all rt.Tlorentzvectors
    #delta eta, phi and R between two muons
    mm_deta, mm_dphi, mm_dr = deltar(leading_muon, subleading_muon)
    
    #delta eta between jets 
    jj_deta = leading_jet.Eta() - subleading_jet.Eta()
    jj_dphi = DeltaPhi(leading_jet.Phi(),subleading_jet.Phi())
   
    
    #create dimuon system 
    mm = leading_muon + subleading_muon

    #create dijet system 
   
    jj = leading_jet+ subleading_jet
    #create dimuon-dijet system 
    
    mmjj = jj+mm
   

    #compute deltaR between all muons and jets
    dr_mjs = []
    for mu in [leading_muon, subleading_muon]:
        for jet in [leading_jet, subleading_jet]:
            _, _, dr_mj = deltar(mu, jet)
            dr_mjs += [dr_mj]
    dr_mj = np.vstack(dr_mjs)
    dRmin_mj = np.min(dr_mj, axis=0) 
    dRmax_mj = np.max(dr_mj, axis=0) 
    
    #compute deltaR between dimuon system and both jets 
    dr_mmjs = []
    for jet in [leading_jet, subleading_jet]:
        _, _, dr_mmj = deltar(mm, jet)
        dr_mmjs += [dr_mmj]
    dr_mmj = np.vstack(dr_mmjs)
    dRmin_mmj = np.min(dr_mmj, axis=0) 
    dRmax_mmj = np.max(dr_mmj, axis=0)

    #Zeppenfeld variable
    Zep = (mm.Eta() - 0.5*(leading_jet.Eta() + subleading_jet.Eta()))

    #Collin-Soper frame variable
    cthetaCS = 2*(leading_muon.Pz() * subleading_muon.E() - leading_muon.E()*subleading_muon.Pz()) / (mm.M() * np.sqrt(np.power(mm.M(), 2) + np.power(mm.Pt(), 2)))
    ret = {
        "dEtamm": mm_deta, "dPhimm": mm_dphi, "dRmm": mm_dr,
        "M_jj": jj.M(), "pt_jj": jj.Pt(), "eta_jj": jj.Eta(), "phi_jj": jj.Phi(),
        "M_mmjj": mmjj.M(), "eta_mmjj": mmjj.Eta(), "phi_mmjj": mmjj.Phi(),
        "dEta_jj": jj_deta,
        "dPhi_jj": jj_dphi,
        "leadingJet_pt": leading_jet.Pt(),
        "subleadingJet_pt": subleading_jet.Pt(),
        "leadingJet_eta": leading_jet.Eta(),
        "subleadingJet_eta": subleading_jet.Eta(),
        "dRmin_mj": dRmin_mj,
        "dRmax_mj": dRmax_mj,
        "dRmin_mmj": dRmin_mmj,
        "dRmax_mmj": dRmax_mmj,
        "Zep": Zep,
        "leadingJet_qgl": j1_qgl,
        "subleadingJet_qgl": j2_qgl, 
        "cthetaCS": cthetaCS,
        "softJet5": nsoft,
        "Higgs_pt": mm.Pt(),
        "Higgs_eta": mm.Eta(),
        "Higgs_mass": mm.M(),
    }


    for k in ret.keys():
        msk = np.isnan(ret[k])
        if np.sum(msk) > 0:
            print("dnn vars nan", k, np.sum(msk))

    return ret

def calc_dnn_pred(j1,j2,m1,m2,n_softJet,j1_qgl,j2_qgl): #j1,j2,m1,m2 are all rt.TLorentzvectors
    dnn_vars = dnn_variables(m1,m2,j1,j2,n_sofjet,j1_qgl,j2_qgl)
    varlist_order = ['softJet5', 'dRmm','dEtamm','M_jj','pt_jj','eta_jj','phi_jj','M_mmjj','eta_mmjj','phi_mmjj','dEta_jj','Zep','dRmin_mj', 'dRmax_mj'
                               ,'dRmin_mmj','dRmax_mmj','dPhimm','leadingJet_pt','subleadingJet_pt',
                               'leadingJet_eta','subleadingJet_eta','leadingJet_qgl','subleadingJet_qgl','cthetaCS','Higgs_pt','Higgs_eta','Higgs_mass']
    dnn_vars_arr = np.vstack([dnn_vars[k] for k in varlist_order]).T

    model_STDwt = np.load('DNN27vars_sig_vbf_bkg_dyvbf_dy105To160_ewk105To160_split_60_40__mod10_190820.npy')
    dnn_mean = model_STDwt[0,:]
    dnn_std = model_STDwt[1,:]

    dnn_vars_arr-=dnn_mean
    dnn_vars_arr/=dnn_std

    m_caltech = load_model('DNN27vars_sig_vbf_bkg_dyvbf_dy105To160_ewk105To160_split_60_40_mod10_190820.h5')
    caltech_dnn_score=m_caltech.predict(dnn_vars_arr).ravel()
    return caltech_dnn_score
