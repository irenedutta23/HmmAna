import os, sys
import ROOT as rt
from root_numpy import root2array, tree2array
#from root_pandas import read_root
import h5py 

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold

from xgboost import XGBClassifier
from xgboost import plot_tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_curve
import keras
import numpy as np
import numpy.lib.recfunctions as nlr
import pandas as pd

from matplotlib import pyplot as plt
import math
from math import pow

from keras.models import Sequential, Model
from keras.layers import Dense, Activation,Dropout, Input, BatchNormalization
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
try:
    import setGPU
except:
    os.environ['KERAS_BACKEND'] = 'tensorflow'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(sys.argv[3])


import advTrain_2018_defFunc

from keras.models import load_model
def convert_Z_data(tree, target=0): 
    feature = tree2array(tree,
                         branches = variables+mass_var+id_variables+wt_variables ,
                         #branches = variables,
                         selection = 'cat_index==7 && (Higgs_mass<120. || Higgs_mass > 130.)')
    if target == 0:
        label = np.zeros(shape = feature.shape, dtype=[('label','f4')])
    else:
        label = np.ones(shape = feature.shape, dtype=[('label','f4')])
        #data = nlr.merge_arrays([label,feature], flatten=True)
        #auxInfo = tree2array(tree,
        #                     branches = ['mass_jj','mass_gg'])
    return feature

def convert_Z(tree, target=0): 
    feature = tree2array(tree,
                         branches = variables+mass_var+id_variables+wt_variables ,
                         #branches = variables,
                         selection = 'cat_index==7')
    if target == 0:
        label = np.zeros(shape = feature.shape, dtype=[('label','f4')])
    else:
        label = np.ones(shape = feature.shape, dtype=[('label','f4')])
        #data = nlr.merge_arrays([label,feature], flatten=True)
        #auxInfo = tree2array(tree,
        #                     branches = ['mass_jj','mass_gg'])
    return feature


def fill_discriminator(oldTree, newTree1, newTree2,newTree3,disc_lookup, disc_lookup_adv1,disc_lookup_adv2,disc_lookup_adv3, final_evt_weight,s,s_adv1,s_adv2,s_adv3,s_wt):
    """
    Args:
        oldTree: tree from the input ROOT file
        newTree: clone of the input tree
        disc_lookup: dictionary of the form {(run, luminosityBlock, event):disc, ...}
        s: struct linked to discriminator branch in newTree
    Returns: None
    Fills newTree with all events from oldTree, with discriminator
    values filled from the lookup table where present.
    """
    num_entries = oldTree.GetEntries()
    for i in range(num_entries):
        oldTree.GetEntry(i)
        if i % 10000 == 0:
            print "Processing event {} of {}".format(i, num_entries)
        s.disc_simpleNN = advTrain_2018_defFunc.get_disc_or_minusone(oldTree, disc_lookup)
        s_adv1.disc_advNN = advTrain_2018_defFunc.get_disc_or_minusone(oldTree, disc_lookup_adv1)
        s_adv2.disc_advNN = advTrain_2018_defFunc.get_disc_or_minusone(oldTree, disc_lookup_adv2)
        s_adv3.disc_advNN = advTrain_2018_defFunc.get_disc_or_minusone(oldTree, disc_lookup_adv3)
        s_wt.evt_wt = advTrain_2018_defFunc.get_disc_or_minusone(oldTree,final_evt_weight )
        newTree1.Fill()
        newTree2.Fill()
        newTree3.Fill()


def output_file_nonTrain(y_pred_cls,y_pred_adv1_cls,y_pred_adv2_cls,y_pred_adv3_cls,evt_weight, x_test_index,z_test,INPUT_FILE,proc_tree,alpha_str1,alpha_str2,alpha_str3):
        # Create a lookup table for discriminator value by event number
        disc_lookup_signal = {}
        disc_lookup_signal_adv1 = {}
        disc_lookup_signal_adv2 = {}
        disc_lookup_signal_adv3 = {}
        proc_evt_weight={}
        sum=0
        h_sig_simple=rt.TH1F("simpleNN_VBF","simpleNN_VBF",50,0.,1.);
        h_sig_adv1=rt.TH1F("advNN_VBF","advNN_VBF",50,0.,1.);
        h_sig_adv2=rt.TH1F("advNN_VBF","advNN_VBF",50,0.,1.);
        h_sig_adv3=rt.TH1F("advNN_VBF","advNN_VBF",50,0.,1.);

        h_sig_simple_Z=rt.TH1F("simpleNN_Z","simpleNN_Z",50,0.,1.);
        h_sig_adv1_Z=rt.TH1F("advNN_Z","advNN_Z",50,0.,1.);
        h_sig_adv2_Z=rt.TH1F("advNN_Z","advNN_Z",50,0.,1.);
        h_sig_adv3_Z=rt.TH1F("advNN_Z","advNN_Z",50,0.,1.);
        for disc_val, disc_val_adv1, disc_val_adv2,disc_val_adv3,evt_wt,mass, run, lumi, event in zip(y_pred_cls, y_pred_adv1_cls, y_pred_adv2_cls,y_pred_adv3_cls,evt_weight,z_test, x_test_index[:,0],x_test_index[:,1], x_test_index[:,2]):
        
            disc_lookup_signal[(run, lumi, event)] = disc_val
            disc_lookup_signal_adv1[(run, lumi, event)] = disc_val_adv1
            disc_lookup_signal_adv2[(run, lumi, event)] = disc_val_adv2
            disc_lookup_signal_adv3[(run, lumi, event)] = disc_val_adv3
            proc_evt_weight[(run, lumi, event)] = evt_wt
            if(mass>110. and mass < 150. and  (mass<120. or mass>130.)):
                h_sig_simple.Fill(disc_val,evt_wt)
                h_sig_adv1.Fill(disc_val_adv1, evt_wt)
                h_sig_adv2.Fill(disc_val_adv2, evt_wt)
                h_sig_adv3.Fill(disc_val_adv3, evt_wt)
            if (mass>76. and mass < 106. ):
                h_sig_simple_Z.Fill(disc_val,evt_wt)
                h_sig_adv1_Z.Fill(disc_val_adv1, evt_wt)
                h_sig_adv2_Z.Fill(disc_val_adv2, evt_wt)
                h_sig_adv3_Z.Fill(disc_val_adv3, evt_wt)
            sum+=evt_wt
        print "Sum of evt_Weight: ",sum
            # We write out the signal and background events to a new ROOT file.
            # For events in the test set, we append the GBM discriminator.
            # For events in the train set, we set the discriminator to -1
            # to avoid accidentally reusing the training set in the analysis.
        print len(disc_lookup_signal), len(disc_lookup_signal_adv1)
        #print INPUT_FILE
        #print os.path.basename(INPUT_FILE)
        out_name = os.path.basename(INPUT_FILE).replace('.root', alpha_str1+'_NNscore.root')
        out_proc1 = rt.TFile(out_name, 'RECREATE')
        out_proc_tree1 = proc_tree.CloneTree(0)
        out_name = os.path.basename(INPUT_FILE).replace('.root', alpha_str2+'_NNscore.root')
        out_proc2 = rt.TFile(out_name, 'RECREATE')
        out_proc_tree2 = proc_tree.CloneTree(0)
        out_name = os.path.basename(INPUT_FILE).replace('.root', alpha_str3+'_NNscore.root')
        out_proc3 = rt.TFile(out_name, 'RECREATE')
        out_proc_tree3 = proc_tree.CloneTree(0)
            
        # This is the boilerplate code for writing something
        # to a ROOT tree from Python.
        
        rt.gROOT.ProcessLine("struct MyStruct2{float disc_advNN;};")
        rt.gROOT.ProcessLine("struct MyStruct3{float disc_simpleNN;};")
        rt.gROOT.ProcessLine("struct MyStruct4{float disc_advNN;};")
        rt.gROOT.ProcessLine("struct MyStruct5{float disc_advNN;};")
        
        rt.gROOT.ProcessLine("struct MyStruct{double evt_wt;};")
        from ROOT import MyStruct3, MyStruct2, MyStruct, MyStruct4, MyStruct5
        s = MyStruct3()
        s_adv1 = MyStruct2()
        s_adv2 = MyStruct4()
        s_adv3 = MyStruct5()
        s_wt = MyStruct()
        disc_branch_proc1 = out_proc_tree1.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
        disc_branch_proc2 = out_proc_tree2.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
        disc_branch_proc3 = out_proc_tree3.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
        disc_branch_sig_proc1 = out_proc_tree1.Branch('disc_advNN', rt.AddressOf(s_adv1, 'disc_advNN'), 'disc_advNN/F');
        disc_branch_sig_proc2 = out_proc_tree2.Branch('disc_advNN', rt.AddressOf(s_adv2, 'disc_advNN'), 'disc_advNN/F');
        disc_branch_sig_proc3 = out_proc_tree3.Branch('disc_advNN', rt.AddressOf(s_adv3, 'disc_advNN'), 'disc_advNN/F');
        evt_weight_branch_proc1 = out_proc_tree1.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
        evt_weight_branch_proc2 = out_proc_tree2.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
        evt_weight_branch_proc3 = out_proc_tree3.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
        print "Writing new ROOT process file with discriminator appended"
        fill_discriminator(proc_tree, out_proc_tree1, out_proc_tree2,out_proc_tree3,disc_lookup_signal,disc_lookup_signal_adv1,disc_lookup_signal_adv2,disc_lookup_signal_adv3, proc_evt_weight,s, s_adv1,s_adv2,s_adv3,s_wt)
        
        # Cristian's code uses GetCurrentFile() for this part.
        # I will do that too just in case (paranoia).
        out_proc1.cd()
        h_sig_simple.Write()
        h_sig_adv1.Write()
        h_sig_simple_Z.Write()
        h_sig_adv1_Z.Write()
        out_proc_tree1.GetCurrentFile().Write()
        #signalNevents.Write()
        out_proc_tree1.GetCurrentFile().Close()


        out_proc2.cd()
        h_sig_simple.Write()
        h_sig_adv2.Write()
        h_sig_simple_Z.Write()
        h_sig_adv2_Z.Write()
        out_proc_tree2.GetCurrentFile().Write()
        #signalNevents.Write()
        out_proc_tree2.GetCurrentFile().Close()

        out_proc3.cd()
        h_sig_simple.Write()
        h_sig_adv3.Write()
        h_sig_simple_Z.Write()
        h_sig_adv3_Z.Write()
        out_proc_tree3.GetCurrentFile().Write()
        #signalNevents.Write()
        out_proc_tree3.GetCurrentFile().Close()
        print 'done writing output'

def data_writer(fileName, x_mean,x_std,dataType='MC'):
    
    print fileName
    Sample_file = rt.TFile.Open(fileName)
    sample_tree = Sample_file.Get("cattree")
    #datanp1 = advTrain_2018_defFunc.convert(sample_tree,1)
    if(dataType=='data'):
        datanp1 = convert_Z_data(sample_tree)
    else:
        datanp1 = convert_Z(sample_tree)
    data_frame1 = pd.DataFrame.from_records(datanp1)
    
    data_frame1.loc[:,variables] = (data_frame1.loc[:,variables]-x_mean)/x_std
    data1 = data_frame1.values
    #print len(data)
    data_reduced1 = data1[:,var_indices]
    data_mass1=data1[:,mass_indices]
    
    data_pred_simple = m.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv = m_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv1 = m1_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv2 = m2_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    #output_file_nonTrain(data_pred_simple,data_pred_adv,data[:,wt_var_indices], data[:,id_var_indices],Sample_file,sample_tree,'0p51em4')
    #print var_indices, mass_indices,wt_var_indices, id_var_indices
    output_file_nonTrain(data_pred_simple,data_pred_adv,data_pred_adv1,data_pred_adv2,data1[:,wt_var_indices], data1[:,id_var_indices],data_mass1,fileName,sample_tree,'1em4','5em5','1em5')


m = load_model('simple_cls_DYnVBF.h5')
m_cls = load_model('advCls_1em4_DYnVBF.h5')
m1_cls = load_model('advCls_1em5_DYnVBF.h5')
m2_cls = load_model('advCls_5em5_DYnVBF.h5')

print(keras.__version__)
variables = ['softJet5','dRmm','dEtamm','dPhimm','M_jj','pt_jj','eta_jj','phi_jj','M_mmjj',#'pt_mmjj',
             'eta_mmjj','phi_mmjj','dEta_jj','Zep','dRmin_mj', 'dRmax_mj'
                                   ,'dRmin_mmj','dRmax_mmj',#'dPhijj',
             'leadingJet_pt','subleadingJet_pt',
                                   'leadingJet_eta','subleadingJet_eta','leadingJet_qgl','subleadingJet_qgl','cthetaCS','Higgs_pt','Higgs_eta']#,'Higgs_mass' ]
mass_var=['Higgs_mass']
id_variables = ['run','lumi','event']
wt_variables = ['genWeight']

SIGNAL_FILE='/bigdata/shared/idutta/VBFHToMuMu_2018.root'
BKG_FILE_2018='/bigdata/shared/idutta/DYJetsToLL_VBFfilter_2018.root'
sigfile = rt.TFile.Open(SIGNAL_FILE)
bkgfile_2018 = rt.TFile.Open(BKG_FILE_2018)
sigtree = sigfile.Get("cattree")
bkgtree3 = bkgfile_2018.Get("cattree")
signp =advTrain_2018_defFunc. convert(sigtree, 1)
bkgnp3 = advTrain_2018_defFunc.convert(bkgtree3, 0)

sig_frame0 = pd.DataFrame.from_records(signp)
bkg_frame3= pd.DataFrame.from_records(bkgnp3)

  
var_indices = [sig_frame0.columns.get_loc(v) for v in variables] # get positions of all the variables set above
mass_indices = [sig_frame0.columns.get_loc(v) for v in mass_var]
id_var_indices = [sig_frame0.columns.get_loc(v) for v in id_variables]
wt_var_indices = [sig_frame0.columns.get_loc(v) for v in wt_variables]

# Standardize only the actual variables not event ID or genweights
x_mean_sig0 = (sig_frame0.loc[:,variables]).mean()
x_std_sig0 = (sig_frame0.loc[:,variables]).std()

x_mean_bkg3 = (bkg_frame3.loc[:,variables]).mean()
x_std_bkg3 = (bkg_frame3.loc[:,variables]).std()

inputfname='/nfshome/idutta/NN_bkg_files_2016_jes.txt'
with open(inputfname) as inputfile:
    for line in inputfile:
        line=line.rstrip()
        print line.rstrip()
        data_writer(line,x_mean_bkg3,x_std_bkg3,'MC')
inputfname='/nfshome/idutta/NN_sig_files_2016_jes.txt'
with open(inputfname) as inputfile:
    for line in inputfile:
        line=line.rstrip()
        print line.rstrip()
        data_writer(line,x_mean_sig0,x_std_sig0,'MC')

data_writer('/bigdata/shared/idutta/2016_files_jes/Data_2016.root',x_mean_bkg3,x_std_bkg3,'data')
