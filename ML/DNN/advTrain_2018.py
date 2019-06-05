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


HLF = ['softJet5','dRmm','dEtamm','dPhimm','M_jj','pt_jj','eta_jj','phi_jj','M_mmjj',#'pt_mmjj',
       'eta_mmjj','phi_mmjj','dEta_jj','Zep','dRmin_mj', 'dRmax_mj','dRmin_mmj','dRmax_mmj',#'dPhijj',
        'leadingJet_pt','subleadingJet_pt','leadingJet_eta','subleadingJet_eta','leadingJet_qgl','subleadingJet_qgl','cthetaCS','Higgs_pt',
        'Higgs_eta','Higgs_mass']
xlabel_HLF= ['Number of EWK jets ($p_{T}$ > 5 GeV)','$\Delta R(\mu \mu)$','$\Delta \eta(\mu \mu)$','$\Delta \phi(\mu \mu)$','$M_{jj} (GeV)$','$p_{T}$ (jj) (GeV)','$\eta_{jj}$','$\phi_{jj}$','$M_{\mu\mu + jj}$ (GeV)','$p_{T} (\mu\mu + jj) $ (Gev)','$\eta_{\mu\mu + jj}$','$\phi_{\mu\mu + jj}$','$\Delta \eta(jj)$','Zeppenfeld Variable','Min $\Delta R(\mu j)$ ', 'Max $\Delta R(\mu j)$','Min $\Delta R(\mu \mu j)$','Max $\Delta R(\mu \mu j)$','$\Delta \phi(jj)$',
             'leading Jet $p_{T}$ (GeV)','sub-leading Jet $p_{T}$ (GeV)','leading Jet $\eta$','sub-leading Jet $\eta$','leading Jet QGL','sub-leading Jet QGL','cos($\\theta_{CS}^{*}$)','Dimuon $p_{T}$ (GeV)',
        'Dimuon $\eta$','Higgs_mass']
for hlf,xlabel_hlf in zip(HLF,xlabel_HLF):
    plt.figure()
    plt.hist(sig_frame0[hlf], bins=40, normed=True, histtype='step', label='VBF signal')
    plt.hist(bkg_frame3[hlf], bins=40, normed=True, histtype='step', label='DY+jets VBF filter')
    plt.xlabel(xlabel_hlf)
    plt.ylabel('Events (Normalized to unity)')
    plt.legend(loc='best')
    plt.savefig("/bigdata/shared/idutta/Higgs_plots/SigBkgCompare_"+hlf+".png")

# Standardize only the actual variables not event ID or genweights
x_mean_sig0 = (sig_frame0.loc[:,variables]).mean()
x_std_sig0 = (sig_frame0.loc[:,variables]).std()

x_mean_bkg3 = (bkg_frame3.loc[:,variables]).mean()
x_std_bkg3 = (bkg_frame3.loc[:,variables]).std()

sig_frame0.loc[:,variables] = (sig_frame0.loc[:,variables]-x_mean_sig0)/x_std_sig0
bkg_frame3.loc[:,variables] = (bkg_frame3.loc[:,variables]-x_mean_bkg3)/x_std_bkg3

#frames_sig = [sig_frame0,sig_frame3,sig_frame1,sig_frame2]
frames_sig = [sig_frame0]
sig_frame = pd.concat(frames_sig)
signal = sig_frame.values
background3 = bkg_frame3.values

#frames = [bkg_frame,bkg_frame0,bkg_frame1,bkg_frame2]
#frames = [bkg_frame3,bkg_frame,bkg_frame0]
frames = [bkg_frame3]
background_full_frame = pd.concat(frames)
background_full = background_full_frame.values

DY_label = np.zeros(len(bkg_frame3.values)) #useful later for evt_wt of test events
#tt2l2v_label = np.ones(len(bkg_frame.values)) #useful later for evt_wt of test events
#ttsl_label = 2*np.ones(len(bkg_frame0.values)) #useful later for evt_wt of test events

bkg_process_label=[]
for i in range(0,len(bkg_frame3.values)):
    bkg_process_label.append(0)
#for i in range(0,len(bkg_frame.values)):
#    bkg_process_label.append(2)
#for i in range(0,len(bkg_frame0.values)):
#    bkg_process_label.append(3)
    
#print "After standardizing"
print "length of Signal :",len(sig_frame.values)
print "length of background DY Jet VBF filter",len(background_full_frame.values)


print "VBF", np.sum(sig_frame0.values[:,wt_var_indices])
scale_VBF=np.sum(sig_frame0.values[:,wt_var_indices]) #useful later for evt_wt o

scale_DY_VBFfilter = np.sum(bkg_frame3.values[:,wt_var_indices])#useful later for evt_wt of test events


randix4 = np.arange(len(background_full))
np.random.shuffle(randix4)


randix2 = np.arange(len(signal))
np.random.shuffle(randix2)
signal_ = signal[randix2]

background_full_ = background_full[randix4]
bkg_process_label = np.array( bkg_process_label)
bkg_process_label_ = bkg_process_label[randix4]

randix5 = np.arange(len(background3))
np.random.shuffle(randix5)
background_VBFfilter_ = background3[randix5]

#background_full_ = background_full_[:len(signal_)]


signal_wt = signal_[:,wt_var_indices]
bkg_wt_full = background_full_[:,wt_var_indices]
bkg3_wt = background3[:,wt_var_indices]
print "sum of weights for signal before scaling: " 
print np.sum(signal_wt)
print "sum of weights for reduced background length before scaling: "
print np.sum(bkg_wt_full)
sig_scale= np.sum(bkg_wt_full)/np.sum(signal_wt)
#sig_scale= 1/np.sum(signal_wt)
signal_[:,wt_var_indices] = np.multiply(signal_wt,sig_scale)
signal_wt = signal_[:,wt_var_indices]
#background_full_[:,wt_var_indices]=1.
bkg_wt_full = background_full_[:,wt_var_indices]

print "sum of weights for signal after scaling: " 
print np.sum(signal_wt)
print "sum of weights for bkg after scaling: " 
print np.sum(bkg_wt_full)
#print "sum of weights for DY bkg: " 
#print np.sum(bkg3_wt)



sig_label = np.ones(len(signal_))

bkg_full_label = np.zeros(len(background_full_))

data_train = np.concatenate((signal_,background_full_))
process_label= np.concatenate((sig_label,bkg_process_label_))
print len(process_label)
#label_train = np.concatenate((sig_label,bkg_VBF_label))
label_train = np.concatenate((sig_label,bkg_full_label))


randix3 = np.arange(len(data_train))
np.random.shuffle(randix3)
data_train= data_train[randix3,...]
process_label=process_label[randix3,...]
label_train = label_train[randix3,...]



print data_train.shape
print label_train.shape



seed = 7
test_size = 0.3
x_train, x_test, y_train, y_test, process_train, process_test = train_test_split(data_train, label_train, process_label, test_size=test_size, random_state=seed)
print y_test.shape
is_sig_test = y_test==1
is_bkg_test = np.invert(is_sig_test)

# For training we ignore the columns with the event ID information
x_train_reduced = x_train[:,var_indices]
x_test_reduced = x_test[:,var_indices]
z_train=x_train[:,mass_indices]
z_test = x_test[:,mass_indices]
   
x_test_index = x_test[:,id_var_indices]
x_train_index = x_train[:,id_var_indices]

x_train_wt = x_train[:,wt_var_indices]
x_test_wt = x_test[:,wt_var_indices]

sig_test_wt=[]
bkg_test_wt_DY=[] #process_label 0

for i in range(0,len(x_test_wt)):
    if y_test[i]==1:
        sig_test_wt.append(x_test_wt[i])
    if y_test[i]==0:
        if(process_test[i]==0):bkg_test_wt_DY.append(x_test_wt[i])
print len(sig_test_wt)
sig_sum = np.sum(sig_test_wt)
DY_sum = np.sum(bkg_test_wt_DY)


print scale_VBF, scale_DY_VBFfilter
sig_test_scale = scale_VBF/sig_sum
DY_test_scale = scale_DY_VBFfilter/DY_sum

evt_weight=[]
print np.sum(y_test)
for i in range(0,len(x_test_wt)):

    if(process_test[i]==0):evt_weight.append(x_test_wt[i]*DY_test_scale)
    elif(process_test[i]==1):evt_weight.append(x_test_wt[i]*sig_test_scale)
   
       
    
print len(evt_weight)
print np.sum(evt_weight)
print x_train_wt
print x_train_reduced.shape, x_test_reduced.shape
print y_train.shape, y_test.shape
print y_train.sum(), y_test.sum()
print x_train_reduced[:,0].shape
'''
f, ax_arr = plt.subplots(5, 6, figsize=(50,50))
i=0
    
for j in range(0,5):
    for k in range(0,6):
        x_train_sig = []
        x_train_bkg =[]
        x_test_sig = []
        x_test_bkg = []
        temp = x_train_reduced[:,i]
        test_tmp = x_test_reduced[:,i]
        for m in range (0,len(y_train)):
            if(y_train[m]==1): x_train_sig.append(temp[m])
            else: x_train_bkg.append(temp[m])
        for m in range (0,len(y_test)): 
            if(y_test[m]==1): x_test_sig.append(test_tmp[m])
            else: x_test_bkg.append(test_tmp[m])    
        ax_arr[j,k].hist(x_train_sig, bins=40, normed=True, histtype='step', label='VBF signal train',color='b')
        ax_arr[j,k].hist(x_test_sig, bins=40, normed=True, histtype='step', linestyle =('dashed'), label='VBF signal test', color = 'b')
        ax_arr[j,k].hist(x_train_bkg, bins=40, normed=True, histtype='step', label='DY Jets VBF filter sample train',color='r')
        ax_arr[j,k].hist(x_test_bkg, bins=40, normed=True, histtype='step', linestyle =('dashed'),label='DY Jets VBF filter sample test', color = 'r')
        ax_arr[j,k].set_xlabel(HLF[i])
        ax_arr[j,k].legend(prop={'size': 10})
        ax_arr[j,k].legend(loc='best')
        ax_arr[j,k].set_yscale("log")
        i+=1
        if i == 27: break
    if i ==27 : break
plt.savefig("/bigdata/shared/idutta/Higgs_plots/TrainTest_variables.png")  


'''


m = advTrain_2018_defFunc.make_cls_model((x_train_reduced.shape[1],),True, 5, 100, 0.1)
m.compile(loss='binary_crossentropy', optimizer='adam',metrics=['accuracy'])
h = advTrain_2018_defFunc.fit_cls_model(m, x_train_reduced,y_train, x_train_wt.ravel(),10000)

advTrain_2018_defFunc.show_losses( [("VBF entropy", h)] )




inp = Input(shape=(x_train_reduced.shape[1],))
print inp
alpha_num=float(sys.argv[1])
alpha_str=str(sys.argv[2])
m_cls, m_adv, m_clsOnly_adv_train, m_cls_advOnly_train = advTrain_2018_defFunc.make_cls_adv_model(inp,True,0.1, 5,3, 100,50, 1e-5,alpha_num , "elu")


losses_train_cls, losses_test_cls= advTrain_2018_defFunc.trainer(m_clsOnly_adv_train, m_cls_advOnly_train,x_train_reduced,x_test_reduced,y_train, z_train,y_test, z_test,x_train_wt,x_test_wt,alpha_num,alpha_str)

m.save('simple_cls_DYnVBF.h5')
m_cls.save('advCls_'+alpha_str+'_DYnVBF.h5')
# In[61]:


#from keras.models import load_model
#m = load_model('simple_cls_DYnVBF.h5')
#m_cls = load_model('advCls_0p51em4_DYnVBF.h5')
pred_simple_cls = m.predict(x_test_reduced, batch_size=10000)[:, 0]
pred_adv_cls = m_cls.predict(x_test_reduced, batch_size=10000)[:, 0]
advTrain_2018_defFunc.output_file(pred_simple_cls,pred_adv_cls,evt_weight, y_test, z_test, x_test_index,SIGNAL_FILE,sigtree,BKG_FILE_2018,bkgtree3,alpha_str)
advTrain_2018_defFunc.cls_adv_score_roc(pred_simple_cls,pred_adv_cls,is_sig_test,is_bkg_test,y_test,z_test,x_test_wt,alpha_str)

advTrain_2018_defFunc.compare_rocs(pred_simple_cls,pred_adv_cls,y_test,z_test,alpha_str)


