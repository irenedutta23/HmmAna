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
from keras.layers import Dense, Activation,Dropout, Input, BatchNormalization, Flatten
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint




# In[2]:


variables = ['softJet5','dRmm','dEtamm','dPhimm','M_jj','pt_jj','eta_jj','phi_jj','M_mmjj',
             #'pt_mmjj',
             'eta_mmjj','phi_mmjj','dEta_jj','Zep','dRmin_mj', 'dRmax_mj'
                                   ,'dRmin_mmj','dRmax_mmj',#'dPhijj',
             'leadingJet_pt','subleadingJet_pt',
                                   'leadingJet_eta','subleadingJet_eta','leadingJet_qgl','subleadingJet_qgl','cthetaCS','Higgs_pt','Higgs_eta']#,'Higgs_mass' ]
mass_var=['Higgs_mass']
id_variables = ['run','lumi','event']
wt_variables = ['genWeight']


# In[3]:


def convert(tree, target=0):
    feature = tree2array(tree,
                        branches = variables+mass_var+id_variables+wt_variables ,
                        #branches = variables,
                         selection = 'cat_index==7 && Higgs_mass>110. && Higgs_mass < 150.')
    if target == 0:
        label = np.zeros(shape = feature.shape, dtype=[('label','f4')])
    else:
        label = np.ones(shape = feature.shape, dtype=[('label','f4')])
    #data = nlr.merge_arrays([label,feature], flatten=True)
    #auxInfo = tree2array(tree,
    #                     branches = ['mass_jj','mass_gg'])
    return feature


# In[4]:


def convert_data(tree):
    feature = tree2array(tree,
                        branches = variables+mass_var+id_variables+wt_variables ,
                        #branches = variables,
                         selection = 'cat_index==7 && Higgs_mass>110. && Higgs_mass < 150. && (Higgs_mass>130. || Higgs_mass<120.)')
    return feature



# In[7]:





# In[9]:


def get_disc_or_minusone(event, disc_lookup):
    """
    Checks the event ID (run/luminosityBlock/event) from the tree (first argument)
    and gets the discriminator value corresponding to the event, 
    if available in the lookup dictionary (second argument).
    Returns -1 if the event is not available.
    """
    #return disc_lookup.get((event.run, event.luminosityBlock, event.event), -1)
    return disc_lookup.get((event.run, event.lumi, event.event), -1)


# In[10]:

def fill_discriminator(oldTree, newTree, disc_lookup, disc_lookup_adv,final_evt_weight,s,s_adv,s_wt):
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
        s.disc_simpleNN = get_disc_or_minusone(oldTree, disc_lookup)
        s_adv.disc_advNN = get_disc_or_minusone(oldTree, disc_lookup_adv)
        s_wt.evt_wt = get_disc_or_minusone(oldTree,final_evt_weight )
        newTree.Fill()


def output_file(y_pred_cls,y_pred_adv_cls,evt_weight, y_test, z_test,x_test_index,SIGNAL_FILE,sigtree,BKG_FILE_2018,bkgtree3,alpha_str):
    # Create a lookup table for discriminator value by event number
    disc_lookup_signal = {}
    disc_lookup_bkg = {}
    disc_lookup_signal_adv = {}
    disc_lookup_bkg_adv = {}
    sig_evt_weight={}
    bkg_evt_weight={}
    ys=[]
    h_sig_simple=rt.TH1F("simpleNN_VBF","simpleNN_VBF",100,0.,1.);
    h_sig_adv=rt.TH1F("advNN_VBF","advNN_VBF",100,0.,1.);
    h_bkg_simple=rt.TH1F("simpleNN_VBF","simpleNN_VBF",100,0.,1.);
    h_bkg_adv=rt.TH1F("advNN_VBF","advNN_VBF",100,0.,1.);
    print len(y_pred_cls), len(y_pred_adv_cls), len(evt_weight), len(y_test), len(x_test_index[:,0]),len(x_test_index[:,1]), len(x_test_index[:,2])
    print y_test.sum()
    
    ctr=0
    for disc_val, disc_val_adv, evt_wt, y_val, mass, run, lumi, event in zip(y_pred_cls, y_pred_adv_cls, evt_weight, y_test, z_test,x_test_index[:,0],x_test_index[:,1], x_test_index[:,2]):
        ctr+=1
        if y_val!=1. and y_val!=0.: print y_val
        if y_val == 1.:
            disc_lookup_signal[(run, lumi, event)] = disc_val
            disc_lookup_signal_adv[(run, lumi, event)] = disc_val_adv
            sig_evt_weight[(run, lumi, event)] = evt_wt
            if (mass<120. or mass>130.):
                h_sig_simple.Fill(disc_val,evt_wt)
                h_sig_adv.Fill(disc_val_adv, evt_wt)
            ys.append(1)
            #print y_val
            #print disc_val
        elif y_val == 0:
            disc_lookup_bkg[(run, lumi, event)] = disc_val
            disc_lookup_bkg_adv[(run, lumi, event)] = disc_val_adv
            bkg_evt_weight[(run, lumi, event)] = evt_wt
            if (mass<120. or mass>130.):
                h_bkg_simple.Fill(disc_val,evt_wt)
                h_bkg_adv.Fill(disc_val_adv, evt_wt)
            #print y_val
            #print disc_val
          
    # We write out the signal and background events to a new ROOT file.
    # For events in the test set, we append the GBM discriminator.  
    # For events in the train set, we set the discriminator to -1
    # to avoid accidentally reusing the training set in the analysis.
    print np.sum(ys)
    print len(ys)
    print ctr
    print len(disc_lookup_signal), len(disc_lookup_signal_adv)
    print len(disc_lookup_bkg), len(disc_lookup_bkg_adv)
    out_signal_name = os.path.basename(SIGNAL_FILE).replace('.root', alpha_str+'_NNscore_new.root')
    out_signal = rt.TFile(out_signal_name, 'RECREATE')
    out_signal_tree = sigtree.CloneTree(0)

    out_bkg_name = os.path.basename(BKG_FILE_2018).replace('.root', alpha_str+'_NNscore_new.root')
    out_bkg = rt.TFile(out_bkg_name, 'RECREATE')
    out_bkg_tree3 = bkgtree3.CloneTree(0)
    
    rt.gROOT.ProcessLine("struct MyStruct2{float disc_advNN;};")
    rt.gROOT.ProcessLine("struct MyStruct3{float disc_simpleNN;};")
    rt.gROOT.ProcessLine("struct MyStruct{double evt_wt;};")
    from ROOT import MyStruct3, MyStruct2, MyStruct
    s = MyStruct3()
    s_adv = MyStruct2()
    s_wt = MyStruct()
    disc_branch_sig = out_signal_tree.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
    disc_branch_bkg3 = out_bkg_tree3.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
   
    disc_branch_sig_adv = out_signal_tree.Branch('disc_advNN', rt.AddressOf(s_adv, 'disc_advNN'), 'disc_advNN/F');
    disc_branch_bkg_adv3 = out_bkg_tree3.Branch('disc_advNN', rt.AddressOf(s_adv, 'disc_advNN'), 'disc_advNN/F');
     
    evt_weight_branch_sig = out_signal_tree.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
    evt_weight_branch_bkg3 = out_bkg_tree3.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
   
    print "Writing new ROOT signal file with discriminator appended"
    fill_discriminator(sigtree, out_signal_tree, disc_lookup_signal,disc_lookup_signal_adv, sig_evt_weight,s, s_adv,s_wt)
    print "Writing new ROOT background file with discriminator appended"
    fill_discriminator(bkgtree3, out_bkg_tree3, disc_lookup_bkg,disc_lookup_bkg_adv, bkg_evt_weight,s, s_adv,s_wt)

    
    out_signal.cd()
    h_sig_simple.Write()
    h_sig_adv.Write()
    out_signal_tree.GetCurrentFile().Write()
    out_signal_tree.GetCurrentFile().Close()

    out_bkg.cd()
    h_bkg_simple.Write()
    h_bkg_adv.Write()
    out_bkg_tree3.GetCurrentFile().Write()
    out_bkg_tree3.GetCurrentFile().Close()
    
    print 'done writing output'


# In[12]:


def make_cls_model(inp, b_norm = False, n_hidden = 1, hidden= 100, do = 0.2):
    #i = Input((26,), name='features')
    i = inp
    di = i
    for h in range(n_hidden):
        d = Dense(hidden, name='hidden_%d'%h, activation='tanh')(di)
        if do:
            d = Dropout(do)(d)
        if b_norm:
            d = BatchNormalization()(d)
        di = d
    o = Dense(1, name='categorization', activation='sigmoid')(d)
    model = Model(i,o,name='classifier')
    model.summary()
    return model


# In[13]:


def fit_cls_model(model,x_train_reduced,y_train, x_train_wt,batch_size = 20, val_split=0.2):
    hist = model.fit(x_train_reduced, y_train, sample_weight= x_train_wt,batch_size=batch_size, nb_epoch=500,verbose=0,validation_split=val_split,
                  callbacks=[EarlyStopping(monitor='val_loss', patience=5, verbose=1, mode='min'),
                            ModelCheckpoint(filepath='model.h5', verbose=0)])
    return hist


# In[14]:


def make_adv_model(b_norm, nlayers_adv=3, nunits_adv=50, activation='elu',do=0.2, input_size=(1,)):

    i = Input(input_size, name='r_input')
    #flatten the input if need be
    if len(input_size)>1: i = Flatten()(i)
    di = i
    for h in range(nlayers_adv):
        d = Dense(nunits_adv, name='hidden_%d'%h, activation=activation)(di)
        if do:
            d = Dropout(do)(d)
        if b_norm:
            d = BatchNormalization()(d)
        di = d
    o = Dense(1,name='mass_predict')(d)
    model = Model(i,o,name='adversarial')
    model.summary()
    return model

# In[15]:


def make_cls_adv_model(inp, b_norm, do, nlayers, nlayers_adv,nunits,nunits_adv, lr, alpha, activation):

    model_cls = make_cls_model(inp,b_norm,nlayers, nunits, do)  #build the classifier model

    sc = model_cls(inp)
    feed_on_layer = -1 ## feed on (1,) the discriminator output
    feed_on_layer = -2 ## feed on (100,) the layer just above
    layer_to_adv = model_cls.layers[feed_on_layer]    
    shape_to_adv = layer_to_adv.output_shape[1:]
    print ("the input shape to the adversary is ",shape_to_adv)
    
    model_adv = make_adv_model(b_norm,nlayers_adv, nunits_adv, activation,do, input_size = shape_to_adv) #build the regression model
   
    #Construct a combined model of input (inp) -> classifier (model_cls) -> classifier output (sc) -> regression (model_cls_adv) -> regression output (out)
    input_for_adv = Model(input = inp,
          output = layer_to_adv.output)(inp)
    out = model_adv(input_for_adv)
    
    opt = keras.optimizers.adam(lr=lr)
    model_cls.compile(loss='binary_crossentropy', optimizer=opt,metrics=['accuracy'])
    #Set the regression model (adversarial part) weights to be fixed
    model_adv.trainable = False
    for l in model_adv.layers:
        l.trainable = False

    print inp, sc, out
    print [sc,out]
    #Construct the combined model (classifier -> regression), where only the classifier is trainable 
    model_clsOnly_adv_train = Model(inp, [sc, out])
    model_clsOnly_adv_train.compile(loss=['binary_crossentropy', 'mse'], optimizer=opt, loss_weights=[1.0, -alpha])
    model_clsOnly_adv_train.summary()

    #Set the regression part to be trainable, classifier to be fixed
    model_cls.trainable = False
    for l in model_cls.layers:
        l.trainable = False
    model_adv.trainable = True
    for l in model_adv.layers:
        l.trainable = True
    
    #Construct the combined model (classifier -> regression), where only the regression part (adversarial) is trainable 
    model_cls_advOnly_train = Model(inp, out)
    model_cls_advOnly_train.compile(loss='mse', optimizer=opt)
    model_cls_advOnly_train.summary()
 
    return model_cls, model_adv, model_clsOnly_adv_train, model_cls_advOnly_train

# In[18]:


def show_losses( histories ):
    plt.figure(figsize=(10,10))
    #plt.ylim(bottom=0)
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training Error by Epoch')
    colors=[]
    do_acc=False
    for label,loss in histories:
        color = tuple(np.random.random(3))
        colors.append(color)
        l = label
        vl= label+" validation"
        if 'acc' in loss.history:
            l+=' (acc %2.4f)'% (loss.history['acc'][-1])
            do_acc = True
        if 'val_acc' in loss.history:
            vl+=' (val acc %2.4f)'% (loss.history['val_acc'][-1])
            do_acc = True
        plt.plot(loss.history['loss'], label=l, color=color)
        if 'val_loss' in loss.history:
            plt.plot(loss.history['val_loss'], lw=2, ls='dashed', label=vl, color=color)


    plt.legend()
    #plt.yscale('log')
    plt.savefig("/bigdata/shared/idutta/Higgs_plots/simpleCLS_loss.png")
  
    if not do_acc: return
    plt.figure(figsize=(10,10))
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    for i,(label,loss) in enumerate(histories):
        color = colors[i]
        if 'acc' in loss.history:
            plt.plot(loss.history['acc'], lw=2, label=label+" accuracy", color=color)
        if 'val_acc' in loss.history:
            plt.plot(loss.history['val_acc'], lw=2, ls='dashed', label=label+" validation accuracy", color=color)
    plt.legend(loc='lower right')
    plt.savefig("/bigdata/shared/idutta/Higgs_plots/simpleCLS_acc.png")


# In[19]:


def DNN_ROC(fpr_keras,tpr_keras,auc_keras,fpr_keras_tmp,tpr_keras_tmp,auc_keras_tmp) :
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='Keras Test (area = {:.3f})'.format(auc_keras))
    plt.plot(fpr_keras_tmp, tpr_keras_tmp, label='Keras Train (area = {:.3f})'.format(auc_keras_tmp))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.show()
    # Zoom in view of the upper left corner.
    plt.figure(2)
    plt.xlim(0, 0.2)
    plt.ylim(0.8, 1) 
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='Keras Test (area = {:.3f})'.format(auc_keras))
    plt.plot(fpr_keras_tmp, tpr_keras_tmp, label='Keras Train (area = {:.3f})'.format(auc_keras_tmp))
    plt.ylabel('True positive rate')
    plt.xlabel('False positive rate')
    plt.title('ROC curve (zoomed in at top left)')
    plt.legend(loc='best')
    plt.show()


# In[20]:


def Combined_ROC(fpr_keras,tpr_keras,auc_keras,fpr_bdt, tpr_bdt,area_bdt):
    plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='DNN (area = {:.3f})'.format(auc_keras),color='b')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    #plt.legend(loc='best')
    #plt.show()
    #plt.figure(figsize=(9,7))
    plt.plot(fpr_bdt,tpr_bdt,label="BDT (area ={:.3f})".format(area_bdt),color='r')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.show()
    #plt.xlabel('Background contamination')
    #plt.ylabel('Signal efficiency')






# In[24]:


from scipy import interp
from scipy.stats import wasserstein_distance
def cls_adv_score_roc(pred_simple_cls,pred_adv_cls,is_sig_test,is_bkg_test,y_test,z_test,x_test_wt,alpha_str):
    

    plt.figure(figsize=(2*10,2*10))
    plt.subplot(4,2,1)

    plt.hist(pred_simple_cls[is_sig_test], bins=np.linspace(0, 1, 100), density=1.0, lw=2.0, histtype="step", label="VBF events")
    plt.hist(pred_simple_cls[is_bkg_test], bins=np.linspace(0, 1, 100), density=1.0, lw=2.0, histtype="step", label="DY events");
    plt.title("Simple classifier")
    plt.xlabel("Classifier output")
    plt.legend(loc="best")
    plt.yscale("log")

    plt.subplot(4,2,2)
    plt.hist(pred_adv_cls[is_sig_test], bins=np.linspace(0, 1, 100), density=1.0, lw=2.0, histtype="step", label="VBF events")
    plt.hist(pred_adv_cls[is_bkg_test], bins=np.linspace(0, 1, 100), density=1.0, lw=2.0, histtype="step", label="DY events");
    plt.title("Adversarially trained classifier")
    plt.xlabel("Adversarially trained classifier output")
    plt.yscale("log")

    plt.subplot(4,2,3)
    bins = np.linspace(100., 150., 25)
    mass_shapes=[]
    for icut in np.linspace(0.0, 0.9, 10):
        cut = pred_simple_cls[:] > icut
        mass_shapes.append(z_test[np.invert(is_sig_test) & cut, 0])
        plt.hist(z_test[np.invert(is_sig_test) & cut, 0], bins=bins, density=1.0, lw=2.0, histtype="step", label="cls>{0:2f}".format(icut));
    plt.legend(fontsize=6, frameon=False)
    #plt.yscale("log")
    plt.xlabel("dilepton mass")
    '''
    print "Simple Classifier"
    # Kolmogorov Test, arrays should be in ascending order
    
    for i in range(0,10):
        mass_shapes[i].sort()
    for i in range(0,10):
        kval=rt.TMath.KolmogorovTest(len(mass_shapes[0]),mass_shapes[0],len(mass_shapes[i]),mass_shapes[i],"X")
        print "Kolmogorov Test value for cls > ",0.1*i, " is ", kval
        wd = wasserstein_distance(mass_shapes[0],mass_shapes[i])
        print "Wasserstein distance for cls > ",0.1*i, " is ", wd
    '''
    plt.subplot(4,2,4)
    mass_shapes_adv=[]
    for icut in np.linspace(0.0, 0.9, 10):
        cut = pred_adv_cls[:] > icut
        mass_shapes_adv.append(z_test[np.invert(is_sig_test) & cut, 0])
        plt.hist(z_test[np.invert(is_sig_test) & cut, 0], bins=bins, density=1.0, lw=2.0, histtype="step", label="cls>{0:2f}".format(icut));
    #plt.yscale("log")
    plt.legend(fontsize=6, frameon=False)
    plt.xlabel("dilepton mass")
    '''
    print "--------------------------------------------------------------------"
    print "Adversarial Training"
    # Kolmogorov Test, arrays should be in ascending order
    
    for i in range(0,10):
        mass_shapes_adv[i].sort()
    for i in range(0,10):
        kval=rt.TMath.KolmogorovTest(len(mass_shapes_adv[0]),mass_shapes_adv[0],len(mass_shapes_adv[i]),mass_shapes_adv[i],"X")
        print "Kolmogorov Test value for cls > ",0.1*i, " is ", kval
        wd = wasserstein_distance(mass_shapes_adv[0],mass_shapes_adv[i])
        print "Wasserstein distance for cls > ",0.1*i, " is ", wd
    '''
    min_mass=[110.,115.,116.5,118, 120.,122,123, 124.,124.5]
    max_mass=[150.,145.,140., 135.,130.,128.,127,126,125.5]   
    mean_fpr = np.linspace(0, 1, 100)
    plt.subplot(4,2,5)
    for i in range(0,len(min_mass)):
        y_pred=[]
        y_test_cls=[]
     
        for k in range(0,len(z_test)):
            if (z_test[k] > min_mass[i] and z_test[k] < max_mass[i]):
                y_pred.append(pred_simple_cls[k])
                y_test_cls.append(y_test[k])
        fpr_keras_tmp, tpr_keras_tmp, thresholds_keras_tmp = roc_curve(y_test_cls,y_pred)
        print thresholds_keras_tmp.shape
        tprs_tmp = interp(mean_fpr, fpr_keras_tmp, tpr_keras_tmp)
        auc_keras_tmp = auc(fpr_keras_tmp, tpr_keras_tmp)
        #plt.plot(mean_fpr, tprs_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]")
        plt.plot(fpr_keras_tmp, tpr_keras_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]") 
        
    plt.plot([0, 1], [0, 1], 'k--')    
    plt.legend(fontsize=6, frameon=False)
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    #plt.xlim(0.,0.01)
    #plt.ylim(0.8,1.0)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc='best')
    
    plt.subplot(4,2,6)
    
    for i in range(0,len(min_mass)):
        y_pred=[]
        y_test_cls=[]
        for k in range(0,len(z_test)):
            if (z_test[k] > min_mass[i] and z_test[k] < max_mass[i] ):
                y_pred.append(pred_adv_cls[k])
                y_test_cls.append(y_test[k])
        fpr_keras_tmp, tpr_keras_tmp, thresholds_keras_tmp = roc_curve(y_test_cls,y_pred)
        tprs_tmp = interp(mean_fpr, fpr_keras_tmp, tpr_keras_tmp)
        auc_keras_tmp = auc(fpr_keras_tmp, tpr_keras_tmp)
        #plt.plot(mean_fpr, tprs_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]")   
        plt.plot(fpr_keras_tmp, tpr_keras_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]") 
    plt.plot([0, 1], [0, 1], 'k--')    
    plt.legend(fontsize=6, frameon=False)
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    #plt.xlim(0.0,.01)
    #plt.ylim(0.8,1.)
    #plt.xscale('log')
    #plt.yscale('log')
    
    plt.legend(loc='best')

    plt.subplot(4,2,7)
    for i in range(0,len(min_mass)):
        y_pred=[]
        y_test_cls=[]
     
        for k in range(0,len(z_test)):
            if (z_test[k] > min_mass[i] and z_test[k] < max_mass[i]):
                y_pred.append(pred_simple_cls[k])
                y_test_cls.append(y_test[k])
        fpr_keras_tmp, tpr_keras_tmp, thresholds_keras_tmp = roc_curve(y_test_cls,y_pred)
        print thresholds_keras_tmp.shape
        tprs_tmp = interp(mean_fpr, fpr_keras_tmp, tpr_keras_tmp)
        auc_keras_tmp = auc(fpr_keras_tmp, tpr_keras_tmp)
        #plt.plot(mean_fpr, tprs_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]")
        plt.plot(fpr_keras_tmp, tpr_keras_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]") 
        
    plt.plot([0, 1], [0, 1], 'k--')    
    plt.legend(fontsize=6, frameon=False)
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.xlim(0.,0.01)
    plt.ylim(0.8,1.0)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc='best')
    
    plt.subplot(4,2,8)
    
    for i in range(0,len(min_mass)):
        y_pred=[]
        y_test_cls=[]
        for k in range(0,len(z_test)):
            if (z_test[k] > min_mass[i] and z_test[k] < max_mass[i] ):
                y_pred.append(pred_adv_cls[k])
                y_test_cls.append(y_test[k])
        fpr_keras_tmp, tpr_keras_tmp, thresholds_keras_tmp = roc_curve(y_test_cls,y_pred)
        tprs_tmp = interp(mean_fpr, fpr_keras_tmp, tpr_keras_tmp)
        auc_keras_tmp = auc(fpr_keras_tmp, tpr_keras_tmp)
        #plt.plot(mean_fpr, tprs_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]")   
        plt.plot(fpr_keras_tmp, tpr_keras_tmp, label="$M_{\mu\mu} \ \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]") 
    plt.plot([0, 1], [0, 1], 'k--')    
    plt.legend(fontsize=6, frameon=False)
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.xlim(0.0,.2)
    plt.ylim(0.1,1.)
    #plt.xscale('log')
    #plt.yscale('log')
    
    plt.legend(loc='best')
    plt.savefig('/bigdata/shared/idutta/Higgs_plots/Model_'+alpha_str+'_performance.png')
    '''
  #-----------------------------------------------------------------------------------------------------# 
    print 'Adversarially trained Classifier'
    print '------------------------------------------------------------------------------------------- '
    histo=[]
    cols=[1,2,3,4,5,6,7,8,9,42]
    for icut in np.linspace(0.0, 0.9, 10):
        cut = pred_adv_cls[:] > icut
        name_bkg='Dilep_mass_bkg_'+str(icut)
        h= rt.TH1F(name_bkg,"",20,110.,150.)
        data=z_test[np.invert(is_sig_test) & cut, 0]
        wt=x_test_wt[np.invert(is_sig_test) & cut, 0]
        h.FillN(len(data),data,wt)
        h.GetXaxis().SetTitle("M(#mu#mu)")
        h.GetYaxis().SetTitle("Entries") 
        h.GetYaxis().SetRangeUser(0.001,5000.)
        h.SetLineColor(cols[int(icut*10)])
        h.SetMarkerColor(cols[int(icut*10)])
        h.Scale(1/h.Integral())
        histo.append(h)
        
        name_sig='Dilep_mass_sig_'+str(icut)
        h_sig= rt.TH1F(name_sig,"",20,110.,150.)
        data_sig=z_test[is_sig_test & cut, 0]
        wt_sig=x_test_wt[is_sig_test & cut, 0]
        h_sig.FillN(len(data_sig),data_sig,wt_sig)
        sens =0
        #for i in range(6,11): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        #print "Sensitivity in [120.,130.] GeV of cls > ", icut, " is : ", math.sqrt(sens)  
        
        sens =0
        #for i in range(7,10): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        #print "Sensitivity in [122.,128.] GeV of cls > ", icut, " is : ", math.sqrt(sens) 
        
        sens =0
        #for i in range(7,9): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        #print "Sensitivity in [124.,126.] GeV of cls > ", icut, " is : ", math.sqrt(sens)
        #plt.hist(z_test[np.invert(is_sig_test) & cut, 0], bins=bins, density=1.0, lw=2.0, histtype="step", label="cls>{0:2f}".format(icut));
    #plt.yscale("log")
    c =rt.TCanvas("c","c",900,600)
    legend = rt.TLegend(0.9, 0.5, 0.99, 0.9)
    legend.SetHeader("Adversarial Training")
    for icut in np.linspace(0.0, 0.9, 10):
        histo[int(icut*10)].Draw("e hist same")
        name="cls > "+str(icut)
        legend.AddEntry(histo[int(icut*10)],name)
    
    legend.Draw("same")  
    rt.gPad.SetLogy()
    c.Draw()
    print 'Simple Classifier'
    print '------------------------------------------------------------------------------------------- '
    histo_simple=[]
    for icut in np.linspace(0.0, 0.9, 10):
        cut = pred_simple_cls[:] > icut
        name_bkg='Dilep_mass_bkg_'+str(icut)
        h= rt.TH1F(name_bkg,"",20,110.,150.)
        data=z_test[np.invert(is_sig_test) & cut, 0]
        wt=x_test_wt[np.invert(is_sig_test) & cut, 0]
        h.FillN(len(data),data,wt)
        h.GetXaxis().SetTitle("M(#mu#mu)")
        h.GetYaxis().SetTitle("Entries") 
        h.GetYaxis().SetRangeUser(0.001,5000.) 
        h.SetLineColor(cols[int(icut*10)])
        h.SetMarkerColor(cols[int(icut*10)])
        h.Scale(1/h.Integral())
        histo_simple.append(h)
        
        
        name_sig='Dilep_mass_sig_'+str(icut)
        h_sig= rt.TH1F(name_sig,"",20,110.,150.)
        data_sig=z_test[is_sig_test & cut, 0]
        wt_sig=x_test_wt[is_sig_test & cut, 0]
        h_sig.FillN(len(data_sig),data_sig,wt_sig)
        sens =0
        #for i in range(6,11): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        print "Sensitivity in [120.,130.] GeV of cls > ", icut, " is : ", math.sqrt(sens)  
        
        sens =0
        #for i in range(7,10): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        print "Sensitivity in [122.,128.] GeV of cls > ", icut, " is : ", math.sqrt(sens) 
        
        sens =0
        #for i in range(7,9): sens+=((h_sig.GetBinContent(i)**2)/(h.GetBinContent(i)*sig_scale**2)) 
        print "Sensitivity in [124.,126.] GeV of cls > ", icut, " is : ", math.sqrt(sens) 
        #plt.hist(z_test[np.invert(is_sig_test) & cut, 0], bins=bins, density=1.0, lw=2.0, histtype="step", label="cls>{0:2f}".format(icut));
    #plt.yscale("log")
    c1 =rt.TCanvas("c1","c1",900,600)
    c1.cd()
    legend1 = rt.TLegend(0.7, 0.5, 0.99, 0.9)
    legend1.SetHeader("Simple Classifier")
    for icut in np.linspace(0.0, 0.9, 10):
        histo_simple[int(icut*10)].Draw("e hist same")
        name="cls > "+str(icut)
        legend1.AddEntry(histo_simple[int(icut*10)],name)
    
    legend1.Draw("same")  
    rt.gPad.SetLogy()
    c1.Draw()
   ''' 


# In[25]:


def compare_rocs(pred_simple_cls,pred_adv_cls,y_test,z_test,alpha_str):
    plt.figure(figsize=(2*11,2*11))
    min_mass=[110.,115.,116.5,118, 120.,122,123, 124.,124.5]
    max_mass= [150.,145.,140., 135.,130.,128.,127,126,125.5]
    mean_fpr = np.linspace(0, 1, 1000)
    for i in range(0,len(min_mass)):
        plt.subplot(3,3,i+1)
        y_pred=[]
        y_pred_adv=[]
        y_test_cls=[]
        for k in range(0,len(z_test)):
            if (z_test[k] > min_mass[i] and z_test[k] < max_mass[i] ):
                y_pred.append(pred_simple_cls[k])
                y_pred_adv.append(pred_adv_cls[k])
                y_test_cls.append(y_test[k])
                
        fpr_keras_tmp, tpr_keras_tmp, thresholds_keras_tmp = roc_curve(y_test_cls,y_pred)
        
        
        fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test_cls,y_pred_adv)
        #print y_test_cls,y_pred,y_pred_adv
        plt.plot(y_pred,y_pred_adv)
        plt.plot(fpr_keras_tmp, tpr_keras_tmp, label="Simple DNN : $M_{\mu\mu} \  \in \ [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]")  
        plt.plot(fpr_keras, tpr_keras, label="Adversarially trained Classifier: $M_{\mu\mu} \  \in \  [$ "+str(min_mass[i])+" , "+str(max_mass[i])+" ]") 
        plt.plot([0, 1], [0, 1], 'k--')    
        plt.legend(fontsize=6, frameon=False)
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.title('ROC curve')
        #plt.xlim(0.0,.1)
        #plt.ylim(0.8,1.)
        #plt.xscale('log')
        #plt.yscale('log')
    
        plt.legend(loc='best')
    plt.savefig('/bigdata/shared/idutta/Higgs_plots/compareRocs_'+alpha_str+'.png')

        
def draw_losses(losses_train_cls,losses_test_cls,alpha_str):
    
    plt.subplot(3,1,1)
    plt.plot([l[2] for l in losses_train_cls], color="red", ls="--", label="Adversarial train loss")
    plt.plot([l[2] for l in losses_test_cls], color="red", label="Adversarial test loss")


    plt.xlabel("epoch")
    #plt.ylim(-1e-5,300.)
    #plt.xlim(10.,100.)
    #plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.tight_layout()

    plt.subplot(3,1,2)
    plt.plot([l[1] for l in losses_train_cls], color="blue", ls="--", label="Classification train loss ")
    plt.plot([l[1]for l in losses_test_cls], color="blue", label="Classification test loss")


    plt.xlabel("epoch")
    #plt.ylim(-1e-5,300.)
    #plt.xlim(10.,100.)
    #plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.tight_layout()
    
    plt.subplot(3,1,3)
    plt.plot([l[0] for l in losses_train_cls], color="orange", ls="--",label="Total train loss")
    plt.plot([l[0] for l in losses_test_cls], color="orange", label="Total test loss")

    plt.xlabel("epoch")
    #plt.ylim(-1e-5,300.)
    #plt.xlim(10.,100.)

    plt.xscale('log')
    plt.legend()
    plt.tight_layout()

    plt.savefig('/bigdata/shared/idutta/Higgs_plots/Losses_'+alpha_str+'.png')
# In[27]:


def trainer(m_clsOnly_adv_train,m_cls_advOnly_train,x_train_reduced,x_test_reduced,y_train, z_train,y_test, z_test,x_train_wt,x_test_wt,alpha, alpha_str):
    print x_train_wt[:,0].shape, x_test_wt[:,0].shape
    batch_size = 10000
    losses_train_cls = []
    losses_test_cls = []
    best_cls_loss=1000.
    best_adv_loss=1000.
    best_cls_epoch = 0
    best_adv_epoch =0
    last_epoch=0
    for epoch in range(300):
        iAdv=1
        iCls=1
        l0 = m_clsOnly_adv_train.evaluate(x_train_reduced, [y_train, z_train[:, 0]], batch_size=batch_size, verbose=False, sample_weight=[x_train_wt[:,0], x_train_wt[:,0]])
        l1 = m_clsOnly_adv_train.evaluate(x_test_reduced, [y_test, z_test[:, 0]], batch_size=batch_size, verbose=False, sample_weight=[x_test_wt[:,0], x_test_wt[:,0]])
        l2 = m_cls_advOnly_train.evaluate(x_train_reduced, z_train[:, 0], batch_size=batch_size, verbose=False, sample_weight=x_train_wt[:,0])

        #train_names = ['train_loss', 'train_cls_loss', 'train_adv_loss']
        #val_names = ['val_loss', 'val_cls_loss', 'val_adv_loss']
        #write_log(tensorboard, train_names, l0, epoch)
        #write_log(tensorboard, val_names, l1, epoch)

        losses_train_cls += [l0]
        losses_test_cls += [l1]
        #totLoss=l0[0]
        #if epoch==0:
            #prev_totLoss=totLoss
        #if (totLoss-prev_totLoss)<0:
            #last_epoch =epoch
        #else:
            #print "turning of adversarial at epoch: ", last_epoch
            #iAdv=0
        if epoch%10==0:
            print epoch, l0, l1, l2
        #if l2>60. : last_epoch =epoch
        #if l2<60. and (epoch-last_epoch)<6:
            #iAdv=0
        #prev_totLoss=totLoss 
        for k in range(iAdv):
            for ibatch in range(0, len(x_train_reduced), batch_size):
                xb = x_train_reduced[ibatch:ibatch+batch_size]
                yb = y_train[ibatch:ibatch+batch_size]
                zb = z_train[ibatch:ibatch+batch_size]
                wb = x_train_wt[ibatch:ibatch+batch_size]
                l = m_cls_advOnly_train.train_on_batch(xb, zb[:, 0], sample_weight=wb[:,0])

        for k in range(iCls):
            for ibatch in range(0, len(x_train_reduced), batch_size):
                xb = x_train_reduced[ibatch:ibatch+batch_size]
                yb = y_train[ibatch:ibatch+batch_size]
                zb = z_train[ibatch:ibatch+batch_size]
                wb = x_train_wt[ibatch:ibatch+batch_size]
                l = m_clsOnly_adv_train.train_on_batch(xb, [yb, zb[:, 0]], sample_weight=[wb[:,0], wb[:,0]])
    
    draw_losses(losses_train_cls, losses_test_cls,alpha_str)
    return losses_train_cls, losses_test_cls


def output_file_nonTrain(y_pred_cls,y_pred_adv_cls,evt_weight, x_test_index,INPUT_FILE,proc_tree,alpha_str):
    # Create a lookup table for discriminator value by event number
    disc_lookup_signal = {}
    disc_lookup_signal_adv = {}
    proc_evt_weight={}
    
    for disc_val, disc_val_adv, evt_wt, run, lumi, event in zip(y_pred_cls, y_pred_adv_cls, evt_weight, x_test_index[:,0],x_test_index[:,1], x_test_index[:,2]):
   
        disc_lookup_signal[(run, lumi, event)] = disc_val
        disc_lookup_signal_adv[(run, lumi, event)] = disc_val_adv
        proc_evt_weight[(run, lumi, event)] = evt_wt
      
          
    # We write out the signal and background events to a new ROOT file.
    # For events in the test set, we append the GBM discriminator.  
    # For events in the train set, we set the discriminator to -1
    # to avoid accidentally reusing the training set in the analysis.
    print len(disc_lookup_signal), len(disc_lookup_signal_adv)
    out_name = os.path.basename(INPUT_FILE).replace('.root', alpha_str+'_NNscore.root')
    out_proc = rt.TFile(out_name, 'RECREATE')
    out_proc_tree = proc_tree.CloneTree(0)


    # This is the boilerplate code for writing something
    # to a ROOT tree from Python.
    
    rt.gROOT.ProcessLine("struct MyStruct2{float disc_advNN;};")
    rt.gROOT.ProcessLine("struct MyStruct3{float disc_simpleNN;};")
    rt.gROOT.ProcessLine("struct MyStruct{double evt_wt;};")
    from ROOT import MyStruct3, MyStruct2, MyStruct
    s = MyStruct3()
    s_adv = MyStruct2()
    s_wt = MyStruct()
    disc_branch_proc = out_proc_tree.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
    disc_branch_sig_proc = out_proc_tree.Branch('disc_advNN', rt.AddressOf(s_adv, 'disc_advNN'), 'disc_advNN/F');
    evt_weight_branch_proc = out_proc_tree.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
    
    print "Writing new ROOT process file with discriminator appended"
    fill_discriminator(proc_tree, out_proc_tree, disc_lookup_signal,disc_lookup_signal_adv, proc_evt_weight,s, s_adv,s_wt)
    
    # Cristian's code uses GetCurrentFile() for this part.
    # I will do that too just in case (paranoia).
    out_proc.cd()
    out_proc_tree.GetCurrentFile().Write()
    #signalNevents.Write()
    out_proc_tree.GetCurrentFile().Close()

    print 'done writing output'
