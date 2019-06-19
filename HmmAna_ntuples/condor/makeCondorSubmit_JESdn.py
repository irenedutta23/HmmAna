import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
datasets = [#["DYJetsToLL_VBFfilter_2017", "DYJetsToLL_VBFfilter_2017"]
            #["DYJetsToLL_VBFfilter_2018","DYJetsToLL_VBFfilter_2018"]
            #["DYJetsToLL_VBFfilter_2016","DYJetsToLL_VBFfilter_2016"]
            #["DYJetsToLL","DYJetsToLL"]
    #["DYJetsToLL_M105To160_incl_2017", "DYJetsToLL_M105To160_incl_2017"]
    #["DYJetsToLL_M105To160_incl_2018", "DYJetsToLL_M105To160_incl_2018"]
    #["DYJetsToLL_M105To160_incl_2016", "DYJetsToLL_M105To160_incl_2016"]
    #["DYJetsToLL_2017", "DYJetsToLL_2017"]                                                                                                                                     
    #["DYJetsToLL_2018","DYJetsToLL_2018"]                                                                                                                                     
    #["DYJetsToLL_2016","DYJetsToLL_2016"] 
    #["ggH_2017","ggH_2017"]
            #["ggH_2018","ggH_2018"]
    #["ggH_2016","ggH_2016"]
    #["VBFHToMuMu_2017","VBFHToMuMu_2017"]
    #["VBFHToMuMu_2016","VBFHToMuMu_2016"]
    #["VBFHToMuMu_2018","VBFHToMuMu_2018"]
    #["ZZ","ZZ"]
            #["ttsl_2018","ttsl_2018"]
            #["ZH_2018","ZH_2018"]
            #["ttH_2018","ttH_2018"]
            #["WminusH_2018","WminusH_2018"]
            #["WplusH_2018","WplusH_2018"]
    #["ZH_2016","ZH_2016"]                                                                                                                                                                          
            #["ttH_2016","ttH_2016"]                                                                                                                                                                       
            #["WminusH_2016","WminusH_2016"]                                                                                                                                                               
            #["WplusH_2016","WplusH_2016"]         
    #["ttTo2l2v_2018","ttTo2l2v_2018"]
    #["ttsl_2017","ttsl_2017"]   
    #["ttTo2l2v_2016","ttTo2l2v_2016"]                                                                                                           
    
    #["ttsl_2016","ttsl_2016"] 
    #["ZH_2017","ZH_2017"]                                                                                               
            #["ttH_2017","ttH_2017"]                                                                                                  
            #["WminusH_2017","WminusH_2017"]
            #["WplusH_2017","WplusH_2017"]                                                      
    #["ttTo2l2v_2017","ttTo2l2v_2017"]
    #["ttJets_DiLept_2017","ttJets_DiLept_2017"]
    #["ttJets_DiLept_2016","ttJets_DiLept_2016"]
    #["ttJets_DiLept_2018","ttJets_DiLept_2018"]
    #["EWK_2016","EWK_2016"]
    #["EWK_2017","EWK_2017"]
    #["EWK_2018","EWK_2018"]
            #["WZTo3LNu_2016","WZTo3LNu_2016"]
            #["WZTo2L2Q_2016","WZTo2L2Q_2016"]
            #["WWTo2L2Nu_2016","WWTo2L2Nu_2016"]
            #["WWToLNuQQ_2016","WWToLNuQQ_2016"]
            #["WWW_4F_2016","WWW_4F_2016"]
    #["WWZ_4F_2016","WWZ_4F_2016"] 
           # ["WZZ_2016","WZZ_2016"]
    #["ZZZ_2016","ZZZ_2016"]
    #["WZTo3LNu_2018","WZTo3LNu_2018"]
            #["WZTo2L2Q_2018","WZTo2L2Q_2018"]
            #["WWTo2L2Nu_2018","WWTo2L2Nu_2018"]
            #["WWToLNuQQ_2018","WWToLNuQQ_2018"]
            #["WWW_4F_2018","WWW_4F_2018"]
    #["WWZ_4F_2018","WWZ_4F_2018"] 
    #["WZZ_2018","WZZ_2018"]
    #["ZZZ_2018","ZZZ_2018"]
    #["WZTo3LNu_2017","WZTo3LNu_2017"]
            #["WZTo2L2Q_2017","WZTo2L2Q_2017"]
            #["WWTo2L2Nu_2017","WWTo2L2Nu_2017"]
            #["WWToLNuQQ_2017","WWToLNuQQ_2017"]
    #["WWW_4F_2017","WWW_4F_2017"]
    #["WWZ_4F_2017","WWZ_4F_2017"] 
       #["WZZ_2017","WZZ_2017"]
    #["ZZZ_2017","ZZZ_2017"]
    #["TTZToLLNuNu_2017","TTZToLLNuNu_2017"]
    #["TTWJetsToLNu_2017","TTWJetsToLNu_2017"]
    #["TTZToLLNuNu_2018","TTZToLLNuNu_2018"]
    #["TTWJetsToLNu_2018","TTWJetsToLNu_2018"]
    #["TTZToLLNuNu_2016","TTZToLLNuNu_2016"]
    #["TTWJetsToLNu_2016","TTWJetsToLNu_2016"]
    #["ZZTo4L_2016","ZZTo4L_2016"]
    #["ZZTo2L2Q_2016","ZZTo2L2Q_2016"]
    ["ZZTo2L2Nu_2016","ZZTo2L2Nu_2016"]
    #["ZZTo4L_2017","ZZTo4L_2017"]
    #["ZZTo2L2Q_2017","ZZTo2L2Q_2017"]
    #["ZZTo2L2Nu_2017","ZZTo2L2Nu_2017"]
    #["ZZTo4L_2018","ZZTo4L_2018"]
    #["ZZTo2L2Q_2018","ZZTo2L2Q_2018"]
    #["ZZTo2L2Nu_2018","ZZTo2L2Nu_2018"]
]

NSections = 10
readFiles = ""

for data in datasets:
    jobidx = 0
    if ( data[0]=="DYJetsToLL_VBFfilter_2017"):
        dataname  = "DYJetsToLL_VBFfilter_2017"
        inputfname = "DYJetsToLL_VBFfilter_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    
    elif ( data[0]=="DYJetsToLL_VBFfilter_2018"):
        dataname = "DYJetsToLL_VBFfilter_2018"
        inputfname = "DYJetsToLL_VBFfilter_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1    

    elif ( data[0]=="DYJetsToLL_VBFfilter_2016"):
        dataname = "DYJetsToLL_VBFfilter_2016"
        inputfname = "DYJetsToLL_VBFfilter_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL_M105To160_incl_2017"):
        dataname = "DYJetsToLL_M105To160_incl_2017"
        inputfname = "DYJetsToLL_M105To160_incl_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL_M105To160_incl_2016"):
        dataname = "DYJetsToLL_M105To160_incl_2016"
        inputfname = "DYJetsToLL_M105To160_incl_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL_M105To160_incl_2018"):
        dataname = "DYJetsToLL_M105To160_incl_2018"
        inputfname = "DYJetsToLL_M105To160_incl_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL_2017"):
        dataname  = "DYJetsToLL_2017"
        inputfname = "DYJetsToLL_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    
    elif ( data[0]=="DYJetsToLL_2018"):
        dataname = "DYJetsToLL_2018"
        inputfname = "DYJetsToLL_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1    

    elif ( data[0]=="DYJetsToLL_2016"):
        dataname = "DYJetsToLL_2016"
        inputfname = "DYJetsToLL_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL"):
        dataname = "DYJetsToLL"
        inputfname = "DYJetsToLL_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ggH_2018"):
        dataname = "ggH_2018"
        inputfname = "ggH_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ggH_2017"):
        dataname = "ggH_2017"
        inputfname = "ggH_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ggH_2016"):
        dataname = "ggH_2016"
        inputfname = "ggH_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="VBFHToMuMu_2017"):
        dataname = "VBFHToMuMu_2017"
        inputfname = "VBFHToMuMu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1
    
    elif ( data[0]=="VBFHToMuMu_2016"):
        dataname = "VBFHToMuMu_2016"
        inputfname = "VBFHToMuMu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1
    elif ( data[0]=="VBFHToMuMu_2018"):
        dataname = "VBFHToMuMu_2018"
        inputfname = "VBFHToMuMu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="WZ"):
        dataname = "WZ"
        inputfname = "WZ_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ZZTo4L_2016"):
        dataname = "ZZTo4L_2016"
        inputfname = "ZZTo4L_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Q_2016"):
        dataname = "ZZTo2L2Q_2016"
        inputfname = "ZZTo2L2Q_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Nu_2016"):
        dataname = "ZZTo2L2Nu_2016"
        inputfname = "ZZTo2L2Nu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo4L_2018"):
        dataname = "ZZTo4L_2018"
        inputfname = "ZZTo4L_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Q_2018"):
        dataname = "ZZTo2L2Q_2018"
        inputfname = "ZZTo2L2Q_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Nu_2018"):
        dataname = "ZZTo2L2Nu_2018"
        inputfname = "ZZTo2L2Nu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    
    elif ( data[0]=="ZZTo4L_2017"):
        dataname = "ZZTo4L_2017"
        inputfname = "ZZTo4L_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Q_2017"):
        dataname = "ZZTo2L2Q_2017"
        inputfname = "ZZTo2L2Q_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZTo2L2Nu_2017"):
        dataname = "ZZTo2L2Nu_2017"
        inputfname = "ZZTo2L2Nu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttsl_2018"):
        dataname = "ttsl_2018"
        inputfname = "ttsl_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ZH_2018"):
        dataname = "ZH_2018"
        inputfname = "ZH_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ttH_2018"):
        dataname = "ttH_2018"
        inputfname = "ttH_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WminusH_2018"):
        dataname = "WminusH_2018"
        inputfname = "WminusH_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WplusH_2018"):
        dataname = "WplusH_2018"
        inputfname = "WplusH_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZH_2016"):
        dataname = "ZH_2016"
        inputfname = "ZH_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ttH_2016"):
        dataname = "ttH_2016"
        inputfname = "ttH_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WminusH_2016"):
        dataname = "WminusH_2016"
        inputfname = "WminusH_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WplusH_2016"):
        dataname = "WplusH_2016"
        inputfname = "WplusH_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttTo2l2v_2018"):
        dataname = "ttTo2l2v_2018"
        inputfname = "ttTo2l2v_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttsl_2017"):
        dataname = "ttsl_2017"
        inputfname = "ttsl_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttTo2l2v_2016"):
        dataname = "ttTo2l2v_2016"
        inputfname = "ttTo2l2v_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttsl_2016"):
        dataname = "ttsl_2016"
        inputfname = "ttsl_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttJets_DiLept_2016"):
        dataname = "ttJets_DiLept_2016"
        inputfname = "ttJets_DiLept_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ttJets_DiLept_2017"):
        dataname = "ttJets_DiLept_2017"
        inputfname = "ttJets_DiLept_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ttJets_DiLept_2018"):
        dataname = "ttJets_DiLept_2018"
        inputfname = "ttJets_DiLept_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZH_2017"):
        dataname = "ZH_2017"
        inputfname = "ZH_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ttH_2017"):
        dataname = "ttH_2017"
        inputfname = "ttH_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WminusH_2017"):
        dataname = "WminusH_2017"
        inputfname = "WminusH_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WplusH_2017"):
        dataname = "WplusH_2017"
        inputfname = "WplusH_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ttTo2l2v_2017"):
        dataname = "ttTo2l2v_2017"
        inputfname = "ttTo2l2v_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="EWK_2017"):
        dataname = "EWK_2017"
        inputfname = "EWK_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="EWK_2018"):
        dataname = "EWK_2018"
        inputfname = "EWK_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="EWK_2016"):
        dataname = "EWK_2016"
        inputfname = "EWK_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WZTo3LNu_2016"):
        dataname = "WZTo3LNu_2016"
        inputfname = "WZTo3LNu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="WZTo2L2Q_2016"):
        dataname = "WZTo2L2Q_2016"
        inputfname = "WZTo2L2Q_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="WWTo2L2Nu_2016"):
        dataname = "WWTo2L2Nu_2016"
        inputfname = "WWTo2L2Nu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWToLNuQQ_2016"):
        dataname = "WWToLNuQQ_2016"
        inputfname = "WWToLNuQQ_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWW_4F_2016"):
        dataname = "WWW_4F_2016"
        inputfname = "WWW_4F_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWZ_4F_2016"):
        dataname = "WWZ_4F_2016"
        inputfname = "WWZ_4F_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="WZZ_2016"):
        dataname = "WZZ_2016"
        inputfname = "WZZ_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZZ_2016"):
        dataname = "ZZZ_2016"
        inputfname = "ZZZ_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WZTo3LNu_2017"):
        dataname = "WZTo3LNu_2017"
        inputfname = "WZTo3LNu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="WZTo2L2Q_2017"):
        dataname = "WZTo2L2Q_2017"
        inputfname = "WZTo2L2Q_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="WWTo2L2Nu_2017"):
        dataname = "WWTo2L2Nu_2017"
        inputfname = "WWTo2L2Nu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWToLNuQQ_2017"):
        dataname = "WWToLNuQQ_2017"
        inputfname = "WWToLNuQQ_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWW_4F_2017"):
        dataname = "WWW_4F_2017"
        inputfname = "WWW_4F_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWZ_4F_2017"):
        dataname = "WWZ_4F_2017"
        inputfname = "WWZ_4F_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="WWZ_4F_2018"):
        dataname = "WWZ_4F_2018"
        inputfname = "WWZ_4F_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
            NSections = 10

    elif ( data[0]=="WZZ_2017"):
        dataname = "WZZ_2017"
        inputfname = "WZZ_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZZ_2017"):
        dataname = "ZZZ_2017"
        inputfname = "ZZZ_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WZTo3LNu_2018"):
        dataname = "WZTo3LNu_2018"
        inputfname = "WZTo3LNu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="WZTo2L2Q_2018"):
        dataname = "WZTo2L2Q_2018"
        inputfname = "WZTo2L2Q_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="WWTo2L2Nu_2018"):
        dataname = "WWTo2L2Nu_2018"
        inputfname = "WWTo2L2Nu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWToLNuQQ_2018"):
        dataname = "WWToLNuQQ_2018"
        inputfname = "WWToLNuQQ_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWW_4F_2018"):
        dataname = "WWW_4F_2018"
        inputfname = "WWW_4F_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WZZ_2018"):
        dataname = "WZZ_2018"
        inputfname = "WZZ_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ZZZ_2018"):
        dataname = "ZZZ_2018"
        inputfname = "ZZZ_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ttJets_DiLept"):
        dataname = "ttJets_DiLept"
        inputfname = "ttJets_DiLept_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="TTZToLLNuNu_2017"):
        dataname = "TTZToLLNuNu_2017"
        inputfname = "TTZToLLNuNu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="TTWJetsToLNu_2017"):
        dataname = "TTWJetsToLNu_2017"
        inputfname = "TTWJetsToLNu_2017_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="TTZToLLNuNu_2018"):
        dataname = "TTZToLLNuNu_2018"
        inputfname = "TTZToLLNuNu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="TTWJetsToLNu_2018"):
        dataname = "TTWJetsToLNu_2018"
        inputfname = "TTWJetsToLNu_2018_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
        
    elif ( data[0]=="TTZToLLNuNu_2016"):
        dataname = "TTZToLLNuNu_2016"
        inputfname = "TTZToLLNuNu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="TTWJetsToLNu_2016"):
        dataname = "TTWJetsToLNu_2016"
        inputfname = "TTWJetsToLNu_2016_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
  
    elif ( data[0]=="MuData1"):
        dataname = "MuData1"
        inputfname = "Data1_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="MuData2"):
        dataname = "MuData2"
        inputfname = "Data2_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="MuData3"):
        dataname = "MuData3"
        inputfname = "Data3_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData4"):
        dataname = "MuData4"
        inputfname = "Data4_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData5"):
        dataname = "MuData5"
        inputfname = "Data5_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData6"):
        dataname = "MuData6"
        inputfname = "Data6_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData7"):
        dataname = "MuData7"
        inputfname = "Data7_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="MuData8"):
        dataname = "MuData8"
        inputfname = "Data8_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="MuData9"):
        dataname = "MuData9"
        inputfname = "Data9_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="MuData10"):
        dataname = "MuData10"
        inputfname = "Data10_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData11"):
        dataname = "MuData11"
        inputfname = "Data11_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData12"):
        dataname = "MuData12"
        inputfname = "Data12_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData13"):
        dataname = "MuData13"
        inputfname = "Data13_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData14"):
        dataname = "MuData14"
        inputfname = "Data14_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="MuData15"):
        dataname = "MuData15"
        inputfname = "Data15_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData16"):
        dataname = "MuData16"
        inputfname = "Data16_Muon_JESdn.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal

    print "Dataset ",  data[0], " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    outDir="/mnt/hadoop/store/user/idutta/HmmAna/HistogramRootFiles/"+data[0]+"_JESdn_18June2019"
    print outDir
    #if not os.path.exists(outDir):
        #mkdirComm="scram unsetenv -sh;gfal-mkdir " +outDir
        #os.system(mkdirComm)
    outdir ="HmmAna/HistogramRootFiles/"+data[0]+"_JESdn_18June2019"
    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/data/idutta/CMSSW_9_4_9/src/HmmAna/master/HmmAna_ntuples/condor/condor_submit/runList_"+data[0]+"_"+str(jobidx)+"_JESdn.txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+data[0]+"_JESdn"+str(jobidx)
        jobName      = "18June2019"+data[0]+"_JESdn_job"+str(jobidx)
        outHistFile = data[0]+"_JESdn_job"+str(jobidx)+".root"
        isData       ="F"
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            line=line.replace("ROOTOUT",outHistFile)
            line=line.replace("DATANAME",dataname)
            line=line.replace("ISDATA",isData)
            line=line.replace("OUTDIR",outdir)
            print line.rstrip()
                        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
