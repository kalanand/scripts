#!/usr/bin/env python


import os,sys
import string, re
from time import gmtime, localtime, strftime





def PrintToTextFile(inputRootFile, txtfilename, mode):
    from ROOT import TTree, TFile, gROOT
    import re, array
    from array import array
    gROOT.Reset()
    f = TFile(inputRootFile,"read")
    tree = f.Get("WJet")

    entries = tree.GetEntriesFast()
    txtfile = file(txtfilename,"w")
    
    # loop over events
    for jentry in xrange( entries ):
        # get the next tree in the chain
        ientry = tree.LoadTree(jentry)
        if ientry < 0:
            break

        # copy next entry into memory and verify
        nb = tree.GetEntry(jentry)
        if nb<=0:
            continue       
        if not (tree.gdevtt==1):
            continue

        if jentry%50000 == 0:
            print "Processing entry = "+str(jentry)

        if mode==0:
            line = '%i' %(tree.event_runNo) + " "+ '%i' %(tree.event_lumi) + " "+ '%i' %(tree.event_evtNo) + " "+ '%10.5f' %(tree.W_muon_pt)+" "+ '%10.5f' %(tree.W_muon_eta) + ' ' + '%10.5f' %(tree.W_muon_phi) + ' ' + '%10.5f' %(tree.event_met_pfmet)  +  ' ' + '%10.5f' %(tree.event_met_pfmetPhi) + ' ' + '%10.5f' %(tree.JetPFCor_Pt[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Eta[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Phi[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Pt[1]) + ' ' + '%10.5f' %(tree.JetPFCor_Eta[1]) + ' ' + '%10.5f' %(tree.JetPFCor_Phi[1]) + '\n'
        else:
            line = '%i' %(tree.event_runNo) + " "+ '%i' %(tree.event_lumi) + " "+ '%i' %(tree.event_evtNo) + " "+ '%10.5f' %(tree.W_electron_pt)+" "+ '%10.5f' %(tree.W_electron_eta) + ' ' + '%10.5f' %(tree.W_electron_phi) + ' ' + '%10.5f' %(tree.event_met_pfmet)  +  ' ' + '%10.5f' %(tree.event_met_pfmetPhi) + ' ' + '%10.5f' %(tree.JetPFCor_Pt[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Eta[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Phi[0]) + ' ' + '%10.5f' %(tree.JetPFCor_Pt[1]) + ' ' + '%10.5f' %(tree.JetPFCor_Eta[1]) + ' ' + '%10.5f' %(tree.JetPFCor_Phi[1]) + '\n'
          
        txtfile.write(line)

    print "flat file "+txtfilename +" has been written."
    f.Close()


############################################################
PrintToTextFile("RD_WmunuJets_DataAll_GoldenJSON_1p1invfb.root", "EventsList_muon.txt", 0);
PrintToTextFile("RD_WenuJets_DataAll_GoldenJSON_1p1invfb.root", "EventsList_electron.txt", 1);
############################################################

