#!/usr/bin/env python

import ROOT
import os
import commands


def col2hist(infile, column=0, name="col", 
             nbin=2000, lobin = -1000., hibin = 1000.):

    file = open(infile,"r") 
    max = -1.e10
    min =  1.e10
    vals = {}

    if name == "col": name = name+str(column)
    
    nline = 0
    line = file.readline()   
    while line:
        sline = line.split()
        value = float(sline[int(column)])
        if value > max : max = value
        if value < min : min = value
        vals[nline] = value
        nline = nline+1
        line = file.readline()

    hname = "h"+name
#    cname = "canv"
    hist = ROOT.TH1F(hname,hname,1000, min-0.05*min, max+0.05*max)
    for value in vals.values():
        hist.Fill(value)

    canv = ROOT.TCanvas(name,name,400,400)
    hist.Draw()
    canv.SaveAs("/afs/fnal.gov/files/home/room2/jhirsch/public_html/temp/"+name+".eps")
    canv.SaveAs("/afs/fnal.gov/files/home/room2/jhirsch/public_html/temp/"+name+".gif")

if __name__ == "__main__":
    import sys
    col2hist(sys.argv[1], sys.argv[2])

    

