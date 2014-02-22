#!/usr/bin/env python

import ROOT
import PhysicsTools.PythonAnalysis as cmstools
import time
import getopt
import sys
import optparse
import commands
import array
import os

outdir = "./"
def draw1D(hdict):
    for hist in hdict.values():
        hname = hist.GetName()
        canv= ROOT.TCanvas(hname, hname, 400, 400)
        pad=canv.GetPad(0)
        pad.SetLogy()
        hist.Draw()
        name = outdir+hname+".gif"
        canv.SaveAs(name)
        name = outdir+hname+".eps"
        canv.SaveAs(name)
    return

def draw2D(hdict):
    for hist in hdict.values():
        hname = hist.GetName()
        canv= ROOT.TCanvas(hname, hname, 400, 400)
        hist.Draw("colz")
        name = outdir+hname+".gif"
        canv.SaveAs(name)
        name = outdir+hname+".eps"
        canv.SaveAs(name)
    return

def clean2D(hdict, min, max):
    for hist in hdict.values(): 
        hist.GetZaxis().SetRangeUser(min,max)
        # Set all bins that are identically 0.0 to -99.0 so that they
        # appear as white in histograms        
        # Start at 1, do not include under/overflow bins
        for j in range(1,nbineta+1):
            for k in range (1,nbinphi+1): 
                if hist.GetBinContent(int(j),int(k)) == 0.0:
                    hist.SetBinContent(int(j),int(k),-99.0)
    return

#######################
# Get options
#######################

if __name__ == "__main__":


    parser = optparse.OptionParser("usage: %prog [options]  file1.txt file2.txt\n")

    parser.add_option ('--t', dest='type', type='string',default = 'peds',
                       help="Type of conditions [qies, peds, qual, gain, emap]")
    parser.add_option ('--i', action="store_true", dest="identical", default=False,
                       help="Use flag if you expect files to be identical.")
    parser.add_option ('--pd', action='store_true', dest="printdiffs", default=False, 
                       help="Use flag with '-i' to print all differing lines")
    parser.add_option ('--d', dest='outdir', type='string',
                       default = './',
                       help="directory for output gifs and eps")
    parser.add_option ('--f', dest='fieldsToComp', default='',
                       help="Fields to compare other than defaults.")
    parser.add_option ('--l', dest='threshold', default='0.0001', type='float',
                       help="Set absolute difference for considering values the same.")
    parser.add_option ('--p', dest='prune', default='-1', type='int',
                       help="Set number of header lines to prune from top of files.")
    options, args = parser.parse_args()

    inputfile1 = args[0]
    inputfile2 = args[1]
    type       = options.type
    identical  = options.identical
    threshold  = options.threshold
    prune      = options.prune
    printdiffs = options.printdiffs
    outdir     = options.outdir+"/"

    #  Move ROOT initialization stuff here or
    #  ROOT --help gets reported instead of my --help!?
    #  ROOT.gSystem.Load("libFWCoreFWLite.so")
    #  ROOT.AutoLibraryLoader.enable()
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(111111)
    ROOT.gStyle.SetPalette(1)
    #####################################################################
    # Setup entries to examine for different conditions types
    #####################################################################

    if type == "peds":
        coord  = [0,1,2] # Location of identifier entries 
        comp   = [4,5,6,7,8,9,10,11] # Location of entries to compare
        sort   = [1,3] # Sort only on eta, phi, depth
        nComp  = len(comp)  # Number of entries to compare
        nCoord = len(coord) # Number of identifier entries (eta, phi, depth, detid)
        if prune < 0: prune  = 2 # Number of lines to remove from top of file
        expectedLength   = 13 # Expected length of each line
    elif type == "emap":
        nCoord = 10
        coord  = [0,1,2,4,5,6,7,9,10,11]
        comp   = [0]
        sort   = [1,3]
        nComp  = 0
        expectedLength   = 12
        if prune < 0: prune  = 1 # As dumped from DB
    elif type == "qies":
        coord  = [0,1,2]
        comp   = range (4, 35 + 1) # + 1 for exclusive limit
        sort   = [1,3]
        nCoord = len(coord)
        nComp  = len(comp)
        expectedLength   = 36
        if prune < 0: prune  = 4 # As dumped from DB
    else:
        print 'Unknown type : %s.  Exiting.' % type
        sys.exit()


    print 'Comparing %s files %s and %s.' % (type, inputfile1, inputfile2)
    print 'For %s files, we compare %s fields.\n' % (type, nComp)


    ##############################################
    # Open and check lengths of files
    ##############################################

    command_wc1 = "wc "+inputfile1
    command_wc2 = "wc "+inputfile2

    wc1 = commands.getoutput(command_wc1).split()
    wc2 = commands.getoutput(command_wc2).split()

    if wc1[0] != wc2[0]:
        print 'Files %s and %s are different lengths.  Exiting.' % (type, inputfile1, inputfile2)
        sys.exit()
    else:
        nLines = int(wc1[0])-prune

    ##############################################
    # Prepare files: sort, prune headers
    ##############################################

    command_touch1 = "touch "+inputfile1+".tmp"
    command_touch2 = "touch "+inputfile2+".tmp"

    command_sort1  = "tail -"+str(nLines)+" "+inputfile1+" | sort -k"\
                     +str(sort[0])+","+str(sort[1])+" "+" > "+inputfile1+".tmp"
    command_sort2  = "tail -"+str(nLines)+" "+inputfile2+" | sort -k"\
                     +str(sort[0])+","+str(sort[1])+" "+" > "+inputfile2+".tmp"

    os.system(command_touch1)
    os.system(command_touch2)
    os.system(command_sort1)
    os.system(command_sort2)

    file1 = open(inputfile1+".tmp", "r")
    file2 = open(inputfile2+".tmp", "r")


    #######################
    # Loop over files
    #######################

    print "Looping over %s lines.\n" % (nLines)

    if identical or type == "qies":
        dline1 = {}
        dline2 = {}
        aline1 = {}
        aline2 = {}
        nDiffs = 0
    if not identical:
        dlist = {}

    for i in range(0,nLines):
        sline1 = file1.readline()
        sline2 = file2.readline()
        line1  = sline1.split()
        line2  = sline2.split() 
        len1   = len(line1)
        len2   = len(line2)

        #######################
        # Check length of lines
        #######################

        if len1 != len2:
            print "Lines",i,"have different lengths."
            print "Exiting."
            sys.exit()
        elif len1 != expectedLength:
            print "Length of line",i,\
                  "differs from expected length",\
                  expectedLength,"for",type,"files."
            print "Length of line",i,"=",len1
            print sline1
            print "Exiting."
            sys.exit()

        ############################
        # Make important comparisons
        ############################

        # First check that each line is for the same channel.
        jlist = range(0,nCoord)
        for j in range(0,nCoord):

            # jlist is (eta, phi, depth)
            jlist[j] = int(line1[coord[j]])

            if line1[coord[j]] != line2[coord[j]]:
                print "Files are not sorted correctly."
                print "Lines",i,"do not match in entry",j,"."
                print "Line from",inputfile1,": ",line1
                print "Line from",inputfile2,": ",line2
                print "Exiting."
                sys.exit()

        # EITHER
        # Compare entries of interest for files that are expected
        # to be identical.  Fill dline dictionaries with with lines
        # that differ.

        if identical or type == "qies":
            # klist is list of values to be compared
            klist1     = range(0,nComp)
            klist2     = range(0,nComp)
            difference = False

            for k in range(0,nComp):
                klist1[k] = float(line1[comp[k]])
                klist2[k] = float(line2[comp[k]])

                diff = float(line1[comp[k]]) - float(line2[comp[k]])

                if abs(diff)>threshold:
                    difference = True
            if difference:
                nDiffs = nDiffs+1
                # dline[(eta, phi, depth)] = (values to be compared)
                dline1[tuple(jlist)] = tuple(klist1)
                dline2[tuple(jlist)] = tuple(klist2)

            ####################################
            # Put ALL lines in these dictionary
            ####################################

            # Why only 8938 lines in aline??
            print i
            aline1[i] = tuple(klist1)
            aline2[i] = tuple(klist2)


        # OR
        # Compare entries of interest for files that are expected
        # to differ.  Here we fill histograms of differences.
        if not identical:
           for k in range(0,nComp):
               # dlist will be used to fill histograms
               # fill dlist with absolute difference
               diff = (float(line1[comp[k]]) - float(line2[comp[k]]))
               chanTuple = (jlist[0],jlist[1],jlist[2],diff)

               # dlist dictionary dlist[nComp,line] = (eta, phi, depth, diff)
               dlist[(k,i)] = chanTuple

    # Remove tmp files.
    command_rm1 = "rm -f "+inputfile1+".tmp"
    command_rm2 = "rm -f "+inputfile2+".tmp"
    os.system(command_rm1)
    os.system(command_rm2)

    ################################################
    # Output information about compiled differences
    ################################################

    if identical:
        print "Found %i differences greater than %f between files.\n" % (nDiffs, threshold)

        if printdiffs:
            if type == qies:
                for key in dline1.keys():
                    print "(eta,phi,depth) =", key
                    print "  File 1:"
                    for entry in range(len(dline1[key])):
                        if (entry+1)%8 == 0 : print "%8.4f" % dline1[key][entry] # Newline
                        else : print "%8.4f" % dline1[key][entry],           # no newline
                    print "  File 2:"
                    for entry in range(len(dline2[key])):
                        if (entry+1)%8 == 0 : print "%8.4f" % dline2[key][entry] # Newline
                        else : print "%8.4f" % dline2[key][entry],           # no newline
                        print "\n"
            else:
                for key in dline1.keys(): 
                    print "(eta,phi,depth) =", key
                    print "  File 1:"
                    for entry in range(len(dline1[key])):
                        print "%8.4f" % dline1[key][entry]
                    print "  File 2:"
                    for entry in range(len(dline2[key])):
                        print "%8.4f" % dline2[key][entry]
                        
    else:

        #######################
        # Plot in ROOT
        #######################
        hibin =  500.
        lobin = -500.
        if type == "peds":

            ################################################
            # Define histograms and Canvases
            ################################################

            # Define 1D histograms/canvases of peds and peds widths
            hDict1_1D = {}
            hDict2_1D = {}
            for cap in range (0,4):
                name1 = "ped_mean_cap"+str(cap)
                name2 = "ped_width_cap"+str(cap)
                hDict1_1D[cap] = ROOT.TH1F(name1, name1, 100000,lobin,hibin)
                hDict2_1D[cap] = ROOT.TH1F(name2, name2, 100000,lobin,hibin)

            t_ped_1D    = "Difference in ped mean (ADC)" 
            t_width_1D  = "Difference in ped width (ADC)"

            # Define 2D histograms and canvases for peds and peds widths
            # in eta,phi
            nbineta = int(83)
            nbinphi = int(72)
            hDict1_2D = {}
            hDict2_2D = {}
            for depth in range(0,4):
                for cap in range (0,4):
                    name1 = "ped_mean_cap"+str(cap)+"_d"+str(depth)
                    name2 = "ped_width_cap"+str(cap)+"_d"+str(depth)
                    hDict1_2D[(cap,depth)] = ROOT.TH2F(name1,name1,
                                                       nbineta, -41.5, 41.5,
                                                       nbinphi, 0.5, 72.5)
                    hDict2_2D[(cap,depth)] = ROOT.TH2F(name2, name2,
                                                       nbineta, -41.5, 41.5,
                                                       nbinphi, 0.5, 72.5)

            t_ped_2D_x   = "Difference in ped mean (ADC) : i#eta"
            t_width_2D_x = "Difference in ped width (ADC) : i#eta"
            t_ped_2D_y   = "i#phi"
            t_width_2D_y = "i#phi"

            # Set titles by looping over dictionaries
            for hist in hDict1_1D.values(): hist.GetXaxis().SetTitle(t_ped_1D)
            for hist in hDict2_1D.values(): hist.GetXaxis().SetTitle(t_width_1D)
            for hist in hDict1_2D.values(): 
                hist.GetXaxis().SetTitle(t_ped_2D_x)
                hist.GetYaxis().SetTitle(t_ped_2D_y)
                hist.GetZaxis().SetRangeUser(-10.0,10.0)
                hist.GetZaxis().SetLabelSize(0.025)
            for hist in hDict2_2D.values(): 
                hist.GetXaxis().SetTitle(t_width_2D_x)
                hist.GetYaxis().SetTitle(t_width_2D_y)
                hist.GetZaxis().SetRangeUser(-10.0,10.0)
                hist.GetZaxis().SetLabelSize(0.025)

            ################################################
            # Fill histograms from dlist dictionary
            ################################################

            # dlist[nComp,line] = (eta, phi, depth, diff)
            MaxP = lobin
            MinP = hibin
            MaxW = lobin
            MinW = hibin
            for i in range(0,nComp): 
                # max/min are x-limits for each 1D histogram
                max = lobin
                min = hibin
                for j in range(0,nLines):

                    # Fill 1D histo for each entry to compare:
                    if max < dlist[i,j][3] :
                        max = dlist[i,j][3]
                    if min > dlist[i,j][3] : min = dlist[i,j][3]

                    if i < 4: 
                        hDict1_1D[i].Fill(dlist[i,j][3])
                    else:
                        # -4 b/c hDict is indexed by cap, dlist by nComp
                        hDict2_1D[i-4].Fill(dlist[i,j][3]) 

                    # Now fill 2D (eta,phi) histo for one of four depths for 
                    # each entry to compare using our dictionary:
                    depth = int(dlist[i,j][2])-1
                    # Only make plots for d=0-4, no ZDC channels
                    if depth >= 0 and depth < 4 : 
                        if i < 4: 
                            hDict1_2D[(i,depth)].Fill(float(dlist[i,j][0]),
                                                      float(dlist[i,j][1]),
                                                      float(dlist[i,j][3]))
                        else:
                            hDict2_2D[(i-4,depth)].Fill(float(dlist[i,j][0]), 
                                                        float(dlist[i,j][1]),
                                                        float(dlist[i,j][3]))
                # Set max/min for each 1D histogram
                if i < 4:
                    print min, max
                    hDict1_1D[i].GetXaxis().SetRangeUser(min-0.05*min,
                                                         max+0.05*max)
                else:
                    print min, max
                    hDict2_1D[i-4].GetXaxis().SetRangeUser(min-0.05*min,
                                                           max+0.05*max)

                # Get global Max and Min for 2D plots (Ped and Width)
                if i < 4:
                    if MaxP < max : MaxP = max
                    if MinP > min : MinP = min
                else :
                    if MaxW < max : MaxW = max
                    if MinW > min : MinW = min



            ################################################
            # After filling, clean up 2D histograms
            ################################################

            clean2D(hDict1_2D, MinP, MaxP)
            clean2D(hDict2_2D, MinW, MaxW)

            ###############################
            # Draw histograms on canvases
            ###############################

            draw1D(hDict1_1D)
            draw1D(hDict2_1D)

            ROOT.gStyle.SetOptStat(0)
            draw2D(hDict1_2D)
            draw2D(hDict2_2D)


        elif type == "emap":
            expectedLength   = 12
        elif type == "qies":
            expectedLength   = 36

            ################################################
            # Define histograms and Canvases
            ################################################

            lobin = -100
            hibin =  100
            # Define 1D histograms/canvases of qie
            hDict1_1D = {}
            hDict2_1D = {}
            for cap in range (0,4):
                for drange in range (0,4):
                    name1 = "qie_slope_cap"+str(cap)+"_r"+str(drange)
                    name2 = "qie_offset_cap"+str(cap)+"_r"+str(drange)
                    hDict1_1D[(cap,drange)] = ROOT.TH1F(name1, name1, 100000,lobin,hibin)
                    hDict2_1D[(cap,drange)] = ROOT.TH1F(name2, name2, 100000,lobin,hibin)

            t_slope_1D   = "Difference in QIE slope" 
            t_offset_1D  = "Difference in QIE offset"

            # Define 2D histograms and canvases for peds and peds widths
            # in eta,phi
            nbineta = int(83)
            nbinphi = int(72)
            hDict1_2D = {}
            hDict2_2D = {}
            for depth in range(0,4):
                for cap in range (0,4):
                    for drange in range (0,4):
                    name1 = "qie_slope__cap"+str(cap)+"_d"+str(depth)+"_r"+str(drange)
                    name2 = "qie_offset_cap"+str(cap)+"_d"+str(depth)+"_r"+str(drange)
                    hDict1_2D[(cap,depth,drange)] = ROOT.TH2F(name1,name1,
                                                       nbineta, -41.5, 41.5,
                                                       nbinphi, 0.5, 72.5)
                    hDict2_2D[(cap,depth,drange)] = ROOT.TH2F(name2, name2,
                                                       nbineta, -41.5, 41.5,
                                                       nbinphi, 0.5, 72.5)

            t_slope_2D_x   = "Difference in QIE slope : i#eta"
            t_offset_2D_x  = "Difference in QIE offset : i#eta"
            t_slope_2D_y   = "i#phi"
            t_offset_2D_y  = "i#phi"

            # Set titles by looping over dictionaries
            for hist in hDict1_1D.values(): hist.GetXaxis().SetTitle(t_slope_1D)
            for hist in hDict2_1D.values(): hist.GetXaxis().SetTitle(t_offset_1D)
            for hist in hDict1_2D.values(): 
                hist.GetXaxis().SetTitle(t_slope_2D_x)
                hist.GetYaxis().SetTitle(t_slope_2D_y)
                hist.GetZaxis().SetRangeUser(-100.0,100.0)
                hist.GetZaxis().SetLabelSize(0.025)
            for hist in hDict2_2D.values(): 
                hist.GetXaxis().SetTitle(t_offset_2D_x)
                hist.GetYaxis().SetTitle(t_offset_2D_y)
                hist.GetZaxis().SetRangeUser(-100.0,100.0)
                hist.GetZaxis().SetLabelSize(0.025)

            ################################################
            # Fill histograms from dlist dictionary
            ################################################

            # dlist[nComp,line] = (eta, phi, depth, diff)
            MaxP = lobin
            MinP = hibin
            MaxW = lobin
            MinW = hibin
            for i in range(0,nComp): 
                # max/min are x-limits for each 1D histogram
                max = lobin
                min = hibin
                for j in range(0,nLines):

                    # Fill 1D histo for each entry to compare:
                    if max < dlist[i,j][3] : max = dlist[i,j][3]
                    if min > dlist[i,j][3] : min = dlist[i,j][3]

                    if i < 16: 
                        hDict1_1D[i].Fill(dlist[i,j][3])
                    else:
                        # -16 b/c hDict is indexed by cap, dlist by nComp
                        hDict2_1D[i-16].Fill(dlist[i,j][3]) 

                    # Now fill 2D (eta,phi) histo for one of four depths for 
                    # each entry to compare using our dictionary:
                    depth = int(dlist[i,j][2])-1
                    # Only make plots for d=0-4, no ZDC channels
                    if depth >= 0 and depth < 4 : 
                        if i < 4: 
                            hDict1_2D[(i,depth)].Fill(float(dlist[i,j][0]),
                                                      float(dlist[i,j][1]),
                                                      float(dlist[i,j][3]))
                        else:
                            hDict2_2D[(i-4,depth)].Fill(float(dlist[i,j][0]), 
                                                        float(dlist[i,j][1]),
                                                        float(dlist[i,j][3]))
                # Set max/min for each 1D histogram
                if i < 4:
                    print min, max
                    hDict1_1D[i].GetXaxis().SetRangeUser(min-0.05*min,
                                                         max+0.05*max)
                else:
                    print min, max
                    hDict2_1D[i-4].GetXaxis().SetRangeUser(min-0.05*min,
                                                           max+0.05*max)

                # Get global Max and Min for 2D plots (Slope and Offset)
                if i < 4:
                    if MaxP < max : MaxP = max
                    if MinP > min : MinP = min
                else :
                    if MaxW < max : MaxW = max
                    if MinW > min : MinW = min

            


            #################################################################
            # For QIES draw 32 plots each with 9094 lines :
            # 4 ranges, 4 caps for each channel for two files
            #################################################################

            # This might not be as useful as we thought, so for now, I skip it.
            drawQieLines = False
            if drawQieLines:

                # Format of line is
                #  (slope0)  (slope1)  (slope2)  (slope3)
                # (offset0) (offset1) (offset2) (offset3)

                # input aline, number of canvases
                nCanvs = 16
                nLinePerCanv = len(aline2)/nCanvs
                nline = 0
                cDict1 = {}
                for ncanv in range (16):
                    if ncanv == 1 : break
                    # Create canvas
                    name1 = "c1_"+str(ncanv)
                    cDict1[ncanv] = ROOT.TCanvas(name1,name1,400,400)
                    # Create dict of functions to plot
                    fDict = {}
                    nfunc = 0
                    # Loop over lines to plot on each canvas
                    for line in range(nLinePerCanv):
                        nline = nline+1
                        if nline%1000==0: print "processing line", line

                        # Loop over 16 slopes/offsets in each line of qie file
                        for i in range (16):
                            slope  = str(aline2[nline][i])
                            offset = str(aline2[nline][i+int(16)])
                            curve = slope+"*x+"+offset
                            name = "f"+str(line)+"_"+str(i)
                            fDict[nfunc] = ROOT.TF1(name,curve,-50000.,50000.)
                            fDict[nfunc].SetLineWidth(1)
                            cDict1[ncanv].cd()
                            cname = cDict1[ncanv].GetName()
                            if nfunc == 1: fDict[nfunc].Draw()
                            else : fDict[nfunc].Draw("same")
                            nfunc = nfunc+1
                    # Draw canvas
                    cDict1[ncanv].SaveAs("/afs/fnal.gov/files/home/room2/jhirsch/public_html/temp/"+cname+".gif")
                
k
