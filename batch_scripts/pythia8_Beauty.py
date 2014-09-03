#!/bin/python

# Check $ROOTSYS/tutorials/pythia/pythia8.C
# for example program
from __future__ import division
import sys
import os
import ROOT as r
import numpy as np
import math
import gc

def loadLibs():
    r.gSystem.Load("$PYTHIA8/lib/libpythia8")
    r.gSystem.Load("libEG")
    r.gSystem.Load("libEGPythia8")
    r.gSystem.Load("$PYTHIA8/rootexamples/pythiaDict.so")

def readSettings(pythia, seed):
    pythia.ReadString("SoftQCD:nonDiffractive = on")
    pythia.ReadString("HardQCD:all = on")
    #pythia.ReadString("SoftQCD:singleDiffractive = on")
    #pythia.ReadString("SoftQCD:doubleDiffractive = on")
    pythia.ReadString("Beams:idA = 2212")
    pythia.ReadString("Beams:idB = 2212")
    pythia.ReadString("Beams:frameType = 2")
    pythia.ReadString("Beams:eA = 400.")
    pythia.ReadString("Beams:eB = 0.")
    pythia.ReadString("Random:setSeed = on")
    pythia.ReadString("Random:seed = %s"%seed)
    # No event record printout.
    pythia.ReadString("Next:numberShowInfo = 0")
    pythia.ReadString("Next:numberShowProcess = 0")
    pythia.ReadString("Next:numberShowEvent = 0")

def makeBTree(howmany, seed):
    loadLibs()
    pythia = r.TPythia8()
    readSettings(pythia, seed)
    pyt = pythia.Pythia8()
    nev = howmany
    pyt.init()

    # Variables for the tree
    BeautyPx  = np.zeros(1, dtype=float)
    BeautyPy  = np.zeros(1, dtype=float)
    BeautyPz  = np.zeros(1, dtype=float)
    BeautyE   = np.zeros(1, dtype=float)
    BeautyPID = np.zeros(1, dtype=int)
    PVx = np.zeros(1, dtype=float)
    PVy = np.zeros(1, dtype=float)
    PVz = np.zeros(1, dtype=float)
    SVx = np.zeros(1, dtype=float)
    SVy = np.zeros(1, dtype=float)
    SVz = np.zeros(1, dtype=float)
    DecayTime = np.zeros(1, dtype=float)
    basepath = "/afs/cern.ch/work/e/egraveri/public/SHiP/batch_scripts/out/PythiaData/"
    filepath = basepath + "btree-%s.root"%seed
    print "Writing tree to %s ..."%filepath

    if os.path.isfile(filepath):
        outfile = r.TFile(filepath,"update")
        tree = outfile.Get("newTree")
        #tree.GetEntry(tree.GetEntries()-1)
        tree.SetBranchAddress("BeautyPx", BeautyPx)
        tree.SetBranchAddress("BeautyPy", BeautyPy)
        tree.SetBranchAddress("BeautyPz", BeautyPz)
        tree.SetBranchAddress("BeautyE", BeautyE)
        tree.SetBranchAddress("BeautyPID", BeautyPID)
        tree.SetBranchAddress("PVx", PVx)
        tree.SetBranchAddress("PVy", PVy)
        tree.SetBranchAddress("PVz", PVz)
        tree.SetBranchAddress("SVx", SVx)
        tree.SetBranchAddress("SVy", SVy)
        tree.SetBranchAddress("SVz", SVz)
        tree.SetBranchAddress("DecayTime", DecayTime)
    else:
        outfile = r.TFile(filepath,"recreate")
        tree = r.TTree("newTree","newTree")
        tree.Branch("BeautyPx", BeautyPx, "BeautyPx/D")
        tree.Branch("BeautyPy", BeautyPy, "BeautyPy/D")
        tree.Branch("BeautyPz", BeautyPz, "BeautyPz/D")
        tree.Branch("BeautyE", BeautyE,    "BeautyE/D")
        tree.Branch("BeautyPID", BeautyPID, "BeautyPID/I")
        tree.Branch("PVx", PVx, "PVx/D")
        tree.Branch("PVy", PVy, "PVy/D")
        tree.Branch("PVz", PVz, "PVz/D")
        tree.Branch("SVx", SVx, "SVx/D")
        tree.Branch("SVy", SVy, "SVy/D")
        tree.Branch("SVz", SVz, "SVz/D")
        tree.Branch("DecayTime", DecayTime, "DecayTime/D")

    nb = 0
    nbb = 0

    # Create the tree branches
    for i in xrange(nev):
        if not (i%500):
            print "Generating event %s..."%i
            gc.collect()
        pythia.GenerateEvent()
        #names = [event[part].name() for part in pyt.event[0].daughterList()]
        for p in xrange(pyt.event.size()):
            part = pyt.event[p]
            #print "ID is ", part.id()
            if part.pz() > 0.:
                if (510 < part.id()%1000 < 550):
                    print "ID is ", part.id()
                    nb += 1
                if (550 < part.id()%1000 < 560):
                    print part.id()
                    nbb += 1
                if (510 < part.id()%1000 < 560):
                    BeautyPx[0] = float(part.px())
                    BeautyPy[0] = float(part.py())
                    BeautyPz[0] = float(part.pz())
                    BeautyE[0] =  float(part.e())
                    BeautyPID[0] = int(part.id())
                    PVx[0] = float(part.vProd().px()/1000.) #in mm -> in m
                    PVy[0] = float(part.vProd().py()/1000.) #in mm -> in m
                    PVz[0] = float(part.vProd().pz()/1000.) #in mm -> in m
                    dFirst = part.daughter1()
                    dLast = part.daughter2()
                    if (dFirst != dLast) and (dFirst>0) and (dLast>0): #this particle has children
                        firstDaughter = pyt.event[part.daughter1()]
                        sv = firstDaughter.vProd()
                        if sv.pz() >0.:
                            SVx[0] = float(sv.px()/1000.) #in mm -> in m
                            SVy[0] = float(sv.py()/1000.) #in mm -> in m
                            SVz[0] = float(sv.pz()/1000.) #in mm -> in m
                            DecayTime[0] = float(sv.e())
                            print BeautyPx[0], BeautyPy[0], BeautyPz[0], BeautyE[0], SVx[0], SVy[0], SVz[0], DecayTime[0]
                            tree.Fill()

    print "%s b meson, %s bb mesons"%(nb, nbb)
    tree.Write("",5) # TObject::kOverwrite
    outfile.Close()

if __name__ == '__main__':
    howmany, seed = int(sys.argv[1]), int(sys.argv[2])
    print "Generating %s POTs with seed %s"%(howmany, seed)
    makeBTree(howmany, seed)
