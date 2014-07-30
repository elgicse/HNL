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
    r.gSystem.Load("$PYTHIA8/PythiaDict/pythiaDict.so")

def readSettings(pythia):
    pythia.ReadString("SoftQCD:nonDiffractive = on")
    #pythia.ReadString("SoftQCD:singleDiffractive = on")
    #pythia.ReadString("SoftQCD:doubleDiffractive = on")
    pythia.ReadString("Beams:idA = 2212")
    pythia.ReadString("Beams:idB = 2212")
    pythia.ReadString("Beams:frameType = 2")
    pythia.ReadString("Beams:eA = 400.")
    pythia.ReadString("Beams:eB = 0.")
    pythia.ReadString("Random:setSeed = on")
    pythia.ReadString("Random:seed = 0")
    # No event record printout.
    pythia.ReadString("Next:numberShowInfo = 0")
    pythia.ReadString("Next:numberShowProcess = 0")
    pythia.ReadString("Next:numberShowEvent = 0")

def makeBTree(howmany):
    loadLibs()
    pythia = r.TPythia8()
    readSettings(pythia)
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


    if os.path.isfile("out/PythiaData/btree.root"):
        outfile = r.TFile("out/PythiaData/btree.root","update")
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
        outfile = r.TFile("out/PythiaData/btree.root","recreate")
        tree = r.TTree("newTree","newTree")
        tree.Branch("BeautyPx", BeautyPx, "BeautyPx/F")
        tree.Branch("BeautyPy", BeautyPy, "BeautyPy/F")
        tree.Branch("BeautyPz", BeautyPz, "BeautyPz/F")
        tree.Branch("BeautyE", BeautyE, "BeautyE/F")
        tree.Branch("BeautyPID", BeautyPID, "BeautyPID/I")
        tree.Branch("PVx", PVx, "PVx/F")
        tree.Branch("PVy", PVy, "PVy/F")
        tree.Branch("PVz", PVz, "PVz/F")
        tree.Branch("SVx", SVx, "SVx/F")
        tree.Branch("SVy", SVy, "SVy/F")
        tree.Branch("SVz", SVz, "SVz/F")
        tree.Branch("DecayTime", DecayTime, "DecayTime/F")

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
            if (510 < part.id()%1000 < 550):
                print "ID is ", part.id()
                nb += 1
            if (550 < part.id()%1000 < 560):
                print part.id()
                nbb += 1
                #mum = pyt.event[part.mother1()]
                #children = mum.daughterList()
                ## Conto un solo fotone se ce n'e' piu' d'uno nello stesso decay!
                #childrenList = []
                #mumId[0] = mum.id()
                #mumMass[0] = mum.m0()
                #childrenId.resize(0)
                #childrenMass.resize(0)
                #for c in xrange(len(children)):
                #    childrenList.append(pyt.event[children[c]])
                #    childrenId.push_back(pyt.event[children[c]].id())
                #    childrenMass.push_back(pyt.event[children[c]].m0())
                #childrenNames = [child.name() for child in childrenList]
                #decayName.Resize(0)
                #decayName.Append(mum.name() + " " + " ".join(childrenNames[:]))
                #tree.Fill()

    print "%s b meson, %s bb mesons"%(nb, nbb)
    tree.Write("",5) # TObject::kOverwrite
    outfile.Close()
