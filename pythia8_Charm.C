// pythia8 basic example
//Author: Andreas Morsch
//
// to run, do
//  root > .x pythia8.C
//
// Note that before executing this script, 
//   -the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
//   -the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
//
#include <iostream>
#include <cmath>
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "/Users/nserra/LHCb/Pythia8/pythia8175/include/Pythia.h"
#include "TPythia8.h"
#include "TRandom3.h"


using namespace Pythia8;
using namespace std;

using namespace ROOT;


void pythia8_Charm(Int_t nev  = 10000, Int_t ndeb = 1, Int_t seed=2210)
{
  char *p8dataenv = (char *) gSystem->Getenv("PYTHIA8DATA"); 
  Int_t numC = 0;
  Float_t c = 3.e11;

  //TRandom3 random = TRandom3(999);

   if (!p8dataenv) {
     char *p8env = (char *)gSystem->Getenv("PYTHIA8"); 
      if (!p8env) {
	cout<<"pythia8.C"<<"Environment variable PYTHIA8 must contain path to pythia directory!"<<endl;
         return;
      }
      TString p8d = p8env;
      p8d += "/xmldoc";
      gSystem->Setenv("PYTHIA8DATA", p8d);
   }
      
   char* path = (char *)gSystem->ExpandPathName("$PYTHIA8DATA");
   if (gSystem->AccessPathName(path)) {
     cout<<"pythia8.C"<<" Environment variable PYTHIA8DATA must contain path to $PYTHIA8/xmldoc directory !"<<endl;
      return;
   }
    
// Load libraries
   gSystem->Load("$PYTHIA8/lib/libpythia8");
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");

   ///// Defining the tree
   TTree *tree = new TTree("newTree", "newTree");
   Float_t CharmPx;
   Float_t CharmPy;
   Float_t CharmPz;
   Float_t CharmE;
   Float_t CharmPID;

   tree->Branch("CharmPx", &CharmPx, "CharmPx/F");
   tree->Branch("CharmPy", &CharmPy, "CharmPy/F");
   tree->Branch("CharmPz", &CharmPz, "CharmPz/F");
   tree->Branch("CharmE", &CharmE, "CharmE/F");
   tree->Branch("CharmPID", &CharmPID, "CharmPID/F");

   //// Production vertex 
   Float_t PVx;
   Float_t PVy;
   Float_t PVz;
   tree->Branch("PVx", &PVx, "PVx/F");
   tree->Branch("PVy", &PVy, "PVy/F");
   tree->Branch("PVz", &PVz, "PVz/F");

   // Decay vertex 
   Float_t SVx;
   Float_t SVy;
   Float_t SVz;
   Float_t DecayTime;
   tree->Branch("SVx", &SVx, "SVx/F");
   tree->Branch("SVy", &SVy, "SVy/F");
   tree->Branch("SVz", &SVz, "SVz/F");
   tree->Branch("DecayTime", &DecayTime, "DecayTime/F");


    
   //// Computing the boost
   Float_t mp = 0.9382;
   TLorentzVector CM_vector = TLorentzVector(0., 0., 400., TMath::Sqrt(400.*400.+mp*mp)+mp);
   TVector3 boostCM = CM_vector.BoostVector();


// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 1000);
// Create pythia8 object
   TPythia8* pythia8 = new TPythia8();
    
// Configure    
   pythia8->ReadString("SoftQCD:minBias = on");
   //pythia8->ReadString("SoftQCD:All = on");
   pythia8->ReadString("SoftQCD:singleDiffractive = on");
   pythia8->ReadString("SoftQCD:doubleDiffractive = on");
   //pythia8->ReadString("HardQCD:all= on");

   /// Here I set the seed
   pythia8->ReadString("Random:setSeed = on");
   char seed_string[50];
   sprintf(seed_string, "Random:seed = %d", seed);
   pythia8->ReadString(seed_string);

// Initialize 
    
   //pythia8->Initialize(2212 /* p */, 2212 /* p */, 450. /* GeV */, 10. /* GeV */);
   pythia8->Initialize(2212 /* p */, 2212 /* p */, 27.54 /* GeV */);

    
// Event loop
   for (Int_t iev = 0; iev < nev; iev++) {
      pythia8->GenerateEvent();
      if (iev < ndeb) pythia8->EventListing();
      pythia8->ImportParticles(particles,"All");
      Int_t np = particles->GetEntriesFast();
// Particle loop
      for (Int_t ip = 0; ip < np; ip++) {
         TParticle* part = (TParticle*) particles->At(ip);
	 
	 Int_t charm = TMath::Abs(part->GetPdgCode());
	 if(charm>410 && charm<440){
	   //cout<<" The particle has charm!!!!"<<endl;
	   numC++;
	   ////////
	   Float_t Px = part->Px();
	   Float_t Py = part->Py();
	   Float_t Pz = part->Pz();
	   Float_t E = part->Energy();
	   
	   
	   TLorentzVector pCharm(Px, Py, Pz, E);
	   pCharm.Boost(boostCM);
	   CharmPx = pCharm.X();
	   CharmPy = pCharm.Y();
	   CharmPz = pCharm.Z();
	   CharmE = pCharm.E();
	   CharmPID = charm;

	   TLorentzVector PV = TLorentzVector();
	   part->ProductionVertex(PV);
	   PVx = PV.X();
	   PVy = PV.Y();
	   PVz = PV.Z();


	   ///// Compute decay time
	   TLorentzVector SV = TLorentzVector(); 
	   Int_t ind = part->GetFirstDaughter();
	   TParticle* part2 = (TParticle*) particles->At(ind);
	   part2->ProductionVertex(SV);
	   Float_t gamma = pCharm.Gamma();
	   SV[2] =  gamma*SV[2];
	   SV[3] =  SV[3];
	   SVx = SV[0];
	   SVy = SV[1];
	   SVz = SV[2];
	   DecayTime = gamma*(SV[3]-PV[0])*(1./c)*1e9;	   	   
	   tree->Fill();
	 }
	 
	 

      }
   }

   pythia8->PrintStatistics();
   cout<<"Number of Charm  "<<numC<<endl;

   char fname[50];
   sprintf(fname, "CharmFixTarget_%d.root", seed);

   TFile *f = new TFile(fname, "RECREATE");
   tree->Write();
   tree->Print();
   f->Save();
   f->Close();

 }
