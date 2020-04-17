//____________________________________________________________________
//
// Generation process: H -> ZZ -> 4 mu
// Using Pythia6 with AliROOT
//
// Author: Héctor Bello Martínez <hbelloma@cern.ch>
// Update: 2018-10-12
//
// Modification by: Sergio ...
// Date:
//
//-------------------------------------------------------------------
// To make an event sample (of size 1000) do
//    shell> aliroot
//    root [0] .L pythia6ProcH2ZZ24MU_HEX.C
//    root [1]  gSystem->Load("libpythia6");
//    root [2]  gSystem->Load("libAliPythia6");
//    root [3] pythiaExample(1000)
// will execute makeEventSample(1000) and showEventSample()
//
//____________________________________________________________________


#ifndef __CINT__
#include "TApplication.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <cstdlib>

#include "TMCParticle.h"
using namespace std;
#endif

//------------- define parametros ---
#define FILENAME    "A0nametest_ncom.root "
#define ECM        7000
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define HISTNAME   "ptSpectraformu"
#define HISTNAMEZ  "ptSpectraforZ"
#define PDGNUMBER 13 


// This function just load the needed libraries if we're executing from
void loadLibraries()
{
#ifdef __CINT__
  // Load the Event Generator abstraction library, Pythia 6
  gSystem->Load("libEG"); 
  gSystem->Load("libEGPythia6"); 
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");

#endif
}

// nEvents is how many events we want.
int makeEventSample(Int_t nEvents)
{
  
  Bool_t isfast=kFALSE;
 if(nEvents<101) isfast=kTRUE;
  loadLibraries();

  TPythia6* pythia = new TPythia6;

  pythia->SetMSEL(0); // full user controll;
//---------------------------------------------------------------------------
//-         Realistic process during the generation (1=on,0=off)            -
//---------------------------------------------------------------------------
/* 
    pythia->SetMSTP(61,0);  // switch off ISR
    pythia->SetMSTP(71,0);  // switch off FSR
    pythia->SetMSTP(81,0);  // switch off multiple interactions
    pythia->SetMSTP(111,0); // Switch off fragmentation
*/


//-------------------------------------------------------------------  
//-               controled subprocess                             --
//-------------------------------------------------------------------
 // B) Hard QCD subprocess 
 // pythia->SetMSEL(2); 
 //  pythia->SetMSUB(11, 1); //qiqj->qiqj
 //  pythia->SetMSUB(12, 1); //qibarqi->qkbarqk
 //  pythia->SetMSUB(13, 1); //qibarqi->gg
 //  pythia->SetMSUB(28, 1); //qig->qig
 //  pythia->SetMSUB(53, 1); //gg->qibarqi
//   pythia->SetMSUB(68, 1); //gg->gg

 // D) little SM Higgs process
 pythia->SetMSUB(102, 1);// g + g -> h0

// -----posible background in z0 prod (uncoment all for crosscheck)
//   pythia->SetMSUB(1,1);   // ff -> single Z/gamma*
//   pythia->SetMSUB(15,1);//fifj->gZ0
//   pythia->SetMSUB(19,1);//fifj-->gammaZ0
   pythia->SetMSUB(22,1);//fibarfi-->2Z0
//   pythia->SetMSUB(23,1);//fibarfi-->W+Z0
//   pythia->SetMSUB(24,1);//fibarfi-->h0Z0
//   pythia->SetMSUB(30,1);//fig-->fiZ0 // 30+22,23,24 producen caida el bk hasta 150
 //  pythia->SetMSUB(35,1);//figamma-->fiZ0
// more see MSUB 40,45,50,70,71,73,74,76, 101,116,117,141 etc

//---------------------------------------------------------
//-         some parameters for the resonance            --
//---------------------------------------------------------

  pythia->SetPMAS(6,1, 172);  // mass of TOP
  pythia->SetPMAS(25, 1, 125);  // mass of Higgs
  pythia->SetCKIN(1,100); //allowed value on first resonance (h0)
  pythia->SetCKIN(2,150);

//--------------------------------------------------------------
//--      setting general switches via SetMSTP() D is default --
//--------------------------------------------------------------
// pythia->SetMSTP(1,3);//max num of gen D=3 set<4
// pythia->SetMSTP(2,1);//calculation alpha_s D=1(0=fixed,1=1stOrder,2=2ndOrder)
// pythia->SetMSTP(3,2);//Lambda in alpha_s D=2 (acording pdf)
// pythia->SetMSTP(4,0);//neutral higgs sector D=0 StandarMod couplings
// pythia->SetMSTP(5,0);//anomalous couplings see pythia6 ec 8.6.5
// pythia->SetMSTP(7,0);// HF subprocess D=0
// pythia->SetMSTP(8,0);// EW param D=0 (a_EM(q2) and fixed sin2theta_W)
// pythia->SetMSTP(9,0);//4th generation inclusion D=0 (no)
// pythia->SetMSTP(11,1);//use electron pdf in e+e- and ep interact.
// pythia->SetMSTP(12,0);// use e´,e+,quark and gluon dist func inside e- D=0 
// pythia->SetMSTP(13,1);//choice of Q2 range for e- to radiate gammas
// pythia->SetMSTP(14,30);//structure of inc. phot beam (D=30)
// pythia->SetMSTP(15,0);//anomalous photon component
// pythia->SetMSTP(16,1);//fract of moment tak by phot radiated of a lept
// pythia->SetMSTP(17,4);//extract factor of virtualphotons
// pythia->SetMSTP(18,3);//choice of pTmin for dir procc
// pythia->SetMSTP(19,4);//partonic cross section in DIS process 99
// pythia->SetMSTP(20,3);//supress of reslvd cross sect to comp overlap
// pythia->SetMSTP(21,1);//nature of fermion fermion scatt
// pythia->SetMSTP(22,0);//spetial overide of norm Q2 se sec 10.4 pythia6
// pythia->SetMSTP(23,1);// DIS (10y83)
// pythia->SetMSTP(25,0);//ang decay corr in Higgs dec to 4fermion(0=scalar dec)
// pythia->SetMSTP(31,1);//param of tot elast and diffrac cross sect 
// pythia->SetMSTP(32,8);//Q2 def in hard scatt
// pythia->SetMSTP(33,0);//K factors for hard cross section
// pythia->SetMSTP(34,1);//interference term in matrix elemnts of QCD procc D=1 inc
// pythia->SetMSTP(35,0);//threshold of HF production
// pythia->SetMSTP(36,2);//reg of alpha_s in Q2->0
// pythia->SetMSTP(37,1);//incl of quark masses in h0 prod
// pythia->SetMSTP(38,5);//q-loop masses in box (5-efective masses)
// pythia->SetMSTP(39,2);//q2 for pd and IS part show
// pythia->SetMSTP(40,0);//coloumb correction in W+W-
// pythia->SetMSTP(41,2);//master switch of all resonance decays
// pythia->SetMSTP(42,1);//mass treatment for finite width FS resonance(D=1 Breit-Wigner)
// ----------- try uncoment this two --------------------
// pythia->SetMSTP(43,2);//Z0/gamma* interference matrix D=3 (D=2 only Z0) 
// pythia->SetMSTP(44,2);//Z'/Z0/gamm* interf matrix elements D=7
// pythia->SetMSTP(45,3);//treatment of WW->WW structure
// pythia->SetMSTP(46,1);//treatment of VV->V´v´ D=1 (all graphs)
// .....
// ....
// pythia->SetMSTP(115,0);//colour rearangement w+w- z0z0
// ...
// ... until 185


//-------------------------------------------------------- 
//--         Setting decay modes by SetMDME()           --
//--------------------------------------------------------
   // Force h0 -> ZZ
   for (Int_t i = 210; i <= 288; ++i)
      pythia->SetMDME(i, 1, 0);
      pythia->SetMDME(225, 1, 1);
  // only muonic decays Z -> mumu
    for (Int_t i = 174; i <= 189; ++i) //188,189
     pythia->SetMDME(i, 1, 0);
     pythia->SetMDME(184, 1, 1);

  pythia->Initialize("cms", "p", "p", ECM);
 
  TFile* file = TFile::Open(FILENAME, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("makeEventSample", "Couldn;t open file %s", FILENAME);
    return 1;
  }
  TTree* tree = new TTree(TREENAME, "Pythia 6 tree");

  int NbinIVplot = 100;
  int Nbinsz = 100;
  int Nbins4mu = 500;  
  int NbinsH =500;
  int Nbinsy =50;
  // Histograms
 TH1D* Hevents = new TH1D("Hevents", "" ,10, 0,10);
 TH1F* ggHcross = new TH1F("ggHcrossec","gg->h^{0} production", 500,0.0,500.0);
 TH1F* hmasainvH = new TH1F("hinvmassH"," h^{0} invariant mass from Z^{0}+Z^{0}", NbinsH,0.0,250.0);
 TH1F* hmasainvH4mu = new TH1F("hinvmassH4mu"," h^{0} invariant mass from Z^{0}+Z^{0}->4#mu", Nbins4mu,0.0,250.0);
 TH1F* hmasainv2mu = new TH1F("hinvmass2mu"," invariant mass from 2#mu", Nbins4mu,0.0,250.0);
 TH1F* hmasainvZa= new TH1F("hinvmassZa","1st Z^{0} invariant mass from #mu^{+}+#mu^{-}",Nbinsz,0.0,140.0);
 TH1F* hmasainvZb= new TH1F("hinvmassZb","2nd Z^{0} invariant mass from #mu^{+}+#mu^{-}",Nbinsz,0.0,140.0);

 TH1F* hdNdyZa= new TH1F("hdNdyZa","Rapidity for Z^{0}",Nbinsy,-5.0,5.0);
 TH1F* hdNdyZb= new TH1F("hdNdyZb","Rapidity for Z^{0}",Nbinsy,-5.0,5.0);

 //other histos to fill (Sergio task)
   TH1F *Cosine1plot = new TH1F("Cosine1graph","Cos(#theta) angle of muon1 in Z1 rest frame",NbinIVplot,-1,1);
   TH1F *Cosine2plot = new TH1F("Cosine2graph","Cos(#theta) angle of muon2 in Z2 rest frame",NbinIVplot,-1,1);  
   TH2F *Correlationplot = new TH2F("Correlationplot","Theta angle muon 1 versus theta angle muon 2",NbinIVplot,-1,1,NbinIVplot,-1,1);
   TH1F *Phiplot1 = new TH1F("Phiplot1","Phi angle of muon 1 in Z1 rest frame",NbinIVplot,0,2*TMath::Pi());
   TH1F *Phiplot2 = new TH1F("Phiplot2","Phi angle of muon 2 in Z2 rest frame",NbinIVplot,0,2*TMath::Pi());

  ggHcross->Sumw2();
  ggHcross->SetYTitle("dN_{h^{0}}/dM_{h^{0}}");
  ggHcross->SetXTitle("M_{h^{0}}[GeV/c^{2}]");
  hmasainvH->Sumw2();
  hmasainvH->SetXTitle("M_{invh^{0}}[GeV/c^{2}]");
  hmasainvH->SetYTitle(Form("Events/%1.1f GeV",250./NbinsH));
  hmasainvH4mu->Sumw2();
  hmasainvH4mu->SetXTitle("M_{inv4#mu}[GeV/c^{2}]");
  hmasainvH4mu->SetYTitle(Form("Events/ %1.1f GeV",250./Nbins4mu));
  hmasainv2mu->Sumw2();
  hmasainv2mu->SetXTitle("M_{inv2#mu}[GeV/c^{2}]");
  hmasainv2mu->SetYTitle(Form("Events/ %1.1f GeV",250./Nbins4mu));
  hmasainvZa->Sumw2();
  hmasainvZa->SetXTitle("M_{invZ^{0}}[GeV/c^{2}]");
  hmasainvZa->SetYTitle(Form("Events/%1.1f GeV",140./Nbinsz));
  hmasainvZb->Sumw2();
  hmasainvZb->SetXTitle("M_{invZ^{0}}[GeV/c^{2}]");
  hmasainvZb->SetYTitle(Form("Events/%1.1f GeV",140./Nbinsz));
//declarando algunos histogramas

 TH1F* hetach = new TH1F("hetach", "Pseudorapidity", 240, -12., 12.);
 TH1F* hych = new TH1F("hych", "rapidity", 240, -12., 12.);
 TH1F* hphich = new TH1F("hphich", "phiangle", 700, 0.0, 7.0);
 TH2F* hyphiwpt= new TH2F("hyphiwpt","#eta vs #phi for Muons",100,-7.0,7.0,100,0.0,7.0);

 TH1F* hetachZ = new TH1F("hetachZ", "PseudorapidityZ", 240, -12., 12.);
 TH1F* hychZ = new TH1F("hychZ", "rapidityZ", 240, -12., 12.);
 TH1F* hphichZ = new TH1F("hphichZ", "phiangleZ", 700, 0.0, 7.0);
 TH2F* hyphiwptZ= new TH2F("hyphiwptZ","#eta vs #phi for Z^0",100,-7.0,7.0,100,0.0,7.0);

 TH1F* hetachH = new TH1F("hetachH", "PseudorapidityH", 240, -12., 12.);
 TH1F* hychH = new TH1F("hychH", "rapidityH", 240, -12., 12.);
 TH1F* hphichH = new TH1F("hphichH", "phiangleH", 700, 0.0, 7.0);
 TH2F* hyphiwptH= new TH2F("hyphiwptH","#eta vs #phi for H^0",100,-7.0,7.0,100,0.0,7.0);



  TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
  tree->Branch(BRANCHNAME, &particles);
   Float_t fy=-999999;  
// para obtener eta y phi antes del loop de particulas

 


   Double_t eta=-999.;
   Double_t y=-999.;
   //Double_t phimv=-999.;
   //Double_t phi=-999.;

   Double_t etach=-999.;
   Double_t ych=-999.;
   Double_t phichmv=-999.;
   Double_t phich=-999.;

   Double_t etachZ=-999.;
   Double_t ychZ=-999.;
   Double_t phichmvZ=-999.;
   Double_t phichZ=-999.;

   Double_t etachH=-999.;
   Double_t ychH=-999.;
   Double_t phichmvH=-999.;
   Double_t phichH=-999.;


   TBranch *newBranch = tree->Branch("fy", &fy, "fy/F");

  pythia->Pystat(4);  // prints a table of kinem cuts CKIN(I) p.135 and p.251
  Int_t nHiggs=0; 
  Int_t eventH=0;

  // Now we make some events (for of events)
  for (Int_t i = 0; i < nEvents; i++) {
    // Show how many events we have.
  if(isfast){
      cout << "Event # " << i << endl;
  }
  else{ if (i % 100 == 0)
      cout << "Event # " << i << endl;
  }

    pythia->GenerateEvent();
    if (i == 10) {
      pythia->Pylist(1);   // list of pythia generated event
     // pythia->Pylist(12);  // to see the list of all decay modes for all defined particles
    }

    // declaration of 4-vectors
    TLorentzVector vechiggs,vecmuon,vecz0,vecza,veczb,vecmum1,vecmum2,vecmup1,vecmup2,veczatomu,veczbtomu,vechtoz,vecto4mu;
    
    Int_t np = particles->GetEntriesFast();
    Int_t nh0=0;
    Int_t nz=0;
    Int_t nmum=0;
    Int_t nmup=0;
    Int_t N_muons =0;
    Hevents->Fill(0);

    // loop over particles (particles for)
    for (Int_t ip=0; ip<np; ip++){
      TMCParticle* part = (TMCParticle*) particles->At(ip);
      //if (part->GetParent()==0) continue;
      TMCParticle* parent=dynamic_cast<TMCParticle*>((*particles)[part->GetParent()]);

      
      //if (parent->GetParent()==0) continue;
      TMCParticle* Grand=dynamic_cast<TMCParticle*>((*particles)[parent->GetParent()]);
      TMCParticle* Grand2=dynamic_cast<TMCParticle*>((*particles)[Grand->GetParent()]);
      Int_t idp = part->GetKF();// pdg number
      Int_t pstat = part->GetKS();//particle status
      Int_t idpparent=parent->GetKF();// pdg number
      Int_t idpGrand=Grand->GetKF();// pdg number
      Int_t idpGrand2=Grand2->GetKF();// pdg number
      //cout << "idp= " << idp << ",  parent= " << parent << ", idpparent= " << idpparent << ",  idpGrand= " << idpGrand <<endl;

       Double_t px=part->GetPx();
       Double_t py=part->GetPy();
       Double_t pz=part->GetPz();
       Double_t e=part->GetEnergy();
       Double_t pt=TMath::Sqrt(px*px+py*py);
   	 Double_t costheta=pz/(TMath::Sqrt(pt*pt+pz*pz));
    	y= TMath::ATan(pz/e);      //momento/energia
    	eta=TMath::ATanH(costheta);
       if (TMath::Abs(eta)>7)continue;
		if (pt<25)continue;
     if (pstat<21)continue; //no iniciales
   
    
      fy=TMath::Sqrt(part->GetPx()*part->GetPx()+part->GetPy()*part->GetPy()); //corregir por la funcion correcta para rapidity (Sergio task)
      newBranch->Fill();

      if (idp==21 && pstat==21) //gluons
      if(isfast) cout<<"gluon is produced"<<endl; 
      
      //lets count the produced Higgs
      if (idp==25 && pstat==21){
      // if(isfast) 
      cout<<"---Bingo! a Higgs is produced"<< "evento= " << i << ",  eta= " << eta<<endl;
       if (eventH==0) {
           eventH=i;
	 cout << "EventH= " << eventH << endl;
         }
        
       Double_t pxh = part->GetPx();
       Double_t pyh = part->GetPy();
       Double_t pzh = part->GetPz();
       Double_t eh = part->GetEnergy();
       vechiggs.SetPxPyPzE(pxh,pyh,pzh,eh);
       Double_t mh=vechiggs.Mag();
       ggHcross->Fill(mh);
       nh0++;
      }
      
      // lets check produced mu-
      if (idp==13 ){
      	nmum++;
        if(nmum % 2 != 0){
	  if(isfast) cout<<"1st  mu minus = "<<endl;
          Double_t px1mum = part->GetPx();
          Double_t py1mum = part->GetPy();
          Double_t pz1mum = part->GetPz();
          Double_t e1mum = part->GetEnergy();
          vecmum1.SetPxPyPzE(px1mum,py1mum,pz1mum,e1mum);
          vecmuon.SetPxPyPzE(px1mum,py1mum,pz1mum,e1mum);
          if(isfast) cout<<"vectormum1-> PX,PY,PZ,EN = "<<px1mum<<", "<<py1mum<<", "<<pz1mum<<", "<<e1mum<<endl;
	}else{
	  if(isfast) cout<<"2nd  mu minus = "<<endl;
          Double_t px2mum = part->GetPx();
          Double_t py2mum = part->GetPy();
          Double_t pz2mum = part->GetPz();
          Double_t e2mum = part->GetEnergy();
          vecmum2.SetPxPyPzE(px2mum,py2mum,pz2mum,e2mum);
          vecmuon.SetPxPyPzE(px2mum,py2mum,pz2mum,e2mum);
          if(isfast) cout<<"vectormum2-> PX,PY,PZ,EN = "<<px2mum<<", "<<py2mum<<", "<<pz2mum<<", "<<e2mum<<endl;
	}
      }
      
      // lets check produced mu+
      if(idp==-13){      
        nmup++;
        if(nmup % 2 != 0){
          if(isfast) cout<<"1st  mu plus = "<<endl;		
          Double_t px1mup = part->GetPx();
          Double_t py1mup = part->GetPy();
          Double_t pz1mup = part->GetPz();
          Double_t e1mup = part->GetEnergy();
          vecmup1.SetPxPyPzE(px1mup,py1mup,pz1mup,e1mup);
          vecmuon.SetPxPyPzE(px1mup,py1mup,pz1mup,e1mup);
          if(isfast) cout<<"vectormup1-> PX,PY,PZ,EN = "<<px1mup<<", "<<py1mup<<", "<<pz1mup<<", "<<e1mup<<endl;
	}
	else{
	  if(isfast) cout<<"2nd  mu plus = "<<endl;
          Double_t px2mup = part->GetPx();
          Double_t py2mup = part->GetPy();
          Double_t pz2mup = part->GetPz();
          Double_t e2mup = part->GetEnergy();
          vecmup2.SetPxPyPzE(px2mup,py2mup,pz2mup,e2mup);
          vecmuon.SetPxPyPzE(px2mup,py2mup,pz2mup,e2mup);
          if(isfast) cout<<"vectormup2-> PX,PY,PZ,EN = "<<px2mup<<", "<<py2mup<<", "<<pz2mup<<", "<<e2mup<<endl;
	}
      } 

      //Lets check the Z0 produced
      if (idp==23 && pstat==21){
        nz++;
        if(nz % 2 != 0){
          if(isfast) cout<<"first Z0"<<endl;       	
          Double_t pxza = part->GetPx();
          Double_t pyza = part->GetPy();
          Double_t pzza = part->GetPz();
          Double_t eza = part->GetEnergy();
          vecza.SetPxPyPzE(pxza,pyza,pzza,eza);
          vecz0.SetPxPyPzE(pxza,pyza,pzza,eza);
          if(isfast) cout<<"4-vector Za-> PX,PY,PZ,EN = "<<pxza<<", "<<pyza<<", "<<pzza<<", "<<eza<<endl;
	  Double_t ya = sin(ip); //task to Sergio definir rapidity como funcion de E y p
	  hdNdyZa->Fill(ya);
	}
	else{ 
          if(isfast) cout<<"second Z0 "<<endl;
          Double_t pxzb = part->GetPx();
          Double_t pyzb = part->GetPy();
          Double_t pzzb = part->GetPz();
          Double_t ezb = part->GetEnergy();
          veczb.SetPxPyPzE(pxzb,pyzb,pzzb,ezb);
          vecz0.SetPxPyPzE(pxzb,pyzb,pzzb,ezb);
          if(isfast) cout<<"4-vector zb-> PX,PY,PZ,EN = "<<pxzb<<", "<<pyzb<<", "<<pzzb<<", "<<ezb<<endl;
        }
     }



     if (TMath::Abs(idp)==13){  // mu+ mu-
    	Double_t pxch=part->GetPx();
   	 Double_t pych=part->GetPy();
   	 Double_t pzch=part->GetPz();
         Double_t ech=part->GetEnergy();
   	 Double_t ptch=TMath::Sqrt(pxch*pxch+pych*pych);
   	 Double_t costhetach=pzch/(TMath::Sqrt(ptch*ptch+pzch*pzch));
    	ych= TMath::ATan(pzch/ech);      //momento/energia
    	etach=TMath::ATanH(costhetach);
   	 if(ptch!=0){ phichmv=TMath::ACos(pxch/ptch);
                  if(pych<0) phich=2*TMath::Pi()-phichmv;
                  else phich=phichmv;
                  }
	// llenando los histogramas dentro del loop de particulas
         hetach->Fill(etach);
         hych->Fill(ych);
         hphich->Fill(phich);
	 if(i==eventH) hyphiwpt->Fill(etach,phich,ptch); //Con Background
         // if(i==eventH && idpGrand2==25) hyphiwpt->Fill(etach,phich,ptch); //Solo Muones de Higgs
     }//end if of particle identified

    if (idp==23 && pstat==21){  // Z
    	 Double_t pxchZ=part->GetPx();
   	 Double_t pychZ=part->GetPy();
   	 Double_t pzchZ=part->GetPz();
         Double_t echZ=part->GetEnergy();
   	 Double_t ptchZ=TMath::Sqrt(pxchZ*pxchZ+pychZ*pychZ);
   	 Double_t costhetachZ=pzchZ/(TMath::Sqrt(ptchZ*ptchZ+pzchZ*pzchZ));
    	ychZ = TMath::ATan(pzchZ/echZ);      //momento/energia
    	etachZ=TMath::ATanH(costhetachZ);
   	 if(ptchZ!=0){ phichmvZ=TMath::ACos(pxchZ/ptchZ);
                  if(pychZ<0) phichZ=2*TMath::Pi()-phichmvZ;
                  else phichZ=phichmvZ;
                  }
	// llenando los histogramas dentro del loop de particulas
         hetachZ->Fill(etachZ);
         hychZ->Fill(ychZ);
         hphichZ->Fill(phichZ);
	 if(i==eventH) hyphiwptZ->Fill(etachZ,phichZ,vecz0.Mag());
         //hyphiwptZ->Fill(etachZ,phichZ);
     }//end if of particle identified

    if (idp==25 && pstat==21){ // H0
    	 Double_t pxchH=part->GetPx();
   	 Double_t pychH=part->GetPy();
   	 Double_t pzchH=part->GetPz();
         Double_t echH=part->GetEnergy();
   	 Double_t ptchH=TMath::Sqrt(pxchH*pxchH+pychH*pychH);
   	 Double_t costhetachH=pzchH/(TMath::Sqrt(ptchH*ptchH+pzchH*pzchH));
    	ychH= TMath::ATan(pzchH/echH);      //momento/energia
    	etachH=TMath::ATanH(costhetachH);
   	 if(ptchH!=0){ phichmvH=TMath::ACos(pxchH/ptchH);
                  if(pychH<0) phichH=2*TMath::Pi()-phichmvH;
                  else phichH=phichmvH;
                  }
	// llenando los histogramas dentro del loop de particulas
         hetachH->Fill(etachH);
         hychH->Fill(ychH);
         hphichH->Fill(phichH);
	 if(i==eventH) hyphiwptH->Fill(etachH,phichH,vechiggs.Mag());
          //hyphiwptH->Fill(etachH,phichH);
     }//end if of particle identified


    }// end particleloop

    nHiggs+=nh0;
    if(isfast) cout<<"..... End of event, Statistics"<<endl;
      if(isfast) cout<<"#######N mum="<<nmum<<endl;
      if(isfast) cout<<"#######N mup="<<nmup<<endl;
      if(isfast) cout<<"#######N Z0="<<nz<<endl;
      if(isfast) cout<<"#######N h0="<<nh0<<endl;

    // Sum of 4-vectors
    veczatomu = vecmum1 + vecmup1; //1st Z0 to mu-mu+
    veczbtomu = vecmum2 + vecmup2; //2nd Z0 to mu-mu+
    vecto4mu = veczatomu+veczbtomu;// 4muons vector
    vechtoz = vecza +veczb;        // 2Z0 vector 
    
    // Invariant mass 
    Double_t masaih=vechtoz.Mag();
    Double_t masaiza=veczatomu.Mag();
    Double_t masaizb=veczbtomu.Mag();
    Double_t masai4mu=vecto4mu.Mag();
     if(isfast){ cout<<"-------------------------------"<<endl;
     cout<<"Masaiza="<<masaiza<<endl;
     cout<<"Masaizb="<<masaizb<<endl;
     cout<<"Masaih="<<masaih<<endl;
     cout<<"Masai4mu="<<masai4mu<<endl;
     cout<<"-------------------------------"<<endl;
     }
     hmasainvZa->Fill(masaiza);
     hmasainvZb->Fill(masaizb);
     hmasainvH->Fill(masaih);
     hmasainvH4mu->Fill(masai4mu);

    tree->Fill();
  }//end loop of events
   if(isfast){
   cout<<"---------------------------------------"<<endl;
   cout<<"--  no of generated Higgs="<<nHiggs<<endl;
   cout<<"---------------------------------------"<<endl;
   }
  tree->Print();
  
  //Scaling histograms due to the binning size 
  hmasainvH->Scale(250./NbinsH);
  hmasainvZa->Scale(140./Nbinsz);
  hmasainvZb->Scale(140./Nbinsz);
  hmasainvH4mu->Scale(250.0/Nbins4mu);

  // some summary transverse momentum  plots:
  TH1D* hist = new TH1D(HISTNAME, "p_{T}  spectrum for  #mu^{+}",
                        1000, 0, 150);
  hist->SetXTitle("p_{T}[GeV/c]");
  hist->SetYTitle("dN/dp_{T}(GeV^{-1}c)");
  char expression[64];
  sprintf(expression,"sqrt(pow(%s.fPx,2)+pow(%s.fPy,2))>>%s",
          BRANCHNAME, BRANCHNAME, HISTNAME);
  char selection[64];
  sprintf(selection,"%s.fKF==%d", BRANCHNAME, PDGNUMBER); //selection of mu-
  tree->Draw(expression,selection);
  hist->Sumw2();

  TH1D* histz = new TH1D(HISTNAMEZ, "p_{T}  spectrum for  Z^{0}",
                        1000, 0, 100);
  histz->SetXTitle("p_{T}[GeV/c]");
  histz->SetYTitle("dN/dp_{T}(GeV^{-1}c)");
  char expressionz[64];
  sprintf(expressionz,"sqrt(pow(%s.fPx,2)+pow(%s.fPy,2))>>%s",
          BRANCHNAME, BRANCHNAME, HISTNAMEZ);
  char selectionz[64];
  sprintf(selectionz,"%s.fKF==%d", BRANCHNAME, 23); //selection of Z0
  tree->Draw(expressionz,selectionz);
  histz->Sumw2();

  file->Write();
  file->Close();

  return 0;
}

// ------ second function to show canvas 
int showEventSample()
{
  loadLibraries();

  TFile* file = TFile::Open(FILENAME, "READ");
  if (!file || !file->IsOpen()) {
    Error("showEventSample", "Couldn;t open file %s", FILENAME);
    return 1;
  }
  TTree* tree = (TTree*)file->Get(TREENAME);
  if (!tree) {
    Error("showEventSample", "couldn't get TTree %s", TREENAME);
    return 2;
  }
  TH1D* Hevents = (TH1D*)file->Get("Hevents");
  TH1D* hist = (TH1D*)file->Get(HISTNAME);
  TH1D* histz = (TH1D*)file->Get(HISTNAMEZ);
  TH1F* masaH = (TH1F*)file->Get("hinvmassH4mu");
  TH1F* masaH2 = (TH1F*)masaH->Clone("hinvmassH4mu");
  TH1F* masaZa = (TH1F*)file->Get("hinvmassZa");
  TH1F* masaZa2 = (TH1F*)masaZa->Clone("hinvmassZa2");
  TH1F* masaZb = (TH1F*)file->Get("hinvmassZb");
  TH1F* masaZb2 = (TH1F*)masaZb->Clone("hinvmassZb2");

TH2F* hyphiwpt2=(TH2F*)file->Get("hyphiwpt");
TH1F* hetach2=(TH1F*)file->Get("hetach");
TH1F* hhych2=(TH1F*)file->Get("hych");
TH1F* hphich2=(TH1F*)file->Get("hphich");

TH2F* hyphiwpt2Z=(TH2F*)file->Get("hyphiwptZ");
TH1F* hetach2Z=(TH1F*)file->Get("hetachZ");
TH1F* hhych2Z=(TH1F*)file->Get("hychZ");
TH1F* hphich2Z=(TH1F*)file->Get("hphichZ");

TH2F* hyphiwpt2H=(TH2F*)file->Get("hyphiwptH");
TH1F* hetach2H=(TH1F*)file->Get("hetachH");
TH1F* hhych2H=(TH1F*)file->Get("hychH");
TH1F* hphich2H=(TH1F*)file->Get("hphichH");

  Int_t nevents=Hevents->GetBinContent(1);
  Float_t deta2=14;

 if (!hist) {
    Error("showEventSample", "couldn't get TH1D %s", HISTNAME);
    return 4;
  }
// EL SIGUIENTE PEDAZO DE CODIGO HACE EL AJUSTE -----------------------------------------------------
//------Sergio Task (hacer mejor ajuste y encontrar chi2, lo mismo para las otras dist)
// Función para aproximar la señal con una gaussiana	
	TF1 *FGaus = new TF1("FGaus","[0]*exp(-0.5*((x-[1])/[2])^2)",70,110);
	FGaus->SetParameters(1300,91,1);
//	FGaus->Draw();

/*
//__________________________________________________________________________________
// Función para ajustar la PseudorapidityMu
	TF1 *FGausR = new TF1("FGausR","[0]*exp(-0.5*((x-[1])/[2])^2)",-10,10);
	FGausR->SetParameters(1,0,1);
//__________________________________________________________________________________
// Función para ajustar la PseudorapidityZ
	TF1 *FGausRZ = new TF1("FGausRZ", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",-8,8);
	FGausRZ->SetParameters(1,-3,1,1,3,1);
//__________________________________________________________________________________
// Función para ajustar la PseudorapidityH
	TF1 *FGausRH = new TF1("FGausRH", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",-8,8);
	FGausRH->SetParameters(1,-3,1,1,3,1);
*/

//__________________________________________________________________________________
// Función para ajustar la masa del Higgs y su background
	TF1 *FGausH = new TF1("FGausH", "[0]*exp(-0.5*((x-[1])/[2])^2)",122,128);
	FGausH->SetParameters(1,125,1);

	TF1 *FGausHbk = new TF1("FGausHbk", "[0]*exp([1]*x)", 115, 135);
	FGausHbk->SetParameters(400,-.01);
// --------------------------------------------------------------------
// Función para aproximar el background con una exponencial
	TF1 *FBack = new TF1("FBack","[0]*exp([1]*x)", 20, 140);
	FBack->SetParameters(400,-.001);
//	FBack->Draw();
// ----------------------------------------------------------------------


   masaZa->Fit("FGaus", "r");
   TF1 *fitGaus = masaZa->GetFunction("FGaus");
     Double_t chi2 = fitGaus->GetChisquare();
     Double_t ndf = fitGaus->GetNDF();
     cout << "Valor de Chi2= " << chi2 << ",  ndf= " << ndf << endl; 

  masaZa2->Fit("FBack", "r"); 
 TF1 *fitBack = masaZa2->GetFunction("FBack");       
     Double_t chi2bk = fitBack->GetChisquare();
     Double_t ndfbk = fitBack->GetNDF();
     cout << "Valor de Chi2bk= " << chi2bk << ",  ndfbk= " << ndfbk << endl; 


   masaZb->Fit("FGaus", "r");
   TF1 *fitGausb = masaZb->GetFunction("FGaus");
     Double_t chi2b = fitGausb->GetChisquare();
     Double_t ndfb = fitGausb->GetNDF();
     cout << "Valor de Chi2= " << chi2b << ",  ndf= " << ndfb << endl; 

  masaZb2->Fit("FBack", "r"); 
 TF1 *fitBackb = masaZb2->GetFunction("FBack");       
     Double_t chi2bkb = fitBack->GetChisquare();
     Double_t ndfbkb = fitBack->GetNDF();
     cout << "Valor de Chi2bk= " << chi2bkb << ",  ndfbk= " << ndfbkb << endl; 



/*
//_______________________________________________________________________________
// Ajuste de la Pseudorapidity
   hetach2->Fit("FGausR", "r");
   TF1 *fitGausR = hetach2->GetFunction("FGausR");
//_______________________________________________________________________________
// Ajuste de la PseudorapidityZ
   hetach2Z->Fit("FGausRZ", "r");
   TF1 *fitGausRZ = hetach2Z->GetFunction("FGausRZ");
//_______________________________________________________________________________
// Ajuste de la PseudorapidityZ
   hetach2H->Fit("FGausRH", "r");
   TF1 *fitGausRH = hetach2H->GetFunction("FGausRH");
//_______________________________________________________________________________
*/

// Ajuste de la masa del Higgs.
   masaH->Fit("FGausH", "r");
   TF1 *fitGausH = masaH->GetFunction("FGausH");
//_______________________________________________________________________________
// Ajuste del Background de la masa del Higgs.
   masaH2->Fit("FGausHbk", "r");
   TF1 *fitGausHbk = masaH2->GetFunction("FGausHbk");





//-------------------------------------------------------------------------------------------------

  // Draw the histogram in a canvas
  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->SetLogy();
  masaH->GetYaxis()->SetRangeUser(1,10000);
  masaH->GetXaxis()->SetRangeUser(100,150);
  fitGausHbk->SetLineColor(kBlack);    
  masaH->Draw("hist");
  fitGausH->Draw("sames");
  fitGausHbk->Draw("sames");
  
/*  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2");
//  canvas2->SetLogy();
   gStyle->SetOptFit(1111);
  masaZa->GetYaxis()->SetRangeUser(1,300);
  masaZa->GetXaxis()->SetRangeUser(20,140);
  masaZa->SetMarkerStyle(23);
  masaZb->SetMarkerStyle(25);
  masaZa->SetMarkerColor(kRed);
  masaZb->SetMarkerColor(kBlue);
   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->SetHeader("Invariant mass","C"); // option "C" allows to center the header
   legend->AddEntry(masaZa,"Z^{0}_{A}","lep");
    legend->AddEntry(masaZb,"Z^{0}_{B}","lep");
    masaZa->Draw();
    fitGaus->Draw("same");
    fitBack->Draw("same");
//    masaZb->Draw("sames");
  legend->Draw();
*/


// CANVAS DE PRUEBA --------------------------------------------------------------------------

  TCanvas* canvasPRUEBA = new TCanvas("canvasPRUEBA", "canvasPRUEBA");                     //
  gStyle->SetOptStat(11);                                                                  //
  gStyle->SetOptFit(111);                                                                  //
  masaZa->GetYaxis()->SetRangeUser(1,300);                                                 //
  masaZa->GetXaxis()->SetRangeUser(20,140);                                                //
  masaZa->SetMarkerStyle(25);                                                              //
  masaZa->SetMarkerColor(kBlue);                                                           //
  masaZa2->SetMarkerStyle(9);                                                              //
  masaZa2->SetMarkerColor(kBlue);                                                          //
  fitBack->SetLineColor(kBlack);                                                           //
                                                                                           //
  TLegend* legendPRUEBA = new TLegend(0.1,0.7,0.48,0.9);                                   //
  legendPRUEBA->SetHeader("Invariant mass","C"); // option "C" allows to center the header //
  legendPRUEBA->AddEntry(masaZa,"Z^{0}_{A}","lep");                                        //
  legendPRUEBA->AddEntry(fitGaus,"signal","lep");                                          //
  legendPRUEBA->AddEntry(fitBack,"background","lep");                                      //
 // legendPRUEBA->AddEntry(masaZa2,"Z^{0}_{a}","lep");                                     //                   
  masaZa->Draw();                                                                          //
 // masaZa2->Draw("same");                                                                 //
  fitBack->Draw("sames");                                                                  //
  fitGaus->Draw("sames");                                                                  //
  legendPRUEBA->Draw();                                                                    //


  TCanvas* canvasPRUEBAb = new TCanvas("canvasPRUEBAb", "canvasPRUEBAb");                     //
  gStyle->SetOptStat(11);                                                                  //
  gStyle->SetOptFit(111);                                                                  //
  masaZb->GetYaxis()->SetRangeUser(1,300);                                                 //
  masaZb->GetXaxis()->SetRangeUser(20,140);                                                //
  masaZb->SetMarkerStyle(25);                                                              //
  masaZb->SetMarkerColor(kBlue);                                                           //
  masaZb2->SetMarkerStyle(9);                                                              //
  masaZb2->SetMarkerColor(kBlue);                                                          //
  fitBackb->SetLineColor(kBlack);                                                           //
                                                                                           //
  TLegend* legendPRUEBAb = new TLegend(0.1,0.7,0.48,0.9);                                   //
  legendPRUEBAb->SetHeader("Invariant mass","C"); // option "C" allows to center the header //
  legendPRUEBAb->AddEntry(masaZb,"Z^{0}_{B}","lep");                                        //
  legendPRUEBAb->AddEntry(fitGaus,"signal","lep");                                          //
  legendPRUEBAb->AddEntry(fitBack,"background","lep");                                      //
 // legendPRUEBA->AddEntry(masaZb2,"Z^{0}_{a}","lep");                                     //                   
  masaZb->Draw();                                                                          //
 // masaZa2->Draw("same");                                                                 //
  fitBackb->Draw("sames");                                                                  //
  fitGaus->Draw("sames");                                                                  //
  legendPRUEBAb->Draw();                                                                    //


//luego en el show event sample   ya para dibujar el plot


 TCanvas* cetaphiwpt=new TCanvas("cetaphiwpt","canvas for y phi and p for muons");
 hyphiwpt2->GetXaxis()->SetTitle("#eta");
 hyphiwpt2->GetYaxis()->SetTitle("#phi (rads)");
 hyphiwpt2->GetZaxis()->SetTitle("|p| (GeV/c)");
 hyphiwpt2->Draw("lego2z");

 TCanvas* cetaphiwptZ=new TCanvas("cetaphiwptZ","canvas for y phi and p for Z^0");
 hyphiwpt2Z->GetXaxis()->SetTitle("#eta");
 hyphiwpt2Z->GetYaxis()->SetTitle("#phi (rads)");
// hyphiwpt2Z->GetZaxis()->SetTitle("|p| (GeV/c)");
 hyphiwpt2Z->GetZaxis()->SetTitle("Multiplicity");
 hyphiwpt2Z->Draw("lego2z");

 TCanvas* cetaphiwptH=new TCanvas("cetaphiwptH","canvas for y phi and p for h^0");
 hyphiwpt2H->GetXaxis()->SetTitle("#eta");
 hyphiwpt2H->GetYaxis()->SetTitle("#phi (rads)");
 //hyphiwpt2H->GetZaxis()->SetTitle("|p| (GeV/c)");
hyphiwpt2H->GetZaxis()->SetTitle("multiplicity");
 hyphiwpt2H->Draw("lego2z");

///* 
//hay que agregar una normalizacion al eta antes de graficar, hay que checar
  TCanvas* ceta=new TCanvas("ceta","canvas eta");
  gStyle->SetOptStat(0);
  //ceta->SetLogy();
  hetach2->GetYaxis()->SetTitle("#LT dN/d#eta #GT");
  hetach2->GetXaxis()->SetTitle("#eta");
  hetach2->GetYaxis()->SetRangeUser(4.6,9.0);
  hetach2->SetMarkerStyle(23);
  hetach2->SetMarkerColor(kBlue);
   TLegend* legeta = new TLegend(0.1,0.7,0.48,0.9);
   legeta->SetHeader("pseudorapidity density","C"); //"C" to center header
   legeta->AddEntry(hetach2,Form("Pythia6 gg->H0->ZZ->4Mu at #sqrt{s}=%d GeV", ECM),"lep");
  hetach2->Scale(1./(deta2*nevents));
  hetach2->Draw();
//  fitGausR->Draw("sames");
  legeta->Draw();

  TCanvas* cetaZ=new TCanvas("cetaZ","canvas eta");
  gStyle->SetOptStat(0);
  //ceta->SetLogy();
  hetach2Z->GetYaxis()->SetTitle("#LT dN/d#eta #GT");
  hetach2Z->GetXaxis()->SetTitle("#eta");
  hetach2Z->GetYaxis()->SetRangeUser(4.6,9.0);
  hetach2Z->SetMarkerStyle(23);
  hetach2Z->SetMarkerColor(kBlue);
   TLegend* legetaZ = new TLegend(0.1,0.7,0.48,0.9);
   legetaZ->SetHeader("pseudorapidity density","C"); //"C" to center header
   legetaZ->AddEntry(hetach2,Form("Pythia6 gg->H0->ZZ at #sqrt{s}=%d GeV", ECM),"lep");
  hetach2->Scale(1./(deta2*nevents));
  hetach2Z->Draw();
  legetaZ->Draw();

  TCanvas* cetaH=new TCanvas("cetaH","canvas eta");
  gStyle->SetOptStat(0);
  //ceta->SetLogy();
  hetach2H->GetYaxis()->SetTitle("#LT dN/d#eta #GT");
  hetach2H->GetXaxis()->SetTitle("#eta");
  hetach2H->GetYaxis()->SetRangeUser(4.6,9.0);
  hetach2H->SetMarkerStyle(23);
  hetach2H->SetMarkerColor(kBlue);
   TLegend* legetaH = new TLegend(0.1,0.7,0.48,0.9);
   legetaH->SetHeader("pseudorapidity density","C"); //"C" to center header
   legetaH->AddEntry(hetach2,Form("Pythia6 gg->H0 at #sqrt{s}=%d GeV", ECM),"lep");
  hetach2->Scale(1./(deta2*nevents));
  hetach2H->Draw();
  legetaH->Draw();


//*/








//--------------------------------------------------------------------------------------------


  TCanvas* canvas3 = new TCanvas("canvas3", "canvasPtforMu");
  canvas3->SetLogy();
  hist->GetYaxis()->SetRangeUser(1,1000000000);
  hist->Draw();

  TCanvas* canvas4 = new TCanvas("canvas4", "canvasPtforZ");
  canvas4->SetLogy();
  histz->GetYaxis()->SetRangeUser(1,1000000000);
  histz->Draw();

  //---Sergio task (crear mas canvas para dibujar Cosine1plot Cosine2plot Phiplot1 Phiplot2 Correlationplot)
  return 0;
}

void pythiaExample(Int_t n=10000) {
   makeEventSample(n);
   showEventSample();
}

#ifndef __CINT__
int main(int argc, char** argv)
{
  TApplication app("app", &argc, argv);
  Int_t n = 100;
  if (argc > 1)
    n = strtol(argv[1], NULL, 0);
  int retVal = 0;
  if (n > 0)
    retVal = makeEventSample(n);
  else {
    retVal = showEventSample();
    app.Run();
  }
  return retVal;
}
#endif

//____________________________________________________________________
// EOF

