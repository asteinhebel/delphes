
// 20 July 2020, Chris Potter. Usage:
// 1) place the this macro in your Delphes installation directory
// 2) modify inputFile, normalization and cme as appropriate
// 3) bash> root -b -q macro1.cc

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include<string.h>

void chainInputs(TChain& chain, const char inFile[]) {

  //read in txt file list of samples, chains for each sample and polarization
  //TChain chain("Delphes");
  FILE* filePointer;
  int bufferLength = 255;
  char buffer[bufferLength];
  filePointer = fopen(inFile, "r");
  while(fgets(buffer, bufferLength, filePointer)) {
    //puts(buffer);
    buffer[strlen(buffer)-1] = '\0';
    if(!gSystem->AccessPathName(buffer)){//if file exists
        chain.Add(buffer);
    } else {
        std::cout << buffer<<" DOES NOT EXIST - skipping" << std::endl;
    }

  }
  fclose(filePointer);
}

double getNorm(int i, bool realHinvBR) {
  //weighting for samples with single channel
  //
  //someday get rid of the hardcoding
  double norm=-99;
  double xsec=1.; //fb
  double luminosity=900.; //fb-1
  double br=1.; 

  if (i==0) {//hinv eLpR
    //file lumi = 1x100/ab
    if (realHinvBR) br=0.001;
    else br=0.1;
    xsec=313.;
    norm=xsec*luminosity*br*0.11/100000;//0.11=BR Z-lep, division for lumi
    //norm=xsec*luminosity*br/100000;//1=BR Z-lep (inclusive sample), division for lumi
  }
  else if (i==1) {//hinv eRpL
    //file lumi = 1x100/ab
    if (realHinvBR) br=0.001;
    else br=0.1;
    xsec=211.;
    norm=xsec*luminosity*br*0.11/100000;//0.11=BR Z-lep, division for lumi
    //norm=xsec*luminosity*br/100000;//1=BR Z-lep (inclusive sample), division for lumi
  }
  else if (i==2||i==3) {//2f 
    //file lumi = 1x1/ab
    norm=luminosity/1000.;
  }
  else if (i==4||i==5||i==6||i==7) { //3f 
    //file lumi = 1x10/ab
    norm=luminosity/10000;
  }
  else if (i==8) { //4f Wev eLpR
    xsec=10200.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==9) { //4f Wev eRpL
    xsec=109.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==10) { //4f WW eLpR
    xsec=37500.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==11) { //4f WW eRpL
    xsec=2580.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==12) { //4f Zee eLpR
    xsec=2510.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==13) { //4f Zee eRpL
    xsec=2630.;
    //norm=xsec*luminosity*br;
    norm=luminosity/1000.;
  }
  else if (i==14) { //4f Zvv eLpR
    //file lumi = 2x10/ab
    //xsec=354.;
    norm=xsec*luminosity*br/10000.;
  }
  else if (i==15) { //4f Zvv eRpL
    //file lumi = 2x10/ab
    //xsec=117.;
    norm=xsec*luminosity*br/10000.;
  }
  else if (i==16) { //4f ZZ eLpR
    //file lumi = 2x10/ab
    //xsec=1800.;
    norm=xsec*luminosity*br/10000.;
  }
  else if (i==17) { //4f ZZ eRpL
    //file lumi = 2x10/ab
    //xsec=827.;
    norm=xsec*luminosity*br/10000.;
  }
  else if (i==18||i==19) { //2f inclusive H
    //file lumi = 2x100/ab
    norm=luminosity/100000.;
  }

  //nevent = xs * int lumi * BR

  return norm;
}

double calculateSignificance (vector<double> vec, bool eLpR) {
  if (eLpR) return vec[0] / sqrt(vec[0]+vec[2]); 
  else return vec[1] / sqrt(vec[1]+vec[3]); 
}

void macro1() {

  bool DEBUG = false;
  bool realHinvBR=false;

  gSystem->Load("libDelphes");

  //setup output tree
  Int_t isEl = 0;
  Double_t m_z;
  Double_t m_rec;
  Double_t eventWeight;
  Double_t normalization = 1.;
  Double_t pt_lep0;
  Double_t pt_lep1;
  Int_t nmbf = -1; //0=signal, 1=SM H, 2=2f, 31=3f ea, 32=3f ap, 41= 4f WW, 42=4f Wev, 43=4f ZZ, 44=4f Zee, 45=4f Zvv
  Int_t iseLpR = 0;
  //AMANDA - store bkg process info
  auto f = TFile::Open("out.root","RECREATE");
  TTree *tree = new TTree("zlepTree","z(lep)");
  tree->Branch("isEl",&isEl,"isEl/I");
  tree->Branch("iseLpR",&iseLpR,"iseLpR/I");
  tree->Branch("nmbf",&nmbf,"nmbf/I");
  tree->Branch("m_z",&m_z,"m_z/D");
  tree->Branch("m_rec",&m_rec,"m_rec/D");
  tree->Branch("eventWeight",&eventWeight,"eventWeight/D");
  tree->Branch("normalization",&normalization,"normalization/D");
  tree->Branch("pt_lep0",&pt_lep0,"pt_lep0/D");
  tree->Branch("pt_lep1",&pt_lep1,"pt_lep1/D");

  //declare histograms
  TH1 *histRecoilMass=new TH1F("RecoilMass", "m_{recoil}", 200, 0.0, 200.0);
  TH1 *histZMass=new TH1F("ZMass", "m_{ee}", 200, 0.0, 200.0);
  //recoil mass in categories
  TH1 *r0= new TH1F("r0",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r1= new TH1F("r1",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r2= new TH1F("r2",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r31= new TH1F("r31",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r32= new TH1F("r32",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r41= new TH1F("r41",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r42= new TH1F("r42",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r43= new TH1F("r43",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r44= new TH1F("r44",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r45= new TH1F("r45",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r0RL= new TH1F("r0RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r1RL= new TH1F("r1RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r2RL= new TH1F("r2RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r31RL= new TH1F("r31RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r32RL= new TH1F("r32RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r41RL= new TH1F("r41RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r42RL= new TH1F("r42RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r43RL= new TH1F("r43RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r44RL= new TH1F("r44RL",";recoil mass [GeV];Entries",60,100.,160);
  TH1 *r45RL= new TH1F("r45RL",";recoil mass [GeV];Entries",60,100.,160);

  //for cutflow
  vector<double> all(6,0);
  vector<double> cut1(6,0);
  vector<double> cut2(6,0);
  vector<double> cut3(6,0);
  vector<double> cut4(6,0);
  vector<double> cut5(6,0);
  vector<double> cut6(6,0);

  const char * inputs[] = { "/home/Amanda/delphes_git/delphes/run/2f1h_inv_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/2f1h_inv_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/2f_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/2f_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/3f_ap_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/3f_ap_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/3f_ea_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/3f_ea_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/4f_Wev_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/4f_Wev_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/4f_WW_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/4f_WW_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/4f_Zee_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/4f_Zee_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/4f_Zvv_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/4f_Zvv_eRpL.txt", "/home/Amanda/delphes_git/delphes/run/4f_ZZ_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/4f_ZZ_eRpL.txt","/home/Amanda/delphes_git/delphes/run/2f1h_eLpR.txt","/home/Amanda/delphes_git/delphes/run/2f1h_eRpL.txt"};
  //const char * inputs[] = { "/home/Amanda/delphes_git/delphes/run/2f1h_inv_eLpR.txt", "/home/Amanda/delphes_git/delphes/run/2f1h_inv_eRpL.txt" };
  int size = sizeof inputs / sizeof inputs[0];

  Double_t cme=250.;
  TLorentzVector beamsP4(0,0,0,cme);

  for (int i=0; i<size; i++) {
    cout<<"Input: "<<inputs[i]<<endl;
    TChain chain("Delphes");
    chainInputs(chain,inputs[i]);
    ExRootTreeReader *treeReader=new ExRootTreeReader(&chain);
    TClonesArray *branchEvent=treeReader->UseBranch("Event");
    TClonesArray *branchParticle=treeReader->UseBranch("Particle");
    TClonesArray *branchElectron=treeReader->UseBranch("Electron");
    TClonesArray *branchMuon=treeReader->UseBranch("Muon");
    TClonesArray *branchMET=treeReader->UseBranch("MissingET");
    //TClonesArray *branchJet=treeReader->UseBranch("Jet");

    //get normalization
    normalization=getNorm(i,realHinvBR);
    
    //store polarization
    if (i%2==0) iseLpR=1;
    else iseLpR=0;

    for(Int_t entry=0; entry<treeReader->GetEntries(); ++entry) {
    //for(Int_t entry=0; entry<200; ++entry) {

      isEl=0; 
      if (entry%100000==0) cout<<"Entry "<<entry<<endl;
 
      treeReader->ReadEntry(entry);
      LHEFEvent* event=(LHEFEvent*) branchEvent->At(0);
      eventWeight=event->Weight;
      m_z=0.0;
      m_rec=0.0;
      Double_t met=-1;
      MissingET* MET;
      /*
      if(branchJet->GetEntries()>0) continue;
      if (i<2) cut1[i]+=normalization;//jet veto signal
      else cut1[(i%2)+2]+=normalization;//cumulate bkg
      */ 

      //truth process ID for inclusive SM H - DO NOT CONSIDER HINV EVENTS FROM THIS FILE
      //index of daughters in particle list
      Int_t d1_i=-1;
      Int_t d2_i=-1;
      Int_t d3_i=-1;//for Z decay
      //daughter PID
      int d1=0;
      int d2=0;
      int d3=0;//for Z decay
      if (i==18||i==19) {//if incluive SM H file
        GenParticle *particleH;
        d1_i=-1;
        d2_i=-1;
        d3_i=-1;
        d1=0;
        d2=0;
        d3=0;
        for (int j=0; j<branchParticle->GetEntries(); j++){
          particleH=(GenParticle*) branchParticle->At(j);
          if (DEBUG) cout<<"entry "<<entry<<", particle "<<j<<", PID "<<particleH->PID<<endl;
          if (particleH->PID==25){
            d1_i=particleH->D1;
            d2_i=particleH->D2;
            if (DEBUG) cout<<"found higgs: particle "<<j<<"\n daughter particles indices "<<d1_i<<" and "<<d2_i<<endl;
          }
          if (j==d1_i) d1=TMath::Abs(particleH->PID);
          if (j==d2_i) d2=TMath::Abs(particleH->PID);
          if (j==d3_i) d3=TMath::Abs(particleH->PID);
          if (d1>0&&d2>0 && d1!=12 && d1!=14 && d1!=16 && d1!=23) {//both daugheters found and identified and allowed for consideration
            if (DEBUG) cout<<"daughter particle PID "<<d1<<" and "<<d2<<"\n Consider this event in cutflow"<<endl;
            break;
          }
          else if (d1>0&&d2>0) {//both daugheters found and identified and still may be Hinv
            if (DEBUG) cout<<"daughter particle PID "<<d1<<" and "<<d2<<"\n Could be Hinv - look for one more daughter"<<endl;
            if (d3_i<0) d3_i=particleH->D1;
	    if (d3==12||d3==14||d3==16) {//identified Hinv event
	      if (DEBUG) cout<<"HINV EVENT (daughter pid "<<d3<<") - SKIPPING"<<endl; 
	      goto new_event;
	    }
	    else if (d3>0) break;//Other Z decay - allow for consideration in cutflow
          }
        }
      }

      //fill counter for all events
      if (i<2) all[i]+=normalization;//all events signal
      else all[(i%2)+2]+=normalization;//cumulate bkg

      MET=(MissingET*) branchMET->At(0);
      met=MET->MET;
      if(met<15) continue;//require MET>15 GeV
      if (i<2) cut1[i]+=normalization;//met signal
      else cut1[(i%2)+2]+=normalization;//cumulate bkg

      Electron *electron1, *electron2;
      if(branchElectron->GetEntries()==2) {
        if (i<2) cut2[i]+=normalization;//cut 2 leptons signal
        else cut2[(i%2)+2]+=normalization;//cumulate bkg
        electron1=(Electron*) branchElectron->At(0);
        electron2=(Electron*) branchElectron->At(1);
	if ((electron1->Charge != electron2->Charge) && electron1->PT>=10 && electron2->PT>=10) {//SFOS, pT>10 GeV
          if (DEBUG) std::cout<<"2 SFOS ele"<<std::endl;
          if (i<2) cut3[i]+=normalization;//cut SFOS signal
          else cut3[(i%2)+2]+=normalization;//cumulate bkg
          TLorentzVector ZP4=electron1->P4()+electron2->P4();
          if(ZP4.M()>=75 && ZP4.M()<=105) {
            if (DEBUG) cout << "Pass M_vis" << endl;
            if (i<2) cut4[i]+=normalization;//cut m_vis signal
            else cut4[(i%2)+2]+=normalization;//cumulate bkg
            if(TMath::Abs(ZP4.Pt())>=20 && TMath::Abs(ZP4.Pt())<=70) {
              if (DEBUG) cout << "Pass Pt_vis" << endl;
              if (i<2) cut5[i]+=normalization;//cut pt_vis signal
              else cut5[(i%2)+2]+=normalization;//cumulate bkg
              TLorentzVector recoilP4=beamsP4-ZP4;
	      if (recoilP4.M()>=110 && recoilP4.M()<=150) {
		if (DEBUG) cout<< "save"<<endl;
		if (DEBUG) cout<<"Higgs decayed to: "<<d1<<" and "<<d2<<endl;
                if (i<2) cut6[i]+=normalization;//cut recoil mass signal
                else cut6[(i%2)+2]+=normalization;//cumulate bkg
                isEl=1;
		//identify if signal or 2,3,4f bkg
		if (i<2) {
		  nmbf=0;
		  if (iseLpR) r0->Fill(recoilP4.M(),eventWeight*normalization);
		  else r0RL->Fill(recoilP4.M(),eventWeight*normalization);
		
		}
		else if (i<4) {
		  nmbf=2;
		  if (iseLpR) r2->Fill(recoilP4.M(),eventWeight*normalization);
		  else r2RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<6) {
		  nmbf=31;
		  if (iseLpR) r31->Fill(recoilP4.M(),eventWeight*normalization);
		  else r31RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<8) {
		  nmbf=32;
		  if (iseLpR) r32->Fill(recoilP4.M(),eventWeight*normalization);
		  else r32RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<10) {
		  nmbf=41;
		  if (iseLpR) r41->Fill(recoilP4.M(),eventWeight*normalization);
		  else r41RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<12) {
		  nmbf=42;
		  if (iseLpR) r42->Fill(recoilP4.M(),eventWeight*normalization);
		  else r42RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<14) {
		  nmbf=43;
		  if (iseLpR) r43->Fill(recoilP4.M(),eventWeight*normalization);
		  else r43RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<16) {
		  nmbf=44;
		  if (iseLpR) r44->Fill(recoilP4.M(),eventWeight*normalization);
		  else r44RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<18) {
		  nmbf=45;
		  if (iseLpR) r45->Fill(recoilP4.M(),eventWeight*normalization);
		  else r45RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else {
		  nmbf=1;
		  if (iseLpR) r1->Fill(recoilP4.M(),eventWeight*normalization);
		  else r1RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
                histZMass->Fill(ZP4.M(), eventWeight);
                histRecoilMass->Fill(recoilP4.M(), eventWeight);
                m_z=ZP4.M();
                m_rec=recoilP4.M();
		pt_lep0=electron1->PT;
		pt_lep1=electron1->PT;
                //tree->Fill();
	      }//end if recoil mass
            }//end if pt vis
          }//end if m vis
  	}//end if opp sign
      }//end electrons
 
      //if electron event accepted, skip muon check
      if (isEl==0) {    
 
      m_z=0.0;
      m_rec=0.0;
      Muon *muon1, *muon2;
      if(branchMuon->GetEntries()==2) {
        if (i<2) cut2[i]+=normalization;//cut 2 leptons signal
        else cut2[(i%2)+2]+=normalization;//cumulate bkg
        muon1=(Muon*) branchMuon->At(0);
        muon2=(Muon*) branchMuon->At(1);
        if ((muon1->Charge != muon2->Charge) && muon1->PT>=10 && muon2->PT>=10) { //SFOS
          if (DEBUG) std::cout<<"2 SFOS mu"<<std::endl;
          if (i<2) cut3[i]+=normalization;//cut SFOS signal
          else cut3[(i%2)+2]+=normalization;//cumulate bkg
          TLorentzVector ZP4=muon1->P4()+muon2->P4();
          if(ZP4.M()>=75 && ZP4.M()<=105) {
            if (DEBUG) cout << "Pass M_vis" << endl;
            if (i<2) cut4[i]+=normalization;//cut m_vis signal
            else cut4[(i%2)+2]+=normalization;//cumulate bkg
            if(TMath::Abs(ZP4.Pt())>20 && TMath::Abs(ZP4.Pt())<70) {
              if (DEBUG) cout << "Pass Pt_vis" << endl;
              if (i<2) cut5[i]+=normalization;//cut pt_vis signal
              else cut5[(i%2)+2]+=normalization;//cumulate bkg
              TLorentzVector recoilP4=beamsP4-ZP4;
	      if (recoilP4.M()>=110 && recoilP4.M()<=150) {
		if (DEBUG) cout<< "save"<<endl;
		if (DEBUG) cout<<"Higgs decayed to: "<<d1<<" and "<<d2<<endl;
                if (i<2) cut6[i]+=normalization;//cut recoil mass signal
                else cut6[(i%2)+2]+=normalization;//cumulate bkg
		//identify if signal or 2,3,4f bkg
		if (i<2) {
		  nmbf=0;
		  if (iseLpR) r0->Fill(recoilP4.M(),eventWeight*normalization);
		  else r0RL->Fill(recoilP4.M(),eventWeight*normalization);
		
		}
		else if (i<4) {
		  nmbf=2;
		  if (iseLpR) r2->Fill(recoilP4.M(),eventWeight*normalization);
		  else r2RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<6) {
		  nmbf=31;
		  if (iseLpR) r31->Fill(recoilP4.M(),eventWeight*normalization);
		  else r31RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<8) {
		  nmbf=32;
		  if (iseLpR) r32->Fill(recoilP4.M(),eventWeight*normalization);
		  else r32RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<10) {
		  nmbf=41;
		  if (iseLpR) r41->Fill(recoilP4.M(),eventWeight*normalization);
		  else r41RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<12) {
		  nmbf=42;
		  if (iseLpR) r42->Fill(recoilP4.M(),eventWeight*normalization);
		  else r42RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<14) {
		  nmbf=43;
		  if (iseLpR) r43->Fill(recoilP4.M(),eventWeight*normalization);
		  else r43RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<16) {
		  nmbf=44;
		  if (iseLpR) r44->Fill(recoilP4.M(),eventWeight*normalization);
		  else r44RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else if (i<18) {
		  nmbf=45;
		  if (iseLpR) r45->Fill(recoilP4.M(),eventWeight*normalization);
		  else r45RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
		else {
		  nmbf=1;
		  if (iseLpR) r1->Fill(recoilP4.M(),eventWeight*normalization);
		  else r1RL->Fill(recoilP4.M(),eventWeight*normalization);
		}
                histZMass->Fill(ZP4.M(), eventWeight);
                histRecoilMass->Fill(recoilP4.M(), eventWeight);
                m_z=ZP4.M();
                m_rec=recoilP4.M();
		pt_lep0=muon1->PT;
		pt_lep1=muon1->PT;
                //tree->Fill();
	      }//end if recoil mass
            }//end if pt vis
          }//end if m vis
        }//end SFOS
      }}//end muons (and if check for isEl)
      tree->Fill();
      new_event:;//skip to here if found Hinv event in 2f1h inclusive higgs sample
    }//end loop over events
    chain.Reset();

  }//end loops over input files

  if(normalization!=1.) {
    histZMass->Scale(normalization/histZMass->GetEntries());
    histRecoilMass->Scale(normalization/histZMass->GetEntries());
  }

  //print cutflow 
  cout<<"LUMINOSITY: 0.9/ab per polarization"<<endl;
  //cutflow saved: signal eLpR, signal eRpL, bkg eLpR, bkg eRpL, signif eLpR, signif eRpL
  all[4]=calculateSignificance(all,true);
  all[5]=calculateSignificance(all,false);
  cut1[4]=calculateSignificance(cut1,true);
  cut1[5]=calculateSignificance(cut1,false);
  cut2[4]=calculateSignificance(cut2,true);
  cut2[5]=calculateSignificance(cut2,false);
  cut3[4]=calculateSignificance(cut3,true);
  cut3[5]=calculateSignificance(cut3,false);
  cut4[4]=calculateSignificance(cut4,true);
  cut4[5]=calculateSignificance(cut4,false);
  cut5[4]=calculateSignificance(cut5,true);
  cut5[5]=calculateSignificance(cut5,false);
  cut6[4]=calculateSignificance(cut6,true);
  cut6[5]=calculateSignificance(cut6,false);
  cout<<"Input: \t\t\t"<<"Signal eLpR\t"<<"Bkgd eLpR\t"<<"Significance\t"<<"Signal eRpL\t"<<"Bkgd eRpL\t"<<"Significance"<<endl;
  cout<<"All events: \t\t"<<all[0]<<"\t\t"<<all[2]<<"\t\t"<<all[4]<<"\t\t"<<all[1]<<"\t\t"<<all[3]<<"\t\t"<<all[5]<<endl;
  cout<<"MET>15 GeV: \t\t"<<cut1[0]<<"\t\t"<<cut1[2]<<"\t\t"<<cut1[4]<<"\t\t"<<cut1[1]<<"\t\t"<<cut1[3]<<"\t\t"<<cut1[5]<<endl;
  cout<<"2 leptons: \t"<<cut2[0]<<"\t\t"<<cut2[2]<<"\t\t"<<cut2[4]<<"\t\t"<<cut2[1]<<"\t\t"<<cut2[3]<<"\t\t"<<cut2[5]<<endl;
  cout<<"SFOS leptons >10 GeV: \t"<<cut3[0]<<"\t\t"<<cut3[2]<<"\t\t"<<cut3[4]<<"\t\t"<<cut3[1]<<"\t\t"<<cut3[3]<<"\t\t"<<cut3[5]<<endl;
  cout<<"75<M_vis<105 GeV: \t"<<cut4[0]<<"\t\t"<<cut4[2]<<"\t\t"<<cut4[4]<<"\t\t"<<cut4[1]<<"\t\t"<<cut4[3]<<"\t\t"<<cut4[5]<<endl;
  cout<<"20<pt_vis<70 GeV: \t"<<cut5[0]<<"\t\t"<<cut5[2]<<"\t\t"<<cut5[4]<<"\t\t"<<cut5[1]<<"\t\t"<<cut5[3]<<"\t\t"<<cut5[5]<<endl;
  cout<<"110<m_recoil<150 GeV: \t"<<cut6[0]<<"\t\t"<<cut6[2]<<"\t\t"<<cut6[4]<<"\t\t"<<cut6[1]<<"\t\t"<<cut6[3]<<"\t\t"<<cut6[5]<<endl;
  
  //print plots
  /*
  TCanvas* c0=new TCanvas("c0", "c0", 500,300);
  c0->cd(1);
  c0->Divide(2,1);

  c0->cd(1);
  histZMass->Draw("HIST");
  c0->cd(2);
  histRecoilMass->Draw("HIST");

  //c0->SaveAs("macro.eps");
  c0->SaveAs("macro_hinvTest.png");
  */

  f->Write();
  f->Close();
}

