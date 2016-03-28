#include <math.h>
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TNtuple.h"
using namespace std;

void res(const char* inFile = "retrack.root", const char* configureFile = "config.txt")
{
  ifstream config(configureFile);
  char COMMENT[256], RADIUS[11], PHI[9], THETA[11], PTOT[10];
  float radius, phi, theta, ptot;
  config.getline(COMMENT,256);
  config.getline(COMMENT,256);
  config.getline(COMMENT,256);
  config.getline(COMMENT,256);
  config.getline(COMMENT,256);
  config>>RADIUS>>radius>>PHI>>phi>>THETA>>theta>>PTOT>>ptot;
  config.close();

  gSystem->Load("/gpfs/mnt/gpfs02/phenix/plhf/plhf3/wenqing/test_fullmap/install/lib/libMyEvent"); // contains MyTrack, MyPair etc.
  gSystem->Load("libReconstruction"); // reconstruction class, contains findIntersection() and findMomentum()

  TFile* f1 = new TFile(inFile,"READ");
  if(!(f1))
  {
    cout<<"can't find input file..."<<endl;
    exit(1);
  }
  TTree* retrack;
  f1->GetObject("retrack",retrack);

  float alpha_e, phi_e, zed_e, px_e, py_e, pz_e;
  float alpha_p, phi_p, zed_p, px_p, py_p, pz_p;

  retrack->SetBranchAddress("alpha_e",&alpha_e);
  retrack->SetBranchAddress("phi_e",&phi_e);
  retrack->SetBranchAddress("zed_e",&zed_e);

  retrack->SetBranchAddress("px_e",&px_e);
  retrack->SetBranchAddress("py_e",&py_e);
  retrack->SetBranchAddress("pz_e",&pz_e);

  retrack->SetBranchAddress("alpha_p",&alpha_p);
  retrack->SetBranchAddress("phi_p",&phi_p);
  retrack->SetBranchAddress("zed_p",&zed_p);

  retrack->SetBranchAddress("px_p",&px_p);
  retrack->SetBranchAddress("py_p",&py_p);
  retrack->SetBranchAddress("pz_p",&pz_p);


  //==============================================================
  //     reconstructed position resolution (spherical coord)
  //==============================================================
  //_s stands for standard reconstruction
  //_r stands for new reconstruction using true conversion radius
  //_r stands for new reconstruction using reconstructed conversion radius
  TH1F* dradius_r = new TH1F("dradius_r","difference from true conversion radius; r_{conv}-r_{conv}^{true} (cm)",1000,-10,10);

  TH1F* dphi_s = new TH1F("dphi_s","difference from true #phi; #phi_{conv}-#phi_{conv}^{true} (rad)",1000,-0.5,0.5);
  TH1F* dphi_r = new TH1F("dphi_r","difference from true #phi; #phi_{conv}-#phi_{conv}^{true} (rad)",1000,-0.5,0.5);

  TH1F* dtheta_s = new TH1F("dtheta_s","difference from true #theta #theta{conv}-#theta{conv}^{true} (rad)",1000,-0.5,0.5);
  TH1F* dtheta_r = new TH1F("dtheta_r","difference from true #theta #theta{conv}-#theta{conv}^{true} (rad)",1000,-0.5,0.5);

  //======================================================
  //                 momentum resolution
  //======================================================
  //transverse plane
  TH1F* DpT_s = new TH1F("DpT_s","relative difference from true p_{T};p_{T}-p_{T}^{true}/p_{T}^{true}",1000,-1,1);
  TH1F* DpT_r = new TH1F("DpT_r","relative difference from true p_{T};p_{T}-p_{T}^{true}/p_{T}^{true}",1000,-1,1);
  TH1F* DpT_rr = new TH1F("DpT_rr","relative difference from true p_{z};p_{z}-p_{z}^{true}/p_{z}^{true}",1000,-1,1);

  //Z direction
  TH1F* Dpz_s = new TH1F("Dpz_s","relative difference from true p_{z};p_{z}-p_{z}^{true}/p_{z}^{true}",1000,-1,1);
  TH1F* Dpz_r = new TH1F("Dpz_r","relative difference from true p_{z};p_{z}-p_{z}^{true}/p_{z}^{true}",1000,-1,1);
  TH1F* Dpz_rr = new TH1F("Dpz_rr","relative difference from true p_{z};p_{z}-p_{z}^{true}/p_{z}^{true}",1000,-1,1);


  int npair = retrack->GetEntries();

  float pT_t, pT_s, pT_r, pT_rr; //_t stands for true value
  float pz_t, pz_s, pz_r, pz_rr; //_t stands for true value
  float radius_t, radius_r;
  float phi_t, phi_s, phi_r; // reconstructed position
  float theta_t, theta_s, theta_r;
  float zvertex = 0;

  for (int i = 0; i < npair; ++i)
  {   
    if(i%1000 == 0) cout << "Pair:  " << i << endl;
    retrack->GetEntry(i);

    //======================================================
    //            true info for the conversions
    //======================================================
    pT_t     = ptot/2*sin(theta); // each single track has half of the total mom of photon
    pz_t     = ptot/2*cos(theta);
    radius_t = radius;
    phi_t    = phi;
    theta_t  = theta;

    //======================================================
    //              standard reconstruction
    //======================================================
    pT_s    = TMath::Hypot(px_p, py_p);
    pz_s    = pz_p;
    phi_s   = TMath::ATan2(py_p, px_p);
    theta_s = TMath::ATan2(pT_s, pz_p);

    //======================================================
    //                new reconstruction
    //======================================================
    MyTrack trk_e, trk_p;

    trk_e.SetPhiDC(phi_e);
    trk_e.SetZDC(zed_e);
    trk_e.SetAlpha(alpha_e);
    trk_e.SetPx(px_e); // standard reconstruction
    trk_e.SetPy(py_e);
    trk_e.SetPz(pz_e);

    trk_p.SetPhiDC(phi_p);
    trk_p.SetZDC(zed_p);
    trk_p.SetAlpha(alpha_p);
    trk_p.SetPx(px_p);
    trk_p.SetPy(py_p);
    trk_p.SetPz(pz_p);

    MyPair pair;
    Reconstruction::findIntersection(&trk_e, &trk_p, &pair, zvertex);
    radius_r = pair.GetRPair();
    phi_r    = (pair.GetPhiElectron()+pair.GetPhiPositron())/2;
    theta_r  = (pair.GetThetaElectron()+pair.GetThetaPositron())/2;

    TVector3 mom_r  = Reconstruction::findMomentum(&trk_p, radius_t, phi_t, theta_t, zvertex);
    TVector3 mom_rr = Reconstruction::findMomentum(&trk_p, radius_r, phi_r, theta_r, zvertex);

    pT_r     = mom_r.Pt();
    pz_r     = mom_r.Pz();
    pT_rr    = mom_rr.Pt();
    pz_rr    = mom_rr.Pz();

    //======================================================
    //           conversion angle resollution
    //======================================================
    dradius_r->Fill(radius_r-radius_t);

    dphi_s->Fill(phi_s-phi_t);
    dphi_r->Fill(phi_r-phi_t);
    dtheta_s->Fill(theta_s-theta_t);
    dtheta_r->Fill(theta_r-theta_t);

    //======================================================
    //                 pT resolution
    //======================================================
    DpT_s->Fill((pT_s-pT_t)/pT_t);
    DpT_r->Fill((pT_r-pT_t)/pT_t);
    DpT_rr->Fill((pT_rr-pT_t)/pT_t);
    Dpz_s->Fill((pz_s-pz_t)/pz_t);
    Dpz_r->Fill((pz_r-pz_t)/pz_t);
    Dpz_rr->Fill((pz_rr-pz_t)/pz_t);

    // rvp->Fill((pT_rr-pT_t)/pT_t, radius_r-radius_t);
  } 
      
  TFile* f2 = new TFile("resolution.root","RECREATE");
  
  dradius_r->Write();
  dphi_s->Write();
  dphi_r->Write();
  dtheta_s->Write();
  dtheta_r->Write();

  DpT_s->Write();
  DpT_r->Write();
  DpT_rr->Write();
  Dpz_s->Write();
  Dpz_r->Write();
  Dpz_rr->Write();

  // rvp->Write();

  f2->Write();

  f1->Close();
  f2->Close();

  delete f1;
  delete f2;
}
