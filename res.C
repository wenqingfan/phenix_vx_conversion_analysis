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
#include "MyPion.h"
// using namespace std;
static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;
static const double c = 299792458;



struct PointVal2D
{
  float x,y;// independent variables
  float z;// dependent variable
  PointVal2D(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
  PointVal2D(){}
};

class PointVal2DSorted
{
  protected:
    float x_lo, x_hi;
    float y_lo, y_hi;

    int level;
    int maxlevel;

    PointVal2DSorted* containers[2][2];

    public:
      PointVal2DSorted( float X_LO, float X_HI, float Y_LO, float Y_HI, int MLEV, int LEV=0 ) : x_lo(X_LO), x_hi(X_HI), y_lo(Y_LO), y_hi(Y_HI), level(LEV), maxlevel(MLEV)
      {
        for(unsigned int i=0;i<2;++i){for(unsigned int j=0;j<2;++j){containers[i][j]=NULL;}}
      }
      virtual ~PointVal2DSorted()
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
              if(containers[i][j]!=NULL)
              {
                delete containers[i][j];
              }}}
      }
      virtual bool insert( PointVal2D const& point );
      virtual void append_list( vector<PointVal2D>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI )
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
              if(containers[i][j]==NULL){continue;}
              if( (containers[i][j]->x_hi<X_LO) || (containers[i][j]->x_lo>X_HI) || (containers[i][j]->y_hi<Y_LO) || (containers[i][j]->y_lo>Y_HI) ){continue;}
              containers[i][j]->append_list( point_list, X_LO, X_HI, Y_LO, Y_HI );}}
      }
};

class PointVal2DSortedEnd : public PointVal2DSorted
{
  public:
    PointVal2DSortedEnd( float X_LO, float X_HI, float Y_LO, float Y_HI, int MLEV, int LEV=0 ) : PointVal2DSorted( X_LO,X_HI,Y_LO,Y_HI,MLEV,LEV ){}
    virtual ~PointVal2DSortedEnd(){}
    virtual bool insert( PointVal2D const& point )
    {
      points.push_back(point);
      return true;
    }
    virtual void append_list( vector<PointVal2D>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI )
    {
      for(unsigned int i=0;i<points.size();++i)
      {

        if( (points[i].x<X_LO) || (points[i].x>X_HI) || (points[i].y<Y_LO) || (points[i].y>Y_HI) ){continue;}
        point_list.push_back( points[i] );
      }
    }
  protected:
    vector<PointVal2D> points;
};

bool PointVal2DSorted::insert( PointVal2D const& point )
{
  if( (point.x < x_lo) || (point.y < y_lo) || (point.x > x_hi) || (point.y > y_hi) )
  {
    return false;
  }

  int x_ind = 0;
  if(point.x > (x_lo + 0.5*(x_hi-x_lo))){x_ind=1;}
  int y_ind = 0;
  if(point.y > (y_lo + 0.5*(y_hi-y_lo))){y_ind=1;}

  if( containers[x_ind][y_ind] == NULL )
  {
    float x_lo_new = x_lo + (float(x_ind))*0.5*(x_hi-x_lo);
    float x_hi_new = x_lo_new + 0.5*(x_hi-x_lo);

    float y_lo_new = y_lo + (float(y_ind))*0.5*(y_hi-y_lo);
    float y_hi_new = y_lo_new + 0.5*(y_hi-y_lo);

    if(level < maxlevel)
    {
      containers[x_ind][y_ind] = new PointVal2DSorted( x_lo_new,x_hi_new, y_lo_new,y_hi_new, maxlevel, level+1 );
    }
    else
    {
      containers[x_ind][y_ind] = new PointVal2DSortedEnd( x_lo_new,x_hi_new, y_lo_new,y_hi_new, maxlevel, level+1 );
    }
  }
  return containers[x_ind][y_ind]->insert( point );
}


PointVal2DSorted alpha_r_phi( 0., 0.7, 0., 0.3, 10 );
PointVal2DSorted alpha_r_p( 0., 0.7, 0., 0.3, 10 );





float fit(float alpha, float r)//r in unit of cm
{
  float r_m = r/100.;
  vector<PointVal2D> points;
  alpha_r_phi.append_list( points, alpha - 0.01, alpha + 0.01, r_m - 0.06, r_m + 0.06 );

  // implementing local polynomial fitting 
  // notation the same as in https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)

  TMatrixD X( points.size(), 6 );
  TMatrixD y( points.size(), 1 );
  for(unsigned int i=0;i<points.size();i+=1)
  {
    y(i,0) = points[i].z;
    // p0 + p1*x + p2*y + p3*x^2 + p4*xy + p5*y^2
    X(i,0) = 1.;
    X(i,1) = points[i].x;
    X(i,2) = points[i].y;
    X(i,3) = points[i].x*points[i].x;
    X(i,4) = points[i].x*points[i].y;
    X(i,5) = points[i].y*points[i].y;
  }
  TMatrixD Xt( 6, points.size() );
  Xt.Transpose(X);

  TMatrixD XtX = Xt * X;
  XtX.Invert();

  TMatrixD beta = XtX * (Xt * y);

  return ( beta(0,0) + beta(1,0)*alpha + beta(2,0)*r_m + beta(3,0)*alpha*alpha + beta(4,0)*alpha*r_m + beta(5,0)*r_m*r_m );
}

float fitp(float alpha, float r)//r in unit of cm
{
    float r_m = r/100.;
    float absalpha;
    absalpha=TMath::Abs(alpha);
    vector<PointVal2D> points;
    alpha_r_p.append_list( points, absalpha - 0.01, absalpha + 0.01, r_m - 0.06, r_m + 0.06 );

    // implementing local polynomial fitting 
    // notation the same as in https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)

    TMatrixD X( points.size(), 6 );
    TMatrixD y( points.size(), 1 );
    for(unsigned int i=0;i<points.size();i+=1)
    {
      y(i,0) = points[i].z;
      // p0 + p1*x + p2*y + p3*x^2 + p4*xy + p5*y^2
      X(i,0) = 1.;
      X(i,1) = points[i].x;
      X(i,2) = points[i].y;
      X(i,3) = points[i].x*points[i].x;
      X(i,4) = points[i].x*points[i].y;
      X(i,5) = points[i].y*points[i].y;
    }
    TMatrixD Xt( 6, points.size() );
    Xt.Transpose(X);

    TMatrixD XtX = Xt * X;
    XtX.Invert();

    TMatrixD beta = XtX * (Xt * y);

    return ( beta(0,0) + beta(1,0)*absalpha + beta(2,0)*r_m + beta(3,0)*absalpha*absalpha + beta(4,0)*absalpha*r_m + beta(5,0)*r_m*r_m );
}


float phir(float alpha, float phiDC, float r)
{
  //DC coord sys: phi -0.5Pi~1.5Pi
    float absalpha, phirDC;
    absalpha=TMath::Abs(alpha);
  if ((r<0)||(r>26))
  {
      return TMath::Pi()/2;
  }

    if (alpha>0)//electron!
    {
        phirDC = phiDC + fit(absalpha, r);
    }
    else//positron!
    {
        phirDC = phiDC - fit(absalpha, r);
    }
    if (phirDC>1.5*TMath::Pi())
    {
        phirDC = phirDC-TMath::Pi();
    }
    if (phirDC<-0.5*TMath::Pi())
    {
      phirDC = phirDC+2*TMath::Pi();
    }
    return phirDC;
}

float Radius(float alpha_e, float alpha_p, float phi_e, float phi_p, float r)
{
    float fi_e, fi_p, dfi;

    fi_e = phir(alpha_e, phi_e, r);
    fi_p = phir(alpha_p, phi_p, r);
    dfi = fi_p-fi_e;

    return dfi;
}

float rootFind(float alpha_e, float alpha_p, float phi_e, float phi_p)
{
  // use Newton's method with finite differences
  float tolerance = 0.001;
  float r_delta = 0.1;
  int max_iter = 10;

  float r = 10.;

  int iter = 0;
  while(iter<max_iter)
  {
    float f0 = Radius( alpha_e, alpha_p, phi_e, phi_p, r );
    if( TMath::Abs(f0)<tolerance ){return r;}
    float f1 = Radius( alpha_e, alpha_p, phi_e, phi_p, r+r_delta );
    float slope = (f1-f0)/r_delta;
    // f0 + slope*d = 0
    // d = -f0/slope;
    float d = -f0/slope;
    if( (r+d)<0. ){ r *= 0.5; }
    else{r += d;}

    iter += 1;
  }

  return -99.;
}

float getPT(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;//pair corresponds to photon if the pair matches

  return pair.Pt();
}

float getPhiv(float px_e, float py_e, float pz_e, float px_p, float py_p, float pz_p)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px_e);
  p1.SetY(py_e);
  p1.SetZ(pz_e);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px_p);
  p2.SetY(py_p);
  p2.SetZ(pz_p);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;//pair corresponds to photon if the pair matches


  TVector3 P1, P2, Photon, z, v, u, w, wc;

  z.SetX(0);
  z.SetY(0);
  z.SetZ(1);//unit vector along z

  P1.SetX(px_e);
  P1.SetY(py_e);
  P1.SetZ(pz_e);

  P2.SetX(px_p);
  P2.SetY(py_p);
  P2.SetZ(pz_p);

  Photon.SetX(pair.Px());
  Photon.SetY(pair.Py());
  Photon.SetZ(pair.Pz());


  v = (P1.Cross(P2)).Unit();//unit vector v corresponds to the unit vector of the cross product of ep pair
  u = Photon.Unit();
  w = (u.Cross(v)).Unit();
  wc = (u.Cross(z)).Unit();

  float phiv =-9999;
  phiv = acos(-w.Dot(wc));


  return phiv;
}

float getMass(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;

 
  return pair.M();
}

float getPi0Mass(float px1, float py1, float pz1, float px2, float py2, float pz2, float px, float py, float pz)
{// px, py, pz correspond to emcal photon cluster 
  TLorentzVector p1, p2;
  TLorentzVector convPhoton, emcPhoton;
  TLorentzVector pi0;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  convPhoton = p1 + p2;

  emcPhoton.SetX(px);
  emcPhoton.SetY(py);
  emcPhoton.SetZ(pz);
  emcPhoton.SetE(emcPhoton.P());

  pi0 = convPhoton+emcPhoton;

  return pi0.M();
}

float getPi0Mass(float px1, float py1, float pz1, float px2, float py2, float pz2)
{// px, py, pz correspond to emcal photon cluster 
  TLorentzVector convPhoton, emcPhoton;
  TLorentzVector pi0;

  convPhoton.SetX(px1);
  convPhoton.SetY(py1);
  convPhoton.SetZ(pz1);
  convPhoton.SetE(convPhoton.P());

  emcPhoton.SetX(px2);
  emcPhoton.SetY(py2);
  emcPhoton.SetZ(pz2);
  emcPhoton.SetE(emcPhoton.P());

  pi0 = convPhoton+emcPhoton;

  return pi0.M();
}


float getEmcPhotonMass(float px, float py, float pz, float E)
{
  TLorentzVector photon;

  photon.SetX(px);
  photon.SetY(py);
  photon.SetZ(pz);
  photon.SetE(E);

  return photon.M();
}

float getConvertedPhotonE(float px1, float py1, float pz1, float px2, float py2, float pz2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector pair;

  p1.SetX(px1);
  p1.SetY(py1);
  p1.SetZ(pz1);
  p1.SetE(sqrt(pow(p1.P(),2) + Me2));

  p2.SetX(px2);
  p2.SetY(py2);
  p2.SetZ(pz2);
  p2.SetE(sqrt(pow(p2.P(),2) + Me2));

  pair = p1 + p2;

  return pair.E();
}

void fill_alpha_r()
{
  TFile f("alpha_p_r_phi.root");
  TNtuple* t=0;
  f.GetObject("alpha_p_r_phi", t);
  float alpha,p,r,phi;
  t->SetBranchAddress("alpha",&alpha);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("p",&p);
  t->SetBranchAddress("phi",&phi);
  for(int i=0;i<t->GetEntries();i+=1)
  {
    t->GetEntry(i);
    PointVal2D point(alpha,r,phi);
    alpha_r_phi.insert(point);
    PointVal2D point2(alpha,r,p);
    alpha_r_p.insert(point2);
  }
  f.Close();
}

void res(const char* inFile = "retrack.root")
{
  gSystem->Load("libMyPion");
  fill_alpha_r();
  TFile* f1 = new TFile(inFile,"READ");
  if(!(f1)){
        cout<<"can't find input file..."<<endl;
        exit(1);
    }
  TTree* retrack = (TTree*)f1->Get("retrack");
  TBranch* br = retrack->GetBranch("MyPion");
  MyPion* pion=0;
  br->SetAddress(&pion);

  TH1F* dconvptT = new TH1F("dconvptT","r resolution;r_{conv}-r_{true} (cm)",1000,-10,10);
  TH1F* dphiT = new TH1F("dphiT","phi resolution;#phi_{conv}-#phi_{true} (rad)",1000,-0.04,0.04);
  TH2F* drdphi  = new TH2F("drdphi","dr versus dphi in transverse plane;dr (cm);dphi (rad)",500,0,10,500,0,0.04);

  // TH1F* phoT_conv = new TH1F("phoT_conv","reconstructed pT_{#gamma};pT_{ee} (GeV)",1000,0,1);
  // TH1F* phoT_emc = new TH1F("phoT_emc","reconstructed pT_{#gamma};pT_{#gamma} (GeV)",1000,0,1);

  const int ptbin=11;
  //======================================================
  //                    p resolution
  //======================================================
  TH1F* dpT_r[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dpT_r[i] = new TH1F(Form("dpT_r[%d]", i),"difference of MCsingle pT_{e^{+}} and reconstructed pT_{e^{+}};dpT_{e^{+}} (GeV)",1000,-1,1);
  }
  TH1F* dphoT_r[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dphoT_r[i] = new TH1F(Form("dphoT_r[%d]", i),"difference of true pT_{#gamma} and reconstructed pT_{#gamma};dpT_{#gamma} (GeV)",1000,-1,1);
  }
  TH1F* dpT_s[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dpT_s[i] = new TH1F(Form("dpT_s[%d]", i),"difference of MCsingle pT_{e^{+}} and std reco pT_{e^{+}};dpT_{e^{+}} (GeV)",1000,-1,1);
  }
  TH1F* dphoT_s[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dphoT_s[i] = new TH1F(Form("dphoT_s[%d]", i),"difference of true pT_{#gamma} and std reco pT_{#gamma};dpT_{#gamma} (GeV)",1000,-1,1);
  }
  TH1F* dpz[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dpz[i] = new TH1F(Form("dpz[%d]", i),"difference of true pz_{e^{+}} and std reco pz_{e^{+}};pz_{e^{+}} (GeV)",1000,-1,1);
  }

  TH1F* DpT_r = new TH1F("DpT_r","difference of MCsingle pT_{#gamma} and reconstructed pT_{#gamma};dpT_{#gamma} ",1000,-1,1);
  TH1F* DphoT_r = new TH1F("DphoT_r","difference of True pT_{#gamma} and reconstructed pT_{#gamma};dpT_{#gamma} (GeV)",1000,-1,1);
  TH1F* DpT_s = new TH1F("DpT_s","difference of MCsingle pT_{e^{+}} and std reco pT_{e^{+}};dpT_{e^{+}}",1000,-1,1);
  TH1F* DphoT_s = new TH1F("DphoT_s","difference of True pT_{#gamma} and std reco pT_{#gamma};dpT_{#gamma}",1000,-1,1);
  TH1F* Dpz = new TH1F("Dpz","difference of True pz_{e^{+}} and std reco pz_{e^{+}};pz_{e^{+}}",1000,-1,1);

  //======================================================
  //                    E resolution
  //======================================================
  TH1F* dEmcE[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dEmcE[i] = new TH1F(Form("dEmcE[%d]", i),"difference of emcal photon E_{#gamma} and true E_{#gamma};dE_{#gamma} (GeV)",1000,-1,1);
  }
  TH1F* dConvE[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    dConvE[i] = new TH1F(Form("dConvE[%d]", i),"difference of converted photon E_{ee} and true E_{ee};dE_{ee} (GeV)",1000,-1,1);
  }

  TH1F* DEmcE = new TH1F("DEmcE","difference of emcal photon E_{#gamma} and true E_{#gamma};dE_{#gamma} (GeV)",1000,-1,1);
  TH1F* DConvE = new TH1F("DConvE","difference of converted photon E_{ee} and true E_{ee};dE_{ee} (GeV)",1000,-1,1);
  

  //======================================================
  //                Asymmetry Distribution
  //======================================================
  TH1F* trueAsy[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    trueAsy[i] = new TH1F(Form("trueAsy[%d]", i),"asymmetry distribution;asymmetry",100,-1,1);
  }
  TH1F* recoAsy[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    recoAsy[i] = new TH1F(Form("recoAsy[%d]", i),"asymmetry distribution;asymmetry",100,-1,1);
  }
  TH1F* TrueAsy = new TH1F("TrueAsy","asymmetry distribution;asymmetry",100,-1,1);
  TH1F* RecoAsy = new TH1F("RecoAsy","asymmetry distribution;asymmetry",100,-1,1);

  //======================================================
  //                Mass Distribution
  //======================================================
  TH1F* trueMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    trueMass[i] = new TH1F(Form("trueMass[%d]", i),"m_{ee} distribution;m_{ee} (GeV)",1000,0,0.04);
  }
  TH1F* recoMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    recoMass[i] = new TH1F(Form("recoMass[%d]", i),"m_{ee} distribution;m_{ee} (GeV)",1000,0,0.04);
  }
  TH1F* stdMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    stdMass[i] = new TH1F(Form("stdMass[%d]", i),"m_{ee} distribution (corrected using #phi_{conv});m_{ee} (GeV)",1000,0,0.04);
  }
  TH1F* recoPionMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    recoPionMass[i] = new TH1F(Form("recoPionMass[%d]", i),"m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }
  TH1F* emcPionMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    emcPionMass[i] = new TH1F(Form("emcPionMass[%d]", i),"m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }
  TH1F* stdPionMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    stdPionMass[i] = new TH1F(Form("stdPionMass[%d]", i),"m_{ee#gamma} distribution (corrected using #phi_{conv});m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }
  TH1F* convPionMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    convPionMass[i] = new TH1F(Form("convPionMass[%d]", i),"m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }
  TH1F* convPionMass_r[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    convPionMass_r[i] = new TH1F(Form("convPionMass_r[%d]", i),"m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }
  TH1F* truePionMass[ptbin] ={0};
  for(int i = 0; i < ptbin; i++)
  {
    truePionMass[i] = new TH1F(Form("truePionMass[%d]", i),"m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0.06,0.22);
  }

  TH1F* TrueMass = new TH1F("TrueMass","m_{ee} distribution;m_{#gamma} (GeV)",1000,0,0.04);
  TH1F* RecoMass = new TH1F("RecoMass","m_{ee} distribution;m_{#gamma} (GeV)",1000,0,0.04);
  TH1F* StdMass = new TH1F("StdMass","m_{ee} distribution;m_{#gamma} (GeV)",1000,0,0.04);

  // TH1F* McPionMass = new TH1F("McPionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);
  TH1F* RecoPionMass = new TH1F("RecoPionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);

  TH1F* EmcPionMass = new TH1F("EmcPionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);
  TH1F* StdPionMass = new TH1F("StdPionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);
  TH1F* ConvPionMass = new TH1F("ConvPionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);
  TH1F* ConvPionMass_r = new TH1F("ConvPionMass_r","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);

  TH1F* TruePionMass = new TH1F("TruePionMass","m_{ee#gamma} distribution;m_{ee#gamma} (GeV)",1000,0,0.5);

  int npion = retrack->GetEntries();
  int npair, nclust;
  float zvtx;

  int arm_p, arm_e;
  float alpha_p, phi_p, zed_p, alpha_e, phi_e, zed_e;//hitting point on DC
  float convptx, convpty, convptz, px0_p, py0_p, pz0_p, px0_e, py0_e, pz0_e, px0, py0, pz0;//true info
  float px_p, py_p, pz_p, px_e, py_e, pz_e;//standard reconstruction
  int emc_arm;
  float emc_ecore, x, y, z, emc_prob;
  float emc_px0, emc_py0, emc_pz0, emc_E0;
  float conv_px0, conv_py0, conv_pz0;

  float pT0_p, pT0_e;
  float pT_pr, pT_er, pT_prr, pT_err;
  float pT_ps, pT_es;
  float convptT, phiT, phoT;
  float convptT_r, phiT_r, phiT_rr;// reconstructed position
  float E1, E1_rr;
  
  float ZEDCUT=4;
  float PHIVCUT=0.1;

  int ibin=0;

  for (int i = 0; i < npion; ++i)
  {
    if(i%10 == 0) cout << "Event:  " << i << endl;
    pion->ClearPion();
    br->GetEntry(i);
    npair  = pion->GetNPair();
    nclust = pion->GetNEmcPhoton();
    zvtx   = pion->GetZVtx();

    if (npair>1) continue;// only keep one converted photon event

    for (int ipair = 0; ipair < npair; ++ipair)
    {
      arm_p   = (pion->GetEntry(ipair)).GetArm_p();
      arm_e   = (pion->GetEntry(ipair)).GetArm_e();

      alpha_p = (pion->GetEntry(ipair)).GetAlpha_p();
      phi_p   = (pion->GetEntry(ipair)).GetPhi_p();
      zed_p   = (pion->GetEntry(ipair)).GetZed_p();
      alpha_e = (pion->GetEntry(ipair)).GetAlpha_e();
      phi_e   = (pion->GetEntry(ipair)).GetPhi_e();
      zed_e   = (pion->GetEntry(ipair)).GetZed_e();

      convptx = (pion->GetEntry(ipair)).GetConvptx();
      convpty = (pion->GetEntry(ipair)).GetConvpty();
      convptz = (pion->GetEntry(ipair)).GetConvptz();
      px0_p   = (pion->GetEntry(ipair)).GetPx0_p();
      py0_p   = (pion->GetEntry(ipair)).GetPy0_p();
      pz0_p   = (pion->GetEntry(ipair)).GetPz0_p();
      px0_e   = (pion->GetEntry(ipair)).GetPx0_e();
      py0_e   = (pion->GetEntry(ipair)).GetPy0_e();
      pz0_e   = (pion->GetEntry(ipair)).GetPz0_e();
      px0     = (pion->GetEntry(ipair)).GetPx0();
      py0     = (pion->GetEntry(ipair)).GetPy0();
      pz0     = (pion->GetEntry(ipair)).GetPz0();

      px_p    = (pion->GetEntry(ipair)).GetPx_p();
      py_p    = (pion->GetEntry(ipair)).GetPy_p();
      pz_p    = (pion->GetEntry(ipair)).GetPz_p();
      px_e    = (pion->GetEntry(ipair)).GetPx_e();
      py_e    = (pion->GetEntry(ipair)).GetPy_e();
      pz_e    = (pion->GetEntry(ipair)).GetPz_e();

      //======================================================
      //             true info for conv photon
      //======================================================
      conv_px0 = px0;
      conv_py0 = py0;
      conv_pz0 = pz0;

      //======================================================
      //             true info for emc photon
      //======================================================
      emc_px0 = (pion->GetTrueEmcPhotonEntry(0)).GetGx();
      emc_py0 = (pion->GetTrueEmcPhotonEntry(0)).GetGy();
      emc_pz0 = (pion->GetTrueEmcPhotonEntry(0)).GetGz();
      emc_E0  = (pion->GetTrueEmcPhotonEntry(0)).GetGe();

      // cout<<"ith pion "<<i<<" emc "<<emc_px0<<" "<<emc_py0<<" "<<emc_pz0<<" "<<emc_E0<<" conv "<<conv_px0<<" "<<conv_py0<<" "<<conv_pz0<<" "<<sqrt(conv_px0*conv_px0+conv_py0*conv_py0+conv_pz0*conv_pz0)<<" pion "<<getPi0Mass(px0, py0, pz0, emc_px0, emc_py0, emc_pz0)<<endl;


      //======================================================
      //           MCsingle info for converted photon
      //======================================================
      convptT = TMath::Hypot(convptx, convpty);
      phiT    = atan2(convpty,convptx);
      pT0_p   = TMath::Hypot(px0_p,py0_p);
      pT0_e   = TMath::Hypot(px0_e,py0_e);
      phoT    = TMath::Hypot(px0,py0); // true photon pT

      if((phoT>=0.6)&&(phoT<0.8))  ibin = 0;
      if((phoT>=0.8)&&(phoT<1.0))  ibin = 1;
      if((phoT>=1.0)&&(phoT<1.2))  ibin = 2;
      if((phoT>=1.2)&&(phoT<1.4))  ibin = 3;
      if((phoT>=1.4)&&(phoT<1.6))  ibin = 4;
      if((phoT>=1.6)&&(phoT<1.8))  ibin = 5;
      if((phoT>=1.8)&&(phoT<2.0))  ibin = 6;
      if((phoT>=2.0)&&(phoT<2.5))  ibin = 7;
      if((phoT>=2.5)&&(phoT<3.0))  ibin = 8;
      if((phoT>=3.0)&&(phoT<3.5))  ibin = 9;
      if((phoT>=3.5)&&(phoT<5))    ibin = 10;

      //======================================================
      //        new reconstruction for converted photon
      //======================================================
      convptT_r = rootFind(alpha_e, alpha_p, phi_e, phi_p);// only transverse plane
      if((convptT_r < 0.) || (convptT_r > 30.)){continue;}
      phiT_r    = phir(alpha_e, phi_e, TMath::Hypot(convptx, convpty));
      phiT_rr   = phir(alpha_e, phi_e, convptT_r);

      pT_pr     = fitp(alpha_p, TMath::Hypot(convptx, convpty));
      pT_prr    = fitp(alpha_p, convptT_r);
      pT_er     = fitp(alpha_e, TMath::Hypot(convptx, convpty));
      pT_err    = fitp(alpha_e, convptT_r);

      dconvptT->Fill(convptT_r-convptT);
      dphiT->Fill(phiT_rr-phiT);
      drdphi->Fill(TMath::Abs(convptT_r-convptT),TMath::Abs(phiT_rr-phiT));

      //======================================================
      //      standard reconstruction for converted photon
      //======================================================
      pT_ps = TMath::Hypot(px_p, py_p);
      pT_es = TMath::Hypot(px_e, py_e);

      //======================================================
      //                 pT resolution
      //======================================================
      DpT_r->Fill((pT_prr-pT0_p)/pT0_p);
      DphoT_r->Fill((getPT(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e)-TMath::Hypot(px0,py0))/TMath::Hypot(px0,py0));

      DpT_s->Fill((pT_ps-pT0_p)/pT0_p);
      DphoT_s->Fill((getPT(px_p, py_p, pz_p, px_e, py_e, pz_e)-TMath::Hypot(px0,py0))/TMath::Hypot(px0,py0));

      Dpz->Fill((pz_p-pz0_p)/pz0_p);// single track

      dpT_r[ibin]->Fill((pT_prr-pT0_p)/pT0_p);
      dphoT_r[ibin]->Fill((getPT(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e)-TMath::Hypot(px0,py0))/TMath::Hypot(px0,py0));

      dpT_s[ibin]->Fill((pT_ps-pT0_p)/pT0_p);
      dphoT_s[ibin]->Fill((getPT(px_p, py_p, pz_p, px_e, py_e, pz_e)-TMath::Hypot(px0,py0))/TMath::Hypot(px0,py0));

      dpz[ibin]->Fill((pz_p-pz0_p)/pz0_p);// single track


      for (int iclust = 0; iclust < nclust; ++iclust)
      {// based on only one converted photon
        emc_arm   = (pion->GetEmcPhotonEntry(iclust)).GetArm();
        emc_ecore = (pion->GetEmcPhotonEntry(iclust)).GetEcore();
        emc_prob  = (pion->GetEmcPhotonEntry(iclust)).GetProb();

        // if (arm_e!=arm_p) continue;
        // if (emc_arm!=arm_p) continue;

        // if ((pion->GetEntry(0)).GetParent()!=7) continue;
        if (emc_ecore<0.4) continue;

        x = (pion->GetEmcPhotonEntry(iclust)).GetX();
        y = (pion->GetEmcPhotonEntry(iclust)).GetY();
        z = (pion->GetEmcPhotonEntry(iclust)).GetZ()-zvtx;


        //======================================================
        //              asymmetry distribution
        //======================================================
        E1    = sqrt(px0*px0+py0*py0+pz0*pz0);
        E1_rr = getConvertedPhotonE(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e);

        // TrueAsy->Fill(TMath::Abs(E1-emc_ecore)/(E1+emc_ecore));
        // RecoAsy->Fill(TMath::Abs(E1_rr-emc_ecore)/(E1_rr+emc_ecore));
        TrueAsy->Fill((E1-emc_E0)/(E1+emc_E0));
        RecoAsy->Fill((E1_rr-emc_ecore)/(E1_rr+emc_ecore));

        trueAsy[ibin]->Fill((E1-emc_E0)/(E1+emc_E0));
        recoAsy[ibin]->Fill((E1_rr-emc_ecore)/(E1_rr+emc_ecore));

        //======================================================
        //                   Energy resolution
        //======================================================
        DEmcE->Fill((emc_ecore-emc_E0)/emc_E0);
        DConvE->Fill((E1_rr-sqrt(px0*px0+py0*py0+pz0*pz0))/sqrt(px0*px0+py0*py0+pz0*pz0));

        dEmcE[ibin]->Fill((emc_ecore-emc_E0)/emc_E0);
        dConvE[ibin]->Fill((E1_rr-sqrt(px0*px0+py0*py0+pz0*pz0))/sqrt(px0*px0+py0*py0+pz0*pz0));

        if (TMath::Abs(emc_ecore-emc_E0)/emc_E0>0.2) continue;

        //======================================================
        //           converted photon mass resolution
        //======================================================
        TrueMass->Fill(getMass(px0_p, py0_p, pz0_p, px0_e, py0_e, pz0_e));
        StdMass->Fill(getMass(pT_ps*cos(phiT_rr), pT_ps*sin(phiT_rr), pz_p, pT_es*cos(phiT_rr), pT_es*sin(phiT_rr), pz_e));// corrected using phiT_rr
        RecoMass->Fill(getMass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e));

        trueMass[ibin]->Fill(getMass(px0_p, py0_p, pz0_p, px0_e, py0_e, pz0_e));
        stdMass[ibin]->Fill(getMass(pT_ps*cos(phiT_rr), pT_ps*sin(phiT_rr), pz_p, pT_es*cos(phiT_rr), pT_es*sin(phiT_rr), pz_e));// corrected using phiT_rr
        recoMass[ibin]->Fill(getMass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e));

        //======================================================
        //               pion mass resolution
        //======================================================
        // McPionMass->Fill(getPi0Mass(px0_p, py0_p, pz0_p, px0_e, py0_e, pz0_e, emc_ecore*x/sqrt(x*x+y*y+z*z), emc_ecore*y/sqrt(x*x+y*y+z*z), emc_ecore*z/sqrt(x*x+y*y+z*z)));
        RecoPionMass->Fill(getPi0Mass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e, emc_ecore*x/sqrt(x*x+y*y+z*z), emc_ecore*y/sqrt(x*x+y*y+z*z), emc_ecore*z/sqrt(x*x+y*y+z*z)));

        EmcPionMass->Fill(getPi0Mass(px0, py0, pz0, emc_ecore*x/sqrt(x*x+y*y+z*z), emc_ecore*y/sqrt(x*x+y*y+z*z), emc_ecore*z/sqrt(x*x+y*y+z*z)));
        StdPionMass->Fill(getPi0Mass(pT_ps*cos(phiT_rr), pT_ps*sin(phiT_rr), pz_p, pT_es*cos(phiT_rr), pT_es*sin(phiT_rr), pz_e, emc_px0, emc_py0, emc_pz0));// corrected using phiT_rr
        ConvPionMass->Fill(getPi0Mass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e, emc_px0, emc_py0, emc_pz0));
        ConvPionMass_r->Fill(getPi0Mass(pT_pr*cos(phiT_r), pT_pr*sin(phiT_r), pz_p, pT_er*cos(phiT_r), pT_er*sin(phiT_r), pz_e, emc_px0, emc_py0, emc_pz0));

        TruePionMass->Fill(getPi0Mass(px0, py0, pz0, emc_px0, emc_py0, emc_pz0));

        recoPionMass[ibin]->Fill(getPi0Mass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e, emc_ecore*x/sqrt(x*x+y*y+z*z), emc_ecore*y/sqrt(x*x+y*y+z*z), emc_ecore*z/sqrt(x*x+y*y+z*z)));

        emcPionMass[ibin]->Fill(getPi0Mass(px0, py0, pz0, emc_ecore*x/sqrt(x*x+y*y+z*z), emc_ecore*y/sqrt(x*x+y*y+z*z), emc_ecore*z/sqrt(x*x+y*y+z*z)));
        stdPionMass[ibin]->Fill(getPi0Mass(pT_ps*cos(phiT_rr), pT_ps*sin(phiT_rr), pz_p, pT_es*cos(phiT_rr), pT_es*sin(phiT_rr), pz_e, emc_px0, emc_py0, emc_pz0));// corrected using phiT_rr
        convPionMass[ibin]->Fill(getPi0Mass(pT_prr*cos(phiT_rr), pT_prr*sin(phiT_rr), pz_p, pT_err*cos(phiT_rr), pT_err*sin(phiT_rr), pz_e, emc_px0, emc_py0, emc_pz0));
        convPionMass_r[ibin]->Fill(getPi0Mass(pT_pr*cos(phiT_r), pT_pr*sin(phiT_r), pz_p, pT_er*cos(phiT_r), pT_er*sin(phiT_r), pz_e, emc_px0, emc_py0, emc_pz0));

        truePionMass[ibin]->Fill(getPi0Mass(px0, py0, pz0, emc_px0, emc_py0, emc_pz0));
      }
    }   
  }
  
      
  TFile* f2 = new TFile("mass.root","RECREATE");

  dconvptT->Write();
  dphiT->Write();

  drdphi->Write();

  DpT_r->Write();
  DphoT_r->Write();
  DpT_s->Write();
  DphoT_s->Write();
  Dpz->Write();

  for(int i = 0; i < ptbin; ++i) dpT_r[i]->Write(Form("dpT_r%d",i));
  for(int i = 0; i < ptbin; ++i) dphoT_r[i]->Write(Form("dphoT_r%d",i));
  for(int i = 0; i < ptbin; ++i) dpT_s[i]->Write(Form("dpT_s%d",i));
  for(int i = 0; i < ptbin; ++i) dphoT_s[i]->Write(Form("dphoT_s%d",i));
  for(int i = 0; i < ptbin; ++i) dpz[i]->Write(Form("dpz%d",i));

  DEmcE->Write();
  DConvE->Write();

  for(int i = 0; i < ptbin; ++i) dEmcE[i]->Write(Form("dEmcE%d",i));
  for(int i = 0; i < ptbin; ++i) dConvE[i]->Write(Form("dConvE%d",i));

  TrueAsy->Write();
  RecoAsy->Write();

  for(int i = 0; i < ptbin; ++i) trueAsy[i]->Write(Form("trueAsy%d",i));
  for(int i = 0; i < ptbin; ++i) recoAsy[i]->Write(Form("recoAsy%d",i));

  TrueMass->Write();
  StdMass->Write();
  RecoMass->Write();

  // McPionMass->Write();
  RecoPionMass->Write();
  EmcPionMass->Write();
  StdPionMass->Write();
  ConvPionMass->Write();
  ConvPionMass_r->Write();
  TruePionMass->Write();

  for(int i = 0; i < ptbin; ++i) trueMass[i]->Write(Form("trueMass%d",i));
  for(int i = 0; i < ptbin; ++i) stdMass[i]->Write(Form("stdMass%d",i));
  for(int i = 0; i < ptbin; ++i) recoMass[i]->Write(Form("recoMass%d",i));
  for(int i = 0; i < ptbin; ++i) recoPionMass[i]->Write(Form("recoPionMass%d",i));
  for(int i = 0; i < ptbin; ++i) emcPionMass[i]->Write(Form("emcPionMass%d",i));
  for(int i = 0; i < ptbin; ++i) stdPionMass[i]->Write(Form("stdPionMass%d",i));
  for(int i = 0; i < ptbin; ++i) convPionMass[i]->Write(Form("convPionMass%d",i));
  for(int i = 0; i < ptbin; ++i) convPionMass_r[i]->Write(Form("convPionMass_r%d",i));
  for(int i = 0; i < ptbin; ++i) truePionMass[i]->Write(Form("truePionMass%d",i));


  f2->Write();

  f1->Close();
  f2->Close();

  delete f1;
  delete f2;
}
