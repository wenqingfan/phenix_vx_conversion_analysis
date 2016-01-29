#ifndef __MYPION_H__
#define __MYPION_H__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyPion                                                               //
//                                                                      //
// Description of the single pion event                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"


class TDirectory;

class MyPair : public TObject
{

private:
   int        parent;

   int        arm_p;
   int        arm_e;

   float      alpha_p;
   float      phi_p;
   float      zed_p;
   float      alpha_e;
   float      phi_e;
   float      zed_e;

   float      convptx;
   float      convpty;
   float      convptz;
   float      px0_p;
   float      py0_p;
   float      pz0_p;
   float      px0_e;
   float      py0_e;
   float      pz0_e;
   float      px0;// photon
   float      py0;
   float      pz0;

   float      px_p;
   float      py_p;
   float      pz_p;
   float      px_e;
   float      py_e;
   float      pz_e;


public:
   MyPair(){  };
   virtual ~MyPair() {  };

   int      GetParent() const { return parent; };

   int      GetArm_p() const { return arm_p; };
   int      GetArm_e() const { return arm_e; };

   float    GetAlpha_p() const { return alpha_p; };
   float    GetPhi_p() const { return phi_p; };
   float    GetZed_p() const { return zed_p; };
   float    GetAlpha_e() const { return alpha_e; };
   float    GetPhi_e() const { return phi_e; };
   float    GetZed_e() const { return zed_e; };

   float    GetConvptx() const { return convptx; };
   float    GetConvpty() const { return convpty; };
   float    GetConvptz() const { return convptz; };
   float    GetPx0_p() const { return px0_p; };
   float    GetPy0_p() const { return py0_p; };
   float    GetPz0_p() const { return pz0_p; };
   float    GetPx0_e() const { return px0_e; };
   float    GetPy0_e() const { return py0_e; };
   float    GetPz0_e() const { return pz0_e; };
   float    GetPx0() const { return px0; };
   float    GetPy0() const { return py0; };
   float    GetPz0() const { return pz0; };

   float    GetPx_p() const { return px_p; };
   float    GetPy_p() const { return py_p; };
   float    GetPz_p() const { return pz_p; };
   float    GetPx_e() const { return px_e; };
   float    GetPy_e() const { return py_e; };
   float    GetPz_e() const { return pz_e; };

   void     SetParent(int sparent) { parent = sparent; };

   void     SetArm_p(int sarm_p) { arm_p = sarm_p; };
   void     SetArm_e(int sarm_e) { arm_e = sarm_e; };

   void    SetAlpha_p(float salpha_p) { alpha_p = salpha_p; };
   void    SetPhi_p(float sphi_p) { phi_p = sphi_p; };
   void    SetZed_p(float szed_p) { zed_p = szed_p; };
   void    SetAlpha_e(float salpha_e) { alpha_e = salpha_e; };
   void    SetPhi_e(float sphi_e) { phi_e = sphi_e; };
   void    SetZed_e(float szed_e) { zed_e = szed_e; };

   void    SetConvptx(float sconvptx) { convptx = sconvptx; };
   void    SetConvpty(float sconvpty) { convpty = sconvpty; };
   void    SetConvptz(float sconvptz) { convptz = sconvptz; };
   void    SetPx0_p(float spx0_p) { px0_p = spx0_p; };
   void    SetPy0_p(float spy0_p) { py0_p = spy0_p; };
   void    SetPz0_p(float spz0_p) { pz0_p = spz0_p; };
   void    SetPx0_e(float spx0_e) { px0_e = spx0_e; };
   void    SetPy0_e(float spy0_e) { py0_e = spy0_e; };
   void    SetPz0_e(float spz0_e) { pz0_e = spz0_e; };
   void    SetPx0(float spx0) { px0 = spx0; };
   void    SetPy0(float spy0) { py0 = spy0; };
   void    SetPz0(float spz0) { pz0 = spz0; };

   void    SetPx_p(float spx_p) { px_p = spx_p; };
   void    SetPy_p(float spy_p) { py_p = spy_p; };
   void    SetPz_p(float spz_p) { pz_p = spz_p; };
   void    SetPx_e(float spx_e) { px_e = spx_e; };
   void    SetPy_e(float spy_e) { py_e = spy_e; };
   void    SetPz_e(float spz_e) { pz_e = spz_e; };


   ClassDef(MyPair,1)  // A track segment
};

class MyEmcPhoton : public TObject
{

private:
   int        arm;
   float      x;           // xyz position of EmcPhoton [cm]
   float      y;           
   float      z;           
   float      ecore;

   float      prob;

public:
   MyEmcPhoton(){  };
   virtual ~MyEmcPhoton() {  };

   int      GetArm() const { return arm; };
   float    GetX() const { return x; };
   float    GetY() const { return y; };
   float    GetZ() const { return z; };
   float    GetEcore() const { return ecore; };
   float    GetProb() const { return prob; };


   void     SetArm(int sarm) {arm = sarm;};
   void     SetX(float sx) { x = sx; };
   void     SetY(float sy) { y = sy; };
   void     SetZ(float sz) { z = sz; };
   void     SetEcore(float secore) { ecore = secore; };
   void     SetProb(float sprob) { prob = sprob; };


   ClassDef(MyEmcPhoton,1)  // A track segment
};

class MyTrueEmcPhoton : public TObject
{

private:

   float      gx;          // true info of emc photon
   float      gy;
   float      gz;
   float      ge;


public:
   MyTrueEmcPhoton(){  };
   virtual ~MyTrueEmcPhoton() {  };

   float    GetGx() const { return gx; };
   float    GetGy() const { return gy; };
   float    GetGz() const { return gz; };
   float    GetGe() const { return ge; };

   void     SetGx(float sgx) { gx = sgx; };
   void     SetGy(float sgy) { gy = sgy; };
   void     SetGz(float sgz) { gz = sgz; };
   void     SetGe(float sge) { ge = sge; };


   ClassDef(MyTrueEmcPhoton,1)  // A track segment
};

class MyPion : public TObject
{

private:
   float        zvtx;

   std::vector<MyPair> fPairsList;            //->array with all Pairs
   std::vector<MyEmcPhoton> fEmcPhotonsList;            //->array with emcCluster
   std::vector<MyTrueEmcPhoton> fTrueEmcPhotonsList;


public:
   MyPion();
   virtual   ~MyPion();

   void       ClearPion();

   void       SetZVtx(float szvtx) { zvtx = szvtx; }
   float      GetZVtx() { return zvtx; }

   void       AddPair(int sparent, int sarm_p, int sarm_e, float salpha_p, float sphi_p, float szed_p, float salpha_e, float sphi_e, float szed_e, float convptx, float convpty, float convptz, float spx0_p, float spy0_p, float spz0_p, float spx0_e, float spy0_e, float spz0_e, float spx0, float spy0, float spz0, float spx_p, float spy_p, float spz_p, float spx_e, float spy_e, float spz_e);

   void       AddEmcPhoton(int sarm, float sx, float sy, float sz, float secore, float sprob);
   void       AddTrueEmcPhoton(float sgx, float sgy, float sgz, float sge);

   int        GetNPair() { return fPairsList.size(); };
   int        GetNEmcPhoton() { return fEmcPhotonsList.size(); };
   int        GetNTrueEmcPhoton() { return fTrueEmcPhotonsList.size(); };

   MyPair&   GetEntry(int i) { return fPairsList[i]; };
   MyEmcPhoton& GetEmcPhotonEntry(int i) { return fEmcPhotonsList[i]; };
   MyTrueEmcPhoton& GetTrueEmcPhotonEntry(int i) { return fTrueEmcPhotonsList[i]; };

   std::vector<MyPair> GetPairs() { return fPairsList; };
   std::vector<MyEmcPhoton> GetEmcPhotons() { return fEmcPhotonsList; };
   std::vector<MyTrueEmcPhoton> GetTrueEmcPhotons() { return fTrueEmcPhotonsList; };

   ClassDef(MyPion,1)  //MyPion structure
};


#endif /*__MYPION_H__*/

