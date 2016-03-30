#ifndef __MYEVENT_H__
#define __MYEVENT_H__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyEvent                                                              //
//                                                                      //
// Description of the event and track parameters                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"


class TDirectory;

namespace PhotonConversionAnalysis
{
   class MyTrack : public TObject
   {

   private:
      int        arm;

      float      px;           // X component of the momentum
      float      py;           // Y component of the momentum
      float      pz;           // Z component of the momentum
      float      phi_DC;
      float      z_DC;
      float      alpha;
      int        charge;

      int        id;         // cluster id
      float      ecore;         // eid cut
      float      disp;
      int        n0;
      float      chi2;
      float      npe0;
      float      dep;
      float      prob;
      float      emcdz;
      float      emcdphi;

      float      emcx;
      float      emcy;
      float      emcz;

      float      crkphi;
      float      crkz;
   public:
      MyTrack()
      {  
         arm = -999;
         px = -999;
         py = -999;
         pz = -999;
         phi_DC = -999;
         z_DC = -999;
         alpha = -999;
         charge = -999;
         id = -999; 
         ecore = -999;
         disp = -999; 
         n0 = -999;
         chi2 = -999; 
         npe0 = -999; 
         dep = -999;
         prob = -999; 
         emcdz = -999;
         emcdphi = -999;
         emcx = -999;
         emcy = -999;
         emcz = -999;
         crkphi = -999;
         crkz = -999;
      };
      virtual ~MyTrack() {  };

      int      GetArm() const { return arm; };
      float    GetPx() const { return px; };
      float    GetPy() const { return py; };
      float    GetPz() const { return pz; };
      float    GetPt() const { return TMath::Sqrt(px*px + py*py); };
      float    GetPhiDC() const { return phi_DC; };
      float    GetZDC() const { return z_DC; };
      float    GetAlpha() const { return alpha; };
      int      GetCharge() const { return charge; };

      int      GetID() const { return id; };
      float    GetEcore() const { return ecore; };
      float    GetDisp() const { return disp; };
      int      GetN0() const { return n0; };
      float    GetChi2() const { return chi2; };
      float    GetNpe0() const { return npe0; };
      float    GetDep() const { return dep; };
      float    GetProb() const { return prob; };
      float    GetEmcdz() const { return emcdz; };
      float    GetEmcdphi() const { return emcdphi; };

      float    GetEmcx() const { return emcx; };
      float    GetEmcy() const { return emcy; };
      float    GetEmcz() const { return emcz; };
    
      float    GetCrkphi() const { return crkphi; };
      float    GetCrkz() const { return crkz; };

      void     SetArm(int sarm) { arm = sarm; };
      void     SetPx(float spx) { px = spx; };
      void     SetPy(float spy) { py = spy; };
      void     SetPz(float spz) { pz = spz; };
      void     SetPhiDC(float sphi_DC) { phi_DC = sphi_DC; };
      void     SetZDC(float sz_DC) { z_DC = sz_DC; };
      void     SetAlpha(float salpha) { alpha = salpha; };
      void     SetCharge(int scharge) { charge = scharge; };

      void     SetID(int sid) { id = sid; };
      void     SetEcore(float secore) { ecore = secore; };
      void     SetDisp(float sdisp) { disp = sdisp; };
      void     SetN0(int sn0) { n0 = sn0; };
      void     SetChi2(float schi2) { chi2 = schi2; };
      void     SetNpe0(float snpe0) { npe0 = snpe0; };
      void     SetDep(float sdep) { dep = sdep; };
      void     SetProb(float sprob) { prob = sprob; };
      void     SetEmcdz(float semcdz) { emcdz = semcdz; };
      void     SetEmcdphi(float semcdphi) { emcdphi = semcdphi; }

      void     SetEmcx(float semcx) { emcx = semcx; };
      void     SetEmcy(float semcy) { emcy = semcy; };
      void     SetEmcz(float semcz) { emcz = semcz; };

      void     SetCrkphi(float scrkphi) { crkphi = scrkphi; };
      void     SetCrkz(float scrkz) { crkz = scrkz; };

      ClassDef(MyTrack,1)  // A track segment
   };

   class MyPair
   {
   private:
      float phi_e;
      float phi_p;
      float theta_e;
      float theta_p;
      float r_pair;

      int id_e;
      int id_p;

   public:
      MyPair()
      {
         phi_e = -999;
         phi_p = -999;
         theta_e = -999;
         theta_p = -999;
         r_pair = -999;
      };
      virtual ~MyPair() {  };

      float GetPhiElectron() const { return phi_e; }
      float GetPhiPositron() const { return phi_p; }
      float GetThetaElectron() const { return theta_e; }
      float GetThetaPositron() const { return theta_p; }
      float GetRPair() const { return r_pair; }
      float GetIDElectron() const { return id_e; }
      float GetIDPositron() const { return id_p; }

      void SetPhiElectron(float sphi_e) { phi_e = sphi_e; }
      void SetPhiPositron(float sphi_p) { phi_p = sphi_p; }
      void SetThetaElectron(float stheta_e) { theta_e = stheta_e; }
      void SetThetaPositron(float stheta_p) { theta_p = stheta_p; }
      void SetRPair(float sr_pair) { r_pair = sr_pair; }
      void SetIDElectron(int sid_e) { id_e = sid_e; }
      void SetIDPositron(int sid_p ) { id_p = sid_p; }

      ClassDef(MyPair,1)  // A e+e- pair
   };

   class MyCluster : public TObject
   {

   private:
      int        arm;
      int        id;
      float      x;           // xyz position of cluster [cm]
      float      y;           
      float      z;           
      float      ecore;

      float      prob;
      float      emctrkdz;
      float      emctrkdphi;

   public:
      MyCluster()
      {  
        arm = -999;
        id = -999;
        x = -999;
        y = -999;
        z = -999;
        ecore = -999; 
        prob = -999;
        emctrkdz = -999; 
        emctrkdphi = -999;
      };
      virtual ~MyCluster() {  };

      int      GetArm() const { return arm; };
      int      GetID() const { return id; };
      float    GetX() const { return x; };
      float    GetY() const { return y; };
      float    GetZ() const { return z; };
      float    GetEcore() const { return ecore; };
      float    GetProb() const { return prob; };
      float    GetEmcdz() const { return emctrkdz; };
      float    GetEmcdphi() const { return emctrkdphi; };


      void     SetArm(int sarm) { arm = sarm; }
      void     SetID(int sid) { id = sid; }
      void     SetX(float sx) { x = sx; };
      void     SetY(float sy) { y = sy; };
      void     SetZ(float sz) { z = sz; };
      void     SetEcore(float secore) { ecore = secore; };
      void     SetProb(float sprob) { prob = sprob; };
      void     SetEmcdz(float semctrkdz) { emctrkdz = semctrkdz; };
      void     SetEmcdphi(float semctrkdphi) { emctrkdphi = semctrkdphi; };


      ClassDef(MyCluster,1)  // A cluster segment
   };

   class MyEvent : public TObject
   {

   private:
      float        BBCcharge;
      float        zvertex;
      float        centrality;

      std::vector<PhotonConversionAnalysis::MyTrack> TrackList; 
      std::vector<PhotonConversionAnalysis::MyPair> PairList;
      std::vector<PhotonConversionAnalysis::MyCluster> ClusterList;

   public:
      MyEvent()
      {
         BBCcharge  = -999;
         zvertex    = -999;
         centrality = -999;
      };
      virtual ~MyEvent();

      void       ClearEvent();
      
      void       SetBBCcharge(float sBBCcharge) { BBCcharge = sBBCcharge; };
      float      GetBBCcharge(){ return BBCcharge; };

      void       SetVtxZ(float szvertex) { zvertex = szvertex; }
      float      GetVtxZ() { return zvertex; }

      void       SetCentrality(float scentrality) { centrality = scentrality; };
      float      GetCentrality(){ return centrality; };

      void       AddTrack(MyTrack newTrack) { TrackList.push_back(newTrack); };
      void       AddPair(MyPair newPair) { PairList.push_back(newPair); };
      void       AddCluster(MyCluster newCluster) { ClusterList.push_back(newCluster); };

      int        GetNtrack() { return TrackList.size(); };
      int        GetNcluster() { return ClusterList.size(); };

      MyTrack&   GetEntry(int i) { return TrackList[i]; };
      MyCluster& GetClusterEntry(int i) { return ClusterList[i]; };

      std::vector<PhotonConversionAnalysis::MyTrack> GetTracks() { return TrackList; };
      std::vector<PhotonConversionAnalysis::MyCluster> GetClusters() { return ClusterList; };

      ClassDef(MyEvent,1)  //MyEvent structure
   };
}

#endif /*__MYEVENT_H__*/

