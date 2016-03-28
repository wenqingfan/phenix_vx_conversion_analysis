#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "MyEvent.h"


ClassImp(MyEvent)
ClassImp(MyTrack)
ClassImp(MyCluster)

// MyTrack newTrack;
// MyCluster newCluster;
   
MyEvent::MyEvent()
{
  BBCcharge = 0;
  zvertex = 0;
  centrality = 0;
}

////////////////////////////////////////////////////////////////////////////////

MyEvent::~MyEvent(){}

////////////////////////////////////////////////////////////////////////////////

void MyEvent::ClearEvent()
{
   TrackList.clear();
   ClusterList.clear();
}


/////////////////////////////////////////////////////////////////////////////////

// void MyEvent::AddTrack(int arm, float px, float py, float pz, float phi_DC, float z_DC, float alpha, int charge, int sid, float secore, float sdisp, int sn0, float schi2, float snpe0, float sdep, float sprob, float semcdz, float semcdphi, float semcx, float semcy, float semcz, float scrkphi, float scrkz)
// {
   
//    newTrack.SetArm(arm);
//    newTrack.SetPx(px);
//    newTrack.SetPy(py);
//    newTrack.SetPz(pz);
//    newTrack.SetPhiDC(phi_DC);
//    newTrack.SetZDC(z_DC);
//    newTrack.SetAlpha(alpha);
//    newTrack.SetCharge(charge);

//    newTrack.SetID(sid);
//    newTrack.SetEcore(secore);
//    newTrack.SetDisp(sdisp);
//    newTrack.SetN0(sn0);
//    newTrack.SetChi2(schi2);
//    newTrack.SetNpe0(snpe0);
//    newTrack.SetDep(sdep);
//    newTrack.SetProb(sprob);
//    newTrack.SetEmcdz(semcdz);
//    newTrack.SetEmcdphi(semcdphi);
//    newTrack.SetCrkphi(scrkphi);
//    newTrack.SetCrkz(scrkz);

//    TrackList.push_back(newTrack);   
// }

// void MyEvent::AddCluster(int arm, int id, float x, float y, float z, float secore, float sprob, float semctrkdz, float semctrkdphi)
// {
// 	newCluster.SetArm(arm);
//    newCluster.SetID(id);
//    newCluster.SetX(x);
//    newCluster.SetY(y);
//    newCluster.SetZ(z);
//    newCluster.SetEcore(secore);
//    newCluster.SetProb(sprob);
//    newCluster.SetEmcdz(semctrkdz);
//    newCluster.SetEmcdphi(semctrkdphi);

//    ClusterList.push_back(newCluster);   
// }
