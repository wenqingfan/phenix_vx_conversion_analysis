#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "MyPion.h"


ClassImp(MyPion)
ClassImp(MyPair)
ClassImp(MyEmcPhoton)
ClassImp(MyTrueEmcPhoton)


MyPion::MyPion()
{
  zvtx = 0;
}

////////////////////////////////////////////////////////////////////////////////

MyPion::~MyPion(){}

////////////////////////////////////////////////////////////////////////////////

void MyPion::ClearPion()
{
   fPairsList.clear();
   fEmcPhotonsList.clear();
   fTrueEmcPhotonsList.clear();
}


/////////////////////////////////////////////////////////////////////////////////

void MyPion::AddPair(int sparent, int sarm_p, int sarm_e, float salpha_p, float sphi_p, float szed_p, float salpha_e, float sphi_e, float szed_e, float sconvptx, float sconvpty, float sconvptz, float spx0_p, float spy0_p, float spz0_p, float spx0_e, float spy0_e, float spz0_e, float spx0, float spy0, float spz0, float spx_p, float spy_p, float spz_p, float spx_e, float spy_e, float spz_e)
{
   MyPair newPair;

   newPair.SetParent(sparent);

   newPair.SetArm_p(sarm_p);
   newPair.SetArm_e(sarm_e);

   newPair.SetAlpha_p(salpha_p);
   newPair.SetPhi_p(sphi_p);
   newPair.SetZed_p(szed_p);
   newPair.SetAlpha_e(salpha_e);
   newPair.SetPhi_e(sphi_e);
   newPair.SetZed_e(szed_e);

   newPair.SetConvptx(sconvptx);
   newPair.SetConvpty(sconvpty);
   newPair.SetConvptz(sconvptz);
   newPair.SetPx0_p(spx0_p);
   newPair.SetPy0_p(spy0_p);
   newPair.SetPz0_p(spz0_p);
   newPair.SetPx0_e(spx0_e);
   newPair.SetPy0_e(spy0_e);
   newPair.SetPz0_e(spz0_e);
   newPair.SetPx0(spx0);
   newPair.SetPy0(spy0);
   newPair.SetPz0(spz0);

   newPair.SetPx_p(spx_p);
   newPair.SetPy_p(spy_p);
   newPair.SetPz_p(spz_p);
   newPair.SetPx_e(spx_e);
   newPair.SetPy_e(spy_e);
   newPair.SetPz_e(spz_e);

   fPairsList.push_back(newPair);   
}

void MyPion::AddEmcPhoton(int sarm, float sx, float sy, float sz, float secore, float sprob)
{
   MyEmcPhoton newEmcPhoton;
   
   newEmcPhoton.SetArm(sarm);
   newEmcPhoton.SetX(sx);
   newEmcPhoton.SetY(sy);
   newEmcPhoton.SetZ(sz);
   newEmcPhoton.SetEcore(secore);
   newEmcPhoton.SetProb(sprob);

   fEmcPhotonsList.push_back(newEmcPhoton);   
}

void MyPion::AddTrueEmcPhoton(float sgx, float sgy, float sgz, float sge)
{
   MyTrueEmcPhoton newTrueEmcPhoton;

   newTrueEmcPhoton.SetGx(sgx);
   newTrueEmcPhoton.SetGy(sgy);
   newTrueEmcPhoton.SetGz(sgz);
   newTrueEmcPhoton.SetGe(sge);

   fTrueEmcPhotonsList.push_back(newTrueEmcPhoton);   
}
