#ifndef __RECONSTRUCTION_H__
#define __RECONSTRUCTION_H__



#include "TVector3.h"

namespace PhotonConversionAnalysis
{

class MyTrack;
class MyPair;

class Reconstruction
{
public:
   Reconstruction( const char* lookup_file_name );
   ~Reconstruction();

   void findIntersection(MyTrack* trk1, MyTrack* trk2, MyPair* pair, float zvertex);

   TVector3 findMomentum(MyTrack* trk, float r_conv, float phi_conv, float theta_conv, float zvertex);

	private:
		class impl;
		impl* pimpl;
		Reconstruction();
		Reconstruction(Reconstruction const&);
 };

}

#endif /*__RECONSTRUCTION_H__*/

