#include "Reconstruction.h"
#include "photon_conversion_analysis/MyEvent.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMatrixD.h"
#include <iostream>

using namespace std;

namespace PhotonConversionAnalysis
{

struct PointValue_3x3
{
	float x,y,z;// independent variables
	float dep[3];// dependent variables
	PointValue_3x3(float X, float Y, float Z, float A, float B, float C) : x(X), y(Y), z(Z)
	{
		dep[0] = A;
		dep[1] = B;
		dep[2] = C;
	}
	PointValue_3x3(){}
};


class PointVal3x3Sorted
{
  protected:
    float x_lo, x_hi;
    float y_lo, y_hi;
    float z_lo, z_hi;

    int level;
    int maxlevel;

    PointVal3x3Sorted* containers[2][2][2];

    public:
      PointVal3x3Sorted( float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV=0 ) : x_lo(X_LO), x_hi(X_HI), y_lo(Y_LO), y_hi(Y_HI), z_lo(Z_LO), z_hi(Z_HI), level(LEV), maxlevel(MLEV)
      {
        for(unsigned int i=0;i<2;++i){for(unsigned int j=0;j<2;++j){for(unsigned int k=0;k<2;++k){containers[i][j][k]=NULL;}}}
      }
      virtual ~PointVal3x3Sorted()
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
          	for(unsigned int k=0;k<2;++k){
              if(containers[i][j][k]!=NULL)
              {
                delete containers[i][j][k];
              }}}}
      }
      virtual bool insert( PointValue_3x3 const& point );
      virtual void append_list( vector<PointValue_3x3>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI )
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
          	for(unsigned int k=0;k<2;++k){
              if(containers[i][j][k]==NULL){continue;}
              if( (containers[i][j][k]->x_hi<X_LO) || (containers[i][j][k]->x_lo>X_HI) || (containers[i][j][k]->y_hi<Y_LO) || (containers[i][j][k]->y_lo>Y_HI) || (containers[i][j][k]->z_hi<Z_LO) || (containers[i][j][k]->z_lo>Z_HI) ){continue;}
              containers[i][j][k]->append_list( point_list, X_LO, X_HI, Y_LO, Y_HI, Z_LO, Z_HI );}}}
      }
};

class PointVal3x3SortedEnd : public PointVal3x3Sorted
{
  public:
    PointVal3x3SortedEnd( float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV=0 ) : PointVal3x3Sorted( X_LO,X_HI,Y_LO,Y_HI,Z_LO,Z_HI,MLEV,LEV ){}
    virtual ~PointVal3x3SortedEnd(){}
    virtual bool insert( PointValue_3x3 const& point )
    {
      points.push_back(point);
      return true;
    }
    virtual void append_list( vector<PointValue_3x3>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI )
    {
      for(unsigned int i=0;i<points.size();++i)
      {
        if( (points[i].x<X_LO) || (points[i].x>X_HI) || (points[i].y<Y_LO) || (points[i].y>Y_HI) || (points[i].z<Z_LO) || (points[i].z>Z_HI) ){continue;}
        point_list.push_back( points[i] );
      }
    }
  protected:
    vector<PointValue_3x3> points;
};

bool PointVal3x3Sorted::insert( PointValue_3x3 const& point )
{
  if( (point.x < x_lo) || (point.y < y_lo) || (point.z < z_lo) || (point.x > x_hi) || (point.y > y_hi) || (point.z > z_hi) )
  {
    return false;
  }

  int x_ind = 0;
  if(point.x > (x_lo + 0.5*(x_hi-x_lo))){x_ind=1;}
  int y_ind = 0;
  if(point.y > (y_lo + 0.5*(y_hi-y_lo))){y_ind=1;}
  int z_ind = 0;
  if(point.z > (z_lo + 0.5*(z_hi-z_lo))){z_ind=1;}

  if( containers[x_ind][y_ind][z_ind] == NULL )
  {
    float x_lo_new = x_lo + (float(x_ind))*0.5*(x_hi-x_lo);
    float x_hi_new = x_lo_new + 0.5*(x_hi-x_lo);

    float y_lo_new = y_lo + (float(y_ind))*0.5*(y_hi-y_lo);
    float y_hi_new = y_lo_new + 0.5*(y_hi-y_lo);

    float z_lo_new = z_lo + (float(z_ind))*0.5*(z_hi-z_lo);
    float z_hi_new = z_lo_new + 0.5*(z_hi-z_lo);

    if(level < maxlevel)
    {
      containers[x_ind][y_ind][z_ind] = new PointVal3x3Sorted( x_lo_new,x_hi_new, y_lo_new,y_hi_new, z_lo_new,z_hi_new, maxlevel, level+1 );
    }
    else
    {
      containers[x_ind][y_ind][z_ind] = new PointVal3x3SortedEnd( x_lo_new,x_hi_new, y_lo_new,y_hi_new, z_lo_new,z_hi_new, maxlevel, level+1 );
    }
  }
  return containers[x_ind][y_ind][z_ind]->insert( point );
}


class Reconstruction::impl
{
	public:
		impl( const char* lookup_file_name ) : alpha_r_z( 0., 0.7, 0., 150., 0., 100., 10 )
		{
			TFile f(lookup_file_name);
			TNtuple* t=0;
			f.GetObject("alpha_p_r_phi_theta_z", t);
			float alpha, p, r, phi, theta, z;
			t->SetBranchAddress("alpha",&alpha);
			t->SetBranchAddress("p",&p);
			t->SetBranchAddress("r",&r);
			t->SetBranchAddress("phi",&phi);
			t->SetBranchAddress("theta",&theta);
			t->SetBranchAddress("z",&z);
			for(int i=0;i<t->GetEntries();i+=1)
			{
				 t->GetEntry(i);
				 PointValue_3x3 point(TMath::Abs(alpha),r,TMath::Abs(z), phi,p,theta);
				 alpha_r_z.insert(point);
			}
			f.Close();
		}

		PointVal3x3Sorted alpha_r_z;
};

Reconstruction::Reconstruction( const char* lookup_file_name )
{
	pimpl = new impl(lookup_file_name);
}

Reconstruction::~Reconstruction()
{
	delete pimpl;
}




static const float max_r = 26.;
static const float lookup_alpha_delta = 0.01;
static const float lookup_r_delta = 2.0;
static const float lookup_z_delta = 4.0;

// 0 <--> phi
// 1 <--> p
// 2 <--> theta
static float lookup_fit( float alpha, float r, float z_DC_m_z_ver, int ind, PointVal3x3Sorted& alpha_r_z )
{
	vector<PointValue_3x3> points;
	alpha_r_z.append_list( points, alpha - lookup_alpha_delta, alpha + lookup_alpha_delta, r - lookup_r_delta, r + lookup_r_delta, z_DC_m_z_ver - lookup_z_delta, z_DC_m_z_ver + lookup_z_delta );

	if( points.size() < 10 )
	{
		cout<<"too few points : "<<points.size()<<endl;
		cout<<"alpha : "<<alpha - lookup_alpha_delta<<" "<<alpha + lookup_alpha_delta<<endl;
		cout<<"r : "<<r - lookup_r_delta<<" "<<r + lookup_r_delta<<endl;
		cout<<"z : "<<z_DC_m_z_ver - lookup_z_delta<<" "<<z_DC_m_z_ver + lookup_z_delta<<endl;
		return -9999.;
	}

	// p0 + p1*x + p2*y + p3*z + p4*x^2 + p5*xy + p6*xz + p7*y^2 + p8*yz + p9*z^2

	TMatrixD X( points.size(), 10 );
	TMatrixD y( points.size(), 1  );
	for(unsigned int i=0;i<points.size();i+=1)
	{
		y(i,0) = points[i].dep[ind];
		X(i,0) = 1.;
		X(i,1) = points[i].x;
		X(i,2) = points[i].y;
		X(i,3) = points[i].z;
		X(i,4) = points[i].x * points[i].x;
		X(i,5) = points[i].x * points[i].y;
		X(i,6) = points[i].x * points[i].z;
		X(i,7) = points[i].y * points[i].y;
		X(i,8) = points[i].y * points[i].z;
		X(i,9) = points[i].z * points[i].z;
	}

	TMatrixD Xt( 10, points.size() );
	Xt.Transpose(X);

	TMatrixD XtX = Xt * X;
	XtX.Invert();

	TMatrixD beta = XtX * (Xt * y);

	float retval = 0.;
	retval += beta(0,0);
	retval += beta(1,0)*alpha;
	retval += beta(2,0)*r;
	retval += beta(3,0)*z_DC_m_z_ver;
	retval += beta(4,0)*alpha*alpha;
	retval += beta(5,0)*alpha*r;
	retval += beta(6,0)*alpha*z_DC_m_z_ver;
	retval += beta(7,0)*r*r;
	retval += beta(8,0)*r*z_DC_m_z_ver;
	retval += beta(9,0)*z_DC_m_z_ver*z_DC_m_z_ver;

	return retval;
}


static float project_phi( float alpha, float phiDC, float r, float z_DC_m_z_ver, PointVal3x3Sorted& alpha_r_z )
{
	//DC coord sys: phi -0.5Pi~1.5Pi
	float absalpha, phirDC;
	absalpha=TMath::Abs(alpha);
	if ((r<0)||(r>max_r))
	{
		return TMath::Pi()/2;
	}

	if (alpha>0)//electron!
	{
		phirDC = phiDC + lookup_fit(absalpha, r, TMath::Abs(z_DC_m_z_ver), 0, alpha_r_z);
	}
	else//positron!
	{
		phirDC = phiDC - lookup_fit(absalpha, r, TMath::Abs(z_DC_m_z_ver), 0, alpha_r_z);
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


static float delta_phi( float alpha_e, float alpha_p, float phi_e, float phi_p, float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p, PointVal3x3Sorted& alpha_r_z )
{
	float fi_e, fi_p, dfi;

	fi_e = project_phi(alpha_e, phi_e, r, z_DC_m_z_ver_e, alpha_r_z);
	fi_p = project_phi(alpha_p, phi_p, r, z_DC_m_z_ver_p, alpha_r_z);
	dfi = fi_p-fi_e;

	return dfi;
}

static float find_intersection( float alpha_e, float alpha_p, float phi_e, float phi_p, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p, PointVal3x3Sorted& alpha_r_z )
{
	float tolerance = 0.001;
	float r_delta = 0.01;
	int max_iter = 20;

	float r = 10.;

	int iter = 0;
	while(iter<max_iter)
	{
		float f0 = delta_phi( alpha_e, alpha_p, phi_e, phi_p, r, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z );
		if( TMath::Abs(f0)<tolerance && iter>0 ){if(r<0){r=0.;}return r;}
		float f1 = delta_phi( alpha_e, alpha_p, phi_e, phi_p, r+r_delta, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z );
		float slope = (f1-f0)/r_delta;
		// f0 + slope*d = 0
		// d = -f0/slope;
		float d = -f0/slope;
		if( (r+d)<0. ){ r *= 0.5; }
		else{r += d;}

		iter += 1;

		if( r < 0. )
		{
			r = 0.;
			if( TMath::Abs(f0)<(3.*tolerance) ){return r;}
		}
	}
	if(iter >= max_iter){r = -r;}
	return r;
}


void Reconstruction::findIntersection(MyTrack const* trk1, MyTrack const* trk2, MyPair* pair, float zvertex)
{
	float r = find_intersection( trk1->GetAlpha(), trk2->GetAlpha(), trk1->GetPhiDC(), trk2->GetPhiDC(), trk1->GetZDC() - zvertex, trk2->GetZDC() - zvertex, pimpl->alpha_r_z );
	pair->SetRPair( r );

	pair->SetPhiElectron( project_phi( trk1->GetAlpha(), trk1->GetPhiDC(), fabs(r), TMath::Abs(trk1->GetZDC() - zvertex), pimpl->alpha_r_z ) );
	pair->SetPhiPositron( project_phi( trk2->GetAlpha(), trk2->GetPhiDC(), fabs(r), TMath::Abs(trk2->GetZDC() - zvertex), pimpl->alpha_r_z ) );

	float theta_e = lookup_fit( TMath::Abs(trk1->GetAlpha()), fabs(r), TMath::Abs(trk1->GetZDC() - zvertex), 2, pimpl->alpha_r_z );
	if( trk1->GetZDC() < zvertex ){theta_e = TMath::Pi() - theta_e; }
	pair->SetThetaElectron( theta_e );

	float theta_p = lookup_fit( TMath::Abs(trk2->GetAlpha()), fabs(r), TMath::Abs(trk2->GetZDC() - zvertex), 2, pimpl->alpha_r_z );
	if( trk2->GetZDC() < zvertex ){theta_p = TMath::Pi() - theta_p; }
	pair->SetThetaPositron( theta_p );
}

TVector3 Reconstruction::findMomentum(MyTrack* trk, float r, float phi_conv, float theta_conv, float zvertex)
{
	TVector3 vec;
	float mom = lookup_fit( TMath::Abs(trk->GetAlpha()), r, TMath::Abs(trk->GetZDC() - zvertex), 1, pimpl->alpha_r_z );
	float phi = project_phi( trk->GetAlpha(), trk->GetPhiDC(), r, TMath::Abs(trk->GetZDC() - zvertex), pimpl->alpha_r_z );
	float theta = lookup_fit( TMath::Abs(trk->GetAlpha()), r, TMath::Abs(trk->GetZDC() - zvertex), 2, pimpl->alpha_r_z );
	if( trk->GetZDC() < zvertex ){theta = TMath::Pi() - theta; }

	float pt = mom*sin(theta);
	vec.SetPtThetaPhi( pt, theta, phi );

	return vec;
}


}

