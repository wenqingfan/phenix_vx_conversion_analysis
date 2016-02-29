#include <map>
#include <vector>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

struct DC_Event
{
	DC_Event() {}

	vector<float> x_global;
	vector<float> y_global;
	vector<int> plane;

	float p_tot;
	float p_phi;

	float r;

	float get_alpha()
	{
		float p_angle;
		if( TMath::Abs( y_global.back() - y_global.front() ) <= TMath::Abs( x_global.back() - x_global.front() ) )
		{
			TGraph graph(plane.size());
			for(unsigned int i=0;i<plane.size();i+=1)
			{
				graph.SetPoint( i, x_global[i], y_global[i] );
			}
			TF1 f1( "f1", "1 ++ x ++ x*x" );
			graph.Fit(&f1,"NQ");

			TGraph graph2(plane.size());
			for(unsigned int i=0;i<plane.size();i+=1)
			{
				graph2.SetPoint( i, sqrt(x_global[i]*x_global[i] + y_global[i]*y_global[i]), x_global[i] );
			}
			TF1 f2( "f2", "1 ++ x ++ x*x" );
			graph2.Fit(&f2,"NQ");
			double eval_x = f2.Eval(220.);

			p_angle = atan2( f1.Derivative(eval_x) * ( x_global.back() - x_global.front() ),  x_global.back() - x_global.front()  );
		}
		else
		{
			TGraph graph(plane.size());
			for(unsigned int i=0;i<plane.size();i+=1)
			{
				graph.SetPoint( i, y_global[i], x_global[i] );
			}
			TF1 f1( "f1", "1 ++ x ++ x*x" );
			graph.Fit(&f1,"NQ");

			TGraph graph2(plane.size());
			for(unsigned int i=0;i<plane.size();i+=1)
			{
				graph2.SetPoint( i, sqrt(x_global[i]*x_global[i] + y_global[i]*y_global[i]), y_global[i] );
			}
			TF1 f2( "f2", "1 ++ x ++ x*x" );
			graph2.Fit(&f2,"NQ");
			double eval_y = f2.Eval(220.);

			p_angle = atan2(  y_global.back() - y_global.front(), f1.Derivative(eval_y) * ( y_global.back() - y_global.front() )  );
		}

		TGraph graph(plane.size());
		for(unsigned int i=0;i<plane.size();i+=1)
		{
			graph.SetPoint( i, sqrt(y_global[i]*y_global[i] + x_global[i]*x_global[i]) , atan2( y_global[i] , x_global[i] ) );
		}
		TF1 f1( "f1", "1 ++ x ++ x*x" );
		graph.Fit(&f1,"NQ");
		float dc_angle = f1.Eval( 220. );


		float diff = p_angle - dc_angle;
		while( diff < -TMath::Pi() ){ diff += 2.*TMath::Pi(); }
		while( diff >  TMath::Pi() ){ diff -= 2.*TMath::Pi(); }
		return diff;
	}
	float get_phi()
	{
		TGraph graph(plane.size());
		for(unsigned int i=0;i<plane.size();i+=1)
		{
			graph.SetPoint( i, sqrt(y_global[i]*y_global[i] + x_global[i]*x_global[i]) , atan2( y_global[i] , x_global[i] ) );
		}
		TF1 f1( "f1", "1 ++ x ++ x*x" );
		graph.Fit(&f1,"NQ");
		float dc_angle = f1.Eval( 220. );
		float p_angle = p_phi * TMath::Pi() / 180. ;


		float diff = p_angle - dc_angle;
		while( diff < -TMath::Pi() ){ diff += 2.*TMath::Pi(); }
		while( diff >  TMath::Pi() ){ diff -= 2.*TMath::Pi(); }
		return diff;
	}
};





void convert_ancdch_to_lookup_table( const char* input_file = "ancdch.root", const char* output_file = "alpha_p_r_phi.root" )
{
	TFile ancfile(input_file);
	TNtuple* anc_ntuple = (TNtuple*)(ancfile.Get("AncDch"));
	float anc_p_tot, anc_p_phi;
	float anc_x_global, anc_y_global;
	float anc_plane;
	float anc_event;
	float anc_r;
	float anc_parent;
	float anc_arm;
	anc_ntuple->SetBranchAddress( "PTOT",     &anc_p_tot );
	anc_ntuple->SetBranchAddress( "PPHI",     &anc_p_phi );
	anc_ntuple->SetBranchAddress( "XINGLOB",  &anc_x_global );
	anc_ntuple->SetBranchAddress( "YINGLOB",  &anc_y_global );
	anc_ntuple->SetBranchAddress( "PLANE",    &anc_plane );
	anc_ntuple->SetBranchAddress( "EVENT",    &anc_event );
	anc_ntuple->SetBranchAddress( "R_VERTEX", &anc_r );
	anc_ntuple->SetBranchAddress( "IDPARENT", &anc_parent );


	map<int, DC_Event> dc_map;

	for(int ent=0;ent<anc_ntuple->GetEntries();ent+=1)
	{
		anc_ntuple->GetEntry(ent);

		if( (int)(anc_parent) != 0 ){continue;}

		int event = (int)(anc_event);
		if( dc_map.find(event) == dc_map.end() ){ dc_map[event] = DC_Event(); }

		DC_Event& dc_event = dc_map[event];

		dc_event.p_tot = anc_p_tot;
		dc_event.p_phi = anc_p_phi;

		dc_event.r = anc_r;

		dc_event.x_global.push_back(anc_x_global);
		dc_event.y_global.push_back(anc_y_global);
		dc_event.plane.push_back((int)(anc_plane));
	}

	TNtuple alpha_p_r_phi("alpha_p_r_phi", "each entry contains alpha, p, r, and phi", "alpha:p:r:phi");
	float alpha, p, r, phi;

	typedef map<int, DC_Event>::iterator it_type;
	for(it_type iterator = dc_map.begin(); iterator != dc_map.end(); iterator++)
	{
		DC_Event& dc_event = iterator->second;

		if( dc_event.plane.size() != 40 ){continue;}

		alpha = dc_event.get_alpha();
		phi   = dc_event.get_phi();
		p     = dc_event.p_tot;
		r     = dc_event.r;

		alpha_p_r_phi.Fill( alpha, p, r, phi );
	}

	TFile outfile(output_file, "recreate");
	outfile.cd();
	alpha_p_r_phi.Write();
	outfile.Close();
}


