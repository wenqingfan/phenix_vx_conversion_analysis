{
	TNtuple alpha_p_r_phi("alpha_p_r_phi", "each entry contains alpha, p, r, and phi", "alpha:p:r:phi");
	ifstream inf("lookup_finegrid.data");
	double alpha, p, r, phi;
	while(true)
	{
		inf>>alpha>>p>>r>>phi;
		if( inf.good() == false ){break;}
		alpha_p_r_phi.Fill( alpha, p, r, phi );
	}
	inf.close();
	TFile outfile("alpha_p_r_phi.root", "recreate");
	outfile.cd();
	alpha_p_r_phi.Write();
	outfile.Close();
}
