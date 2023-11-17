void macros_pixels_2(){

	double fposition[112][3];

	int tmp = 0;
	char cc[100];
	ifstream readIn("../../DetectorStructure/Parameter_UP/PositionUP.txt");
	for(int i=0;i<112;i++){
		readIn
			>>tmp
			>>fposition[i][0]
			>>fposition[i][1]
			>>fposition[i][2];
		cout
			<<"  i = "<<i
			<<"  x = "<<fposition[i][0]
			<<"  y = "<<fposition[i][1]
			<<"  z = "<<fposition[i][2]
			<<endl;
	}

	string path="../macros/Pixel/Pixel_";
	for(int i=0;i<112;i++){
		sprintf(cc,"%i",i);
		string name = path + cc + ".mac";
		
		double X = fposition[i][0];
		double Y = fposition[i][1];
		double Z = fposition[i][2];
		double R = TMath::Sqrt(X*X + Y*Y + Z*Z);

		double d1 = -X/R;
		double d2 = -Y/R;
		double d3 = -Z/R;

		double A=30.;
		double c1 = -A*d1;
		double c2 = -A*d2;
		double c3 = -A*d3;

		d1 = d1/(-d3);
		d2 = d2/(-d3);
		d3 = d3/(-d3);

		double r_theta = d1;
		double r_phi   = d2/(-r_theta);

		//cout
		//	<<"  Pixel_i = "<<i
		//	<<"  c1 = "<<c1<<"  c2 = "<<c2<<"  c3 = "<<c3
		//	<<"  d1 = "<<d1<<"  d2 = "<<d2<<"  d3 = "<<d3
		//	<<"  r_theta = "<<r_theta<<"  r_phi = "<<r_phi
		//	<<endl;

		
		ofstream output(name.c_str());
		output<<"/run/verbose 1"<<endl;
		output<<"/run/initialize"<<endl;

		output<<"/gun/sourceGen false"<<endl;
		output<<"/gps/particle gamma"<<endl;

		output<<"/gps/ene/type Mono "<<endl;

		output<<"/gps/pos/type   Plane"<<endl; 
		output<<"/gps/pos/shape  Circle"<<endl;
		output<<"/gps/pos/radius 0.001 cm"<<endl;


		output<<"/gps/pos/rot1 1 0 "<<r_theta<<endl;
		output<<"/gps/pos/rot2 "<<r_phi<<" 1 0"<<endl;
		output<<"/gps/direction "<<d1<<" "<<d2<<" "<<d3<<endl;
		output<<"/gps/pos/centre "<<c1<<" "<<c2<<" "<<c3<<endl;

		output.close();
	

	}

}
