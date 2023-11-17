void macros_500kev(){
//Note: This scrip is to set up the macro file for different angle(theta and phi)

	string path = "../macros/aleksei/500kev/";

	int i_phi = 0;
	//for(int i_phi=0; i_phi<180; i_phi+=5)
	for(int i_theta=0; i_theta<=90; i_theta+=5){

			
		double theta = i_theta;
		if(i_theta==90)
			theta = 89.99;
		double phi   = i_phi;
		theta = (theta)/180.*TMath::Pi();
		phi = (phi)/180.*TMath::Pi();

		double d1 = -TMath::Sin(theta)*TMath::Cos(phi);
		double d2 = -TMath::Sin(theta)*TMath::Sin(phi);
		double d3 = -TMath::Cos(theta);
		cout<<"i_theta = "<<i_theta<<"  d3 = "<<d3<<endl;

		double A = 30;
		double c1 = -A*d1;
		double c2 = -A*d2;
		double c3 = -A*d3;

		d1 = d1/(-d3);
		d2 = d2/(-d3);
		d3 = d3/(-d3);

		double r_theta = d1;
		double r_phi   = 0 ;
		if(r_theta != 0 )
			r_phi   = d2/(-r_theta);
		
		string name = "gamma_circle_";
		char cc[50];
		sprintf(cc, "%i_%i",i_theta,i_phi);
		//sprintf(cc,"%.*f",1,i_theta);

		name = path + name + cc + ".mac";
		ofstream output(name.c_str());

		output<<"/run/verbose 1"<<endl;
		output<<"/run/initialize"<<endl;

		output<<"/gun/sourceGen false"<<endl;
		output<<"/gps/particle gamma"<<endl;

		output<<"/gps/ene/type Mono "<<endl;
    output<<"/gps/ene/mono 500 keV"<<endl;
		output<<"/gps/pos/type   Plane"<<endl; 
		output<<"/gps/pos/shape  Circle"<<endl;
		output<<"/gps/pos/radius 20 cm"<<endl;
		output<<"/gps/pos/rot1 1 0 "<<r_theta<<endl;
		output<<"/gps/pos/rot2 "<<r_phi<<" 1 0"<<endl;
		output<<"/gps/pos/centre "<<c1<<" "<<c2<<" "<<c3<<endl;
		output<<"/gps/direction "<<d1<<" "<<d2<<" "<<d3<<endl;

	  output<<"/run/beamOn 10000"<<endl;
		output.close();

	}





}
