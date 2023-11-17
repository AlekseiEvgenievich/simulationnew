void write_macros(){
//NOTE:
//		 rot1(1,     0, r_theta)
//		 rot2(r_phi, 1, 0)
//     D = -(rot1*rot2) = (r_theta, -r_theta*r_phi, -1)
		
		int tmp = 0;
		string path = "./macros/";
		for(double i=0;i<=90;i+=5){

				tmp++;
				cout<<"Number: "<<tmp<<"  theta = "<<i<<endl;
				double theta = i;
				double phi = 0;
				theta = (theta)/180*TMath::Pi();
				phi = (phi)/180*TMath::Pi();

				double r_phi = -TMath::Tan(phi);
				double r_theta = (1./(TMath::Cos(theta)*TMath::Cos(theta)) - 1) / (1+r_phi*r_phi);
				r_theta = TMath::Sqrt(theta);


				double d1 = r_theta;
				double d2 = -r_theta*r_phi;
				double d3 = -1;

				double R=30.;
				double A = R/TMath::Sqrt(d1*d1+d2*d2+d3*d3);

				double c1 = -A*d1;
				double c2 = -A*d2;
				double c3 = -A*d3;

				string name = "gamma_circle_";
				char cc[50];
				sprintf(cc,"%i",i);
				name = path + name + cc + ".mac";
				ofstream output(name.c_str());

				output<<"/run/verbose 1"<<endl;
				output<<"/run/initialize"<<endl;

				output<<"/gun/sourceGen false"<<endl;
				output<<"/gps/particle gamma"<<endl;

				output<<"/gps/ene/type Mono "<<endl;

				output<<"/gps/pos/type   Plane"<<endl; 
				output<<"/gps/pos/shape  Circle"<<endl;
				output<<"/gps/pos/radius 20 cm"<<endl;


				output<<"/gps/pos/rot1 1 0 "<<r_theta<<endl;
				output<<"/gps/pos/rot2 "<<r_phi<<" 1 0"<<endl;
				output<<"/gps/direction "<<d1<<" "<<d2<<" "<<d3<<endl;
				output<<"/gps/pos/centre "<<c1<<" "<<c2<<" "<<c3<<endl;

				output.close();

		}
}
