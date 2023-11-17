void write_macros_one(){
	string name = "../macros/Pixel_109.mac";

	double X = -3.50007e-05;
	double Y = -17.5466;
	double Z = 138.896;

	//double theta = TMath::ACos(Z/D_xyz);
	//double phi   = TMath::ASin(X/D_xy);

	//cout<<"theta = "<<theta/TMath::Pi()*180.<<"  phi = "<<phi/TMath::Pi()*180.<<endl;
	//cout<<"Path: "<<name.c_str()<<endl;

	//double r_phi = -TMath::Tan(phi);
	//double r_theta = (1./(TMath::Cos(theta)*TMath::Cos(theta)) - 1) / (1+r_phi*r_phi);
	//r_theta = TMath::Sqrt(theta);

	//double d1 = r_theta;
	//double d2 = -r_theta*r_phi;
	//double d3 = -1;

	double d1 = -X;
	double d2 = -Y;
	double d3 = -Z;

	double r_theta = X/Z;
	double r_phi   = Y/Z/r_theta;

	double R=30.;
	double A = R/TMath::Sqrt(d1*d1+d2*d2+d3*d3);

	double c1 = -A*d1;
	double c2 = -A*d2;
	double c3 = -A*d3;

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
