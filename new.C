
double Func_DiffuseGamma(double *e, double* par){

	//cm^-2*s^1*sr^-1*keV^-1
	double f;
	double eng = e[0];
	if(eng<=60.0){
		f = 7.877*TMath::Power(eng,-1.29)*TMath::Exp(-eng/41.13);
	}
	else if(eng>60.0){
		f = 4.32e-4*TMath::Power(eng/60., -6.5)
			+ 8.4e-3*TMath::Power(eng/60.,-2.58)
			+ 4.8e-4*TMath::Power(eng/60.,-2.05);
	}

	return f;
}

double Func_AlbedoGamma(double* e, double* par){

	////cm^-2*s^1*sr^-1*keV^-1
	double f;
	double eng = e[0];

	if(eng<=200.0){
		f = 1.87e-2/( TMath::Power(eng/33.7,-5.) + TMath::Power(eng/33.7, 1.72) );
	}
	else if(eng>200. && eng<=20000.){
		f = 1.01e-4*TMath::Power(eng/1000., -1.34);
	}
	else if(eng>20000.){
		f = 7.29e-4*TMath::Power(eng/1000.,-2.0);
	}
	return f;
}

double Func_PrimaryProton(double* e, double* par){

	//cm^-2*s^1*sr^-1*keV^-1
	double f;
	double eng = e[0];
	double F = 1.23e-8;
	double A = 2.39e-6;
	double Z = 1.;//atomic number
	double ee = 1.602e-19;
	double a = 0.155;
	double b = 2.83;
	double E1cut = 5.1e5;//keV
	double phi = 6.5e5;//kV
	double Mc2 = 938272.0813;//Proton Mass in units of keV/c2
	//double E2cut = 12.25e6;//keV
	double E2cut = 1.13e7;//keV


	f = F*TMath::Power(eng/1e6,-a)*TMath::Exp(-TMath::Power(eng/E1cut, -a+1) ) 
		+ A*TMath::Power((eng+Z*ee*phi)/1.e6, -b)
		* (TMath::Power(eng+Mc2,2.)-Mc2*Mc2)/(TMath::Power(eng+Mc2+Z*ee*phi,2.)-Mc2*Mc2)
		*1./(1 + TMath::Power(eng/E2cut, -12.));

	return f;
}

double Func_SecondaryProton(double* e, double* par){

	//cm^-2*s^1*sr^-1*keV^-1
	double f;
	double eng = e[0];
	double a = 0.155;
	double Ecut = 5.1e5;//keV

	f = 1.23e-8*TMath::Power(eng/1.e6,-a)
		* TMath::Exp(-TMath::Power(eng/Ecut,-a+1));

	return f;
}

double Func_AlbedoNeutron(double* e, double* par){

	//cm^-2*s^1*sr^-1*keV^-1
	double f;
	double eng = e[0];
	if(eng>= 10. && eng<1.e3){
		f = 9.98e-8*TMath::Power(eng/1.e6, -0.5);
	}
	else if(eng>=1.e3 && eng<1.e5){
		f = 3.16e-9*TMath::Power(eng/1.e6, -1.);
	}
	else if(eng>=1.e5 && eng<1.e8){
		f = 3.16e-10*TMath::Power(eng/1.e6, -2.);
	}

	return f;
}
