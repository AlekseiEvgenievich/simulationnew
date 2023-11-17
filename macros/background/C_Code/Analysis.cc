int Analysis(){

	//string fname = "../Counts_cosmic_100s_1e1keV_1e5keV_50pts.txt";
	string fname = "../Counts_Albedo_100s_1e1keV_1e5keV_50pts.txt";
	ifstream fopen( fname.c_str() );

	double Energy[50];
	double Entries[50];

	double Weight[3];
	memset(Weight, 0, sizeof(Weight));

	for(int i=0; i<50; i++){
		fopen>>Energy[i]>>Entries[i];

		cout<<Energy[i]<<"  "<<Entries[i]<<endl;
		if(i<16)  
			Weight[0] = Weight[0] + Entries[i];
		else if(i>=16 && i<32)
			Weight[1] = Weight[1] + Entries[i];
		else 
			Weight[2] = Weight[2] + Entries[i];
	}

	cout<<Weight[0]<<"  "<<Weight[1]<<"  "<<Weight[2]<<endl;



	return 1;

}
