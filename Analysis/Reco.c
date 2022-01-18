#include <math.h>
#include <stddef.h>


double momentum(double B,double L, double theta) {
        return (0.3*B*L)/(theta);
}


void Reco() {
	TFile *f = new TFile("rootfiles/jack_1.root");
	f->ls();
	TTree *t = (TTree*)f->Get("B5");
	Int_t nentries = t->GetEntries();

	TBranch *D1Hits = t->GetBranch("Dc1Hits");
	TBranch *D2Hits = t->GetBranch("Dc2Hits");
	//TBranch *BD1_x = t->GetBranch("Dc1HitsVector_x");
	
	std::vector<double> *PD1_x = 0;
	TBranch *BD1_x = 0;
	std::vector<double> *PD1_z = 0;
        TBranch *BD1_z = 0;
	std::vector<double> *PD1_y = 0;
        TBranch *BD1_y = 0;
	t->SetBranchAddress("Dc1HitsVector_x",&PD1_x, &BD1_x);
	t->SetBranchAddress("Dc1HitsVector_z",&PD1_z, &BD1_z);
	t->SetBranchAddress("Dc1HitsVector_y",&PD1_y, &BD1_y);	

	std::vector<double> *PD2_x = 0;
        TBranch *BD2_x = 0;
        std::vector<double> *PD2_z = 0;
        TBranch *BD2_z = 0;
	std::vector<double> *PD2_y = 0;
        TBranch *BD2_y = 0;
        t->SetBranchAddress("Dc2HitsVector_x",&PD2_x, &BD2_x);
        t->SetBranchAddress("Dc2HitsVector_z",&PD2_z, &BD2_z);
	t->SetBranchAddress("Dc2HitsVector_y",&PD2_y, &BD2_y);
	
	TH2F *h2 = new TH2F("h2","PLOT",100,-0.1,4.1,1000,-100,1);
        TH2F *h3 = new TH2F("h3","PLOT",100,-0.1,4.1,1000,-100,1);

	double fVD1_x[nentries][5];
        Int_t fVD1_z[nentries][5];
       	double fVD1_y[nentries][5];
	
	double fVD2_x[nentries][5];
        Int_t fVD2_z[nentries][5];
        double fVD2_y[nentries][5];
	
	//Int_t nHits[nentries];
	Int_t nD1Hits[nentries];
	Int_t nD2Hits[nentries];

	for (Int_t i=0; i<nentries; i++) {
		Long64_t tentry = t->LoadTree(i);
		nD1Hits[i] = D1Hits->GetEntry(i);
		BD1_x->GetEntry(tentry);
		BD1_z->GetEntry(tentry);
		BD1_y->GetEntry(tentry);
	

		//out << PD1_x->size() << "," << PD1_z->size() << endl;
		for (UInt_t t=0; t<PD1_x->size(); t++) {
			fVD1_x[i][t] = PD1_x->at(t);
			fVD1_z[i][t] = PD1_z->at(t);
			fVD1_y[i][t] = PD1_y->at(t);

			//cout << VD1_z[t] << endl;
			h2->Fill(fVD1_z[i][t],fVD1_x[i][t]);
		}

		nD2Hits[i] = D2Hits->GetEntry(i);
		BD2_x->GetEntry(tentry);
                BD2_z->GetEntry(tentry);
		BD2_y->GetEntry(tentry);
                
		//nHits = nD1Hits + nD2Hits; 
                //double VD2_x[nD2Hits];
                //double VD2_z[nD2Hits];
		//double VD2_y[nD2Hits];

		for (UInt_t t=0; t<PD2_x->size(); t++) {
                        fVD2_x[i][t] = PD2_x->at(t);
                        fVD2_z[i][t] = PD2_z->at(t);
			fVD2_y[i][t] = PD2_y->at(t);
                        h3->Fill(fVD2_z[i][t],fVD2_x[i][t]);
                }

		
		
	}
	cout << "Saving Raw output pdf" <<endl;
        TCanvas canvas("canvas");
        canvas.Print("rawoutput.pdf[");
        canvas.Clear();
	h2->Draw();
	canvas.Print("rawoutput.pdf[");
        canvas.Clear();
	h3->Draw();
	canvas.Print("rawoutput.pdf");
        canvas.Print("rawoutput.pdf]");
//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//	\/											\/
//	\/				Start of Reconstruction					\/
//	\/											\/
//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
	
	double fzArm1[5] , fzArm2[5];
	//std::vector<double> vec[10]; 
	
	fzArm1[0] = - 6. ;// * m
	fzArm1[1] = - 5.5;// * m
	fzArm1[2] = - 5. ;// * m
	fzArm1[3] = - 4.5;// * m 
	fzArm1[4] = - 4. ;// * m

	fzArm2[0] =   4.0;// * m
	fzArm2[1] =   4.5;// * m
	fzArm2[2] =   5. ;// * m
	fzArm2[3] =   5.5;// * m
	fzArm2[4] =   6. ;// * m  

	TH1F *h10 = new TH1F("h10","Momentum",100,0,10);
		
	for(Int_t i=0; i<nentries; i++) {
		double vec[10][3];
		Int_t nhitsA1c0,nhitsA1c1,nhitsA1c2,nhitsA1c3,nhitsA1c4 = 0;
		Int_t nhitsA2c0,nhitsA2c1,nhitsA2c2,nhitsA2c3,nhitsA2c4 = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Counting how many Hits there are in each chambers 
//
//

		for (Int_t t=0; t<nD1Hits[i]; t++) {
			switch(fVD1_z[i][t]){
				case 0:
					nhitsA1c0++;	
				case 1:
					nhitsA1c1++;
				case 2: 
					nhitsA1c2++;
				case 3:
					nhitsA1c3++;
				case 4:
					nhitsA1c4++;
			}
		}

		for (Int_t t=0; t<nD2Hits[i]; t++) {
                        switch(fVD2_z[i][t]){
                                case 0:
                                        nhitsA2c0++;
                                case 1:
                                        nhitsA2c1++;
                                case 2:
                                        nhitsA2c2++;
                                case 3:
                                        nhitsA2c3++;
                                case 4:
                                        nhitsA2c4++;
                        }
                }


		// Initialise vectors containing positions of each hit in each chamber
		double vecA1c0[nhitsA1c0][2] , vecA1c1[nhitsA1c1][2], vecA1c2[nhitsA1c2][2] , vecA1c3[nhitsA1c3][2],
	       		vecA1c4[nhitsA1c4][2];
		double vecA2c0[nhitsA2c0][2] , vecA2c1[nhitsA2c1][2], vecA2c2[nhitsA2c2][2] , vecA2c3[nhitsA2c3][2],
                        vecA2c4[nhitsA2c4][2];
		//Loop over all hits in Arm one and save to vector 		
		for (Int_t t=0; t<nD1Hits[i]; t++) {
			switch(fVD1_z[i][t]){
                                case 0:
					// If there is only one hit then that hit is saved
                                        if(nhitsA1c0==1) {
						vecA1c0[0][0], vecA1c0[0][1] = fVD1_x[i][t], fVD1_y[i][t];
					} else {
						// If there is more than one hit
						Int_t w = 0;
						// First check if the next hit recorded is in the same chamber if so 
						// check if it is within 50 micrometer
						if (fVD1_z[i][t] == fVD1_z[i][t+1] && 
								abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
							vecA1c0[w][0] , vecA1c0[w][1] = fVD1_x[i][t], fVD1_y[i][t];
							t++; // If so skip the next hit as it is discarded
							w++; // Next hit index of chamber
						}
						//If no hit is smaller than 50 micrometer just add like normal
						vecA1c0[w][0] , vecA1c0[w][1] = fVD1_x[i][t], fVD1_y[i][t];
						w++;
					}
                                case 1:
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c0==1) {
                                                vecA1c1[0][0], vecA1c1[0][1] = fVD1_x[i][t], fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                Int_t v = 0;
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c1[v][0] , vecA1c1[v][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        v++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c1[v][0] , vecA1c1[v][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                v++;
                                        }
                                case 2:
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c2==1) {
                                                vecA1c2[0][0], vecA1c2[0][1] = fVD1_x[i][t], fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                Int_t y = 0;
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c2[y][0] , vecA1c2[y][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        y++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c2[y][0] , vecA1c2[y][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                y++;
                                        }
                                case 3:
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c3==1) {
                                                vecA1c3[0][0], vecA1c3[0][1] = fVD1_x[i][t], fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                Int_t x = 0;
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c3[x][0] , vecA1c3[x][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        x++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c3[x][0] , vecA1c3[x][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                x++;
                                        }
                                case 4:
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c4==1) {
                                                vecA1c4[0][0], vecA1c4[0][1] = fVD1_x[i][t], fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                Int_t z = 0;
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c4[z][0] , vecA1c4[z][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        z++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c4[z][0] , vecA1c4[z][1] = fVD1_x[i][t], fVD1_y[i][t];
                                                z++;
                                        }
                        }
		}
/*
		for (Int_t t=0; t<nhits; t++) {
			vec[t][0] = fVD1_x[i][t];
			vec[t][1] = fVD1_y[i][t];
			vec[t][2] = fzChmbr1[t];
			vec[t+5][0] = fVD2_x[i][t];
			vec[t+5][1] = fVD2_y[i][t];
			vec[t+5][2] = fzChmbr2[t] ; 
			//cout << fzChmbr2[t] << endl;
		}
		double xdiff = abs(vec[4][0]-vec[9][0]);
		//cout <<vec[0][0] <<vec[9][0] << endl;
		double zdiff = abs(vec[4][2] -vec[9][2]);
		double theta = atan(xdiff/zdiff);
		double length = sqrt(pow(xdiff,2) + pow(zdiff,2));
		double mom = momentum(0.5,length,theta);
		cout << mom  << endl;
*/
	}

}

