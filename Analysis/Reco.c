#include <math.h>
#include "TGraph.h"
#include <stddef.h>
#include "TRandom.h"

struct Fitvals {
    double grad;
    double inter;
    double chi2;
};

double Getmomentum(double B,double L, double grad) {
        double angle = atan(grad);
        double p = (0.3*B*L)/(sin(angle));
	return p;
}

struct Fitvals Getgradient(double xval[5], double zval[5]){
	struct Fitvals v;
	double xerrors[5] {0.1,0.1,0.1,0.1,0.1};
        double zerrors[5] {0  ,  0,  0,  0,  0};
	TGraphErrors* line = new TGraphErrors(5, zval, xval, zerrors, xerrors);
       	line->Fit("pol1", "q");
        double gradient = line->GetFunction("pol1")->GetParameter(1);
	double chi2 = line->GetFunction("pol1")->GetChisquare();
	double c = line->GetFunction("pol1")->GetParameter(0);
	v.grad = gradient/1000;
	v.chi2 = chi2;
	v.inter= c/1000;

	return v;
}

double straightline(double zval, double m, double c) {
	double x;
		x = m*zval + c;
	return x;
}

void Reco() {
	TFile *f = new TFile("rootfiles/newSim/LeadBefore.root");
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
	
	TH2F *h2 = new TH2F("h2","PLOT",100,-0.1,4.1,1000,-1,1);
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

	TRandom* rand = new TRandom();

	for (Int_t i=0; i<nentries; i++) {
		Long64_t tentry = t->LoadTree(i);
		//nD1Hits[i] = D1Hits->GetEntry(i);
		BD1_x->GetEntry(tentry);
		BD1_z->GetEntry(tentry);
		BD1_y->GetEntry(tentry);
	
		nD1Hits[i] = PD1_x->size();
		
		//out << PD1_x->size() << "," << PD1_z->size() << endl;
		for (UInt_t t=0; t<PD1_x->size(); t++) {
			fVD1_x[i][t] = PD1_x->at(t);// +  rand->Gaus(-0.1,0.1);
			fVD1_z[i][t] = PD1_z->at(t);
			fVD1_y[i][t] = PD1_y->at(t);

			h2->Fill(fVD1_z[i][t],fVD1_x[i][t]);
		}

		nD2Hits[i] = PD2_x->size();
		BD2_x->GetEntry(tentry);
                BD2_z->GetEntry(tentry);
		BD2_y->GetEntry(tentry);
                
		//nHits = nD1Hits + nD2Hits; 
                //double VD2_x[nD2Hits];
                //double VD2_z[nD2Hits];
		//double VD2_y[nD2Hits];

		for (UInt_t t=0; t<PD2_x->size(); t++) {
                        fVD2_x[i][t] = PD2_x->at(t);// + rand->Gaus(-0.1,0.1);
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
	
	fzArm1[0] =  - 6. ;// * m
	fzArm1[1] =  - 5.5;// * m
	fzArm1[2] =  - 5. ;// * m
	fzArm1[3] =  - 4.5;// * m 
	fzArm1[4] =  - 4. ;// * m

	fzArm2[0] =   4.0;// * m
	fzArm2[1] =   4.5;// * m
	fzArm2[2] =   5. ;// * m
	fzArm2[3] =   5.5;// * m
	fzArm2[4] =   6. ;// * m  

	TH1F *h10 = new TH1F("h10","Momentum",100,60,140);
	TMultiGraph *mg = new TMultiGraph();	
	TMultiGraph *mg2 = new TMultiGraph();
		 //nentries
	for(Int_t i=0; i<nentries; i++) {
		Int_t nhitsA1c0 = 0;
		Int_t nhitsA1c1 = 0;
		Int_t nhitsA1c2 = 0;
		Int_t nhitsA1c3 = 0;
		Int_t nhitsA1c4 = 0;
		Int_t nhitsA2c0 = 0;
		Int_t nhitsA2c1 = 0;
		Int_t nhitsA2c2 = 0;
		Int_t nhitsA2c3 = 0;
		Int_t nhitsA2c4 = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Counting how many Hits there are in each chambers 
//
		//i=334;
		//std:: cout << nD1Hits[i] << "," << nD2Hits[i] << std::endl;
		for (Int_t t=0; t<nD1Hits[i]; t++) {
			if (fVD1_z[i][t] == 0){
					nhitsA1c0++;	
			}else if (fVD1_z[i][t] == 1) {
					nhitsA1c1++;
			}else if (fVD1_z[i][t] == 2) { 
					nhitsA1c2++;
			}else if (fVD1_z[i][t] == 3) {
					nhitsA1c3++;
			}else if (fVD1_z[i][t] == 4) {
					nhitsA1c4++;
			}
		}
		for (Int_t t=0; t<nD2Hits[i]; t++) {
                        if (fVD2_z[i][t] == 0){
                                        nhitsA2c0++;
                        } else if( fVD2_z[i][t] == 1) {
                                        nhitsA2c1++;
                        } else if (fVD2_z[i][t] == 2) {
                                        nhitsA2c2++;
                        } else if (fVD2_z[i][t] == 3) { 
                                        nhitsA2c3++;
                        } else if (fVD2_z[i][t] == 4) {
                                        nhitsA2c4++;
                        }
                }
		//std::cout << nhitsA1c0 << "," << nhitsA1c1 << "," << nhitsA1c2 << "," << nhitsA1c3 << "," << nhitsA1c4 << std::endl;
		// Initialise vectors containing positions of each hit in each chamber
		double vecA1c0[nhitsA1c0][2];
		double vecA1c1[nhitsA1c1][2];
		double vecA1c2[nhitsA1c2][2];
		double vecA1c3[nhitsA1c3][2];
	       	double vecA1c4[nhitsA1c4][2];
		double vecA2c0[nhitsA2c0][2]; 
		double vecA2c1[nhitsA2c1][2];
	       	double vecA2c2[nhitsA2c2][2];
	       	double vecA2c3[nhitsA2c3][2];
                double vecA2c4[nhitsA2c4][2];
		//Loop over all hits in Arm one and save to vector 		
		Int_t w{0};
		Int_t v{0};
		Int_t y{0};
		Int_t z{0};
		Int_t x{0};
		for (Int_t t=0; t<nD1Hits[i]+1; t++) {
			if (fVD1_z[i][t] == 0){
					// If there is only one hit then that hit is saved
                                        if(nhitsA1c0==1) {
						vecA1c0[0][0] = fVD1_x[i][t];
						vecA1c0[0][1] = fVD1_y[i][t];
					} else {
						// If there is more than one hit
						// First check if the next hit recorded is in the same chamber if so 
						// check if it is within 50 micrometer
						if (fVD1_z[i][t] == fVD1_z[i][t+1] && 
								abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
							vecA1c0[w][0] = fVD1_x[i][t];
							vecA1c0[w][1] = fVD1_y[i][t];
							t++; // If so skip the next hit as it is discarded
							w++; // Next hit index of chamber
						}
						//If no hit is smaller than 50 micrometer just add like normal
						vecA1c0[w][0] = fVD1_x[i][t];
						vecA1c0[w][1] = fVD1_y[i][t];
						w++;
					}
			} 
			if (fVD1_z[i][t] ==1){
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c0==1) {
                                                vecA1c1[0][0] = fVD1_x[i][t];
						vecA1c1[0][1] = fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c1[v][0] = fVD1_x[i][t];
							vecA1c1[v][1] = fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        v++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c1[v][0] = fVD1_x[i][t];
						vecA1c1[v][1] = fVD1_y[i][t];
                                                v++;
                                        }
                        }
			if (fVD1_z[i][t] == 2) {
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c2==1) {
                                                vecA1c2[0][0] = fVD1_x[i][t];
						vecA1c2[0][1] = fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c2[y][0] = fVD1_x[i][t];
							vecA1c2[y][1] = fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        y++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c2[y][0] = fVD1_x[i][t];
						vecA1c2[y][1] = fVD1_y[i][t];
                                                y++;
                                        }
                        }
		       	if (fVD1_z[i][t] == 3) {
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c3==1) {
                                                vecA1c3[0][0] = fVD1_x[i][t];
						vecA1c3[0][1] = fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c3[x][0] = fVD1_x[i][t];
							vecA1c3[x][1] = fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        x++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c3[x][0] = fVD1_x[i][t];
						vecA1c3[x][1] = fVD1_y[i][t];
                                                x++;
                                        }
                       }
			if (fVD1_z[i][t] == 4) {
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA1c4==1) {
                                                vecA1c4[0][0] = fVD1_x[i][t];
						vecA1c4[0][1] = fVD1_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD1_z[i][t] == fVD1_z[i][t+1] &&
                                                                abs(fVD1_x[i][t]-fVD1_x[i][t+1])< pow(5,-5)){
                                                        vecA1c4[z][0] = fVD1_x[i][t];
							vecA1c4[z][1] = fVD1_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        z++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA1c4[z][0] = fVD1_x[i][t];
						vecA1c4[z][1] = fVD1_y[i][t];
                                                z++;
                                        }
		 		}
                        }	
		Int_t w2{0};
                Int_t v2{0};
                Int_t y2{0};
                Int_t z2{0};
		Int_t x2{0};
		//Now save all hits for Arm two
		for (Int_t t=0; t<nD2Hits[i]; t++) {
			if (fVD2_z[i][t] == 0) {
					// If there is only one hit then that hit is saved
                                        if(nhitsA2c0==1) {
						vecA2c0[0][0] = fVD2_x[i][t];
						vecA2c0[0][1] = fVD2_y[i][t];
					} else {
						// If there is more than one hit
						// First check if the next hit recorded is in the same chamber if so 
						// check if it is within 50 micrometer
						if (fVD2_z[i][t] == fVD2_z[i][t+1] && 
								abs(fVD2_x[i][t]-fVD2_x[i][t+1])< pow(5,-5)){
							vecA2c0[w2][0] = fVD2_x[i][t];
							vecA2c0[w2][1] = fVD2_y[i][t];
							t++; // If so skip the next hit as it is discarded
							w2++; // Next hit index of chamber
						}
						//If no hit is smaller than 50 micrometer just add like normal
						vecA2c0[w2][0] = fVD2_x[i][t];
						vecA2c0[w2][1] = fVD2_y[i][t];
						w2++;
					}
                         } else if (fVD2_z[i][t] == 1){
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA2c0==1) {
                                                vecA2c1[0][0] = fVD2_x[i][t];
						vecA2c1[0][1] = fVD2_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD2_z[i][t] == fVD2_z[i][t+1] &&
                                                                abs(fVD2_x[i][t]-fVD2_x[i][t+1])< pow(5,-5)){
                                                        vecA2c1[v2][0] = fVD2_x[i][t];
							vecA2c1[v2][1] = fVD2_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        v2++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA2c1[v][0] = fVD2_x[i][t];
						vecA2c1[v][1] = fVD2_y[i][t];
                                                v2++;
                                        }
			} else if (fVD2_z[i][t] == 2) {
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA2c2==1) {
                                                vecA2c2[0][0] = fVD2_x[i][t];
						vecA2c2[0][1] = fVD2_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD2_z[i][t] == fVD2_z[i][t+1] &&
                                                                abs(fVD2_x[i][t]-fVD2_x[i][t+1])< pow(5,-5)){
                                                        vecA2c2[y2][0] = fVD2_x[i][t];
							vecA2c2[y2][1] = fVD2_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        y2++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA2c2[y2][0] = fVD2_x[i][t];
						vecA2c2[y2][1] = fVD2_y[i][t];
                                                y2++;
                                        }
			} else if (fVD2_z[i][t] == 3){
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA2c3==1) {
                                                vecA2c3[0][0] = fVD2_x[i][t];
						vecA2c3[0][1] = fVD2_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD2_z[i][t] == fVD2_z[i][t+1] &&
                                                                abs(fVD2_x[i][t]-fVD2_x[i][t+1])< pow(5,-5)){
                                                        vecA2c3[x2][0] = fVD2_x[i][t];
							vecA2c3[x2][1] = fVD2_y[i][t];
                                                        t++; // If so skip the next hit as it is discarded
                                                        x2++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA2c3[x2][0] = fVD2_x[i][t]; 
						vecA2c3[x2][1] = fVD2_y[i][t];
                                                x2++;
                                        }
			} else if (fVD2_z[i][t] == 4){
                                        // If there is only one hit then that hit is saved
                                        if(nhitsA2c4==1) {
                                                vecA2c4[0][0] = fVD2_x[i][t];
						vecA2c4[0][1] = fVD2_y[i][t];
                                        } else {
                                                // If there is more than one hit
                                                // First check if the next hit recorded is in the same chamber if so 
                                                // check if it is within 50 micrometer
                                                if (fVD2_z[i][t] == fVD2_z[i][t+1] &&
                                                                abs(fVD2_x[i][t]-fVD2_x[i][t+1])< pow(5,-5)){
                                                        vecA2c4[z2][0] = fVD2_x[i][4];
							vecA2c4[z2][1] = fVD2_y[i][4];
                                                        t++; // If so skip the next hit as it is discarded
                                                        z2++; // Next hit index of chamber
                                                }
                                                //If no hit is smaller than 50 micrometer just add like normal
                                                vecA2c4[z2][0] = fVD2_x[i][t];
						vecA2c4[z2][1] = fVD2_y[i][t];
                                                z2++;
                                        }
				}
                        }
		double zval1[5] {fzArm1[0], fzArm1[1], fzArm1[2], fzArm1[3], fzArm1[4]};
		double zval2[5] {fzArm2[0], fzArm2[1], fzArm2[2], fzArm2[3], fzArm2[4]};
		double grad_before{0};
		double grad_after{0};
		double inter_before{0};
		double inter_after{0};
		
		if (false) {
                double xval1[5] {vecA1c0[0][0], vecA1c1[0][0], vecA1c2[0][0], vecA1c3[0][0], vecA1c4[0][0]};
                struct Fitvals vals_pre;
		vals_pre = Getgradient(xval1,zval1);
		grad_before = vals_pre.grad;
		double chi2_before = vals_pre.chi2;
		
		double xval2[5] {vecA2c0[0][0], vecA2c1[0][0], vecA2c2[0][0], vecA2c3[0][0], vecA2c4[0][0]};
		struct Fitvals vals_post;
		vals_post = Getgradient(xval2,zval2);
		grad_after = vals_post.grad;
		double chi2_after = vals_post.chi2;
		}
	//	double grad = (grad_before+grad_after)/(1+grad_before*grad_after);
	//	double p = Getmomentum(0.5,2,grad);
	//	cout << p << endl;
	//	h10->Fill(p);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		This allows for multi hit detection. It loops over all hits and find the best fit for the 
// 		tragectory. 
// 		As this is a nested loop it takes a while to run so just taking normal hits is fine 
// 		for the moment.
//
	if (true) {
		double xval1op[5];
		struct Fitvals optimise;
		double grad_before_op;
		double inter_before_op;
		double chi2_min {100};

		for (Int_t p=0; p < nhitsA1c0; p++) {
			for (Int_t o=0; o < nhitsA1c1; o++) {
				for (Int_t l=0; l<nhitsA1c2; l++) {
					for (Int_t k=0; k<nhitsA1c3; k++) {
                                		for (Int_t j=0; j<nhitsA1c3; j++) {
							xval1op[0] = vecA1c0[p][0];
							xval1op[1] = vecA1c1[o][0];
							xval1op[2] = vecA1c2[l][0];
							xval1op[3] = vecA1c3[k][0];
							xval1op[4] = vecA1c4[j][0];
							optimise = Getgradient(xval1op,zval1);

							if (optimise.chi2 < chi2_min) {
								grad_before_op = optimise.grad;
								inter_before_op = optimise.inter;
								chi2_min = optimise.chi2;
							}
										
                                        	}
                                	}
				}
			}
		}
		double xval2op[5];
                struct Fitvals optimise2;
                double grad_after_op{0};
		double inter_after_op{0};
                double chi2_min2 {100};

                for (Int_t p=0; p < nhitsA2c0; p++) {
                        for (Int_t o=0; o < nhitsA2c1; o++) {
                                for (Int_t l=0; l<nhitsA2c2; l++) {
                                        for (Int_t k=0; k<nhitsA2c3; k++) {
                                                for (Int_t j=0; j<nhitsA1c3; j++) {
                                                        xval2op[0] = vecA2c0[p][0];
                                                        xval2op[1] = vecA2c1[o][0];
                                                        xval2op[2] = vecA2c2[l][0];
                                                        xval2op[3] = vecA2c3[k][0];
                                                        xval2op[4] = vecA2c4[j][0];
                                                        optimise2 = Getgradient(xval2op,zval2);

                                                        if (optimise2.chi2 < chi2_min2) {
                                                                grad_after_op = optimise2.grad;
								inter_after_op = optimise2.inter;
								chi2_min2 = optimise2.chi2;
                                                        }

                                                }
                                        }
                                }
                        }
                }
		grad_before = grad_before_op;
		inter_before = inter_before_op;
		grad_after  = grad_after_op;
		inter_after = inter_after_op;
	}
	double grad = (grad_before+grad_after)/(1+grad_before*grad_after);
        double p = abs(Getmomentum(0.5,2,grad));
        //cout << p << endl;
        h10->Fill(p);
	std::cout << "Event = " << i << std::endl;
	//cout<<endl;
	double fitval1[5];
	for (int p = 0; p<5; p++) {	
		fitval1[p] = straightline(zval1[p], grad_before, inter_before) * 1000;
	}
	TGraph* gr = new TGraph(5,zval1,fitval1);
	double fitval2[5];
        for (int l = 0; l<5; l++) {
                fitval2[l] = straightline(zval2[l], grad_after, inter_after)*1000;
	}
        TGraph* gr2 = new TGraph(5,zval2,fitval2);
	
	if (i ==20) {
		double ECx = straightline( 7, grad_after,inter_after) *1000;
		std::cout << ECx << std::endl;
	}

	mg->Add(gr);
	mg2->Add(gr2); 
	}

	h10->Fit("gaus");

	cout << "Saving output pdf" <<endl;
	TCanvas canvas2("canvas");
	canvas2.Print("output.pdf[");
	h10->Draw();
	canvas2.Print("output.pdf[");
	canvas2.Clear();
	mg->Draw("AC");
	canvas2.Print("output.pdf[");
        canvas2.Clear();
        mg2->Draw("AC");
	canvas2.Print("output.pdf");
        canvas2.Print("output.pdf]");

}

