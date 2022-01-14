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
        double fVD1_z[nentries][5];
       	double fVD1_y[nentries][5];
	
	double fVD2_x[nentries][5];
        double fVD2_z[nentries][5];
        double fVD2_y[nentries][5];
	
	for (Int_t i=0; i<nentries; i++) {
		Long64_t tentry = t->LoadTree(i);
		Int_t nD1Hits = D1Hits->GetEntry(i);
		BD1_x->GetEntry(tentry);
		BD1_z->GetEntry(tentry);
		BD1_y->GetEntry(tentry);
	
		//double VD1_x[nD1Hits];
		//double VD1_z[nD1Hits];
		//double VD1_y[nD1Hits];

		//out << PD1_x->size() << "," << PD1_z->size() << endl;
		for (UInt_t t=0; t<PD1_x->size(); t++) {
			fVD1_x[i][t] = PD1_x->at(t);
			fVD1_z[i][t] = PD1_z->at(t);
			fVD1_y[i][t] = PD1_y->at(t);

			//cout << VD1_z[t] << endl;
			h2->Fill(fVD1_z[i][t],fVD1_x[i][t]);
		}

		Int_t nD2Hits = D2Hits->GetEntry(i);
		BD2_x->GetEntry(tentry);
                BD2_z->GetEntry(tentry);
		BD2_y->GetEntry(tentry);
        
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
	cout << "Saving pdf" <<endl;
        TCanvas canvas("canvas");
        canvas.Print("output.pdf[");
        canvas.Clear();
	h2->Draw();
	canvas.Print("output.pdf[");
        canvas.Clear();
	h3->Draw();
	canvas.Print("output.pdf");
        canvas.Print("output.pdf]");
//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//	\/											\/
//	\/				Start of Reconstruction					\/
//	\/											\/
//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
	
	double fzChmbr1[5] , fzChmbr2[5];
	//std::vector<double> vec[10]; 
	
	fzChmbr1[0] = - 6. ;// * m
	fzChmbr1[1] = - 5.5;// * m
	fzChmbr1[2] = - 5. ;// * m
	fzChmbr1[3] = - 4.5;// * m 
	fzChmbr1[4] = - 4. ;// * m

	fzChmbr2[0] =   4.0;// * m
	fzChmbr2[1] =   4.5;// * m
	fzChmbr2[2] =   5. ;// * m
	fzChmbr2[3] =   5.5;// * m
	fzChmbr2[4] =   6. ;// * m  
	
	for(Int_t i=0; i<nentries; i++) {
		std::vector<double> vec[10];    
		for (Int_t t=0; t<5; t++) {
			vec[i] = {fVD1_x[i][t], fVD1_y[i][t] , fzChmbr1[t] };
			vec[i+5] = {fVD2_x[i][t], fVD2_y[i][t] , fzChmbr2[t] }; 
		}
	}

}
