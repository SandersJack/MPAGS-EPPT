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

	
	for (Int_t i=0; i<nentries; i++) {
		Long64_t tentry = t->LoadTree(i);
		Int_t nD1Hits = D1Hits->GetEntry(i);
		BD1_x->GetEntry(tentry);
		BD1_z->GetEntry(tentry);
		BD1_y->GetEntry(tentry);
	
		double VD1_x[nD1Hits];
		double VD1_z[nD1Hits];
		double VD1_y[nD1Hits];

		//out << PD1_x->size() << "," << PD1_z->size() << endl;
		for (UInt_t t=0; t<PD1_x->size(); t++) {
			VD1_x[t] = PD1_x->at(t);
			VD1_z[t] = PD1_z->at(t);
			VD1_y[t] = PD1_y->at(t);

			//cout << VD1_z[t] << endl;
			h2->Fill(VD1_z[t],VD1_x[t]);
		}

		Int_t nD2Hits = D2Hits->GetEntry(i);
		BD2_x->GetEntry(tentry);
                BD2_z->GetEntry(tentry);
		BD2_y->GetEntry(tentry);
        
                double VD2_x[nD2Hits];
                double VD2_z[nD2Hits];
		double VD2_y[nD2Hits];

		for (UInt_t t=0; t<PD2_x->size(); t++) {
                        VD2_x[t] = PD2_x->at(t);
                        VD2_z[t] = PD2_z->at(t);
			VD2_y[t] = PD2_y->at(t);
                        h3->Fill(VD2_z[t],VD2_x[t]);
                }

		
		
	}
	h2->Draw();
	h3->Draw();

//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//	\/											\/
//	\/				Start of Reconstruction					\/
//	\/											\/
//	\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

	double fzChmbr0 = - 6. ;// * m
	double fzChmbr1 = - 5.5;// * m
	double fzChmbr2 = - 5. ;// * m
	double fzChmbr3 = - 4.5;// * m 
	double fzChmbr4 = - 4. ;// * m

	double fzChmbr5 =   4.0;// * m
	double fzChmbr6 =   4.5;// * m
	double fzChmbr7 =   5. ;// * m
	double fzChmbr8 =   5.5;// * m
	double fzChmbr9 =   6. ;// * m  
}
