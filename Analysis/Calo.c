#include <stddef.h>

void Calo() {
	TFile *f = new TFile("rootfiles/100GeV0_5T.root");
        f->ls();
        TTree *t = (TTree*)f->Get("B5");
        Int_t nentries = t->GetEntries();
	t->Print();

	std::vector<double> *ECEnergy = 0;
	TBranch *EC = 0;

	std::vector<double> *HCEnergy = 0;
	TBranch *HC = 0;

	double V_EC[nentries][80];
	double V_HC[nentries][20];

	TH1F *h1 = new TH1F("h1","EC",301,-0.1,300);
	TH2F *h2 = new TH2F("h2","EC",81,-0.1,80.1,301,-0.1,300);
	auto *h3 = new TH3F("h3","EC Energy in Cells",20, -0.1, 19.1, 4, -0.1,3.1, 301, -0.1, 300.1); 
	TH1F *h4 = new TH1F("h1","HC",101,-0.1,100);
        TH2F *h5 = new TH2F("h2","HC",21,-0.1,20.1,101,-0.1,100);
        auto *h6 = new TH3F("h3","HC Energy in Cells",10, -0.1, 9.1, 2, -0.1,1.1, 301, -0.1, 300.1);


	t->SetBranchAddress("ECEnergyVector",&ECEnergy, &EC);
//	std::cout << ECEnergy << std::endl;
	for (Int_t i=0; i<nentries; i++) {
                Long64_t tentry = t->LoadTree(i);
		EC->GetEntry(tentry);
		//std::cout << ECEnergy->size() << std::endl;
		for (UInt_t t=0; t<ECEnergy->size(); t++) {
			V_EC[i][t] = ECEnergy->at(t);
			if (V_EC[i][t] != 0) {
				h1->Fill(V_EC[i][t]);
				h2->Fill(t,V_EC[i][t]);
			}
		}
		
	}

	t->SetBranchAddress("HCEnergyVector",&HCEnergy, &HC);
//      std::cout << ECEnergy << std::endl;
        for (Int_t i=0; i<nentries; i++) {
                Long64_t tentry = t->LoadTree(i);
                HC->GetEntry(tentry);
                for (UInt_t t=0; t<HCEnergy->size(); t++) {
                        V_HC[i][t] = HCEnergy->at(t);
                        if (V_HC[i][t] != 0) {
                                h4->Fill(V_HC[i][t]);
                                h5->Fill(t,V_HC[i][t]);
                        }
                }

        }
	

	for (Int_t t=0; t<80; t++) {
		double col = floor(t/4);
		int row = t-(4*col);
		std::cout << col << "," << row <<std::endl;
		h3->Fill(col,row,V_EC[10][t]); 
	}

	for (Int_t t=0; t<20; t++) {
                double col = floor(t/2);
                int row = t-(2*col);
                std::cout << col << "," << row <<std::endl;
                h6->Fill(col,row,V_HC[10][t]);
		//std:cout << V_HC[10][t] << std::endl;
        }
	//h6->Fill(4.,0,0);
	cout << "Saving Raw Calo output pdf" <<endl;
        TCanvas canvas("canvas");
        canvas.Print("rawCalooutput.pdf[");
        canvas.Clear();
        h1->Draw();
        canvas.Print("rawCalooutput.pdf[");
        canvas.Clear();
        h2->Draw();
        canvas.Print("rawCalooutput.pdf");
        canvas.Clear();
        h3->Draw("LEGO");
        canvas.Print("rawCalooutput.pdf");
        canvas.Clear();
        h4->Draw();
        canvas.Print("rawCalooutput.pdf[");
        canvas.Clear();
        h5->Draw();
        canvas.Print("rawCalooutput.pdf");
        canvas.Clear();
        h6->Draw("LEGO");
        canvas.Print("rawCalooutput.pdf");
	canvas.Print("rawCalooutput.pdf]");


}
