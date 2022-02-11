#include <stddef.h>

void Calo() {
	TFile *f = new TFile("rootfiles/1000evnt.root");
        f->ls();
        TTree *t = (TTree*)f->Get("B5");
        Int_t nentries = t->GetEntries();
	t->Print();

	std::vector<double> *ECEnergy = 0;
	TBranch *EC = 0;

	double V_EC[nentries][80];

	TH1F *h1 = new TH1F("h1","PLOT",301,-0.1,300);
	TH2F *h2 = new TH2F("h2","PLOT",81,-0.1,80.1,1000,-0.1,300);

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
			//std::cout << V_EC[i][t] << std::endl;
		}
		
	}
	cout << "Saving Raw Calo output pdf" <<endl;
        TCanvas canvas("canvas");
        canvas.Print("rawCalooutput.pdf[");
        canvas.Clear();
        h1->Draw();
        canvas.Print("rawCalooutput.pdf[");
        canvas.Clear();
        h2->Draw();
        canvas.Print("rawCalooutput.pdf");
        canvas.Print("rawCalooutput.pdf]");


}
