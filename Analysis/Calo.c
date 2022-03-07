#include <stddef.h>

void Calo() {
	TFile *f = new TFile("rootfiles/newSim/protons.root");
        f->ls();
        TTree *t = (TTree*)f->Get("B5");
        Int_t nentries = t->GetEntries();
	t->Print();

	std::vector<double> *ECEnergy = 0;
	TBranch *EC = 0;

	std::vector<double> *HCEnergy = 0;
	TBranch *HC = 0;

	double *trueEC = 0;
        TBranch *tEC = 0;

	double V_EC[nentries][80];
	double V_HC[nentries][20];

	TH1F *h11 = new TH1F("h1","EC",111,97.1,103.1);
        h11->GetXaxis()->SetTitle("E [GeV]");
	TH1F *h1 = new TH1F("h1","EC",111,0,110);
	h1->GetXaxis()->SetTitle("E [GeV]");
	TH2F *h2 = new TH2F("h2","EC",81,-0.1,80.1,111,-0.1,61);
	h2->GetXaxis()->SetTitle("Channel ID");
	h2->GetYaxis()->SetTitle("E [GeV]");	
	auto *h3 = new TH3F("h3","EC Energy in Cells",20, -0.1, 19.1, 4, -0.1,3.1, 111, -0.1, 110.1); 
	h3->GetXaxis()->SetTitle("Column");
	h3->GetYaxis()->SetTitle("Row");
	h3->GetZaxis()->SetTitle("E [GeV]");
	TH1F *h4 = new TH1F("h4","HC",51,-0.1,10);
	h4->GetXaxis()->SetTitle("E [GeV]");
        TH2F *h5 = new TH2F("h5","HC",21,-0.1,20.1,101,-0.1,10);
        h5->GetXaxis()->SetTitle("Channel ID");
        h5->GetYaxis()->SetTitle("E [GeV]");
	auto *h6 = new TH3F("h6","HC Energy in Cells",10, -0.1, 9.1, 2, -0.1,1.1, 111, -0.1, 10.1);
	h6->GetXaxis()->SetTitle("Column");
        h6->GetYaxis()->SetTitle("Row");
        h6->GetZaxis()->SetTitle("E [GeV]");
	TH1F *h45 = new TH1F("h45","Total Energy",111,50,151.1);
	TH1F *h44 = new TH1F("h44","Beta Const",111,0,110.1);

	t->SetBranchAddress("ECEnergyVector",&ECEnergy, &EC);
	t->SetBranchAddress("ECEnergy",&trueEC, &tEC);
//	std::cout << ECEnergy << std::endl;
	double sum_EC[nentries];

	for (Int_t i=0; i<nentries; i++) {
                Long64_t tentry = t->LoadTree(i);
		EC->GetEntry(tentry);
		tEC->GetEntry(tentry);
		sum_EC[i] = 0;
		//std::cout << ECEnergy->size() << std::endl;
		for (UInt_t t=0; t<ECEnergy->size(); t++) {
			V_EC[i][t] = ECEnergy->at(t);
			if (V_EC[i][t] != 0) {
				sum_EC[i] += V_EC[i][t];
				//h1->Fill(V_EC[i][t]);
				h2->Fill(t,V_EC[i][t]/1000);
			}
		}
		//cout << sum_EC << endl;
		h1->Fill(sum_EC[i]/1000);
		h11->Fill(tEC->GetEntry(tentry));
		
	}
	//h1->Fit("gaus");
	double sum_HC[nentries];
	t->SetBranchAddress("HCEnergyVector",&HCEnergy, &HC);
        //std::cout << ECEnergy << std::endl;
        for (Int_t i=0; i<nentries; i++) {
                Long64_t tentry = t->LoadTree(i);
                HC->GetEntry(tentry);
		sum_HC[i] = 0;
		for (UInt_t t=0; t<HCEnergy->size(); t++) {
                        V_HC[i][t] = HCEnergy->at(t);
                        if (V_HC[i][t] != 0) {
				sum_HC[i] += V_HC[i][t];
                                //h4->Fill(V_HC[i][t]);
                                h5->Fill(t,V_HC[i][t]/1000);
                        }
                }
		h4->Fill(sum_HC[i]/1000);
		std::cout << sum_HC[i] << std::endl;
        }
	//h4->Fit("gaus");
	//

	double beta[nentries];
	double betasum = 0;
	for (Int_t f=0; f<nentries; f++) {
		beta[f] = 25; //Starting to help
		for (Int_t d=0; d<10000; d++) {
			double E = sum_EC[f]*1.005/1000 + sum_HC[f]*beta[f]/1000;
			double chi = 100.004399-E;
			if ( chi < .01 && chi > -0.1) {
				break;
			} else if ( chi >= .01) {
				beta[f] += 0.1;
			} else if ( chi <= -0.1) {
				beta[f] -= 0.1;
			}
			//cout << E << endl;
		}
		betasum += beta[f];
		h44->Fill(beta[f]);
	}

	double betaavg = 30.87;
	double sumTot[nentries];
	for (Int_t v=0; v<nentries; v++) {	
		sumTot[v] = 0;
		sumTot[v] = sum_EC[v] + sum_HC[v]*betaavg;
		cout << sum_HC[v] << endl;
		h45->Fill(sumTot[v]/1000);
	}

	for (Int_t t=0; t<80; t++) {
		double col = floor(t/4);
		int row = t-(4*col);
		//std::cout << col << "," << row <<std::endl;
		h3->Fill(col,row,V_EC[20][t]/1000); 
	}

	for (Int_t t=0; t<20; t++) {
                double col = floor(t/2);
                int row = t-(2*col);
                std::cout << col << "," << row <<std::endl;
                h6->Fill(col,row,V_HC[20][t]/1000);
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
        canvas.Clear();
        h45->Draw();
        canvas.Print("rawCalooutput.pdf");
	canvas.Clear();
        h44->Draw();
        canvas.Print("rawCalooutput.pdf");
	canvas.Clear();
        h11->Draw();
        canvas.Print("rawCalooutput.pdf");
	canvas.Print("rawCalooutput.pdf]");


}
