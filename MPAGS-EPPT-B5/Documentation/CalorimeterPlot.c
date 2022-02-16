#include "Reader.h"

int CalorimeterPlot()
{
  TH2F*h=new TH2F("h",":x:y",20,0,20,4,0,4);
  TH1F*h0=new TH1F("h0","",80,0,80);
  TH2F*hcal=new TH2F("hcal",":x:y",10,0,10,2,0,2);
  TH1F*hcal0=new TH1F("hcal0","",20,0,20);

  Reader reader;
  reader.GetEntry(0);

  double totalEnergy= 0;
  for(int i=0;i<80;i++)
    {
      std::cout<< i<< " "<<i/4+0.1 << " "<<i%4+0.1 << " "<<reader.ECEnergyVector->at(i)<<std::endl;
      h->SetBinContent(i/4+1, i%4+1, reader.ECEnergyVector->at(i));
      double energy = reader.ECEnergyVector->at(i);
      h0->SetBinContent(i+1, energy);
      totalEnergy+= energy;
    }
  TCanvas*c=new TCanvas("c","",600,600);
  
  h->Draw("colz");
  c->Print("test.pdf");
  std::cout<< h0->GetEntries()<< " "<< totalEnergy << " " << reader.ECEnergy<<std::endl;
  double totalEnergy2= 0;
  for(int i=0;i<20;i++)
    {
      std::cout<< i<< " "<<i/2+0.1 << " "<<i%2+0.1 << " "<<reader.HCEnergyVector->at(i)<<std::endl;
      hcal->SetBinContent(i/2+1,i%2+1, reader.HCEnergyVector->at(i));
      hcal0->SetBinContent(i+1, reader.HCEnergyVector->at(i));
      totalEnergy2+= reader.HCEnergyVector->at(i);
    }
  std::cout<< hcal0->GetEntries()<< " "<< totalEnergy2 << " " << reader.HCEnergy<<std::endl;
  TCanvas*cc=new TCanvas("cc","",600,600);
  hcal->Draw("colz");

  return 0;
}
