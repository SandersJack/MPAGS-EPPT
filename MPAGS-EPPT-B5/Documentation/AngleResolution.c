int AngleResolution()
{
  TRandom3*rndm=new TRandom3(0);
  TH1F*h1=new TH1F("h1","",1000,-1.e-3,1.e-3);
  TH1F*h2=new TH1F("h2","",1000,-1.e-3,1.e-3);
  TH1F*h3=new TH1F("h3","",1000,-1.e-3,1.e-3);
  TH1F*h4=new TH1F("h4","",1000,-1.e-3,1.e-3);
  for(int iExp=0;iExp<100000;iExp++)
    {
      double x[5]={0.,0.5,1.,1.5,2.};
      double y[5]={0.};
      double ex[5]={0.};
      double ey[5]={0.1e-3,0.1e-3,0.1e-3,0.1e-3,0.1e-3};
      for(int i=0;i<5;i++)
	y[i]+=rndm->Gaus(0.,0.1e-3);// 100mum resolution

      double x1[4]={0.,0.5,1.5,2.};
      double y1[4]={y[0],y[1],y[3],y[4]};
      double ex1[4]={0.};
      double ey1[4]={0.1e-3,0.1e-3,0.1e-3,0.1e-3};

      TGraphErrors*gr=new TGraphErrors(5,x,y,ex,ey);
      gr->Fit("pol1","q");
      TGraphErrors*gr1=new TGraphErrors(4,x1,y1,ex1,ey1);
      gr1->Fit("pol1","q");
      h1->Fill(gr->GetFunction("pol1")->GetParameter(1));
      h2->Fill((y[4]-y[0])/(x[4]-x[0]));
      h3->Fill((y[3]-y[1])/(x[3]-x[1]));
      h4->Fill(gr1->GetFunction("pol1")->GetParameter(1));
      
    }
  TCanvas*c1=new TCanvas("c1","",600,600);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h1->Fit("gaus");
  h2->Fit("gaus","","same");
  h3->Fit("gaus","","same");
  h4->Fit("gaus","","same");
  return 0;
}
