const float IntLUMI = 19.3;
RooRealVar *mjj_;
using namespace RooFit;


void fitWWa(){

  TFile* fin 
    = new TFile("/uscms_data/d2/jdamgov/NTuples_53X_MORIOND13/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test/Utilities/fracFit/mu_a0W_WWA_PhotonEt.root ");


  TH1D* data_obs = (TH1D*) fin->Get("data_obs");
  TH1D* th1fkdata = (TH1D*) fin->Get("th1fkdata");
  TH1D* th1wwa = (TH1D*) fin->Get("th1wwa");
  TH1D* th1wza = (TH1D*) fin->Get("th1wza");
  TH1D* th1zz = (TH1D*) fin->Get("th1zz");
  TH1D* th1wajets = (TH1D*) fin->Get("th1wajets");
  TH1D* th1zajets = (TH1D*) fin->Get("th1zajets");
  TH1D* th1Top = (TH1D*) fin->Get("th1Top");

  int NbinsX = data_obs->GetNbinsX();
  double xmin = data_obs->GetXaxis()->GetBinLowEdge(1);
  double xmax = data_obs->GetXaxis()->GetBinLowEdge(NbinsX+1);

  data_obs->GetXaxis()->SetRangeUser(xmin, xmax);
  th1fkdata->GetXaxis()->SetRangeUser(xmin, xmax);
  th1wwa->GetXaxis()->SetRangeUser(xmin, xmax);
  th1wza->GetXaxis()->SetRangeUser(xmin, xmax);
  th1zz->GetXaxis()->SetRangeUser(xmin, xmax);
  th1wajets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1zajets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1Top->GetXaxis()->SetRangeUser(xmin, xmax);


  mjj_ = new RooRealVar( "Mjj", "m_{jj}", xmin, xmax, "GeV");
  RooRealVar Mass = *mjj_;
  RooDataHist* data = new RooDataHist("data","data", *mjj_, data_obs);

  RooHistPdf* pdffk = makePdf(th1fkdata, "pdffk");
  RooHistPdf* pdfwwa = makePdf(th1wwa, "pdfwwa");
  RooHistPdf* pdfzz = makePdf(th1zz, "pdfzz");
  RooHistPdf* pdfwajets = makePdf(th1wajets, "pdfwajets");
  RooHistPdf* pdfzajets = makePdf(th1zajets, "pdfzajets");
  RooHistPdf* pdfTop = makePdf(th1Top, "pdfTop");

  double fkNorm = th1fkdata->Integral();
  double wwaNorm = th1wwa->Integral() + th1wza->Integral();
  double zzNorm = th1zz->Integral();
  double wajetsNorm = th1wajets->Integral();
  double zajetsNorm = th1zajets->Integral();
  double TopNorm = th1Top->Integral();


  RooRealVar nfk("nfk","nfk",                 fkNorm,     0.0,   1000.);
  RooRealVar nwwa("nwwa","nwwa",              wwaNorm);
  RooRealVar nzz("nzz","nzz",                 zzNorm);
  RooRealVar nwajets("nwajets","nwajets",     400.0,     0.0,   10000.);
  RooRealVar nzajets("nzajets","nzajets",     zajetsNorm);
  RooRealVar nTop("nTop","nTop",              TopNorm);


  RooArgList* components 
    = new RooArgList(*pdffk, *pdfwwa, *pdfzz, *pdfwajets, *pdfzajets, *pdfTop);
  RooArgList* yields = new RooArgList(nfk, nwwa, nzz, nwajets, nzajets, nTop);


  RooAddPdf totalPdf("totalPdf","extended sum pdf", *components, *yields);

  RooGaussian consNfk("consNfk","", nfk, RooConst(fkNorm),RooConst(0.2*fkNorm)) ;
  RooGaussian consNwwa("consNwwa","", nwwa, RooConst(wwaNorm),RooConst(0.28*wwaNorm)) ;
  RooGaussian consNzz("consNzz","", nzz, RooConst(zzNorm),RooConst(0.1*zzNorm)) ;
  RooGaussian consNzajets("consNzajets","", nzajets, RooConst(zajetsNorm),RooConst(0.22*zajetsNorm)) ;
  RooGaussian consNTop("consNTop","", nTop, RooConst(TopNorm),RooConst(0.21*TopNorm)) ;


  RooFitResult *fitResult 
    = totalPdf.fitTo(*data, Save(true), 
		     ExternalConstraints(consNfk),
		     //ExternalConstraints(consNwwa),
		     //ExternalConstraints(consNzz),
		     //ExternalConstraints(consNzajets),
		     //ExternalConstraints(consNTop),
		     RooFit::Extended(true), 
		     //RooFit::Minos(true), 
		     //RooFit::Hesse(false),
		     //PrintEvalErrors(-1),
		     // RooFit::Range(rangeString),
		     Warnings(false) 
		     );

   fitResult->Print("v");


   std::cout << "===================== Wa+jets k-factor = " << 
     nwajets.getVal() / wajetsNorm << "  +- "  << nwajets.getError() / wajetsNorm << std::endl;


   // ********** Make and save Canvas for the plots ********** //
   gROOT->ProcessLine(".L ~kalanand/tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.19);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);
  RooAbsData::ErrorType errorType = RooAbsData::SumW2;



   TCanvas* c = new TCanvas("fit","",500,500);
   RooPlot* frame1 = Mass.frame();
   data->plotOn(frame1,RooFit::DataError(errorType), Name("h_data"));
   totalPdf.plotOn(frame1,ProjWData(*data), Name("h_total"));
   totalPdf.plotOn(frame1,ProjWData(*data),Components("pdfwajets"),
		     LineColor(kRed), LineStyle(2), Name("h_wajets"));

   totalPdf.plotOn(frame1,ProjWData(*data),Components("pdffk,pdfwwa, pdfzz,pdfzajets,pdfTop"),
		     LineColor(kBlack), LineStyle(2), Name("h_others"));
   totalPdf.plotOn(frame1,ProjWData(*data));

   frame1->SetMinimum(0);
   frame1->SetMaximum(1.35* frame1->GetMaximum());
   frame1->Draw("e0");


   std::cout << "===================== chi2/ dof = " << frame1->chiSquare() << std::endl;

   TPaveText *plotlabel4 = new TPaveText(0.25,0.66,0.5,0.81,"NDC");
   plotlabel4->SetTextColor(kBlack);
   plotlabel4->SetFillColor(kWhite);
   plotlabel4->SetBorderSize(0);
   plotlabel4->SetTextAlign(12);
   plotlabel4->SetTextSize(0.04);
   char temp[50];
   sprintf(temp, "#chi^{2} / dof = %.2f", frame1->chiSquare());
   plotlabel4->AddText(temp);
   plotlabel4->Draw();

   cmsPrelim2();

   TLegend* legend = new TLegend(0.55,0.72,0.88,0.91);
   RooHist* datahist = frame1->getHist("h_data");
   RooCurve* totalhist = frame1->getCurve("h_total");
   RooCurve* wjetshist = frame1->getCurve("h_wajets");
   RooCurve* otherhist = frame1->getCurve("h_others");

   legend->AddEntry( datahist, "Data", "PE");  
   legend->AddEntry( totalhist, "Fit", "L");
   legend->AddEntry( wjetshist, "W#gamma+jets", "L");
   legend->AddEntry( otherhist, "Other processes", "L");
   legend->SetFillColor(0);
   legend->Draw();
   c->SaveAs( "WVa_WjetsKfactorFit.png");
   c->SaveAs( "WVa_WjetsKfactorFit.pdf");

}




RooHistPdf*  makePdf(TH1D* hist, char* name) {

  RooDataHist* rdh = new RooDataHist("rdh","", *mjj_, hist);
  RooHistPdf* pdf = new RooHistPdf(name,name,*mjj_,*rdh);
  return pdf;
}




////CMS Preliminary label and lumu -- upper left corner
void cmsPrelim2()
{
   const float LUMINOSITY = IntLUMI;
   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);

   latex.SetTextAlign(31); // align right
   latex.DrawLatex(0.90,0.96,"#sqrt{s} = 8 TeV");
   if (LUMINOSITY > 0.) {
      latex.SetTextAlign(11); // align left
      latex.DrawLatex(0.21,0.85,Form("#int #font[12]{L} dt = %.1f fb^{-1}", LUMINOSITY));
   }
   latex.SetTextAlign(11); // align left
   latex.DrawLatex(0.18,0.96,"CMS preliminary");
}

