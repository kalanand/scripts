// RooMCStudy. Generate and fit 750 epxperiments. 
// Done only for signal pdf events here.


{
gROOT->Reset();
bool doNorm = false;

//calculate the dalitz boundary
Double_t md0 = 1.8645;
Double_t mpi0=0.1349766;
Double_t mpi = 0.13957;

Double_t upper=pow((md0-mpi),2);
Double_t lower=pow((mpi0+mpi),2);



RooRealVar S23("S23","S23",lower,upper); //This is pi+ pi0  rho+
RooRealVar S31("S31","S31",lower,upper); //This is pi- pi0  rho-

//Construct signal pdf
BdkPdfDDalitz mypdf("mypdf","mypdf",S23,S31,1);  //using pi0pi+ , pi0pi- as dalitz variable
if(doNorm==true) {               
  mypdf.parameters(); 
  BdkDDalitzAmp::normalizeAll(); 
  mypdf.parameters().writeToFile("pars.txt");  
}
else { mypdf.parameters().readFromFile("pars.txt"); }


RooAbsPdf* gen = mypdf.getPdf();


                       


RooArgSet allPars = mypdf.parameters();
RooRealVar* var1 = allPars.find("mypdf.pdf.dalitzAmp.Nonres_amp");
RooRealVar* var2 = allPars.find("mypdf.pdf.dalitzAmp.Nonres_phase");
RooRealVar* var3 = allPars.find("mypdf.pdf.dalitzAmp.Rho+_amp");
RooRealVar* var4 = allPars.find("mypdf.pdf.dalitzAmp.Rho+_phase");
RooRealVar* var5 = allPars.find("mypdf.pdf.dalitzAmp.Rho-_amp");
RooRealVar* var6 = allPars.find("mypdf.pdf.dalitzAmp.Rho-_phase");
RooRealVar* var7 = allPars.find("mypdf.pdf.dalitzAmp.Rho0_amp");
RooRealVar* var8 = allPars.find("mypdf.pdf.dalitzAmp.Rho0_phase");
var1->setConstant(false);
var2->setConstant(false);
var5->setConstant(false);
var6->setConstant(false);
var7->setConstant(false);
var8->setConstant(false);
var3->setVal(1.0);
var3->setConstant(true);
var4->setVal(0.0);
var4->setConstant(true);

var1->setRange("var1Range", 0., 10.0);
var2->setRange("var2Range", 0., 150.); 
var5->setRange("var5Range", 0.0, 10.0);
var6->setRange("var6Range", -150., 150.);
var7->setRange("var3Range", 0., 10.0);
var8->setRange("var4Range", -150., 150.);


RooAbsPdf* fit = mypdf.getPdf();



RooMCStudy toymc(*gen,*fit,RooArgSet(S23,S31),"","m");
//toymc.generate(1000,20000,true,"toymc/toymc_%03d.dat");

toymc.fit(750,"toymc/toymc_%03d.dat");






// toymc.generateAndFit(1000,20000);


RooPlot* pullplot1 = toymc.plotPull(*var1,-4,4,80,true);
RooPlot* pullplot2 = toymc.plotPull(*var2,-4,4,80,true);
RooPlot* pullplot3 = toymc.plotPull(*var7,-4,4,80,true);
RooPlot* pullplot4 = toymc.plotPull(*var8,-4,4,80,true);
RooPlot* pullplot5 = toymc.plotPull(*var5,-4,4,80,true);
RooPlot* pullplot6 = toymc.plotPull(*var6,-4,4,80,true);

pullplot1->SetTitle("Nonres amp");
pullplot2->SetTitle("Nonres phase");
pullplot3->SetTitle("Rho0 amp");
pullplot4->SetTitle("Rho0 phase");
pullplot5->SetTitle("Rho- amp");
pullplot6->SetTitle("Rho- phase");
pullplot1->GetXaxis()->SetTitle("Nonres amp Pull");
pullplot2->GetXaxis()->SetTitle("Nonres phase Pull");
pullplot3->GetXaxis()->SetTitle("Rho0 amp Pull");
pullplot4->GetXaxis()->SetTitle("Rho0 phase Pull");
pullplot5->GetXaxis()->SetTitle("Rho- amp Pull");
pullplot6->GetXaxis()->SetTitle("Rho- phase Pull");


TCanvas c1("c1","Pull distribution of fit parameters",1200,600);
c1.Divide(3,2);
c1.cd(1);pullplot1->Draw();
c1.cd(2);pullplot2->Draw();
c1.cd(3);pullplot3->Draw();
c1.cd(4);pullplot4->Draw();
c1.cd(5);pullplot5->Draw();
c1.cd(6);pullplot6->Draw();
c1.Update();
c1.SaveAs("toymc/mcStudyPull.ps");
c1.SaveAs("toymc/mcStudyPull.eps");
c1.Close();

RooPlot* errorplot1 = toymc.plotError(*var1,0.,0.1);
RooPlot* errorplot2 = toymc.plotError(*var2,0.,3.);
RooPlot* errorplot3 = toymc.plotError(*var7,0.,0.05);
RooPlot* errorplot4 = toymc.plotError(*var8,0.,3.);
RooPlot* errorplot5 = toymc.plotError(*var5,0.,0.05);
RooPlot* errorplot6 = toymc.plotError(*var6,0.,3.);

errorplot1->SetTitle("Nonres amp");
errorplot2->SetTitle("Nonres phase");
errorplot3->SetTitle("Rho0 amp");
errorplot4->SetTitle("Rho0 phase");
errorplot5->SetTitle("Rho- amp");
errorplot6->SetTitle("Rho- phase");
errorplot1->GetXaxis()->SetTitle("Nonres amp Error");
errorplot2->GetXaxis()->SetTitle("Nonres phase Error");
errorplot3->GetXaxis()->SetTitle("Rho0 amp Error");
errorplot4->GetXaxis()->SetTitle("Rho0 phase Error");
errorplot5->GetXaxis()->SetTitle("Rho- amp Error");
errorplot6->GetXaxis()->SetTitle("Rho- phase Error");


TCanvas c2("c2","Error distribution of fit parameters",1200,600);
c2.Divide(3,2);
c2.cd(1);errorplot1->Draw();
c2.cd(2);errorplot2->Draw();
c2.cd(3);errorplot3->Draw();
c2.cd(4);errorplot4->Draw();
c2.cd(5);errorplot5->Draw();
c2.cd(6);errorplot6->Draw();
c2.Update();
c2.SaveAs("toymc/mcStudyError.ps");
c2.SaveAs("toymc/mcStudyError.eps");
c2.Close();


RooPlot* meanplot1 = var1->frame(0.6,1.2);
toymc.plotParamOn(meanplot1);
RooPlot* meanplot2 = var2->frame(70,90);
toymc.plotParamOn(meanplot2);
RooPlot* meanplot3 = var7->frame(0.4,0.7);
toymc.plotParamOn(meanplot3);
RooPlot* meanplot4 = var8->frame(0.0,20.0);
toymc.plotParamOn(meanplot4);
RooPlot* meanplot5 = var5->frame(0.4,0.8);
toymc.plotParamOn(meanplot5);
RooPlot* meanplot6 = var6->frame(-6.,-2.);
toymc.plotParamOn(meanplot6);

meanplot1->SetTitle("Nonres amp");
meanplot2->SetTitle("Nonres phase");
meanplot3->SetTitle("Rho0 amp");
meanplot4->SetTitle("Rho0 phase");
meanplot5->SetTitle("Rho- amp");
meanplot6->SetTitle("Rho- phase");
meanplot1->GetXaxis()->SetTitle("Nonres amp Mean");
meanplot2->GetXaxis()->SetTitle("Nonres phase Mean");
meanplot3->GetXaxis()->SetTitle("Rho0 amp Mean");
meanplot4->GetXaxis()->SetTitle("Rho0 phase Mean");
meanplot5->GetXaxis()->SetTitle("Rho- amp Mean");
meanplot6->GetXaxis()->SetTitle("Rho- phase Mean");

TCanvas c3("c3","Distribution of mean-values of fit parameters",1200,600);
c3.Divide(3,2);
c3.cd(1);meanplot1->Draw();
c3.cd(2);meanplot2->Draw();
c3.cd(3);meanplot3->Draw();
c3.cd(4);meanplot4->Draw();
c3.cd(5);meanplot5->Draw();
c3.cd(6);meanplot6->Draw();
c3.Update();
c3.SaveAs("toymc/mcStudyMean.ps");
c3.SaveAs("toymc/mcStudyMean.eps");
c3.Close();


TCanvas canvas("canvas", "", 880, 680); 
pullplot1->Draw();
canvas.SaveAs("toymc/mcStudyPull-1.ps");
canvas.SaveAs("toymc/mcStudyPull-1.eps");
canvas.Clear();
pullplot2->Draw();
canvas.SaveAs("toymc/mcStudyPull-2.ps");
canvas.SaveAs("toymc/mcStudyPull-2.eps");
canvas.Clear();
pullplot3->Draw();
canvas.SaveAs("toymc/mcStudyPull-3.ps");
canvas.SaveAs("toymc/mcStudyPull-3.eps");
canvas.Clear();
pullplot4->Draw();
canvas.SaveAs("toymc/mcStudyPull-4.ps");
canvas.SaveAs("toymc/mcStudyPull-4.eps");
canvas.Clear();
pullplot5->Draw();
canvas.SaveAs("toymc/mcStudyPull-5.ps");
canvas.SaveAs("toymc/mcStudyPull-5.eps");
canvas.Clear();
pullplot6->Draw();
canvas.SaveAs("toymc/mcStudyPull-6.ps");
canvas.SaveAs("toymc/mcStudyPull-6.eps");
canvas.Clear();

errorplot1->Draw();
canvas.SaveAs("toymc/mcStudyError-1.ps");
canvas.SaveAs("toymc/mcStudyError-1.eps");
canvas.Clear();
errorplot2->Draw();
canvas.SaveAs("toymc/mcStudyError-2.ps");
canvas.SaveAs("toymc/mcStudyError-2.eps");
canvas.Clear();
errorplot3->Draw();
canvas.SaveAs("toymc/mcStudyError-3.ps");
canvas.SaveAs("toymc/mcStudyError-3.eps");
canvas.Clear();
errorplot4->Draw();
canvas.SaveAs("toymc/mcStudyError-4.ps");
canvas.SaveAs("toymc/mcStudyError-4.eps");
canvas.Clear();
errorplot5->Draw();
canvas.SaveAs("toymc/mcStudyError-5.ps");
canvas.SaveAs("toymc/mcStudyError-5.eps");
canvas.Clear();
errorplot6->Draw();
canvas.SaveAs("toymc/mcStudyError-6.ps");
canvas.SaveAs("toymc/mcStudyError-6.eps");
canvas.Clear();

meanplot1->Draw();
canvas.SaveAs("toymc/mcStudyMean-1.ps");
canvas.SaveAs("toymc/mcStudyMean-1.eps");
canvas.Clear();
meanplot2->Draw();
canvas.SaveAs("toymc/mcStudyMean-2.ps");
canvas.SaveAs("toymc/mcStudyMean-2.eps");
canvas.Clear();
meanplot3->Draw();
canvas.SaveAs("toymc/mcStudyMean-3.ps");
canvas.SaveAs("toymc/mcStudyMean-3.eps");
canvas.Clear();
meanplot4->Draw();
canvas.SaveAs("toymc/mcStudyMean-4.ps");
canvas.SaveAs("toymc/mcStudyMean-4.eps");
canvas.Clear();
meanplot5->Draw();
canvas.SaveAs("toymc/mcStudyMean-5.ps");
canvas.SaveAs("toymc/mcStudyMean-5.eps");
canvas.Clear();
meanplot6->Draw();
canvas.SaveAs("toymc/mcStudyMean-6.ps");
canvas.SaveAs("toymc/mcStudyMean-6.eps");
canvas.Clear();
canvas.Close();

}  //end the macro






















