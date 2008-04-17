/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: 
 *    File: $Id: BdkBatchMCStudy.cc,v 1.20 2007/08/22 07:59:04 fwinkl Exp $
 * Authors:
 *   Frank Winklmeier
 * Description:
 *   Extension of RooMCStudy that supports batch queues
 *
 * Copyright (C) 2006 Colorado State University
 *****************************************************************************/

#include <fstream>
#include <sys/time.h>

#include "BToDKTo3piK/BdkBatchMCStudy.hh"
#include "BToDKTo3piK/BdkPdfDKDalitz.hh"
#include "BToDKTo3piK/BdkMath.hh"

#include "RooFitCore/RooAbsArg.hh"
#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooSimGenContext.hh"
#include "RooFitCore/RooGenContext.hh"
#include "RooFitCore/RooSimultaneous.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooFitResult.hh"

#include "TKey.h"
#include "TList.h"
#include "TMatrix.h"

using namespace std;

ClassImp(BdkBatchMCStudy)


BdkBatchMCStudy::BdkBatchMCStudy(const RooAbsPdf& model, const RooArgSet& observables,
                                 RooCmdArg arg1, RooCmdArg arg2,
                                 RooCmdArg arg3,RooCmdArg arg4,RooCmdArg arg5,
                                 RooCmdArg arg6,RooCmdArg arg7,RooCmdArg arg8) :
  RooMCStudy(model, observables, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
{
  cout << "This constructor is currently not supported!"<<endl;
  assert(0);
}

BdkBatchMCStudy::BdkBatchMCStudy(const RooAbsPdf& genModel, const RooAbsPdf& fitModel, 
                                 const RooArgSet& dependents, const char* genOptions, 
                                 const char* fitOptions, const RooDataSet* genProtoData, 
                                 const RooArgSet& projDeps) :
  RooMCStudy(genModel, fitModel, dependents, genOptions, fitOptions, genProtoData, projDeps),
  _onResPdf(0),
  _fitResults(0),
  _submitOption(""),
  _useRandomCP(kFALSE),
  _seed(0)
{
  _genOptions = genOptions;   // this is not stored in base class
}

// Use "p" in fitOptions to only use the yields penalty term for fitting
// Use "s" in fitOptions to only use the shape/Dalitz term for fitting
// Use "x" in fitOptions to save additional information, see saveExtInfo()
// Use "m" in fitOptions to skip the MINOS step
BdkBatchMCStudy::BdkBatchMCStudy(const BdkPdfOnRes& genModel, const BdkPdfOnRes& fitModel, 
                                 const RooArgSet& dependents, const char* genOptions, 
                                 const char* fitOptions, const RooDataSet* genProtoData, 
                                 const RooArgSet& projDeps) :
  RooMCStudy(*((BdkPdfOnRes&)genModel).getPdf(), 
             *((BdkPdfOnRes&)fitModel).getPdf(), 
             dependents, genOptions, fitOptions, genProtoData, projDeps),
  _fitResults(0),
  _submitOption(""),
  _useRandomCP(kFALSE),
  _seed(0)
{
  _genOptions = genOptions;
  _onResPdf = &((BdkPdfOnRes&)fitModel);

  // Delete objects created by RooMCStudy constructor
  delete _fitParams;
  delete _fitInitParams;
  delete _fitParData;

  // Initialize them correctly from floating parameters only
  _fitParams = new RooArgSet(((BdkPdfOnRes&)fitModel).parametersFree());
  _fitInitParams = (RooArgSet*) _fitParams->snapshot(kTRUE) ;

  // Add some additional variables
  _xyNllVar = new RooRealVar("xyNLL","-log(Likelihood)",0);

  // Put all in fit parameter dataset
  RooArgSet tmp(*_fitParams);
  tmp.add(*_nllVar);
  tmp.add(*_xyNllVar);
  tmp.setAttribAll("StoreError",kTRUE);
  tmp.setAttribAll("StoreAsymError",kTRUE) ;
  _fitParData = new RooDataSet("fitParData","Fit Parameters DataSet",tmp);
  tmp.setAttribAll("StoreError",kFALSE) ;
  tmp.setAttribAll("StoreAsymError",kFALSE) ;  
}



BdkBatchMCStudy::~BdkBatchMCStudy()
{
}


void BdkBatchMCStudy::generateBatch(const char* mcFile, Int_t nSamples, Int_t nEvtPerSample)
{
  
  //  generate(nSamples, nEvtPerSample, kTRUE);

  TFile f(mcFile,"RECREATE");

  // Save the parameters used to generate the MC to the file  
  //  _genModel->getParameters(RooArgSet())->Write();
  
  // Save toy MC data to file
  
  RooAbsGenContext* genContext;
  if (_genModel->InheritsFrom(RooSimultaneous::Class()))
    genContext = new RooSimGenContext(*((RooSimultaneous*)_genModel),_dependents,_genProtoData);
  else
    genContext = new RooGenContext(*_genModel,_dependents,_genProtoData);
  
  for (int i=0; i<nSamples; i++) {
    //    genData(i)->Write();
    cout << "Generating MC sample "<<i<<endl;
    //RooDataSet *data = _genModel->generate(_dependents,*_genProtoData,nEvtPerSample);
    RooDataSet *data = genContext->generate(nEvtPerSample);

    data->Write();
    delete data;
  }
  f.Close();
  delete genContext;
}


Bool_t BdkBatchMCStudy::init(const char* mcFile)
{
  _mcFile = mcFile;
  if (!mcFile) return kFALSE;

  TFile f(mcFile);
  if (f.IsZombie()) return kFALSE;

  // Read parameters used to generate the samples
  RooArgSet *genParams = (RooArgSet*)f.Get("parameters");
  if (!genParams) return kFALSE;

  // This is important for correct pulls
  _genParams = (RooArgSet*)genParams->snapshot(kFALSE);

  f.Close();
  return kTRUE;
}


TList* getDataSetKeys(TFile &f)
{
  TList *mcList = new TList();

  // Make a list of all RooDataSet in the file
  TList *keyList = f.GetListOfKeys();
  TIter nextKey(keyList);
  while (TKey *key = (TKey*)nextKey()) {
    TString className(key->GetClassName());
    if (className!="RooDataSet") continue;   // no other objects
    mcList->Add(key->Clone());     // add the key to the list
  }

  return mcList;
}

Bool_t BdkBatchMCStudy::fitBatch(const char* fitFile, Int_t firstSample, Int_t nSamples)
{

  if (_onResPdf) {
    cout << "BdkBatchMCStudy::fitBatch() not supported for BdkPdfOnRes."<<endl;
    return kFALSE;
  }
                  
  TFile f(_mcFile);
  if (f.IsZombie()) return kFALSE;

  TList *mcList = getDataSetKeys(f);

  if (nSamples<0) nSamples = mcList->GetEntries();

  // Now read the required RooDataSets from the file
  TList fitList;
  for (Int_t i=firstSample; i<firstSample+nSamples; i++) {
    fitList.Add(((TKey*)mcList->At(i))->ReadObj()->Clone());
  }
  f.Close();   // object are cloned so we don't need the file anymore


  // Fit the selected samples
  fit(nSamples, fitList);

  // Write fit results to file
  TFile f2(fitFile,"RECREATE");
  if (f2.IsZombie()) return kFALSE;

  fitParDataSet().Write();
  f2.Close();

  return kTRUE;
}


Bool_t BdkBatchMCStudy::generateAndFitBatch(const char* fitFile, Int_t nEvents, Int_t nSamples)
{
  if (_seed==0) {
    // Set random seed to system time
    struct timeval tv;
    gettimeofday(&tv, 0);
    _seed = tv.tv_usec;
  }
  RooRandom::randomGenerator()->SetSeed(_seed);
  cout <<"Setting random seed to "<<_seed<<endl;

  // for RooAbsPdf we can use the base class function
  if (_onResPdf==0) generateAndFit(nSamples,nEvents,kFALSE);

  // for BdkPdfOnRes we need our own setup
  else {
    onResGenerateAndFit(nEvents, nSamples);
  }

  // Write fit results to file
  TFile f2(fitFile,"RECREATE");
  if (f2.IsZombie()) return kFALSE;

  fitParDataSet().Write();
  if (_fitResults) {
    for (int i=0; i<_fitResults->GetEntries(); i++)
      _fitResults->At(i)->Write();
  }
    
  f2.Close();

  return kTRUE;
}


// mostly copied from run() and fitSample() in RooMCStudy.cc
// If nEvtPerSample<0 call BdkPdfOnRes::totalNumEvts() to get sample size
Bool_t BdkBatchMCStudy::onResGenerateAndFit(Int_t nEvtPerSample, Int_t nSamples)
{
  while(nSamples--) {
    
    cout << "BdkBatchMCStudy::run: "
         << "Generating and fitting sample " << nSamples << endl ;
    
    //
    // Generate sample
    //
    RooDataSet* genSample = 0;
    
    Int_t nEvt(nEvtPerSample) ;
        
    if (_randProto && _genProtoData && _genProtoData->numEntries()!=nEvt) {
      cout << "BdkBatchMCStudy: (Re)randomizing event order not supported." << endl;
      return kFALSE;
    }

    // Reset all fit parameters to their initial values  
    *_fitParams = *_fitInitParams ;

    // If requested, randomize CP parameters for each experiment
    // Pulls will be wrong in this case !!!
    if (_useRandomCP) setRandomCP();
    
    // Now, we can get the total number of events to generate if not specified
    if (nEvt<0) nEvt = (Int_t)round(_onResPdf->totalNumEvts());
    
    // For some reason we can't use the generator context.
    // Probably because the PDF gets re-created after useXXX() calls.
    //    genSample = _genContext->generate(nEvt) ;

    // Make sure we generate data for all variables
    _onResPdf->useDalitz();
    _onResPdf->useNnComb();
    _onResPdf->useNnCont();
    _onResPdf->useDE();
    _onResPdf->getPdf();

    // Use the PDF itself to generate
    gROOT->cd();
    //    _onResPdf->setNsigAsymFromXY();

    genSample = _onResPdf->generate(nEvt,_extendedGen);

    //    genSample->tree().SetScanField(0);
    //    genSample->Scan();
    
    //
    // Fit sample
    //
    
    if (_fitOptions.Contains("m")) _onResPdf->useMinos(false);
    else _onResPdf->useMinos();

    Bool_t usePenalty = kTRUE;
    Bool_t usePenaltyOnly = kFALSE;
    if (_fitOptions.Contains("p")) usePenaltyOnly = kTRUE;
    if (_fitOptions.Contains("s")) usePenalty = kFALSE;

    Bool_t fitOK = _onResPdf->fit(*genSample,usePenalty,usePenaltyOnly);
    if (fitOK) {
      if (_onResPdf->yieldFitResult()->status()==0 &&
          _onResPdf->xyFitResult()->status()==0) {

        // Store NLL's and fit parameters
        _nllVar->setVal(_onResPdf->yieldFitResult()->minNll());
        _xyNllVar->setVal(_onResPdf->xyFitResult()->minNll());

        RooArgSet tmp(*_fitParams);
        tmp.add(*_nllVar);
        tmp.add(*_xyNllVar);

        // Save additional information if requested
        if (_fitOptions.Contains("x")) saveExtInfo(tmp);
        else _fitParData->add(tmp);
      }
          
    }
    else cout << "Fit failed."<<endl;

    delete genSample ;
  }

  _canAddFitResults = kFALSE ;
  calcPulls() ;
  return kFALSE ;
}


void BdkBatchMCStudy::setRandomCP()
{  
  if (_onResPdf->rhoPlus()) {    
    // Generate rB, gamma and delta
    _rB = RooRandom::randomGenerator()->Uniform(0, 0.4);
    _gamma = RooRandom::randomGenerator()->Uniform(0, 180);
    _delta = RooRandom::randomGenerator()->Uniform(-180, 180);
  
    // Set rho/theta parameters (gamma phase flip is done automatically)
    ((BdkPdfDKDalitz*)_onResPdf->sigGoodD0N().getDalitzPdf())->setCPparams(_rB, _gamma, _delta);
    ((BdkPdfDKDalitz*)_onResPdf->sigGoodD0P().getDalitzPdf())->setCPparams(_rB, _gamma, _delta);
  
    // Save these values for later storage
    _rhoPgen = _onResPdf->rhoPlus()->getVal();
    _rhoNgen = _onResPdf->rhoMinus()->getVal();
    _thetaPgen = _onResPdf->thetaPlus()->getVal();
    _thetaNgen = _onResPdf->thetaMinus()->getVal();
  }
  else {
    _xPgen = RooRandom::randomGenerator()->Uniform(-0.4, 0.4);
    _xNgen = RooRandom::randomGenerator()->Uniform(-0.4, 0.4);
    _yPgen = RooRandom::randomGenerator()->Uniform(-0.4, 0.4);
    _yNgen = RooRandom::randomGenerator()->Uniform(-0.4, 0.4);
                
    _onResPdf->xPlus()->setVal(_xPgen);
    _onResPdf->xMinus()->setVal(_xNgen);
    _onResPdf->yPlus()->setVal(_yPgen);
    _onResPdf->yMinus()->setVal(_yNgen);
    cout << "BdkBatchMCStudy: Setting cartesian coordinates: "
         << "x+=" << _xPgen << ", "
         << "x-=" << _xNgen << ", "
         << "y+=" << _yPgen << ", "
         << "y-=" << _yNgen << endl;
      
  }
  _onResPdf->setNsigAsymFromXY();
}


// Save initial floating parameters and correlation matrix of x/y fit
void BdkBatchMCStudy::saveExtInfo(const RooArgSet& parSet)
{
  RooFitResult* r = _onResPdf->xyFitResult();

  // Store CP parameters used to generate data
  RooArgList* pars = new RooArgList(); // (RooArgList*)r->floatParsInit().snapshot(kTRUE);

  if (_onResPdf->rhoPlus()) {
    pars->addOwned(*(new RooRealVar("rhoPgen","rhoPgen",_rhoPgen)));
    pars->addOwned(*(new RooRealVar("rhoNgen","rhoNgen",_rhoNgen)));
    pars->addOwned(*(new RooRealVar("thetaPgen","thetaPgen",_thetaPgen)));
    pars->addOwned(*(new RooRealVar("thetaNgen","thetaNgen",_thetaNgen)));
  }
  else {
    pars->addOwned(*(new RooRealVar("xPgen","xPgen",_xPgen)));
    pars->addOwned(*(new RooRealVar("xNgen","xNgen",_xNgen)));
    pars->addOwned(*(new RooRealVar("yPgen","yPgen",_yPgen)));
    pars->addOwned(*(new RooRealVar("yNgen","yNgen",_yNgen)));
  }
  
  pars->addOwned(*(new RooRealVar("rB","rB",_rB)));
  pars->addOwned(*(new RooRealVar("gamma","gamma",_gamma)));
  pars->addOwned(*(new RooRealVar("delta","delta",_delta)));
    
  // Add covariance matrix
  TMatrix cov, cor;
  BdkMath::getCovCorMatrix(r, cov, cor);

  for (int i=0; i<cov.GetNcols(); i++) {
    for (int j=0; j<=i; j++) {
      TString name("cov_");
      name += i;
      name += j;
      RooRealVar* mij = new RooRealVar(name,name,cov(i,j));
      pars->addOwned(*mij);
    }
  }
  
  pars->setAttribAll("StoreError",kFALSE);
  pars->setAttribAll("StoreAsymError",kFALSE);
  
  // Add new columns if this is the first entry
  if (_fitParData->numEntries()==0) _fitParData->addColumns(*pars);

  // Add our parameters to existing ones
  RooArgSet tmp(parSet);
  tmp.add(*pars);

  // Save all parameters in dataset
  _fitParData->add(tmp);

  // Cleanup all RooAbsArgs created in this function
  delete pars;
}


// Convert a RooArgSet into its constructing string
TString BdkBatchMCStudy::argSetToString(RooArgSet &set)
{
  TString s = "RooArgSet(";
  TIterator *iter = set.createIterator();
  Bool_t firstArg = kTRUE;
  while (RooAbsArg *arg = (RooAbsArg*)iter->Next()) {
    if (!firstArg) s += ",";
    firstArg = kFALSE;
    s += arg->GetName();
  }
  s += ")";
  return s;
}


// create one fit script
void BdkBatchMCStudy::writeScript(const char* mcFile, const char* fitFile, 
                                  const char* scriptFile,
                                  const char* setupScript,
                                  Int_t firstSample, Int_t nSamples, Int_t nEvents)
{
  ofstream of;
  of.open(scriptFile);

  of <<"{"<<endl;

  // run init script
  of <<"gROOT->ProcessLine(\".x "<<setupScript<<"\");"<<endl;

  // create BdkBatchMCStudy object
  of <<"BdkBatchMCStudy *toyMC = new BdkBatchMCStudy("
     <<_genModel->GetName()<<","
     <<_fitModel->GetName()<<","
     <<argSetToString(_dependents);

  // generator, fit options
  of <<",\""<<_genOptions<<"\",\""<<_fitOptions<<"\",";

  // prototype dataset
  if (_genProtoData) of << _genProtoData->GetName();
  else of << "0";

  // RooArgSet of projDeps
  of << ","<<argSetToString(_projDeps);

  of <<");"<<endl;

  if (mcFile) {
    of << "toyMC->init(\""<<mcFile<<"\");"<<endl;
    of << "toyMC->fitBatch(\""<<fitFile<<"\","
       << firstSample <<","<<nSamples<<");"<<endl;
  }
  else {
    of << "toyMC->generateAndFitBatch(\""<<fitFile<<"\","
       << nEvents <<","<<nSamples<<");"<<endl;
  }
  of << "return 0;"<<endl;
  of << "}"<<endl;

  of.close();
}

// Create the batch scripts
// Returns list of root files where the fit output will be stored
// This list also gets written to "<fitFileBase>.files"
vector<string> BdkBatchMCStudy::createScripts(const char* mcFile, 
                                              const char* fitFileBase,
                                              const char* setupScript,
                                              Int_t fitsPerJob,
                                              Int_t nSamples, Int_t nEvtPerSample)
{
  vector<string> fileList;

  if (mcFile) {
    // open file with toy MC
    TFile f(mcFile);
    if (f.IsZombie()) return fileList;

    // get number of datasets in toy MC file
    TList *mcList = getDataSetKeys(f);
    f.Close();
    nSamples = mcList->GetEntries();
  }

  Int_t first = 0;
  Int_t job = 1;
  while (first<nSamples) {

    // File name of fit result
    TString fitFile = fitFileBase;
    fitFile += "-";
    fitFile += job;

    // File name of root script
    TString scriptFile = fitFile;
    scriptFile += ".cc";
    fitFile += ".root";

    // Write script for one job
    Int_t N = fitsPerJob;
    if (first+fitsPerJob>=nSamples) N = nSamples-first;
    cout << "Creating fit script "<<scriptFile<<endl;
    writeScript(mcFile,fitFile,scriptFile,setupScript,first,N,nEvtPerSample);
    
    // Put file in list and move next
    string s(fitFile);
    fileList.push_back(s);
    job++;
    first += fitsPerJob;
  }

  // Write the submit script
  ofstream submit;
  TString submitScript(fitFileBase);
  submitScript += "-sub";
  cout << "Creating submit script "<<submitScript<<endl;
  submit.open(submitScript);

  // Store the list of root files
  ofstream ofiles;
  TString ofilesPath(fitFileBase);
  ofilesPath += ".files";
  cout << "Creating file list "<<ofilesPath<<endl;
  ofiles.open(ofilesPath);

  for (Int_t i=0; i<fileList.size(); i++) {
    ofiles << fileList[i] << endl;
    TString logFile(fileList[i]);
    logFile.ReplaceAll(".root",".log");
    TString script(fileList[i]);
    script.ReplaceAll(".root",".cc");
    submit << "bsub "<<_submitOption<<" -o "<<logFile<<" bbrroot -q -b "<<script<<endl;
  }
  submit.close();
  ofiles.close();

  return fileList;
}


  

