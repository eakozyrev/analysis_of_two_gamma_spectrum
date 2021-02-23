#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TText.h"
#include "TMatrixTSym.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <fstream>
#include "RooChi2Var.h"
#include <string>
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

//#include "../mygenerator/Cmd3Generator2pi2pi0_ke.cc"


#define chi25Ccut 25.

#include <stdlib.h>

using namespace std;


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;


int nbinsgg = 70;
int nbinsggf = 200;

double signalN = 0.;
double dsignalN = 0.;
double effmc = 0.;
double deffmc = 0.;
double energy = 0.;
double emeas = 0;
int YEAR = 0;

double paramm[100];

TH1D *M2g3pi = new TH1D("M2g3pi","M2g3pi",nbinsgg,50,250);
TH1D *M2gkskl = new TH1D("M2gkskl","M2gkskl",nbinsgg,50,250);
TH1D *M2getag = new TH1D("M2getag","M2getag",nbinsgg,50,250);
TH1D *M2gbg = new TH1D("M2gbg","M2gbg",nbinsgg,50,250);
TH1D *M2gbck = new TH1D("M2gbck","M2gbck",nbinsgg,50,250);
TH1D *M2g = new TH1D("M2g","M2g",nbinsgg,50,250);
TH1D *M2gmc = new TH1D("M2gmc","M2gmc",nbinsgg,50,250);

TH1D *M2g3pif = new TH1D("M2g3pif","M2g3pif",nbinsggf,50,250);
TH1D *M2gksklf = new TH1D("M2gksklf","M2gksklf",nbinsggf,50,250);
TH1D *M2getagf = new TH1D("M2getagf","M2getagf",nbinsggf,50,250);
TH1D *M2gbgf = new TH1D("M2gbgf","M2gbgf",nbinsggf,50,250);


TF1 *flat = new TF1("flat","1",0,10000.);



void fill(){
    paramm[1] = 1.35779e+02;
    paramm[2] = 4;
    paramm[3] = 2.22365e-01;
    paramm[4] = 0;
    paramm[5] = 8.35561e+00;
    paramm[6] = 0;
    paramm[7] = 0;
    paramm[8] = 0;
    paramm[9] = 0;
    paramm[10] = 1;
}

inline bool exists_file(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


string exp_mc(double energy, double lowbin){
    int i = 0;
    string s = "";
    string name = "../4pisel/histograms/4pi_mc_";

    for(double i = 0; i < 10; i = i + 0.1){
	for(double j = 0; j < 800; j++){
                string btv = Form("../4pisel/histograms/mc_%g_%g.root",energy + i, lowbin + j);
		if(exists_file(btv) == true){
            	s = btv;
                return s;
		}
                btv = Form("../4pisel/histograms/mc_%g_%g.root",energy - i, lowbin + j);
		if(exists_file(btv) == true){
            	s = btv;
                return s;
		}
                btv = Form("../4pisel/histograms/mc_%g_%g.root",energy + i, lowbin - j);
		if(exists_file(btv) == true){
            	s = btv;
                return s;
		}
                btv = Form("../4pisel/histograms/mc_%g_%g.root",energy - i, lowbin - j);
		if(exists_file(btv) == true){
            	s = btv;
                return s;
		}
        }
    }

    return s;
}


double fit2g1(string filee, string filemc, string filebkgr){

  // Observable
  double startt = 50.;//startt, endd
  double endd = 250.;
  RooRealVar x("mpi0","m_{#pi^{0}}, MeV/c^{2}",startt, endd);
  TH1F *hexp = new TH1F("hexp","dfsdf",100,startt,endd);
  RooRealVar frac("frac","frac",0.9,0.40,1.0);

  TFile *newfile = TFile::Open(("../4pisel/histograms/"+filee).c_str());
  TH1D *Emeas = (TH1D*)newfile->Get("Emeas");
  emeas = Emeas->GetMean();
  Emeas = (TH1D*)newfile->Get("Year");
  double year = Emeas->GetMean();

  TTree* tree = (TTree*)newfile->Get("Tree");
   
  // Generate a data sample of 1000 events in x from model

   //RooRealVar x("mpi0","mpi0",50,250);
   x.setBins(100); 
   RooRealVar E_poln("E_poln","E_poln",-400,400); 
   RooRealVar P_poln("P_poln","P_poln",-400,400); 
   RooRealVar chi25C("chi25C","chi25C",-400,40000);
   RooRealVar tthg1("tthg1","tthg1",-5,5); 
   RooRealVar tthg2("tthg2","tthg2",-5,5); 
   RooRealVar tthg3("tthg3","tthg3",-5,5); 
   RooRealVar tthg4("tthg4","tthg4",-5,5);
   RooRealVar tthpic1("tthpic1","tthpic1",-5,5);
   RooRealVar tthpic2("tthpic2","tthpic2",-5,5);
   RooRealVar helicity1("helicity1","helicity1",-5,5);
   RooRealVar helicity2("helicity2","helicity2",-5,5);
   RooRealVar momg1("momg[0]","momg[0]",-5,5000);
   RooRealVar momg2("momg[1]","momg[1]",-5,5000);
   RooRealVar momg3("momg[2]","momg[2]",-5,5000);
   RooRealVar momg4("momg[3]","momg[3]",-5,5000);
   double cut_polar_angle = 0.95; 
   double cut_E = 500.;
   string selection = Form("helicity1 < 1. && helicity2 < 1.  && chi25C < 600 && mpi0 > %g && mpi0 < %g && abs(tthg1 - 1.571) < %g && abs(tthg2 - 1.571) < %g && abs(tthg3 - 1.571) < %g && abs(tthg4 - 1.571) < %g && (tthpic1 - 1.571) < %g && (tthpic2 - 1.571) < %g",startt, endd, cut_polar_angle, cut_polar_angle, cut_polar_angle, cut_polar_angle, cut_polar_angle, cut_polar_angle);
   TCut cutt = TCut(selection.c_str());
   cout << selection << endl;
   tree->Draw("mpi0 >> hexp",cutt,"goff");
   if((double)tree->GetEntries(cutt) < 1){signalN = 0.;return 0;}
   TH1F *myh = (TH1F *)gDirectory->Get("hexp"); 
   RooDataSet *data = new RooDataSet("data","data",RooArgSet(RooArgSet(x,helicity1,helicity2,chi25C,tthg1, tthg2, tthg3, tthg4),RooArgSet(tthpic1,tthpic2,momg1,momg2,momg3,momg4)),Import(*tree),Cut(selection.c_str()));
   RooKeysPdf kest_exp("kest_exp","kest_exp",x,*data,RooKeysPdf::MirrorBoth);
   
  //*************************************************************************
  //===========================MC============================================

  TFile *newfilemc = TFile::Open(("../4pisel/histograms/"+filemc).c_str());
  TH1D *fhEff = (TH1D*)newfilemc->Get("hEff");
  TTree* treemc = (TTree*)newfilemc->Get("Tree");
  effmc = (double)treemc->GetEntries(cutt)/100000.;
  deffmc = sqrt(effmc*(1.-effmc)/100000.);
  RooDataSet datamc("datamc","datamc",RooArgSet(RooArgSet(x,helicity1,helicity2,chi25C,tthg1, tthg2, tthg3, tthg4),RooArgSet(tthpic1,tthpic2,momg1,momg2,momg3,momg4)),Import(*treemc),Cut(selection.c_str()));
  RooKeysPdf kest1("kest1","kest1",x,datamc,RooKeysPdf::MirrorBoth) ;
  RooRealVar mg("mg","mg",-0.8,-3,3); 
  RooRealVar sg("sg","sg",1,0.01,10); 
  RooGaussian gauss("gauss","gauss",x,mg,sg);
  cout << "===========================================================================9" << endl;
  RooFFTConvPdf resol("lxg","landau (X) gauss",x,kest1,gauss);
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++9" << endl;
  //**************************************************************************
  //========================== BCKGR ===========================================

  TFile *newfileBCKGR = TFile::Open(("../4pisel/histograms/"+filebkgr).c_str());
  TTree* treeBCKGR = (TTree*)newfileBCKGR->Get("Tree");
  RooDataSet dataBCKGR("dataBCKGR","dataBCKGR",RooArgSet(RooArgSet(x,helicity1,helicity2,chi25C,tthg1, tthg2, tthg3, tthg4),RooArgSet(tthpic1,tthpic2,momg1,momg2,momg3,momg4)),Import(*treeBCKGR),Cut(selection.c_str()));
  RooKeysPdf kestBCKGR("kestBCKGR","kestBCKGR",x,dataBCKGR,RooKeysPdf::MirrorBoth);
  
  //**************************************************************************
  //==========================MODEL===========================================
  //RooRealVar c0("c0","c0",1.,0.,30.); 
  //RooRealVar c1("c1","c1",0.);//,-10.,10.); 
  //RooRealVar c2("c2","c2",0.,-1,0.); 
  //RooGenericPdf P("P","abs(c0 + c1*(mpi0-135) + c2*(mpi0-135)*(mpi0-135))",RooArgSet(x,c0,c1,c2)); 
  //RooGenericPdf P("P","abs(mpi0*0.0000001 + 0.0001 + c0 + c2*(mpi0-150)*(mpi0-150))",RooArgSet(x,c0,c2));
  RooAddPdf model("model","model",RooArgList(resol,kestBCKGR),frac) ;
  //model.plotOn(xframe1);
  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(*data) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  //RooMinuit m(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  //m.setVerbose(kTRUE) ;
  //m.setVerbose(kFALSE) ;

  // Call MIGRAD to minimize the likelihood
  RooFitResult* r = model.fitTo(*data,Save());
  //m.migrad();
  cout << "-log(L) at minimum = " << r->minNll() << endl ;
  cout << " (double)frac.getVal() = " <<  (double)frac.getVal() << endl;
  cout << "tree->GetEntries(cutt) = "  << tree->GetEntries(cutt) << endl;
  signalN = (double)frac.getVal()*tree->GetEntries(cutt);
  dsignalN = frac.getError()*tree->GetEntries(cutt);
  // Print values of all parameters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  //model.getParameters(x)->Print("s") ;
  RooPlot* xframe1 = x.frame(Title(Form("E_{beam} = %g MeV",emeas)));
  data->plotOn(xframe1);
  model.plotOn(xframe1); 
  model.plotOn(xframe1,Components(kestBCKGR),LineStyle(kDashed));
  paramm[0] = (double)sg.getVal();
  paramm[1] = (double)sg.getError();
  paramm[2] = (double)mg.getVal();
  paramm[3] = (double)mg.getError();  
  TCanvas *s = new TCanvas();
  TPad*    upperPad = new TPad("upperPad", "upperPad", 0.0,0.3,1.0,1.0);
  TPad*    lowerPad = new TPad("lowerPad", "lowerPad", 0.0,0.1,1.0,0.3);
  upperPad->Draw();
  lowerPad->Draw();
  upperPad->cd();
  //gPad->SetLeftMargin(0.15) ; 
  xframe1->GetYaxis()->SetTitleOffset(0.6);
  xframe1->SetNdivisions(8,"Y");
  //xframe1->GetYaxis()->SetLabelOffset(1.5); 
  xframe1->Draw();
  
  
  TH1F *hmc = new TH1F("hmc","dfssddf",myh->GetNbinsX(),myh->GetBinLowEdge(1),myh->GetBinCenter(myh->GetNbinsX())+myh->GetBinWidth(1)/2.);
  model.fillHistogram(hmc,RooArgList(x),tree->GetEntries(cutt));
  myh->Divide(hmc);
  myh->SetYTitle("data/fit");
  myh->SetTitle("");
  myh->SetAxisRange(0.5,1.5,"Y");
  myh->SetAxisRange(startt,endd,"X");
  myh->SetXTitle("m_{#pi^{0}}, MeV/c^{2}");
  
  
  lowerPad->cd();
  gStyle->SetOptStat(0);
  myh->SetLineWidth(3.);
  myh->SetLineColor(1);
  gStyle->SetTextSize(1.8);
  myh->SetLabelSize(0.19,"xy");
  myh->SetTitleSize(0.19,"xy");
  myh->SetNdivisions(2,"Y");
  myh->SetTitleOffset(1.2,"x");
  myh->SetTitleOffset(0.2,"y");
  lowerPad->SetLeftMargin(0.15) ;	
  myh->Draw();
  s->SaveAs(("plots/"+filee+".png").c_str());
  s->SaveAs(("plots/"+filee+".root").c_str());
  s->Close(); 
  newfile->Close();
  newfilemc->Close();
  newfileBCKGR->Close();
  cout << "===============================================" << endl;
  cout << "                     end                       " << endl;
  cout << "===============================================" << endl;
  return 1.;

}


void resolution_study(){
    //string name = "data.dat";
    string name = "data1920.dat";
    ifstream streamm(name.c_str());
    TH1D *hChi2_NDF = new TH1D("hChi2_NDF","hChi2_NDF",100,-1.,10.);
    int i = 0;
    double number[10000];
    double width1[10000],width2[10000],mean1[100000],mean2[10000],frac1[10000],factor[10000];
    double dwidth1[10000],dwidth2[10000],dmean1[100000],dmean2[10000],dfrac1[10000],dfactor[10000];
    
    while(i < 10000){
      string ss,ss1,ss2,ss3,ssb;
      double a;
      streamm >> ss >> ss1 >> ssb >> a >> a;
      if(ss=="stop" || streamm.eof()==1){break;}
      fit2g1(("data/" + ss).c_str(),("2pi2pi0_v2/"+ss1).c_str(),("multih/" + ssb).c_str());
      width1[i] = paramm[0];dwidth1[i] = paramm[1];
      width2[i] = paramm[2];dwidth2[i] = paramm[3];
      mean1[i] = paramm[1]; dmean1[i] = paramm[21];
      mean2[i] = paramm[4]; dmean2[i] = paramm[24];
      frac1[i] = paramm[3]; dfrac1[i] = paramm[23];
      factor[i] = paramm[50]; dfactor[i] = paramm[51];
      hChi2_NDF->Fill(paramm[0]);
      number[i] = emeas;
      i++;
    }
	
   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(0-10,-100.,2000+10,10000.);
   frd->SetXTitle("E_{beam}, MeV");
   frd->SetYTitle("MeV");

   TGraphErrors *Crossu  = new TGraphErrors(i,number,factor,0,dfactor);
   Crossu->SetMarkerColor(2);
   Crossu->SetMarkerStyle(20);
   Crossu->SetLineColor(2);
   Crossu->SetLineWidth(2.);
   Crossu->SetTitle("width1");
   Crossu->Draw("P");

   TGraphErrors *Cross  = new TGraphErrors(i,number,width1,0,dwidth1);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(20);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(2.);
   Cross->SetTitle("additional smearing");
   Cross->Draw("P");

   TGraphErrors *Cross1  = new TGraphErrors(i,number,width2,0,dwidth2);
   Cross1->SetMarkerColor(1);
   Cross1->SetMarkerStyle(20);
   Cross1->SetLineColor(1);
   Cross1->SetLineWidth(2.);
   Cross1->SetTitle("additional shift");
   Cross1->Draw("P");

   TGraphErrors *Cross2  = new TGraphErrors(i,number,mean1,0,dmean1);
   Cross2->SetMarkerColor(4);
   Cross2->SetMarkerStyle(20);
   Cross2->SetLineColor(4);
   Cross2->SetLineWidth(2.);
   Cross2->SetTitle("mean1");
   Cross2->Draw("P");

   TGraphErrors *Cross3  = new TGraphErrors(i,number,mean2,0,dmean2);
   Cross3->SetMarkerColor(3);
   Cross3->SetMarkerStyle(20);
   Cross3->SetLineColor(3);
   Cross3->SetLineWidth(2.);
   Cross3->SetTitle("mean2");
   Cross3->Draw("P");

   TGraphErrors *Cross4  = new TGraphErrors(i,number,frac1,0,dfrac1);
   Cross4->SetMarkerColor(1);
   Cross4->SetMarkerStyle(20);
   Cross4->SetLineColor(1);
   Cross4->SetLineWidth(2.);
   Cross4->SetTitle("frac1");
   Cross4->Draw("P");

   TCanvas *s1 = new TCanvas();
   hChi2_NDF->Draw();

}

int fill_new_cross_table(){

  TSQLResult *res;
  TSQLRow *row;
  TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
  string iwant; 
  res = db->Query(iwant.c_str());
  int i = 0;
  double step = 0.02;
  for(double en = 0.8; en < 2.03; en = en + step){
    iwant = Form("select emeas, demeas, sigma, dsigma, lum  from epoints where comment = 'v1' and emeas > %g and emeas < %g  and Lum IS NOT NULL and year < 2017",(en-step)*500,(en)*500.);
    //cout << iwant << endl;
    if(en > 1.005 && en < 1.04)step = 0.002;
    else step = 0.02;
    res = db->Query(iwant.c_str());
    i = 0;
    double totalen = 0;
    double dtotalen = 0;
    double totalcr=0.;
    double dtotalcr = 0;
    double totallum = 0;
    double norm = 0.;
    while((row = res->Next())){
      TString Lum = row->GetField(4);
      totallum = totallum+(double)Lum.Atof();
      TString En = row->GetField(0);
      TString dEn = row->GetField(1);
      TString Cr = row->GetField(2);
      TString dCr = row->GetField(3);
      double lum = (double)Lum.Atof();
      totalen = totalen + lum*(double)En.Atof();
      dtotalen = dtotalen + lum*lum*pow((double)dEn.Atof(),2.);
      totalcr = totalcr + (double)Cr.Atof()/dCr.Atof();
      dtotalcr = dtotalcr + 1.;
      norm = norm + 1./(double)dCr.Atof();
      i++;
    }
    if(i==0)continue;
    totalen = totalen/totallum;
    dtotalen = sqrt(dtotalen)/totallum;
    totalcr = totalcr/norm;
    dtotalcr = sqrt(dtotalcr)/norm;

    cout << totalen << " +/- " << dtotalen << "     "  << totalcr << " +/- " << dtotalcr << endl;
    iwant = Form("INSERT into cross_sec (energy, denergy, sigma, dsigma) VALUES  (%g, %g, %g, %g)",totalen,dtotalen,totalcr, dtotalcr);
    res = db->Query(iwant.c_str());
  }
  
  return 1;
  
}

double efficiency(double en, int year){

    TSQLResult *res;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    string iwant = Form("select emeas, eff4pi  from epoints where comment = 'v1' and emeas > %g and emeas < %g and year = %i",en-30.,en+30,year); 
    res = db->Query(iwant.c_str());
    int i = 0;
    double resul = 0;
    double Ntotal = 0;
    double dNtotal = 0;
    while((row = res->Next())){
          TString En = row->GetField(0);
          TString N = row->GetField(1);
	  i++;
	  resul = resul + (double)N.Atof();
    }

    db->Close();
    return resul/(double)i;
}


void fill_eff(){

    TSQLResult *res;
    TSQLResult *res1;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    char sql[4096];
    const char *ins = "select emeas, name, year from epoints where comment = 'v1' and Lum IS NOT NULL";
    res = db->Query(ins);
    char ins1[10000];
    while((row = res->Next())){
          TString En = row->GetField(0);
          string Name =  row->GetField(1);
          TString Year =  row->GetField(2);
          string insert = "UPDATE epoints SET eff = %lf WHERE (comment = 'v1' AND name = '"+Name+"')";
          strcpy(ins1, insert.c_str());
          sprintf(sql,ins1,efficiency((double)En.Atof(),(int)Year.Atof()));
          cout << sql << endl;
          res1 = db->Query(sql);
          delete res1;
    }

}

void fill_effictr(){

    TSQLResult *res;
    TSQLResult *res1;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    char sql[4096];
    const char *ins = "select emeas, name, year from epoints where comment = 'v1' and Lum IS NOT NULL";
    res = db->Query(ins);
    char ins1[10000];
    while((row = res->Next())){
          TString En = row->GetField(0);
          string Name =  row->GetField(1);
          TString Year =  row->GetField(2);
          string insert = "UPDATE epoints SET effictr1exp = 1,deffictr1exp=0,effictr1mc = 1,deffictr1mc = 0 WHERE (comment = 'v1' AND name = '"+Name+"')";
          strcpy(ins1, insert.c_str());
          sprintf(sql,ins1);
          cout << sql << endl;
          res1 = db->Query(sql);
          delete res1;
    }

}


void fill_emeas(){

    TSQLResult *res;
    TSQLResult *res1;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    char sql[4096];
    const char *ins = "select emeas, name, year from epoints where comment = 'v1' and Lum IS NOT NULL";
    res = db->Query(ins);
    char ins1[10000];

    string insert = "UPDATE epoints SET emeas = %lf WHERE name = '509.1_18019'";
    strcpy(ins1, insert.c_str());
    sprintf(sql,ins1,509.7);
    cout << sql << endl;
    res1 = db->Query(sql);


    delete res1;

}


void fill_emeas1(){

  ifstream stream("../data/names.dat");
  TSQLResult *res;
  TSQLResult *res1;
  TSQLRow *row;
  TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
  char sql[4096];
  while(stream.eof()==0){
    string A,B;
    stream >> A >> B;
    if(stream.eof()==1)break;
    ifstream stream1("../data/run_data_fit.dat");
    string CC,D;
    cout << A << endl;
    while(A.compare(CC)!=0){stream1 >> CC;}
    stream1 >> CC >> CC;
    double lum;
    stream1 >> lum;
    double energ;
    stream1 >> energ;
    string ins = Form(("UPDATE epoints set emeas = %g, ebeam = %i, lum = %g where year > 2018 and name = '" +B + "' ").c_str(),energ,int(energ),lum);
    cout << ins << endl;
    strcpy(sql,ins.c_str());
    res = db->Query(sql);
  }

}



/*
 *
 * ebeam emeas demeas sigma dsigma Lum dLum rad drad trig dtrig N4pi dN4pi eff4pi N4piwc dN4piwc eff4piwc N4piw0 dN4piw0 eff4piw0 eff1 eff2 eff3 eff4 comment
 *
 * */
void runnew(string namefile){
    bool write_update = false;
    string name = namefile;//"runnew_b.dat";
	//string name = "runnew1.dat";
    ifstream streamm(name.c_str());
    ofstream stream("../data/fitN.dat");
    int i = 0;
    char sql[4096];
    TSQLResult *res;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    string names;char insert[10000];
    string ins;// = "INSERT INTO epoints (name,ebeam, emeas, demeas,N4pi,dN4pi,eff4pi,year,comment) VALUES (" + names + "%lf, %lf, %lf, %lf, %lf, %lf, %i, 'v1')";
    //const char *ins = "UPDATE epoints SET N4pi = %lf, dN4pi = %lf, eff4pi = %lf, year = %li WHERE (comment = 'v1' AND ebeam = %lf)";
    while(i < 10000){
        string ss,ss1,ss2,ss3,ssb;
        double a;
        streamm >> ss >> ss1 >> ssb >> a >> a;
        if(ss=="stop" || streamm.eof()==1){break;}
        fit2g1(("data/" + ss).c_str(),("2pi2pi0_v2/"+ss1).c_str(),("multih/" + ssb).c_str());
        ifstream stream_name("../data/names.dat");
	while(stream_name.eof()==0){
        	stream_name >> ss2 >> ss3;
        	if(ss2 == ss){names = ss3;break;}
	}
	cout << "HERE 1 " << endl;
        if(write_update)ins = "INSERT INTO epoints (name,ebeam, emeas, demeas,N4pi,dN4pi,eff4pi,deff4pi,year,comment) VALUES ('" + names + "',%lf, %lf, %lf, %lf, %lf, %lf,%lf, %i, 'v1')";
	cout << "HERE 2 " << endl;
        if(!write_update)ins = "UPDATE epoints SET N4pi = %lf, dN4pi = %lf, eff4pi = %lf WHERE (comment = 'v1' AND name = '" + names + "')";
        strcpy(insert, ins.c_str());
	cout << "HERE 3 " << endl;
        energy = ((int)(energy*10.))/10.;
        stream  << energy << " " << signalN << " " << dsignalN << " " << effmc << " " << deffmc << endl;
	cout << "HERE 4 " << endl;
        if(write_update)sprintf(sql,insert,energy,emeas,0.,signalN,dsignalN,effmc,deffmc,YEAR);
		cout << "HERE 5 " << endl;
        if(!write_update)sprintf(sql,insert,signalN,dsignalN,effmc);
		cout << "HERE 6 " << endl;
        cout << sql << endl;
		cout << "HERE 7 " << endl;
        res = db->Query(sql);
		cout << "HERE 8 " << endl;
        delete res;
		cout << "HERE 9 " << endl;
        i++;
    }
}

// select emeas, Lum, rad, trig, N4pi from epoints where comment = 'v1' and Lum IS NOT NULL order by 'emeas';

void fill_cross(){

    TSQLResult *res;
    TSQLResult *res1;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    char sql[4096];
    const char *ins = "select ebeam, Lum, rad, trig, N4pi, dN4pi, eff4pi, deff4pi, dLum, name, effictr1exp,deffictr1exp,effictr1mc,deffictr1mc, year, efficpi0, defficpi0, efficpi0mc, defficpi0mc, Npi01, dNpi01, Npi01mc, dNpi01mc from epoints where comment = 'v1' and Lum IS NOT NULL";
    res = db->Query(ins);
    char ins1[10000];
    while((row = res->Next())){
          TString En = row->GetField(0);
          TString N = row->GetField(4);
          TString dN = row->GetField(5);
          TString Lum = row->GetField(1);
          TString Rad = row->GetField(2);
          TString Trig = row->GetField(3);
          TString Eff = row->GetField(6);
          TString dEff = row->GetField(7);
	  TString dLum = row->GetField(8);
          string Name =  row->GetField(9);
          TString Year =  row->GetField(14);
          TString effictr1exp =  row->GetField(10);
          TString deffictr1exp =  row->GetField(11);
          TString effictr1mc =  row->GetField(12);
          TString deffictr1mc =  row->GetField(13);
	  TString efficpi0  =  row->GetField(15);
	  TString defficpi0 =   row->GetField(16);
	  TString efficpi0mc =   row->GetField(17);
	  TString defficpi0mc =   row->GetField(18);
	  TString Npi01  =   row->GetField(19);
	  TString dNpi01  =   row->GetField(20);
	  TString Npi01mc  =   row->GetField(21);
	  TString dNpi01mc  =   row->GetField(22);          
          double effic = Eff.Atof();
	  double factorpi0 = 1;
          if(effic < 0.0001){
		cout << "too low efficiency " << endl;
		effic = 0.001;
	  }
          double factor = 1;
          double dfactor = 0;
          if(Year.Atof()==2011){

	    //	factor = pow(effictr1exp.Atof()/effictr1mc.Atof(),2.);
	    //	dfactor = pow(deffictr1exp.Atof()/effictr1exp.Atof(),4.);
	    //	if(fabs(factor - 1) > 0.4)factor = 1;
	  }
	  if(/*En.Atof() > 699 && */efficpi0.Atof() > 0.3){
	      //factorpi0 = efficpi0.Atof()*efficpi0.Atof()/efficpi0mc.Atof()/efficpi0mc.Atof();
	  }
          double cross = N.Atof()/Lum.Atof()/Rad.Atof()/Trig.Atof()/effic/factor/factorpi0;
          double dcross = 0;
	  //cross = Npi01.Atof()/Lum.Atof()/Rad.Atof()/Trig.Atof()/Npi01mc.Atof()*100000.;
	  cout << Npi01.Atof() << " " << Lum.Atof()  << " " << Rad.Atof()  << " " << Trig.Atof()  << " " <<Npi01mc.Atof()  << " " << effic << endl;
	  double rell = dN.Atof()/N.Atof();
	  //rell = dNpi01.Atof()/Npi01.Atof();
          if(cross > 0)dcross = sqrt(pow(rell,2.)/* + pow(dEff.Atof()/Eff.Atof(),2.)*/ + pow(dLum.Atof()/Lum.Atof(),2.) +  dfactor)*cross;
	  
          cout << cross  << " +/- " << dNpi01.Atof()/Npi01.Atof() << " +/- " << dLum.Atof()/Lum.Atof() << endl;
          string insert = "UPDATE epoints SET sigma = %lf, dsigma = %lf WHERE (comment = 'v1' AND name = '"+Name+"')";
          strcpy(ins1, insert.c_str());
          sprintf(sql,ins1,cross,dcross,En.Atof());
          cout << sql << endl;
          res1 = db->Query(sql);
          delete res1;
    }

}



void luminosity_fill(){
    ifstream stream11("../data/lum2011.dat");
    TSQLResult *res;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    const char *ins = "UPDATE epoints SET Lum = %lf, dLum = %lf  WHERE (year = %i AND comment = 'v1' AND ebeam = %lf)";
    char sql[4096];char insert[4096];
    double a,b,c;
    while(stream11.eof()==0){
        stream11 >> a >> b >>c;
        sprintf(sql,ins,b,c,2011,a);
        cout << sql << endl;
        res = db->Query(sql);
        delete res;
    }
    ifstream stream12("../data/lum2012.dat");
    while(stream12.eof()==0){
        stream12 >> a >> b >>c;
        sprintf(sql,ins,b,c,2012,a);
        cout << sql << endl;
        res = db->Query(sql);
        delete res;
    }
    ifstream stream13("../data/lum2013.dat");
    while(stream13.eof()==0){
        stream13 >> a >> b >>c;
        sprintf(sql,ins,b,c,2013,a);
        cout << sql << endl;
        res = db->Query(sql);
        delete res;
    }
    ifstream stream14("../data/lum2014.dat");
    while(stream14.eof()==0){
        string ss;
        stream14 >> a >> b >>c >> ss;
        string comand = "UPDATE epoints SET Lum = %lf, dLum = %lf  WHERE ( name = '"+ss+"');";
        strcpy(insert, comand.c_str());
        sprintf(sql,insert,b,c);
        cout << sql << endl;
        res = db->Query(sql);
        delete res;
    }
    ifstream stream17("../data/lum2017.dat");
    while(stream17.eof()==0){
        string ss;
        stream17 >> a >> b >>c >> ss;
        string comand = "UPDATE epoints SET Lum = %lf, dLum = %lf  WHERE ( name = '"+ss+"');";
        strcpy(insert, comand.c_str());
        sprintf(sql,insert,b,c);
        cout << sql << endl;
        res = db->Query(sql);
        delete res;
    }
}



void plot(string name,string dname, double years, bool drawcanvas, int color){
    TSQLResult *res;
    TSQLRow *row;
    TSQLServer *db = TSQLServer::Connect("mysql://localhost/2pi2pi0","root","");
    string btv = "select emeas, " + name + ',' +dname + Form(" from epoints where comment = 'v1' and Lum IS NOT NULL and year = %f",years);
    const char *ins = btv.c_str();    
    res = db->Query(ins);
    if(drawcanvas == true){
		TCanvas *s0 = new TCanvas();
    		TH1F *fr0  = s0->DrawFrame(0.7,1.,2.070,50.);
    		fr0->SetXTitle("E_{c.m.}, MeV");
    		fr0->SetYTitle("#sigma (e^{+}e^{-} #rightarrow 2#pi^{#pm}2#pi^{0}) ,nb");
    }


    TGraphErrors *Cross0 = new TGraphErrors();
    Cross0->SetTitle(Form("SCAN%g",years));
    int i = 0;
    double Ntotal = 0;
    double dNtotal = 0;
    double lum = 0.;
    while((row = res->Next())){
          TString En = row->GetField(0);
	  cout << En;
          int k = 1;int j = 2;
          TString draw = row->GetField(k);
          TString ddraw = row->GetField(j);
          TString title = row->GetField(12);
          double magnitude = (double)draw.Atof();
          double dmagnitude = (double)ddraw.Atof();
          Cross0->SetPoint(i,(double)2.*En.Atof()/1000.,magnitude);
          Cross0->SetPointError(i,0,(double)dmagnitude);
	  //         Ntotal = Ntotal + (double)N.Atof();
	  // dNtotal = dNtotal + (double)dN.Atof()*(double)dN.Atof();
	  // lum = lum + Lum.Atof();
          i++;
    }
    //cout << Ntotal << " +/- " << sqrt(dNtotal) << endl;
    // cout << " Luminosity = " << lum << endl;
    Cross0->SetMarkerColor(color);
    Cross0->SetMarkerStyle(20);
    Cross0->SetLineColor(color);
    Cross0->SetLineWidth(2.);
    Cross0->Draw("P");
}





//Cmd3Generator2pi2pi0_ke gen;
/*
void draw_model(){

    gen.isISR = 0;
    gen.omegapi0_ = 1.0;
    gen.omegapi0phase_ = 0.;

    gen.a1pi_rhopi_c_phase = 2.31409e+00;
    gen.a1pi_sigmapi_c_phase = -2.52700e-01;

    gen.a1pi_rhopi_ =  gen.a1pi_rhopi_c;
    gen.a1pi_rhopi_phase =  gen.a1pi_rhopi_c_phase;
    gen.a1pi_sigmapi_ = gen.a1pi_sigmapi_c;
    gen.a1pi_sigmapi_phase = gen.a1pi_sigmapi_c_phase;

    gen.rhoprhom_ = 0.;
    gen.rhof0_ = 0.;
    gen.rhof0phase_ = 0.;
    gen.a2pi_ = 0.;
    int n = 15;
    double enfix = 1.1;
    double en[1000], cr[1000];
    for(int i = 0; i < n; i++){
        en[i] = 700. + 30*i;
        gen.Energy = en[i]/1000.;
        double e = gen.Energy;
        gen.a1pi_rhopi_c = 0.1*sqrt(gen.cross_a1_rhopi_pi(&e,&e)/gen.cross_a1_rhopi_pi(&enfix,&enfix));
        gen.omegapi0_ = sqrt(gen.cross_omega_pi0(&e,&e)/gen.cross_omega_pi0(&enfix,&enfix));
        gen.a1pi_sigmapi_c = 0.01*sqrt(gen.cross_a1_sigmapi_pi(&e,&e)/gen.cross_a1_sigmapi_pi(&enfix,&enfix));
        gen.a1pi_rhopi_ =  gen.a1pi_rhopi_c;
        gen.a1pi_rhopi_phase =  gen.a1pi_rhopi_c_phase;
        gen.a1pi_sigmapi_ = gen.a1pi_sigmapi_c;
        gen.a1pi_sigmapi_phase = gen.a1pi_sigmapi_c_phase;


        gen.fill_sample();
        gen.fill_sample_c();
        cr[i] = gen.cross_section()/gen.cross_section_c();
        cout << e << " " << cr[i] << endl;
    }


    TGraphErrors *Cross2d   = new TGraphErrors(n,en,cr,0,0);
    Cross2d->SetMarkerColor(6);
    Cross2d->SetMarkerStyle(24);
    Cross2d->SetLineColor(6);
    Cross2d->SetLineWidth(3.);
    Cross2d->Draw("P");

}
*/
void draw_babar(){

    ifstream stream("../data/2pi2pi0babar.dat");
    double en[1000], den[1000], cr[1000], dcr[1000];
    int i = 0;
    while(stream.eof()==0){
        double A, B,C,D,E;
        stream >> A >> B >> C >> D >> E;
        en[i] = (A/2.+B/2.);
        den[i] = (B/2. - A/2.);
        cr[i] = C;
        dcr[i] = D;
        i++;
    }


    TGraphErrors *Cross2   = new TGraphErrors(i-1,en,cr,den,dcr);
    Cross2->SetMarkerColor(4);
    Cross2->SetMarkerStyle(20);
    Cross2->SetLineColor(4);
    Cross2->SetLineWidth(2.);
    Cross2->Draw("P");

}

void draw_SND(){

    ifstream stream("/home/eakozyrev/diskD/4pi_new/data/2pi2pi0SND.dat");
    double en[1000], den[1000], cr[1000], dcr[1000];
    int i = 0;
    while(stream.eof()==0){
        double A, B,C,D,E;
        stream >> A >> B >> C >> D;
        en[i] =A;
        den[i] = B;
        cr[i] = C;
        dcr[i] = D;
        i++;
    }


    TGraphErrors *Cross2f   = new TGraphErrors(i-1,en,cr,den,dcr);
    Cross2f->SetMarkerColor(6);
    Cross2f->SetMarkerStyle(20);
    Cross2f->SetLineColor(6);
    Cross2f->SetLineWidth(2.);
    Cross2f->Draw("P");

}




void compare_hist(string hist, string file1, string hfile2){

    TFile *file = TFile::Open(("" + file1).c_str());
    TH1D *hist1 = (TH1D*)file->Get(hist.c_str());

    TFile *file2 = TFile::Open(("" + hfile2).c_str());
    TH1D *hist2 = (TH1D*)file2->Get(hist.c_str());
     
    hist2->SetNormFactor(hist1->GetEntries());
    //hist2->Rebin(6);
    //hist1->Rebin(6);
    hist1->SetXTitle(hist.c_str());
    hist2->SetXTitle(hist.c_str());
    hist1->Draw();
    hist1->SetLineWidth(2.);
    hist2->SetLineColor(4);
    hist2->SetFillColor(4);
    hist2->SetFillStyle(3002);
    hist2->Draw("same");
    
    
}


int script_compare(){

	TCanvas s("dsds","dsdsd",1000,500);
        s.Divide(2,1);
	ifstream stream("../data/run_data_fit.dat");
	string filee,filemc,filebkg;
	double btv;
        while(stream.eof() == 0){

		stream >> filee >> filemc >> filebkg >> btv >> btv;
		if(btv < 530){continue;}
		if(stream.eof() == 1)return 1;
		s.cd(1);
		compare_hist("Ephoton",("../4pisel/histograms/"+filee).c_str(),("../4pisel/histograms/"+filemc).c_str());
		s.cd(2);
		compare_hist("Ephoton_kin",("../4pisel/histograms/"+filee).c_str(),("../4pisel/histograms/"+filemc).c_str());
		TText tex(0.4,0.8,filee.c_str());
		tex.SetTextSize(0.05);
		tex.Draw();
		s.SaveAs(("plots/E_photon/"+filee+".png").c_str());
	}
return 1;

}





void draw_ratio(){

  ifstream streamneutr("/home/eakozyrev/4pi/cross_sections/myresult.dat");

  ifstream streamchar("/home/eakozyrev/4pi/cross_sections/2pip2pim.dat");

  double crneutr[100],crch[100],dcrneutr[100],dcrch[100],ratio[100],dratio[100],en[100];

  int i = 0;

  while(streamneutr.eof()==0){
      streamneutr >> en[i] >> crneutr[i] >> dcrneutr[i];
      i++;
  }
  i = 0;
  while(streamchar.eof()==0){
      streamchar >> en[i] >> crch[i] >> dcrch[i];
      ratio[i] = crneutr[i]/crch[i];
      dratio[i] = ratio[i]*sqrt(pow(dcrneutr[i]/crneutr[i],2.) + pow(dcrch[i]/crch[i],2.));
      cout << i << "  "  << en[i] << "  " << crneutr[i]<< " " << crch[i] << " " << ratio[i] << endl;
      i++;
  }

  i--;

  TCanvas *s = new TCanvas();
  TH1F *fr  = s->DrawFrame(900,1.,1070.,11.);
  fr->SetXTitle("E_{c.m.}, MeV");
  fr->SetYTitle("#sigma (e^{+}e^{-} #rightarrow 2#pi^{#pm}2#pi^{0}) ,nb");
  TGraphErrors *Cross   = new TGraphErrors(i,en,ratio,0,dratio);
  Cross->SetMarkerColor(2);
  Cross->SetMarkerStyle(20);
  Cross->SetLineColor(2);
  Cross->SetLineWidth(2.);
  Cross->Draw("P");




}





void resolution(string name1, double limit){

    ifstream stream("2files_exp.dat");
    TH1D *Trig;
    TFile *file;
    double res_exp[1000],res_mc[1000];
    double dres_exp[1000], dres_mc[1000];
    double en[1000],en_mc[1000];
    int i = 0;
    int nexp;
    while(stream.eof()==0){

        string ss;
        stream >> ss;
        if(ss=="stop"){break;}
        file = TFile::Open(("../histograms/" + ss).c_str());
        Trig = (TH1D*)file->Get("Emeas");
        en[i] = Trig->GetMean()*2.;
        Trig = (TH1D*)file->Get(name1.c_str());
        Trig->SetAxisRange(-limit,limit,"X");
        res_exp[i] = Trig->GetRMS();
        dres_exp[i] = Trig->GetRMSError();
        cout << ss << " " << res_exp[i]<< " " << dres_exp[i] << endl;
        i++;
    }
    nexp = i;
    ifstream stream_mc("5files_2pi2pi0_phsp_mc.dat");
    i = 0;
    while(stream_mc.eof()==0){

        string ss;
        stream_mc >> ss;
        if(ss=="stop"){break;}
        file = TFile::Open(("../histograms/" + ss).c_str());
        Trig = (TH1D*)file->Get("Emeas");
        en_mc[i] = Trig->GetMean()*2.;
        Trig = (TH1D*)file->Get(name1.c_str());
        Trig->SetAxisRange(-limit,limit,"X");
        res_mc[i] = Trig->GetRMS();
        dres_mc[i] = Trig->GetRMSError();
        cout << ss << " " << res_mc[i]<< " " << dres_mc[i] << endl;
        i++;
    }


    TGraphErrors *Crossbabar   = new TGraphErrors(nexp,en,res_exp,0,dres_exp);
    Crossbabar->SetMarkerColor(1);
    Crossbabar->SetMarkerStyle(20);
    Crossbabar->SetLineColor(1);
    Crossbabar->SetLineWidth(2.);
    Crossbabar->SetTitle("data");
    Crossbabar->Draw("P");

    TGraphErrors *Crossbabar1   = new TGraphErrors(i,en_mc,res_mc,0,dres_mc);
    Crossbabar1->SetMarkerColor(2);
    Crossbabar1->SetMarkerStyle(20);
    Crossbabar1->SetLineColor(2);
    Crossbabar1->SetLineWidth(2.);
    Crossbabar1->SetTitle("MC");
    Crossbabar1->Draw("P");
}


void resol(){

    TCanvas *s = new TCanvas();
    s->Divide(2,2);
    TH1F *fr  = s->cd(1)->DrawFrame(500,0.,2170.,0.02);
    fr->SetXTitle("E_{c.m.}, MeV");
    fr->SetYTitle("#sigma(#theta_{rec}^{#pi^{+/-}} - #theta_{det}^{#pi^{+/-}}), rad");
    resolution("hDelta_th",0.03);

    TH1F *fr1  = s->cd(2)->DrawFrame(500,0.,2170.,20.);
    fr1->SetXTitle("E_{c.m.}, MeV");
    fr1->SetYTitle("#sigma(p_{rec}^{#pi^{+/-}} - p_{det}^{#pi^{+/-}}), MeV/c");
    resolution("hDelta_P",40);

    TH1F *fr2  = s->cd(3)->DrawFrame(500,0.,2170.,0.02);
    fr2->SetXTitle("E_{c.m.}, MeV");
    fr2->SetYTitle("#sigma(#theta_{rec}^{#gamma} - #theta_{det}^{#gamma}), rad");
    resolution("hDelta_th_ph",0.03);

    TH1F *fr3  = s->cd(4)->DrawFrame(500,0.,2170.,20);
    fr3->SetXTitle("E_{c.m.}, MeV");
    fr3->SetYTitle("#sigma(p_{rec}^{#gamma} - p_{det}^{#gamma}), MeV/c");
    resolution("hDelta_P_ph",40);


}



double fit2g2(string filee,string filemc){
    fill();
    bool expmc=true;
    TF1 *fit = new TF1("fit","[0]*([3]/(sqrt(2.*3.1415)*sqrt([2]*[2]+[9]))*exp(-(x-[1])*(x-[1])/2./([2]*[2]+[9]))+(1-[3])/(sqrt(2.*3.1415)*sqrt([5]*[5]+[9]))*exp(-(x-[1]-[4])*(x-[1]-[4])/2./([5]*[5]+[9]))) +([6]+[7]*(x-145.)+[8]*(x-145.)*(x-145.))*[10]",0,1000.);
    fit->SetParName(0,"yields");
    fit->SetParName(1,"mean");
    fit->SetParameters(paramm);
    fit->SetParLimits(0,0,10000.);
    fit->SetParLimits(1,133.,144.);
    fit->SetParLimits(2,1,10.);
    fit->SetParLimits(3,0.,1.);
    fit->SetParLimits(4,-6,6.);
    fit->SetParLimits(5,1,40.);
    fit->SetParLimits(9,1,50.);
    fit->FixParameter(6,0);
    fit->FixParameter(7,0);
    fit->FixParameter(8,0);
    fit->FixParameter(9,0);
    fit->FixParameter(10,1);
    TCanvas *s = new TCanvas();
    TFile *file0;TH1D *Mpi0;

    file0 = TFile::Open(filemc.c_str());
    TH1D *fhEff = (TH1D*)file0->Get("hEff");
    Mpi0 = (TH1D*)file0->Get("Mpi0f");
    cout << "Mpi0->GetRMS() = " << Mpi0->GetRMS() << endl;
    fit->SetParameter(0,Mpi0->GetEntries());
    fit->SetParLimits(0,0.8*Mpi0->GetEntries(),Mpi0->GetEntries()+100.);
    fit->SetParLimits(5,2,Mpi0->GetRMS()*0.9);
    fit->SetParLimits(2,1,Mpi0->GetRMS()*0.25);
    Mpi0->Fit("fit","RL","",115,155); 
    fit->ReleaseParameter(6);
    fit->ReleaseParameter(7);
    Mpi0->Fit("fit","RL","RL",105,165);
    fit->SetParLimits(8,-1,0.);
    Mpi0->Fit("fit","RL");Mpi0->Fit("fit","RL");Mpi0->Fit("fit","RL");Mpi0->Fit("fit","RL");
    for(int i = 0; i < fit->GetNpar(); i++){paramm[i]=fit->GetParameter(i);paramm[i+20]=fit->GetParError(i);}
    paramm[0]= (double)fit->GetChisquare()/(double)fit->GetNDF();
    effmc = fhEff->GetMean()*(double)fit->GetParameter(0)/(double)Mpi0->GetEntries();
    paramm[50] = (double)fit->GetParameter(0)/(double)Mpi0->GetEntries();
    paramm[51] = (double)fit->GetParError(0)/(double)Mpi0->GetEntries();
    deffmc = sqrt(effmc*(1.-effmc)/30000.);
    fit->SetParameter(0,0.);
    fit->SetLineColor(1);
    fit->Draw("same");
    s->SaveAs(("plots/mc"+filee+".png").c_str());
    fit->SetLineColor(2);
    
    for(int i = 1; i <=5; i++){fit->FixParameter(i,fit->GetParameter(i));}

    // ================================== DATA ==============================
    // **********************************************************************
    fit->ReleaseParameter(9);
    file0 = TFile::Open(("../4pisel/histograms/"+filee).c_str());
    TH1D *Ebeeem = (TH1D*)file0->Get("Ebeeem");
    energy = Ebeeem->GetMean();
    if(energy < 100){
	cout << "Some errors are observed here" << endl;
	return 0.;
    }
    double lowbin = 2;
    TH1D *hRun = (TH1D*)file0->Get("Runnumber");
    for(int i = 1000; i < 100000; i++){if(hRun->GetBinContent(i) > 1){lowbin = i;break;}}
    cout << "lowbin = " << lowbin << endl;
    TH1D *Emeas = (TH1D*)file0->Get("Emeas");
    emeas = Emeas->GetMean();
    s->Close();return 1.;
    TH1D *hYear = (TH1D*)file0->Get("Year");
    YEAR =  hYear->GetMean();
    TH1D *Nsim = (TH1D*)file0->Get("Nsim");
    if(Nsim->GetMean() > 1)expmc= false;
    Mpi0 = (TH1D*)file0->Get("Mpi0f");
    if(Mpi0->GetEntries() < 2){signalN = 0; dsignalN = 0;return 0.;}
    fit->SetParLimits(0,0.2*Mpi0->GetEntries(),Mpi0->GetEntries()+100.);
    fit->SetParLimits(9,-7,7);
    fit->SetParameter(0,Mpi0->GetEntries());
    Mpi0->Sumw2(true);
    //Mpi0->Fit("fit","RL","",120,160);                     // FIRST FIT ==============================
    //fit->ReleaseParameter(6);
    //if(Mpi0->GetEntries() < 400)fit->FixParameter(3,1);
    //if(Mpi0->GetEntries() > 100){fit->ReleaseParameter(7);}
    Mpi0->Fit("fit","RL");Mpi0->Fit("fit","RL");Mpi0->Fit("fit","RL");          // Second FIT ==================
    //cout << "Mpi0->GetEntries() = " << Mpi0->GetEntries() << endl;
    //if(Mpi0->GetEntries() > 400){
    //    fit->SetParLimits(8,-1,0.);
    //}
    //for(int i = 1; i <=5; i++){fit->FixParameter(i,fit->GetParameter(i));}
    //Mpi0->Fit("fit","RL");                   // THIRD FIT ==============================
    
    signalN = fit->GetParameter(0);
    dsignalN = fit->GetParError(0);
    Mpi0->SetTitle(Form("E = %g,   Nsignal = %g",emeas,signalN));
    fit->SetParameter(0,0.);
    fit->SetLineColor(1);
    fit->Draw("same");
    cout << "energy*2. = " << energy*2. << endl;
    s->SaveAs(("plots/"+filee+".png").c_str());
    s->Close();
    return 1.;
}





