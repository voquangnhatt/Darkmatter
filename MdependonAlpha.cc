#include <iostream>
using namespace std;
#include <string.h>
#include <math.h>
#include <TROOT.h>
#include <TF1.h>
#include <TCanvas.h>

//    double a = x[0];
    double mh2 = 10;
//    double gva = p[0];
    double v1 = 0.165;
    double v = 246;
    double mv = 0.05;
    double mx = 0.1;
    double mn = 0.9395;
    double gx = mv/v1;
    double yx = (pow(2, 0.5)*mx)/v1;
    double ghnn = 0.283;
//    double GF = 1.16637E-5;
    double mh1 = 125;
    double bex = 0.6;
    double ben = 0.6;
    double gax = 1/(pow(1- bex*bex, 0.5));
    double gan = 1/(pow(1- ben*ben, 0.5));
    double px = gax*mx*bex;
    double pn = gan*mn*ben;
    double Ex = pow(mx*mx + px*px, 0.5);
    double En = pow(mn*mn + pn*pn, 0.5);
    double s = pow(Ex + En, 2);
    double px1 = pow(0.5*(s - mx*mx - mn*mn), 0.5);
    double pn1 = pow(0.5*(s - mx*mx - mn*mn), 0.5);
    double Ex1 = pow(mx*mx + px1*px1, 0.5);
    double En1 = pow(mn*mn + pn1*pn1, 0.5);
    double the = 30;
    double t = 2*mx*mx - 2*(Ex*Ex1 - px*px1*cos(30));
    double u = mx*mx + mn*mn - 2*(Ex*Ex1 + px*px1*cos(30));
    double e = pow((4*3.141592654)/137, 0.5); //TMath::Pi()
    double D = (gx*e)/(t - mv*mv);


//    return ( (1./32.)*pow(sin(2*a), 2)*yx*yx*ghnn*ghnn*pow(h*h - mh1*mh1, 2)*(5*mx*mx -2*t)*(5*mn*mn -2*t) ) / ( pow(t - h*h, 2)*pow(t - mh1*mh1, 2) )  ;

    double EQ(double *x, double *p)
    {
      double val = ( (1./32.)*pow(sin(2*x[0]), 2)*yx*yx*ghnn*ghnn*pow(mh2*mh2 - mh1*mh1, 2)*(5*mx*mx -2*t)*(5*mn*mn -2*t) ) / ( pow(t - mh2*mh2, 2)*pow(t - mh1*mh1, 2) );
      	return val;
    }
    TF1 *EQ_1;


    double EP(double *x, double *p)
    {
    double val = 365.6;
    	return val;
    }
    TF1 *EP_1;


    double ET(double *x, double *p )
    {
    double val = 1.47;
      return val;
    }
    TF1 *ET_1;


void MdependonAlpha()
{
  EQ_1 = new TF1("", EQ, 0.0, 90.0, 1);
  EP_1 = new TF1("", EP, 0.0, 90.0, 1);
  ET_1 = new TF1("", ET, 0.0, 90.0, 1);
  TCanvas *c1 = new TCanvas("c1","Matrix element depend on alpha", 900, 700);
  gPad->SetTicks(1,1);
  c1->SetLogy();

	EQ_1->SetNpx(100000);
  EQ_1->GetXaxis()->CenterTitle();
  EQ_1->GetYaxis()->CenterTitle();
  EQ_1->SetTitle("; #alpha; Matrix element (#LT|M|^{2}#GT)");
  EQ_1->SetLineColor(2);
  EQ_1->SetLineWidth(3);
  EQ_1->SetMaximum(1E4);
  EQ_1->SetMinimum(1E-20);

  EP_1->GetXaxis()->CenterTitle();
  EP_1->GetYaxis()->CenterTitle();
  EP_1->SetTitle("; #alpha; Matrix element (#LT|M|^{2}#GT)");
  EP_1->SetLineColor(4);
  EP_1->SetLineWidth(3);

  ET_1->GetXaxis()->CenterTitle();
  ET_1->GetYaxis()->CenterTitle();
  ET_1->SetTitle("; #alpha; Matrix element (#LT|M|^{2}#GT)");
  ET_1->SetLineColor(1);
  ET_1->SetLineWidth(3);


  EQ_1->Draw();
  ET_1->Draw("same");
  EP_1->Draw("same");
//  EQ_1->Draw("same");

  TLegend *legend = new TLegend( 0.65, 0.59, 0.85, 0.75);
  legend->SetHeader("m_{h_{2}} = 10 (GeV)","c");
  legend->AddEntry(EQ_1, "scalar propagator");
  legend->AddEntry(EP_1,"vector propagator");
  legend->AddEntry(ET_1,"vector propagator, c_{A} = 0");
//  legend->AddEntry(EQ_1,"KLM");
  legend->SetBorderSize(0);
  legend->Draw();
}
