#include <iostream>
using namespace std;
#include <string.h>
#include <math.h>
#include <TROOT.h>
#include <TF1.h>
#include <TCanvas.h>

    double gv = 1;
    double cv = 1;
    double ca = 1;
    double ga = 1;
    double gc = gv*gv + ga*ga;
    double gt = gv*gv - ga*ga;
    double cc = cv*cv + ca*ca;
    double ct = cv*cv - ca*ca;
    double al = 40;
    double v1 = 0.2;
    double v = 246;
    double mv = 0.05;
    double mx = 0.1;
    double mn = 0.9395;
    double gx = mv/v1;
    double ghnn = 0.283;
    double mh1 = 125;
    double bex = 0.6;
    double ben = 0.6;
    double the = 30;
    double e = pow((4*3.141592654)/137, 0.5); //TMath::Pi()
    double gax = 1/(pow(1- bex*bex, 0.5));
    double gan = 1/(pow(1- ben*ben, 0.5));
    double pn = gan*mn*ben;
    double En = pow(mn*mn + pn*pn, 0.5);


    double px(double *x, double *p)
    {
      double val = gax*x[0]*bex;
      return val;
    }


    double Ex(double *x, double *p)
    {
      double val = pow(x[0]*x[0] + px(x,p)*px(x,p), 0.5);
      return val;
    }


    double s(double *x, double *p)
    {
      double val = pow(Ex(x,p) + En, 2);
      return val;
    }

    double px1(double *x, double *p)
    {
      double val = pow(0.5*(s(x,p) - x[0]*x[0] - mn*mn), 0.5);
      return val;
    }

    double pn1(double *x, double *p)
    {
      double val = pow(0.5*(s(x,p) - x[0]*x[0] - mn*mn), 0.5);
      return val;
    }

    double Ex1(double *x, double *p)
    {
      double val = pow(x[0]*x[0] + px1(x,p)*px1(x,p), 0.5);
      return val;
    }


    double En1(double *x, double *p)
    {
      double val = pow(mn*mn + pn1(x,p)*pn1(x,p), 0.5);
      return val;
    }


    double t(double *x, double *p)
    {
      double val = 2*x[0]*x[0] - 2*(Ex(x,p)*Ex1(x,p) - px(x,p)*px1(x,p)*cos(30));
      return val;
    }


    double u(double *x, double *p)
    {
      double val = x[0]*x[0] + mn*mn - 2*(Ex(x,p)*Ex1(x,p) + px(x,p)*px1(x,p)*cos(30));
      return val;
    }

//    double px = gax*x[0]*bex;
//    double Ex = pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5);
//    double s = pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2);
//    double px1 = pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5);
//    double pn1 = pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5);
//    double Ex1 = pow(x[0]*x[0] + pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5), 0.5);
//    double En1 = pow(mn*mn + pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5), 0.5);
//    double t = 2*x[0]*x[0] - 2*(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5)*pow(x[0]*x[0] + pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5), 0.5) - gax*x[0]*bex*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*cos(30));
//    double u = x[0]*x[0] + mn*mn - 2*(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5)*pow(x[0]*x[0] + pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5), 0.5) + gax*x[0]*bex*pow(0.5*(pow(pow(x[0]*x[0] + gax*x[0]*bex*gax*x[0]*bex, 0.5) + En, 2) - x[0]*x[0] - mn*mn), 0.5)*cos(30));

    double GR(double *x, double *p)
    {
      double val = ( (px(x,p)*4*gx*gx*e*e)/(16*3.141592654*s(x,p)*px1(x,p)*pow(t(x,p) - mv*mv, 2)) )*(4*x[0]*x[0]*mn*mn*gt*ct + x[0]*x[0]*gt*cc*(t(x,p) - 2*mn*mn) + mn*mn*gc*ct*(t(x,p) - 2*x[0]*x[0] ) )
                 + ( (px(x,p)*2*gx*gx*e*e)/(16*3.141592654*s(x,p)*px1(x,p)*pow(t(x,p) - mv*mv, 2)) )*gc*cc*( pow( s(x,p) - x[0]*x[0] - mn*mn, 2) + pow(u(x,p) - x[0]*x[0] - mn*mn, 2) )
                 + ( (px(x,p)*16*gx*gx*e*e)/(16*3.141592654*s(x,p)*px1(x,p)*pow(t(x,p) - mv*mv, 2)) )*ga*ca*(6*gv*cv*pow(s(x,p) - x[0]*x[0] - mn*mn, 2) + (( 2*ga*ca*x[0]*x[0]*mn*mn*t(x,p) )/(mv*mv)) + (( ga*ca*x[0]*x[0]*mn*mn*t(x,p)*t(x,p) )/(mv*mv*mv*mv))  );
        return val;
    }
    TF1 *GR_1;


    double NGR(double *x, double *p)
    {
      double val = ( (px(x,p)*4*gx*gx*e*e)/(16*3.141592654*s(x,p)*px1(x,p)*pow(t(x,p) - mv*mv, 2)) )*(4*x[0]*x[0]*mn*mn*gt*ct + x[0]*x[0]*gt*cc*(t(x,p) - 2*mn*mn) + mn*mn*gc*ct*(t(x,p) - 2*x[0]*x[0] ) )
                 + ( (px(x,p)*2*gx*gx*e*e)/(16*3.141592654*s(x,p)*px1(x,p)*pow(t(x,p) - mv*mv, 2)) )*gc*cc*( pow( s(x,p) - x[0]*x[0] - mn*mn, 2) + pow(u(x,p) - x[0]*x[0] - mn*mn, 2) );
        return val;
    }
    TF1 *NGR_1;


void Crosssectiondependonmx()
{
  GR_1 = new TF1("", GR, 0, 400, 1);
  NGR_1 = new TF1("", NGR, 0, 400, 1);
  TCanvas *c1 = new TCanvas("c1","Cross section depenon mass of DM", 900, 700);
  gPad->SetTicks(1,1);
  c1->SetLogy();

//	GR_1->SetNpx(100000);
  GR_1->GetXaxis()->CenterTitle();
  GR_1->GetYaxis()->CenterTitle();
  GR_1->SetTitle("; m_{#chi} (GeV); #sigma");
  GR_1->SetLineColor(2);
  GR_1->SetLineWidth(3);
  GR_1->SetMaximum(5E3);
  GR_1->SetMinimum(1E-7);

  NGR_1->GetXaxis()->CenterTitle();
  NGR_1->GetYaxis()->CenterTitle();
  NGR_1->SetTitle("; m_{#chi} (GeV); #sigma");
  NGR_1->SetLineColor(4);
  NGR_1->SetLineWidth(3);


  GR_1->Draw();
  NGR_1->Draw("same");
//  EQ_1->Draw("same");

  TLegend *legend = new TLegend( 0.65, 0.55, 0.85, 0.79);
  legend->SetHeader("m_{V} = 50 MeV, v' = 200 MeV","l");
  legend->AddEntry(GR_1, "g_{A}, c_{A} #neq 0");
  legend->AddEntry(NGR_1,"g_{A}, c_{A} = 0");
//  legend->AddEntry(EQ_1,"KLM");
  legend->SetBorderSize(0);
  legend->Draw();
}
