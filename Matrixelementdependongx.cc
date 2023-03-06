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
double v1 = 0.001;
double v = 246;
//double mv = 0.05;
double mx = 0.1;
double mn = 0.9395;
//double gx = mv/v1;//from 0 to 1, log scale. Set gx first
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

double px = gax*mx*bex;
double Ex = pow(mx*mx + px*px, 0.5);
double s = pow(Ex + En, 2);
double px1 = pow(0.5*(s - mx*mx - mn*mn), 0.5);
double pn1 = pow(0.5*(s - mx*mx - mn*mn), 0.5);
double Ex1 = pow(mx*mx + px1*px1, 0.5);
double En1 = pow(mn*mn + pn1*pn1, 0.5);
double t = 2*mx*mx - 2*(Ex*Ex1 - px*px1*cos(30));
double u = mx*mx + mn*mn - 2*(Ex*Ex1 + px*px1*cos(30));


double mv(double *x, double *p)
{
  double val = x[0]*v1;
  return val;
}

double GX(double *x, double *p)
{
  double val = ( (4*x[0]*x[0]*e*e)/(pow(t - mv(x,p)*mv(x,p), 2)) )*(4*mx*mx*mn*mn*gt*ct + mx*mx*gt*cc*(t - 2*mn*mn) + mn*mn*gc*ct*(t - 2*mx*mx ) )
             + ( (2*x[0]*x[0]*e*e)/(pow(t - mv(x,p)*mv(x,p), 2)) )*gc*cc*( pow( s - mx*mx - mn*mn, 2) + pow(u - mx*mx - mn*mn, 2) )
             + ( (16*x[0]*x[0]*e*e)/(pow(t - mv(x,p)*mv(x,p), 2)) )*ga*ca*(6*gv*cv*pow(s - mx*mx - mn*mn, 2) + (( 2*ga*ca*mx*mx*mn*mn*t )/( mv(x,p)*mv(x,p) ) ) + (( ga*ca*mx*mx*mn*mn*t*t )/( mv(x,p)*mv(x,p)*mv(x,p)*mv(x,p) ))  );
    return val;
}
TF1 *GX_1 ;

double NGX(double *x, double *p)
{
  double val = ( (4*x[0]*x[0]*e*e)/(pow(t - mv(x,p)*mv(x,p), 2)) )*(4*mx*mx*mn*mn*gt*ct + mx*mx*gt*cc*(t - 2*mn*mn) + mn*mn*gc*ct*(t - 2*mx*mx ) )
             + ( (2*x[0]*x[0]*e*e)/(pow(t - mv(x,p)*mv(x,p), 2)) )*gc*cc*( pow( s - mx*mx - mn*mn, 2) + pow(u - mx*mx - mn*mn, 2) );
    return val;
}
TF1 *NGX_1;


void Matrixelementdependongx()
{
GX_1 = new TF1("", GX, 0, 1, 1);
NGX_1 = new TF1("", NGX, 0, 1, 1);
TCanvas *c1 = new TCanvas("c1","Matrix element depend on gx of DM", 900, 700);
gPad->SetTicks(1,1);
c1->SetLogy();

//	GR_1->SetNpx(100000);
GX_1->GetXaxis()->CenterTitle();
GX_1->GetYaxis()->CenterTitle();
GX_1->SetTitle("; g_{#chi};  #LT|M|^{2}#GT");
GX_1->SetLineColor(2);
GX_1->SetLineWidth(3);


GX_1->SetMaximum(1E19);
GX_1->SetMinimum(8E-6);

NGX_1->GetXaxis()->CenterTitle();
NGX_1->GetYaxis()->CenterTitle();
NGX_1->SetTitle("; g_{#chi}; #LT|M|^{2}#GT");
NGX_1->SetLineColor(4);
NGX_1->SetLineWidth(3);


GX_1->Draw();
NGX_1->Draw("same");
//  EQ_1->Draw("same");

TLegend *legend = new TLegend( 0.62, 0.62, 0.88, 0.89);
legend->SetHeader("m_{V} = (0 ; 1) MeV, v' = 1 MeV","l");
legend->AddEntry(GX_1, "g_{A}, c_{A} #neq 0");
legend->AddEntry(NGX_1,"g_{A}, c_{A} = 0");
//  legend->AddEntry(EQ_1,"KLM");
legend->SetBorderSize(0);
legend->Draw();
}
