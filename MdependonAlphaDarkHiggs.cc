#include <iostream>
using namespace std;
#include <string.h>
#include <math.h>
#include <TROOT.h>
#include <TF1.h>
#include <TCanvas.h>


double Matrix(double *x, double *p)
{
    long double a = x[0];
    double h = x[1];
    double v1 = 0.165;
    double v = 246;
    double mv = 0.05;
    double mx = 0.1;
    double mn = 0.9395;
    double gx = mv/v1;
    double yx = (pow(2, 0.5)*mx)/v1;
    double ghnn = 0.283;
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
    double t = 2*mx*mx - 2*(Ex*Ex1 - px*px1*cos(30));
    double u = mx*mx + mn*mn - 2*(Ex*Ex1 + px*px1*cos(30));
    double e = pow((4*3.141592654)/137, 0.5); //TMath::Pi()

    return ( (1./32.)*pow(sin(2*a), 2)*yx*yx*ghnn*ghnn*pow(h*h - mh1*mh1, 2)*(5*mx*mx -2*t)*(5*mn*mn -2*t) ) / ( pow(t - h*h, 2)*pow(t - mh1*mh1, 2) );
}

void MdependonAlphaDarkHiggs()
{
    TF2 *f2 = new TF2("", Matrix, 0, 90.0, 0.0, 200.0, 1);
    f2->SetTitle("; #alpha ; m_{h_{2}} (GeV); Matrix element (#LT|M|^{2}#GT)");
    f2->GetXaxis()->CenterTitle();
    f2->GetYaxis()->CenterTitle();
    f2->GetZaxis()->CenterTitle();
    f2->Draw("surf2Z");
}
