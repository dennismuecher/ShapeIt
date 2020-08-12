#include "../Include/ShapeRho.h"

ShapeRho::ShapeRho(ShapeSetting* t_setting, ShapeMatrix* t_matrix) {
    m_sett = t_setting;
    m_matrix = t_matrix;
    Read();
}

void ShapeRho::Read() {
    rhoGraph = new TGraphErrors("../OsloFiles/leveldensity_76Ge_EB_ShapeIt.txt","%lg %lg %lg");
    rhoGraph->SetPoint(rhoGraph->GetN(), e_snail, rho_snail);
    rhoGraph->SetPointError(rhoGraph->GetN()-1, 0, drho_snail);
}

void ShapeRho::Draw() {
    rhoGraph->SetMarkerStyle(4);
    rhoGraph->SetMarkerColor(kRed);
    rhoGraph->SetTitle("level density; energy (MeV); level density (1/MeV)");
    rhoGraph->SetFillColorAlpha(4,0.5);
    rhoGraph->SetFillStyle(3010);
    rhoGraph->Draw("P3*");
}

void ShapeRho::Draw(double alpha, double alpha_l, double alpha_h) {

TGraphErrors* graph_t_mid = Transform(1,alpha);
TGraphErrors* graph_t_low =Transform(1,alpha_l);
TGraphErrors* graph_t_high = Transform(1,alpha_h);

TGraphAsymmErrors* graph_t = new TGraphAsymmErrors();

    for (int i=0; i < rhoGraph->GetN(); i++) {
        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i], graph_t_mid->GetY()[i]);

        double EY_l =TMath::Abs(graph_t_low->GetY()[i] -graph_t_mid->GetY()[i]);
        double EY_h =TMath::Abs(graph_t_high->GetY()[i] -graph_t_mid->GetY()[i]);
        if (graph_t_mid->GetX()[i] < e_snail) {
            EY_l = TMath::Sqrt(TMath::Power(EY_l,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
            EY_h = TMath::Sqrt(TMath::Power(EY_h,2) + TMath::Power(graph_t_mid->GetEY()[i],2));
        }
        graph_t->SetPointError(graph_t->GetN()-1,0,0,EY_l,EY_h);
    }
    graph_t->SetMarkerStyle(4);
    graph_t->SetMarkerColor(kBlack);
    graph_t->SetTitle("level density; energy (MeV); level density (1/MeV)");
    graph_t->SetFillColorAlpha(kRed,0.5);
    graph_t->SetFillStyle(3010);
    graph_t->Draw("aP3*");
}


TGraphErrors* ShapeRho::Transform(double A, double alpha) {
    TGraphErrors* graph_t = new TGraphErrors();
    //normalize at 2MeV
    double scale = 1 / TMath::Exp(alpha * 2);
    for (int i=0; i < rhoGraph->GetN(); i++) {
        double Y = scale * TMath::Exp(alpha * rhoGraph->GetX()[i]) * rhoGraph->GetY()[i];
        double EY = scale * TMath::Exp(alpha * rhoGraph->GetX()[i]) * rhoGraph->GetEY()[i];

        graph_t->SetPoint( graph_t->GetN(), rhoGraph->GetX()[i],Y);
        graph_t->SetPointError( graph_t->GetN()-1, 0, EY);

    }

    return graph_t;
}


double ShapeRho::Eval(double ene) {
    return rhoGraph->Eval(ene);
}


