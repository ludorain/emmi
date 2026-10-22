//#include "/home/preghenella/cernbox/EIC/sipm4eic-2023-characterisation/utils/graphutils.C"
//#include "/home/preghenella/EIC/sipm4eic-laser-analysis/macros/qGaus.C"

#include "makeiv.C"

namespace utils {
  TGraphErrors *log(TGraphErrors *gin);
  TGraphErrors *derivate(TGraphErrors *gin);
  std::array<double, 2>  max(TGraphErrors *gin);
  void gstyle(TGraph *g, int color);
  void fstyle(TF1 *f, int color);
}

namespace qstat {
  TF1 *qtof(double norm, double mean, double sigma, double qR, double qL = 1., double min = -1000., double max = 1000.);
}

void ivscan_vbreak(std::string board, std::string channel, std::vector<int> temps, float vmin = 48., float vmax = 55.);

void
ivscan_vbreak(std::string board, std::string channel, std::vector<int> temps, float vmin = 48., float vmax = 55.)
{

  gStyle->SetOptFit();

  auto gout = new TGraphErrors;
  
  for (auto temp : temps) {

    auto tempK = int(temp + 273);
    std::string fname_ivscan = "HAMA3_sn" + board + "_" + std::to_string(tempK) + "K_" + channel + "_zoom.ivscan.csv";
    std::string fname_zero   = "HAMA3_sn" + board + "_" + std::to_string(tempK) + "K_" + channel + ".zero.csv";
    
    auto gin = makeiv(fname_ivscan, fname_zero);
    auto gld = utils::derivate(utils::log(gin));
    auto max = utils::max(gld);
    auto fun = qstat::qtof(max[1], max[0], 0.15, 3., 1., vmin, vmax);
    gld->Fit(fun, "0qIMRE", "", max[0] - 1., max[0] + 1.);
    auto vbreak = fun->GetParameter(1);
    auto vbreak_err = fun->GetParError(1);

    auto c = new TCanvas(fname_ivscan.c_str(), fname_ivscan.c_str(), 800, 800);
    c->SetMargin(0.15, 0.15, 0.15, 0.15);
    auto h = c->DrawFrame(vmin, -2., vmax, max[1] * 1.1, Form("T = %d #circC;bias voltage (V); LD", temp));
    h->GetXaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetTitleOffset(1.5);
    
    utils::gstyle(gld, kAzure-3);
    gld->Draw("samep");
    utils::fstyle(fun, kBlack);
    fun->Draw("same");
    c->SaveAs((fname_ivscan + ".ldfit.png").c_str());

    auto n = gout->GetN();
    gout->SetPoint(n, temp, vbreak);
    gout->SetPointError(n, 0., vbreak_err);

  }

  auto fpol1 = (TF1 *)gROOT->GetFunction("pol1");  
  utils::fstyle(fpol1, kBlack);
  TFitResultPtr fitr = gout->Fit(fpol1, "S");
  TMatrixDSym cov = fitr->GetCovarianceMatrix();
  auto vbmin = fpol1->Eval(10.);
  auto vbmax = fpol1->Eval(30.);

  auto c = new TCanvas("ivscan_vbreak", "ivscan_vbreak", 800, 800);
  c->SetMargin(0.15, 0.15, 0.15, 0.15);
  auto h = c->DrawFrame(10., vbmin, 30., vbmax, ";temperature (#circC);breakdown voltage (V);");
  h->GetXaxis()->SetTitleOffset(1.5);
  h->GetYaxis()->SetTitleOffset(1.5);
  utils::gstyle(gout, kAzure-3);
  gout->Draw("samep");
  c->SaveAs("ivscan_vbreak.png");

  std::ofstream fout;
  fout.open("ivscan_vbreak.txt");
  fout << "#temp vbreak vbreak_err " << std::endl;
  for (int i = 0; i < gout->GetN(); ++i) {
    fout << gout->GetX()[i] << " "
	 << gout->GetY()[i] << " "
	 << gout->GetEY()[i]
	 << std::endl;
  }
  fout.close();

  fout.open("ivscan_vbreak_params.txt");
  fout << "#p0 p1 cov00 cov11 cov01 chi2 ndf" << std::endl;
  fout << fpol1->GetParameter(0) << " "
       << fpol1->GetParameter(1) << " "
       << cov(0,0) << " "
       << cov(1,1) << " "
       << cov(0,1) << " "
       << fpol1->GetChisquare() << " "
       << fpol1->GetNDF() << " "
       << std::endl;
  fout.close();
  
}

namespace utils {
  
  TGraphErrors *
  log(TGraphErrors *gin)
  {
    auto g = new TGraphErrors;
    for (int i = 0; i < gin->GetN(); ++i) {
      auto x  = gin->GetX()[i];
      auto y  = gin->GetY()[i];
      auto ex = gin->GetEX()[i];
      auto ey = gin->GetEY()[i];
      g->SetPoint(i, x, std::log(y));
      g->SetPointError(i, ex, ey / y);
    }
    return g;
  }
  
  TGraphErrors *
  derivate(TGraphErrors *gin)
  {
    auto g = new TGraphErrors;
    for (int i = 0; i < gin->GetN() - 1; ++i) {
      auto x0  = gin->GetX()[i];
      auto y0  = gin->GetY()[i];
      auto ey0 = gin->GetEY()[i];
      auto x1  = gin->GetX()[i + 1];
      auto y1  = gin->GetY()[i + 1];
      auto ey1 = gin->GetEY()[i + 1];
      auto dx  = x1 - x0;
      auto dy  = y1 - y0;
      auto x   = 0.5 * (x0 + x1);
      auto ex  = 0.5 * dx;    
      auto y   = (y1 - y0) / dx;  
      auto ey  = std::hypot(ey0, ey1) / dx;
      g->SetPoint(i, x, y);
      g->SetPointError(i, ex, ey);
    }
    return g;
  }

  std::array<double, 2>
  max(TGraphErrors *gin)
  {
    bool init = false;
    std::array<double, 2> point;
      for (int i = 0; i < gin->GetN(); ++i) {
	auto x = gin->GetX()[i];
	auto y = gin->GetY()[i];
	if (!init || y > point[1]) {
	  init = true;
	  point = {x, y};
	}
      }
    return point;
  }
  
  
  void 
  gstyle(TGraph *g, int color)
  {
    g->SetLineWidth(2);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
  }
  
  void
  fstyle(TF1 *f, int color)
  {
    f->SetLineWidth(2);
    f->SetLineStyle(2);
    f->SetLineColor(color);
    f->SetNpx(1000);
  }

}

namespace qstat {
  
  double
  qexp(double x, double q)
  {
    if (std::fabs(1. - q) < 1.e-6) return std::exp(x);
    if ((1. + (1. - q) * x) > 0.) return std::pow(1. + (1. - q) * x, 1. / (1. - q));
    return std::pow(0., 1. / (1. - q));
  }

  double_t
  qtof(double *_x, double *_par)
  {
    auto norm  = _par[0];
    auto mean  = _par[1];
    auto sigma = _par[2];
    auto qR    = _par[3];
    auto qL    = _par[4];
    auto x     = _x[0] - mean;
    auto q     = x < 0. ? qL : qR;
    auto beta  = 1. / (2. * sigma * sigma);
    return norm * qexp(-beta * x * x, q);
  }

  TF1 *
  qtof(double norm, double mean, double sigma, double qR, double qL = 1., double min = -1000., double max = 1000.)
  { 
    TF1 *f = new TF1("qtof", qtof, min, max, 5);
    f->SetParameter(0, norm);
    f->SetParameter(1, mean);
    f->SetParameter(2, sigma);
    f->SetParameter(3, qR);
    f->SetParameter(4, qL);
    f->SetParLimits(3, 0., 100.);
    f->SetParLimits(4, 0., 100.);
    return f;
  }

}
