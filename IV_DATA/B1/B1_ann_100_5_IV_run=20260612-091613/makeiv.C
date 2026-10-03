#pragma once

#include "graphutils.C"

TGraphErrors *
makeiv(const std::string fnivscan, const std::string fnzero, bool invertX = true, bool invertY = true, bool write = true)
{

  /** compute zero volt level **/
  auto zero = 0.;
  auto ezero = 0.;
  if (!fnzero.empty()) {
    std::cout << " --- open zero volt level file: " << fnzero << std::endl;
    auto gzero = new TGraph(fnzero.c_str(), "%lg %*lg %lg");
    gzero = invertY ? graphutils::invertY(gzero) : gzero;
    auto pzero = new TProfile("pzero", "", 1, 0., 1.);
    for (int i = 0; i < gzero->GetN(); ++i)
      pzero->Fill(0., gzero->GetY()[i]);
    zero = pzero->GetBinContent(1);
    ezero = pzero->GetBinError(1);
  }
  std::cout << " --- zero volt level: " << zero * 1.e12 << " +- " << ezero * 1.e12 << " pA" << std::endl;

  /** IV curve from average of all measurements **/
  std::cout << " --- open iv scan file: " << fnivscan << std::endl;
  auto givscan = new TGraph(fnivscan.c_str(), "%*lg %lg %lg");
  givscan = invertX ? graphutils::invertX(givscan) : givscan;
  givscan = invertY ? graphutils::invertY(givscan) : givscan;
  auto pivscan = new TProfile("pivscan", "", 100001, -0.0005, 100.0005);
  for (int i = 0; i < givscan->GetN(); ++i)
    pivscan->Fill(givscan->GetX()[i], givscan->GetY()[i]);

  /** make TGraphErrors **/
  auto hivscan = new TH1F("hivscan", "", 100001, -0.0005, 100.0005);
  auto gout = new TGraphErrors;
  gout->SetName(fnivscan.c_str());
  auto npoints = 0;
  for (int i = 0; i < pivscan->GetNbinsX(); ++i) {
    if (pivscan->GetBinError(i + 1) == 0.) continue;
    auto cen = pivscan->GetBinCenter(i + 1);
    auto val = pivscan->GetBinContent(i + 1);
    auto eval = pivscan->GetBinError(i + 1);
    val = val - zero;
    eval = std::sqrt(eval * eval + ezero * ezero);
    gout->SetPoint(npoints, cen, val);
    gout->SetPointError(npoints, 0., eval);
    hivscan->SetBinContent(i + 1, val);
    hivscan->SetBinError(i + 1, eval);
    npoints++;
  }

  if (!write) return gout;

  /** write to file **/
  auto lastindex = fnivscan.find_last_of("."); 
  auto fnout = fnivscan.substr(0, lastindex) + ".makeiv.root";
  std::cout << " --- writing output file: " << fnout << std::endl;
  auto fout = TFile::Open(fnout.c_str(), "RECREATE");
  gout->Write("ivscan");
  hivscan->Write("hivscan");
  fout->Close();
  
  return gout;
}

