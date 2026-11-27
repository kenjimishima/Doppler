#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>
#include <fstream>
#include "PlotVoigtAll.C"

Double_t A40Ca = 0.969;
Double_t A42Ca = 0.00647;
Double_t A44Ca = 0.0209;
Double_t A48Ca = 0.00187;

Double_t Mu40Ca = 422.79179;
Double_t Mu42Ca = 422.79156;
Double_t Mu44Ca = 422.79133;
Double_t Mu48Ca = 422.79089;

//Wikipedia...
Double_t mass_40Ca = 39.96259098;
Double_t mass_42Ca = 41.95861801;
Double_t mass_44Ca = 43.9554818;
Double_t mass_48Ca = 47.952534;
const Int_t num_ca = 4;
Double_t ca_mass[num_ca] = {mass_40Ca, mass_42Ca, mass_44Ca, mass_48Ca};

//Use mass correction
//Bool_t bMassCorrection = kTRUE;
Bool_t bMassCorrection = kFALSE;

std::string StripExtAndDir(const std::string& s)
{
  // ディレクトリ区切りの位置
  size_t p = s.find_last_of("/\\");
  std::string base = (p == std::string::npos) ? s : s.substr(p + 1);
  // 拡張子の位置
  size_t dot = base.find_last_of('.');
  if (dot == std::string::npos) return base;  
  return base.substr(0, dot);
}

Double_t WavelengthNmToFreqMHz(double lambda_nm)
{
  // TMath::C() : speed of light [m/s]
  double lambda_m = lambda_nm * 1e-9;       // nm → m
  double freq_Hz  = TMath::C() / lambda_m;  // Hz
  return freq_Hz * 1e-6;                    // Hz → MHz
}

Double_t FreqMHzToWavelengthNm(double freq_MHz)
{
  double freq_Hz = freq_MHz * 1e6;          // MHz → Hz
  double lambda_m = TMath::C() / freq_Hz;   // m
  return lambda_m * 1e9;                    // m → nm
}

Double_t VoigtFunc(Double_t *x, Double_t *par)
{
  Double_t A      = par[0];
  Double_t mu     = par[1];
  Double_t sigma  = par[2];
  Double_t gamma  = par[3];
  Double_t bg     = par[4];
  Double_t xx     = x[0] - mu;
  Double_t voigt  = TMath::Voigt(xx, par[2], par[3]);
  return par[0] * voigt + par[4];
}

// domain-conversion bugs.
// par layout: [0]=sigma, [1]=gamma, [2]=bg, [3..6]=A40..A48, [7..10]=mu40..mu48
Double_t Voigt4(Double_t *x, Double_t *par)
{
  Double_t sum = 0.;
  const int num = 4;
  Double_t sigma = par[0];
  Double_t gamma = par[1];
  Double_t bg    = par[2];
  for (int i = 0; i < 4; i++) {
    Double_t A  = par[3+i];
    Double_t mu = par[3+num+i];
    Double_t xx = x[0] - mu;
    if(bMassCorrection) sigma = par[0]*TMath::Sqrt(ca_mass[0]/ca_mass[i]); //Mass correction
    sum += A * TMath::Voigt(xx, par[0], par[1]);
  }
  return sum + bg;
}

void EditFitCaVoigt(const char* filename = "RUN81_500_OD2.7_E0.5s_BL0.25_Xsm3_nolens_BGC.txt")
{
  std::cout << std::fixed << std::setprecision(2);

  std::ifstream infile(filename);
  if (!infile.is_open()) {
    Error("FitCaVoigt", "Cannot open input file: %s", filename);
    return;
  }

  std::string basename = StripExtAndDir(filename);

  // Save directory = same name as basename (without extension)
  std::string rundir = basename;

  std::vector<double> x, y, ex, ey; //wavelength
  std::vector<double> fx, efx; //frequency
  double xv, yv;
  while (infile >> xv >> yv) {
    x.push_back(xv);
    y.push_back(yv);
    ex.push_back(0.0);          // x誤差なし
    //    ey.push_back(0.01 * yv);    // y誤差 1%
    ey.push_back(TMath::Sqrt(yv) );   // sqrt 誤差
    fx.push_back( WavelengthNmToFreqMHz(xv) );
    efx.push_back(0.0); 
  }

  auto minmax = std::minmax_element(y.begin(), y.end());
  double ymin = *minmax.first;
  double ymax = *minmax.second;
  
  TGraphErrors *gr = new TGraphErrors(x.size(),&x[0], &y[0],&ex[0], &ey[0]);
  gr->SetTitle(Form("Calcium fluorescence spectrum / %s; Wavelength [nm] ; Count / 0.5 s",filename));
  gr->SetMarkerStyle(20);

  // --- キャンバス ---
  auto c0 = new TCanvas("c0", "Voigt Fit - Ca Isotopes", 1000, 700);
  gr->Draw("AP");
  gPad->SetLogy();
  gPad->SetGrid();
  // --- 軸ラベルサイズをさらに小さく ---
  gr->GetXaxis()->SetTitleSize(0.035);
  gr->GetXaxis()->SetLabelSize(0.025);
  gr->GetYaxis()->SetTitleSize(0.035);
  gr->GetYaxis()->SetLabelSize(0.035);
  gr->GetYaxis()->SetRangeUser(ymin*0.3, ymax*3.);

  // パラメータ初期値
  Double_t peakshift = 0.;
  double sigma0 = 3.e-5;
  double gamma0 = 3.e-5;
  // cutoff/erf_sig: not used in fit, retained as display variables
  double cutoff0  = 0.;
  double erf_sig0 = 0.;
  Double_t tolerance = 6.0e-8; //in ratio
  double peak = ymax*sigma0*TMath::Sqrt(TMath::TwoPi())*TMath::Pi();
  double bg0 = y[0];

  // Single-peak pre-fit: VoigtFunc (5 par: A, mean, sigma, gamma, bg)
  TF1 *fit_single = new TF1("fit_single", VoigtFunc, Mu40Ca*(1-8*tolerance), Mu40Ca*(1+50*tolerance), 5);
  fit_single->SetParNames("A", "mean", "sigma", "gamma", "background");
  fit_single->SetParameter(0, peak);
  fit_single->SetParameter(1, Mu40Ca);
  fit_single->SetParameter(2, sigma0);
  fit_single->SetParameter(3, gamma0);
  fit_single->SetParameter(4, bg0);
  fit_single->SetParLimits(0, 0., 1.e10);
  fit_single->SetParLimits(1, Mu40Ca*(1-10*tolerance), Mu40Ca*(1+10*tolerance));
  fit_single->SetParLimits(2, sigma0*0.1, sigma0*10.);
  fit_single->SetParLimits(3, gamma0*0.1, gamma0*10.);
  fit_single->SetParLimits(4, bg0*0.3, bg0*3.);
  fit_single->SetNpx(1000);
  //  fit_single->Draw();
  //  return;
  gr->Fit(fit_single, "R");
  peak      = fit_single->GetParameter(0);
  peakshift = fit_single->GetParameter(1) - Mu40Ca;
  sigma0    = fit_single->GetParameter(2);
  gamma0    = fit_single->GetParameter(3);
  bg0       = fit_single->GetParameter(4);
  cout << "40Ca peak shift = " << peakshift << endl;
  //  return;

  auto c1 = new TCanvas("c1", "Voigt Fit - Ca Isotopes", 1000, 700);
  gr->Draw("AP");
  gPad->SetLogy();
  gPad->SetGrid();

  // Voigt4: 11 par (sigma, gamma, bg, A40-48, mu40-48)
  TF1 *fit = new TF1("fit", Voigt4, x.front(), x.back(), 11);
  Double_t par[11] = {
    sigma0,
    gamma0,
    bg0,   // bg: constant background
    peak*A40Ca,
    peak*A42Ca,
    peak*A44Ca,
    peak*A48Ca,
    Mu40Ca+peakshift,
    Mu42Ca+peakshift,
    Mu44Ca+peakshift,
    Mu48Ca+peakshift
  };
  fit->SetParameters(par);

  fit->SetParName(0,  "sigma");
  fit->SetParName(1,  "gamma");
  fit->SetParName(2,  "background");
  fit->SetParName(3,  "A40");
  fit->SetParName(4,  "A42");
  fit->SetParName(5,  "A44");
  fit->SetParName(6,  "A48");
  fit->SetParName(7,  "mu40");
  fit->SetParName(8,  "mu42");
  fit->SetParName(9,  "mu44");
  fit->SetParName(10, "mu48");
  /////////
  fit->SetParLimits(0, sigma0*0.8, sigma0*1.2);
  fit->SetParLimits(1, gamma0*0.8, gamma0*1.2);
  fit->SetParLimits(2, bg0*0.5, bg0*2.);          // background
  fit->SetParLimits(3,  peak*A40Ca*0.1, peak*A40Ca*5);
  fit->SetParLimits(4,  peak*A42Ca*0.1, peak*A42Ca*5);
  fit->SetParLimits(5,  peak*A44Ca*0.1, peak*A44Ca*5);
  fit->SetParLimits(6,  peak*A48Ca*0.1, peak*A48Ca*5);
  fit->SetParLimits(7,  (Mu40Ca+peakshift)*(1 - tolerance), (Mu40Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(8,  (Mu42Ca+peakshift)*(1 - tolerance), (Mu42Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(9,  (Mu44Ca+peakshift)*(1 - tolerance), (Mu44Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(10, (Mu48Ca+peakshift)*(1 - tolerance), (Mu48Ca+peakshift)*(1 + tolerance));
  fit->SetLineColor(kRed);
  fit->SetNpx(1000);

  // フィット実行
  gr->Fit(fit, "R");
  //  fit->Draw("");
  Double_t p40 = fit->GetParameter(7);
  Double_t p42 = fit->GetParameter(8);
  Double_t p44 = fit->GetParameter(9);
  Double_t p48 = fit->GetParameter(10);
  //  printf("%f %f %f %f\n", p40, p42, p44, p48);

  // 共通パラメータ
  double sigma  = fit->GetParameter(0);
  double gamma  = fit->GetParameter(1);
  double cutoff  = 0.; // no cutoff
  double erf_sig = 0.; // no cutoff
  Double_t sigmaFr  = WavelengthNmToFreqMHz(Mu40Ca) * sigma   / Mu40Ca;
  Double_t gammaFr  = WavelengthNmToFreqMHz(Mu40Ca) * gamma   / Mu40Ca;
  Double_t cutoffFr  = 0.; // no cutoff
  Double_t erf_sigFr = 0.; // no cutoff
  double background = fit->GetParameter(2);

  // 共通パラメータ
  cout << "sigma  = " << sigma  << " [nm] = " << sigmaFr  << " [MHz]" << endl;
  cout << "gamma  = " << gamma  << " [nm] = " << gammaFr  << " [MHz]" << endl;
  cout << "cutoff = " << cutoff << " [nm] = " << cutoffFr << " [MHz]  (no cutoff)" << endl;
  cout << "erf_sigma = " << erf_sig << " [nm] = " << erf_sigFr << " [MHz]  (no cutoff)" << endl;

  TLegend *leg = new TLegend(0.15, 0.65, 0.38, 0.88);  // (x1, y1, x2, y2)
  leg->AddEntry(gr, "Data", "p");
  leg->AddEntry(fit, "Voigt fit", "l");

  const int num = 4;
  int colors[num] = {kBlue, kGreen+2, kMagenta, kOrange+1};
  const char* names[num] = {"^{40}Ca", "^{42}Ca", "^{44}Ca", "^{48}Ca"};
  
  for (int i = 0; i < num; i++) {    
    Double_t A  = fit->GetParameter(3 + i);       // A40,42,44,48
    Double_t mu = fit->GetParameter(3 + num + i); // mu40,42,44,48    
    TF1 *f = new TF1(
		     Form("f_%s", names[i]),
		     VoigtFunc,
		     x.front(), x.back(),
		     5
		     );
    f->SetNpx(1000);
    Double_t sigma_corr = sigma;
    if(bMassCorrection) sigma_corr = sigma*TMath::Sqrt(ca_mass[0]/ca_mass[i]); //Mass correction
    f->SetParameters(A, mu, sigma_corr, gamma, 0.);
    f->SetLineColor(colors[i]);
    f->SetLineStyle(2);
    f->SetLineWidth(2);
    f->Draw("same");
    leg->AddEntry(f, names[i], "l");
  }
  leg->Draw();

  if (gSystem->AccessPathName(rundir.c_str())) gSystem->mkdir(rundir.c_str(), true);
  c1->Update();
  c1->SaveAs(Form("./%s/Fit_%s.root", rundir.c_str(), basename.c_str()));
  c1->SaveAs(Form("./%s/Fit_%s.pdf",  rundir.c_str(), basename.c_str()));
  c1->SaveAs(Form("./%s/Fit_%s.png",  rundir.c_str(), basename.c_str()));
  //  return;

  // Frequency
  auto c2 = new TCanvas("c2", "Voigt Fit - Ca Isotopes ", 1000, 700);
  Double_t offset = WavelengthNmToFreqMHz(p40);
  int n = fx.size();
  for (int i = 0; i < n; i++) fx[i] -= offset;
  TGraphErrors *grf = new TGraphErrors(n, &fx[0], &y[0], &efx[0], &ey[0]);
  grf->SetTitle(Form("Calcium fluorescence spectrum / %s; Frequency [MHz] ; Count / 0.5 s",filename));
  grf->SetMarkerStyle(20);

  grf->Draw("AP");
  gPad->SetLogy();
  gPad->SetGrid();
  //  leg->Draw();
  grf->GetYaxis()->SetRangeUser(ymin*0.3, ymax*3.);

  Double_t f40 = WavelengthNmToFreqMHz(p40);
  Double_t f42 = WavelengthNmToFreqMHz(p42);
  Double_t f44 = WavelengthNmToFreqMHz(p44);
  Double_t f48 = WavelengthNmToFreqMHz(p48);

  peak *= f40/Mu40Ca; // convert units  
  // fx is built from wavelength, which is monotonically increasing while
  // frequency is inversely related, so fx.front() > fx.back() (range must
  // be passed in the correct min,max order for TF1).
  TF1 *fitFr = new TF1("fitFr", Voigt4, fx.back(), fx.front(), 11);
  // [FIX] background is a constant, so it carries over to the frequency
  // domain unchanged - no slope conversion needed (this is what removes
  // the source of the earlier mismatch).
  Double_t parFr[11] = {
    sigmaFr,
    gammaFr,
    background,    // background (unchanged - constant, domain-independent)
    peak*A40Ca,
    peak*A42Ca,
    peak*A44Ca,
    peak*A48Ca,
    f40-f40,
    f42-f40,
    f44-f40,
    f48-f40
  };
  for (int i = 0; i < 11; i++) cout << i << " " << parFr[i] << endl;
  fitFr->SetParameters(parFr);
  /////////
  fitFr->SetParName(0,  "sigma");
  fitFr->SetParName(1,  "gamma");
  fitFr->SetParName(2,  "background");
  fitFr->SetParName(3,  "A40");
  fitFr->SetParName(4,  "A42");
  fitFr->SetParName(5,  "A44");
  fitFr->SetParName(6,  "A48");
  fitFr->SetParName(7,  "Fr40");
  fitFr->SetParName(8,  "Fr42");
  fitFr->SetParName(9,  "Fr44");
  fitFr->SetParName(10, "Fr48");
  fitFr->SetParLimits(0, sigmaFr*0.8,  sigmaFr*1.2);
  fitFr->SetParLimits(1, gammaFr*0.8,  gammaFr*1.2);
  fitFr->SetParLimits(2, background*0.5, background*2.); // background
  fitFr->SetParLimits(3, peak*A40Ca*0.1, peak*A40Ca*5);
  fitFr->SetParLimits(4, peak*A42Ca*0.1, peak*A42Ca*5);
  fitFr->SetParLimits(5, peak*A44Ca*0.1, peak*A44Ca*5);
  fitFr->SetParLimits(6, peak*A48Ca*0.1, peak*A48Ca*5);
  // Peak positions (Fr40..Fr48) constrained the same way as mu40..mu48 in
  // the wavelength-domain fit, scaling the relative `tolerance` into MHz.
  Double_t toleranceFr = tolerance * f40;
  fitFr->SetParLimits(7,  (f40-f40) - toleranceFr, (f40-f40) + toleranceFr);
  fitFr->SetParLimits(8,  (f42-f40) - toleranceFr, (f42-f40) + toleranceFr);
  fitFr->SetParLimits(9,  (f44-f40) - toleranceFr, (f44-f40) + toleranceFr);
  fitFr->SetParLimits(10, (f48-f40) - toleranceFr, (f48-f40) + toleranceFr);
  fitFr->SetLineColor(kRed);
  fitFr->SetNpx(1000);

  // フィット実行
  grf->Fit(fitFr, "R");
  Double_t p40Fr = fitFr->GetParameter(3);
  Double_t p42Fr = fitFr->GetParameter(4);
  Double_t p44Fr = fitFr->GetParameter(5);
  Double_t p48Fr = fitFr->GetParameter(6);
  //  printf("%f %f %f %f\n", p40, p42, p44, p48);

  // 共通パラメータ
  sigmaFr   = fitFr->GetParameter(0);
  gammaFr   = fitFr->GetParameter(1);
  cutoffFr  = 0.; // no cutoff
  erf_sigFr = 0.; // no cutoff

  cout << "ACa40/ACa40 = " << p40Fr/p40Fr*100. << " [%]" << endl;
  cout << "ACa42/ACa40 = " << p42Fr/p40Fr*100. << " [%]" << endl;
  cout << "ACa44/ACa40 = " << p44Fr/p40Fr*100. << " [%]" << endl;
  cout << "ACa48/ACa40 = " << p48Fr/p40Fr*100. << " [%]" << endl;
  
  cout << "sigma = "     << sigmaFr   << " [MHz]" << endl;
  cout << "gamma = "     << gammaFr   << " [MHz]" << endl;
  cout << "cutoff = "    << cutoffFr  << " [MHz]  (no cutoff)" << endl;
  cout << "erf_sigma = " << erf_sigFr << " [MHz]  (no cutoff)" << endl;

  TLegend *legFr = new TLegend(0.65, 0.65, 0.88, 0.88);  // (x1, y1, x2, y2)
  legFr->AddEntry(grf, "Data", "p");
  legFr->AddEntry(fitFr, "Voigt fit", "l");

  for (int i = 0; i < num; i++) {    
    TF1 *f = new TF1(
		     Form("fr_%s", names[i]),
		     VoigtFunc,
		     fx.back(), fx.front(),
		     5
		     );
    double A  = fitFr->GetParameter(3 + i);       // A40,42,44,48
    double mu = fitFr->GetParameter(3 + num + i); // Fr40,42,44,48    
    f->SetNpx(1000);
    Double_t sigmaFr_corr = sigmaFr;
    Double_t gammaFr_corr = gammaFr;
    if(bMassCorrection) sigmaFr_corr = sigmaFr*TMath::Sqrt(ca_mass[0]/ca_mass[i]); //Mass correction
    f->SetParameters(A, mu, sigmaFr_corr, gammaFr_corr, 0.);
    f->SetLineColor(colors[i]);
    f->SetLineStyle(2);
    f->SetLineWidth(2);    
    f->Draw("same");
    legFr->AddEntry(f, names[i], "l");
  }
  legFr->Draw();
  c2->Update();
  c2->SaveAs(Form("./%s/FitFr_%s.root", rundir.c_str(), basename.c_str()));
  c2->SaveAs(Form("./%s/FitFr_%s.pdf",  rundir.c_str(), basename.c_str()));
  c2->SaveAs(Form("./%s/FitFr_%s.png",  rundir.c_str(), basename.c_str()));

  // --- Save fit parameters to text file ---
  std::ofstream fout(Form("./%s/Par_%s.txt", rundir.c_str(), basename.c_str()));
  fout << std::fixed << std::setprecision(6);
  fout << "# Fit parameters for: " << filename << "\n";
  fout << "# -------------------------------------------------------\n";
  fout << "# [Wavelength fit - nm]\n";
  fout << "# Chi2/NDF = " << fit->GetChisquare() << " / " << fit->GetNDF()
       << " = " << (fit->GetNDF() > 0 ? fit->GetChisquare()/fit->GetNDF() : 0.) << "\n";
  fout << "# Name               Value               Error\n";
  for (int i = 0; i < fit->GetNpar(); i++) {
    fout << Form("  %-18s  %20.6e  %20.6e\n",
                 fit->GetParName(i),
                 fit->GetParameter(i),
                 fit->GetParError(i));
  }
  fout << "#\n";
  fout << "# [Frequency fit - MHz]\n";
  fout << "# Chi2/NDF = " << fitFr->GetChisquare() << " / " << fitFr->GetNDF()
       << " = " << (fitFr->GetNDF() > 0 ? fitFr->GetChisquare()/fitFr->GetNDF() : 0.) << "\n";
  fout << "# Name               Value               Error\n";
  for (int i = 0; i < fitFr->GetNpar(); i++) {
    fout << Form("  %-18s  %20.6e  %20.6e\n",
                 fitFr->GetParName(i),
                 fitFr->GetParameter(i),
                 fitFr->GetParError(i));
  }
  fout << "#\n";
  fout << "# [Isotope ratios - from frequency fit]\n";
  fout << Form("  ACa40/ACa40  %20.6f  [%%]\n", p40Fr/p40Fr*100.);
  fout << Form("  ACa42/ACa40  %20.6f  [%%]\n", p42Fr/p40Fr*100.);
  fout << Form("  ACa44/ACa40  %20.6f  [%%]\n", p44Fr/p40Fr*100.);
  fout << Form("  ACa48/ACa40  %20.6f  [%%]\n", p48Fr/p40Fr*100.);
  fout.close();
  cout << "Parameter file saved: ./" << rundir << "/Par_" << basename << ".txt" << endl;

  return;
}