#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>
#include <fstream>

Double_t A40Ca = 0.969;
Double_t A42Ca = 0.00647;
Double_t A44Ca = 0.0209;
Double_t A48Ca = 0.00187;

Double_t Mu40Ca = 422.79179;
Double_t Mu42Ca = 422.79156;
Double_t Mu44Ca = 422.79133;
Double_t Mu48Ca = 422.79089;

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
  // par[0] : amplitude
  // par[1] : mean (center)
  // par[2] : sigma (Gaussian sigma)
  // par[3] : gamma (Lorentz HWHM)    Double_t xx = x[0] - par[1];
  // par[4] : background
  Double_t xx = x[0] - par[1];
  Double_t voigt = TMath::Voigt(xx, par[2], par[3]);
  return par[0] * voigt + par[4];
}

Double_t Voigt4(Double_t *x, Double_t *par)
{
  Double_t sum = 0; // background offset
    const int num = 4;
    Double_t bg    = par[2+2*num];
    for (int i = 0; i < 4; i++) {
        Double_t sigma = par[0];
        Double_t gamma = par[1];
        Double_t A     = par[2+i];
        Double_t mu    = par[2+num+i];
        Double_t xx = x[0] - mu;
        sum += A * TMath::Voigt(xx, sigma, gamma);
    }
    return sum + bg;
}

void FitCaVoigt(const char* filename = "RUN31_frequencyscan.txt")
//void FitCaVoigt(const char* filename = "RUN51_frequencyscan_nm.txt")
{
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    Error("FitCaVoigt", "Cannot open input file: %s", filename);
    return;
  }

  std::string basename = StripExtAndDir(filename);  
  std::vector<double> x, y, ex, ey; //wavelength
  std::vector<double> fx, efx; //frequency
  double xv, yv;
  while (infile >> xv >> yv) {
    x.push_back(xv);
    y.push_back(yv);
    ex.push_back(0.0);          // x誤差なし
    ey.push_back(0.01 * yv);    // y誤差 1%
    fx.push_back( WavelengthNmToFreqMHz(xv) );
    efx.push_back(0.0); 
  }

  auto minmax = std::minmax_element(y.begin(), y.end());
  double ymin = *minmax.first;
  double ymax = *minmax.second;
  
  TGraphErrors *gr = new TGraphErrors(x.size(),&x[0], &y[0],&ex[0], &ey[0]);
  gr->SetTitle("Calcium Isotope Spectrum; Wavelength [nm] ; Count");
  gr->SetMarkerStyle(20);

  // --- キャンバス ---
  auto c1 = new TCanvas("c1", "Voigt Fit - Ca Isotopes", 1000, 700);
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
  double gamma0 = 2.5e-5;
  //  Double_t tolerance = 3.0e-5; //in wavelength
  Double_t tolerance = 6.0e-8; //in ratio
  double peak = ymax*sigma0*TMath::Sqrt(TMath::TwoPi())*TMath::Pi();
  double bg0 = peak * 1.e-3;
  
  TF1 *fit_single = new TF1("fit_single", VoigtFunc, Mu40Ca*(1-5*tolerance), Mu40Ca*(1+30*tolerance), 5);
  fit_single->SetParameter(0,peak);
  fit_single->SetParameter(1,Mu40Ca);
  fit_single->SetParameter(2,sigma0);
  fit_single->SetParameter(3,gamma0);
  fit_single->SetParameter(4,gamma0);
  fit_single->SetParLimits(0, 0,1.e10);
  fit_single->SetParLimits(1, Mu40Ca*(1-5*tolerance), Mu40Ca*(1+5*tolerance));
  fit_single->SetParLimits(2, 1.e-6,5.e-5);
  fit_single->SetParLimits(3, 1.e-6,5.e-5);
  gr->Fit(fit_single, "R");
  sigma0 = fit_single->GetParameter(2);
  gamma0 = fit_single->GetParameter(3);
  peakshift = fit_single->GetParameter(1) - Mu40Ca; 
  cout << "40Ca peak shift = " << peakshift<<endl; 
  //  return;
  
  // Fit関数の定義
  TF1 *fit = new TF1("fit", Voigt4, x.front(), x.back(), 11);
  fit->SetParameters(sigma0,
		     gamma0,
		     peak*A40Ca,
		     peak*A42Ca,
		     peak*A44Ca,
		     peak*A48Ca,
		     Mu40Ca-peakshift,
		     Mu42Ca-peakshift,
		     Mu44Ca-peakshift,
		     Mu48Ca-peakshift,
		     bg0);
  fit->SetParNames("sigma", "gamma", "A40", "A42", "A44", "A48", "mu40", "mu42", "mu44", "mu48", "background");
  fit->SetParLimits(0, sigma0*0.8,sigma0*1.2);
  fit->SetParLimits(1, gamma0*0.8,gamma0*1.2);
  fit->SetParLimits(2, 0.,1.e15);
  fit->SetParLimits(3, 0.,1.e15);
  fit->SetParLimits(4, 0.,1.e15);
  fit->SetParLimits(5, 0.,1.e15);
  fit->SetParLimits(6, (Mu40Ca+peakshift)*(1 - tolerance), (Mu40Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(7, (Mu42Ca+peakshift)*(1 - tolerance), (Mu42Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(8, (Mu44Ca+peakshift)*(1 - tolerance), (Mu44Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(9, (Mu48Ca+peakshift)*(1 - tolerance), (Mu48Ca+peakshift)*(1 + tolerance));
  fit->SetParLimits(10, 0., peak);
  fit->SetLineColor(kRed);
  fit->SetNpx(1000);

  // フィット実行
  gr->Fit(fit, "R");
  //  fit->Draw("");
  Double_t p40 = fit->GetParameter(6);
  Double_t p42 = fit->GetParameter(7);
  Double_t p44 = fit->GetParameter(8);
  Double_t p48 = fit->GetParameter(9);
  //  printf("%f %f %f %f\n", p40, p42, p44, p48);

  // 共通パラメータ
  double sigma = fit->GetParameter(0);
  double gamma = fit->GetParameter(1);
  Double_t sigmaFr =  WavelengthNmToFreqMHz(Mu40Ca) * sigma/Mu40Ca;
  Double_t gammaFr =  WavelengthNmToFreqMHz(Mu40Ca) * gamma/Mu40Ca;

  // 共通パラメータ
  cout << "simga = "<< sigma << " [nm] = "<< sigmaFr<<" [MHz]" << endl;
  cout << "gamma = "<< gamma << " [nm] = "<< gammaFr<<" [MHz]" << endl;

  TLegend *leg = new TLegend(0.15, 0.65, 0.38, 0.88);  // (x1, y1, x2, y2)
  leg->AddEntry(gr, "Data", "p");
  leg->AddEntry(fit, "Voigt", "l");

  const int num = 4;
  int colors[num] = {kBlue, kGreen+2, kMagenta, kOrange+1};
  const char* names[num] = {"40Ca", "42Ca", "44Ca", "48Ca"};
  
  for (int i = 0; i < num; i++) {    
    TF1 *f = new TF1(
		     Form("f_%s", names[i]),
		     VoigtFunc,
		     x.front(), x.back(),
		     5
		     );
    double A  = fit->GetParameter(2 + i); // A40,42,44,48
    double mu = fit->GetParameter(2 + num  + i); // mu40,42,44,48    
    f->SetNpx(1000);
    f->SetParameters(A, mu, sigma, gamma, 0.);
    f->SetLineColor(colors[i]);
    f->SetLineStyle(2);
    f->SetLineWidth(2);    
    f->Draw("same");
    leg->AddEntry(f, names[i], "l");
  }
  leg->Draw();

  if (gSystem->AccessPathName("fig")) gSystem->mkdir("fig", true);


  c1->SaveAs(Form("./fig/Fit_%s.root",basename.c_str()));
  c1->SaveAs(Form("./fig/Fit_%s.pdf",basename.c_str()));
  c1->SaveAs(Form("./fig/Fit_%s.png",basename.c_str()));

  // --- キャンバス ---
  auto c2 = new TCanvas("c2", "Voigt Fit - Ca Isotopes ", 1000, 700);
  Double_t offset = WavelengthNmToFreqMHz(p40);
  int n = fx.size();
  for (int i = 0; i < n; i++) fx[i] -= offset;
  TGraphErrors *grf = new TGraphErrors(n, &fx[0], &y[0], &efx[0], &ey[0]);
  grf->SetTitle("Calcium Isotope Spectrum; Frequency [MHz] ; Count");
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

  // FirFr関数の定義
  TF1 *fitFr = new TF1("fitFr", Voigt4, fx.front(), fx.back(), 11);
  fitFr->SetParameters(sigmaFr,
		       gammaFr,
		       peak*A40Ca,
		       peak*A42Ca,
		       peak*A44Ca,
		       peak*A48Ca,
		       f40,
		       f42,
		       f44,
		       f48,
		       bg0);
  fitFr->SetParNames("sigma", "gamma", "A40", "A42", "A44", "A48", "fr40", "fr42", "fr44", "fr48", "background");
  fitFr->SetParLimits(0, sigmaFr*0.8,sigmaFr*1.2);
  fitFr->SetParLimits(1, gammaFr*0.8,gammaFr*1.2);
  fitFr->SetParLimits(2, 0.,1.e15);
  fitFr->SetParLimits(3, 0.,1.e15);
  fitFr->SetParLimits(4, 0.,1.e15);
  fitFr->SetParLimits(5, 0.,1.e15);
  fitFr->SetParLimits(6,  -tolerance,  +tolerance);
  fitFr->SetParLimits(7, (f42-f40)*(1 - tolerance), (f42-f40)*(1 + tolerance));
  fitFr->SetParLimits(8, (f44-f40)*(1 - tolerance), (f44-f40)*(1 + tolerance));
  fitFr->SetParLimits(9, (f48-f40)*(1 - tolerance), (f48-f40)*(1 + tolerance));
  fitFr->SetParLimits(10, 0., peak);
  fitFr->SetLineColor(kRed);
  fitFr->SetNpx(1000);
  
  // フィット実行
  grf->Fit(fitFr, "R");
  //  fit->Draw("");
  Double_t p40Fr = fitFr->GetParameter(6);
  Double_t p42Fr = fitFr->GetParameter(7);
  Double_t p44Fr = fitFr->GetParameter(8);
  Double_t p48Fr = fitFr->GetParameter(9);
  //  printf("%f %f %f %f\n", p40, p42, p44, p48);

  // 共通パラメータ
  sigmaFr = fitFr->GetParameter(0);
  gammaFr = fitFr->GetParameter(1);

  cout << "simga = "<< sigma << " [nm] = "<< sigmaFr<<" [MHz]" << endl;
  cout << "gamma = "<< gamma << " [nm] = "<< gammaFr<<" [MHz]" << endl;

  TLegend *legFr = new TLegend(0.65, 0.65, 0.88, 0.88);  // (x1, y1, x2, y2)
  legFr->AddEntry(grf, "Data", "p");
  legFr->AddEntry("fitFr", "Voigt", "l");

  for (int i = 0; i < num; i++) {    
    TF1 *f = new TF1(
		     Form("fr_%s", names[i]),
		     VoigtFunc,
		     fx.front(), fx.back(),
		     5
		     );
    double A  = fitFr->GetParameter(2 + i); // A40,42,44,48
    double mu = fitFr->GetParameter(2 + num  + i); // mu40,42,44,48    
    f->SetNpx(1000);
    f->SetParameters(A, mu, sigmaFr, gammaFr, 0.);
    f->SetLineColor(colors[i]);
    f->SetLineStyle(2);
    f->SetLineWidth(2);    
    f->Draw("same");
    legFr->AddEntry(f, names[i], "l");
  }
  legFr->Draw();
  c2->Update();
  c2->SaveAs(Form("./fig/FitFr_%s.root",basename.c_str()));
  c2->SaveAs(Form("./fig/FitFr_%s.pdf",basename.c_str()));
  c2->SaveAs(Form("./fig/FitFr_%s.png",basename.c_str()));
  return;
}
