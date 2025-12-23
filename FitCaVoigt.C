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

// ==== メインマクロ ====
void FitCaVoigt() {
  //  const char* filename = "RUN31_frequencyscan.txt";
  const char* filename = "RUN51_frequencyscan_nm.txt";
  std::ifstream infile(filename);
  
  if (!infile.is_open()) {
    Error("FitCaVoigt",
	  "Cannot open input file: %s", filename);
    return;
  }
    
  std::vector<double> x, y, ex, ey;
  double xv, yv;
  while (infile >> xv >> yv) {
    x.push_back(xv);
    y.push_back(yv);
    ex.push_back(0.0);          // x誤差なし
    ey.push_back(0.01 * yv);    // y誤差 1%
  }
  TGraphErrors *gr = new TGraphErrors(x.size(),&x[0], &y[0],&ex[0], &ey[0]);
  gr->SetTitle("Calcium Isotope Spectrum;Wavelength (nm);Counts");
  gr->SetMarkerStyle(20);

  // --- キャンバス ---
  auto c1 = new TCanvas("c1", "Voigt Fit - Ca Isotopes", 1000, 700);
  gr->SetTitle("Ca Isotope Fluorescence;Wavelength [nm];Counts");
  gr->SetMarkerStyle(20);
  gr->Draw("AP");
  gPad->SetLogy();
  gPad->SetGrid();
  // --- 軸ラベルサイズをさらに小さく ---
  gr->GetXaxis()->SetTitleSize(0.035);
  gr->GetXaxis()->SetLabelSize(0.025); // ← ここを半分に
  gr->GetYaxis()->SetTitleSize(0.035);
  gr->GetYaxis()->SetLabelSize(0.035);

  // パラメータ初期値
  Double_t peakshift = 0.;
  double sigma0 = 3.e-5;
  double gamma0 = 2.5e-5;
  Double_t tolerance = 3.0e-5;
  double peak = 1.e9*sigma0*TMath::Sqrt(TMath::TwoPi())*TMath::Pi();
  double bg0 = peak * 1.e-3;
  
  TF1 *fit_single = new TF1("fit_single", VoigtFunc, Mu40Ca-5*tolerance, Mu40Ca+30*tolerance, 5);
  fit_single->SetParameter(0,peak);
  fit_single->SetParameter(1,Mu40Ca);
  fit_single->SetParameter(2,sigma0);
  fit_single->SetParameter(3,gamma0);
  fit_single->SetParameter(4,gamma0);
  fit_single->SetParLimits(0, 0,1.e10);
  fit_single->SetParLimits(1, Mu40Ca-5*tolerance, Mu40Ca+5*tolerance);
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
  fit->SetParameters(sigma0, gamma0, peak*A40Ca, peak*A42Ca, peak*A44Ca, peak*A48Ca,
		     Mu40Ca-peakshift, Mu42Ca-peakshift, Mu44Ca-peakshift, Mu48Ca-peakshift, bg0);
  fit->SetParNames("sigma", "gamma", "A40", "A42", "A44", "A48", "mu40", "mu42", "mu44", "mu48", "background");
  fit->SetParLimits(0, sigma0*0.8,sigma0*1.2);
  fit->SetParLimits(1, gamma0*0.8,gamma0*1.2);
  fit->SetParLimits(2, 0.,1.e10);
  fit->SetParLimits(3, 0.,1.e10);
  fit->SetParLimits(4, 0.,1.e10);
  fit->SetParLimits(5, 0.,1.e10);
  fit->SetParLimits(6, Mu40Ca+peakshift - tolerance, Mu40Ca+peakshift + tolerance);
  fit->SetParLimits(7, Mu42Ca+peakshift - tolerance, Mu42Ca+peakshift + tolerance);
  fit->SetParLimits(8, Mu44Ca+peakshift - tolerance, Mu44Ca+peakshift + tolerance);
  fit->SetParLimits(9, Mu48Ca+peakshift - tolerance, Mu48Ca+peakshift + tolerance);
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

  // 共通パラメータ
  cout << "simga = "<< sigma << " [nm] = "<< TMath::C()/(Mu40Ca*1.e-9) * sigma/Mu40Ca *1.e-6<<" [MHz]" << endl;
  cout << "gamma = "<< gamma << " [nm] = "<< TMath::C()/(Mu40Ca*1.e-9) * gamma/Mu40Ca *1.e-6<<" [MHz]" << endl;

  TLegend *leg = new TLegend(0.15, 0.65, 0.38, 0.88);  // (x1, y1, x2, y2)
  leg->AddEntry(gr, "Data", "p");
  leg->AddEntry(fit, "Voigt", "l");

  const int num = 4;
  // 色と名前（40,42,44,48Ca）
  int colors[num] = {kBlue, kGreen+2, kMagenta, kOrange+1};
  const char* names[num] = {"40Ca", "42Ca", "44Ca", "48Ca"};
  std::vector<TF1*> f_iso;
  
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
    f_iso.push_back(f);
    leg->AddEntry(f, names[i], "l");
  }
  
  leg->Draw();


  
  c1->SaveAs("FitCaVoigt.root");
  c1->SaveAs("FitCaVoigt.pdf");
  c1->SaveAs("FitCaVoigt.png");

}
