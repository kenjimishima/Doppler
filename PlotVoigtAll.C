
// PlotVoigtAll.C
// ROOT macro to compare:
// 1) Pure Voigt
// 2) Voigt with power-law cutoff
// 3) Voigt with erf cutoff
//
// Usage:
//   root -l PlotVoigtAll.C
//   root -l -q PlotVoigtAll.C
//
// ------------------------------------------------------------

#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TAxis.h"
#include "TStyle.h"
#include <cmath>

Double_t bg = 2.e-6;

// ============================================================
// Core functions
// ============================================================

// --- Pure Voigt wrapper ---
static inline double PureVoigt(double dnu,
                               double sigma,
                               double gamma)
{
  return TMath::Voigt(dnu, sigma, gamma);
}

// --- Voigt with power-law cutoff ---
static inline double VoigtCutoffPL(double dnu,
                                   double sigma,
                                   double gamma,
                                   double nuc,
                                   double q)
{
  const double V = TMath::Voigt(dnu, sigma, gamma);
  if (nuc <= 0.0) return V;

  const double a = std::abs(dnu) / nuc;
  const double C = 1.0 / (1.0 + std::pow(a, q));
  return V * C;
}

// --- Voigt with error-function cutoff ---
static inline double VoigtCutoffErf(double dnu,
                                    double sigma,
                                    double gamma,
                                    double nuc,
                                    double sigc)
{
  const double V = TMath::Voigt(dnu, sigma, gamma);
  if (nuc <= 0.0) return V;

  const double arg =
    (std::abs(dnu) - nuc) / (TMath::Sqrt2() * sigc);

  const double C = 0.5 * TMath::Erfc(arg);
  return V * C;
}

// ============================================================
// TF1 wrappers
// ============================================================

// par[0]=A, [1]=mu, [2]=sigma, [3]=gamma
double PureVoigt_TF1(double *x, double *par)
{
  const double dnu = x[0] - par[1];
  return par[0] * PureVoigt(dnu, par[2], par[3]) + bg;
}

// par[0]=A, [1]=mu, [2]=sigma, [3]=gamma, [4]=nuc, [5]=q
double VoigtCutoffPL_TF1(double *x, double *par)
{
  const double dnu = x[0] - par[1];
  return par[0] *
         VoigtCutoffPL(dnu, par[2], par[3], par[4], par[5]) + bg;
}

// par[0]=A, [1]=mu, [2]=sigma, [3]=gamma, [4]=nuc, [5]=sigc
double VoigtCutoffErf_TF1(double *x, double *par)
{
  const double dnu = x[0] - par[1];
  return par[0] *
         VoigtCutoffErf(dnu, par[2], par[3], par[4], par[5]) + bg;
}

// ============================================================
// Main plot macro
// ============================================================
void PlotVoigtAll()
{
  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c", "Voigt cutoff comparison", 900, 700);

  const double xmin = -2000.0;
  const double xmax =  2000.0;

  // ---- Common Voigt parameters ----
  const double A     = 1.0;
  const double mu    = 0.0;
  const double sigma = 30.0;
  const double gamma = 80.0;

  // ---- Cutoff parameters ----
  const double nuc  = 600.0;
  const double q    = 6.0;
  const double sigc = 140.0;

  // ---- Pure Voigt ----
  TF1 *fV = new TF1("fV", PureVoigt_TF1, xmin, xmax, 4);
  fV->SetParameters(A, mu, sigma, gamma);
  fV->SetLineColor(kBlack);
  fV->SetLineStyle(2);
  fV->SetLineWidth(2);

  // ---- Power-law cutoff ----
  TF1 *fPL = new TF1("fPL", VoigtCutoffPL_TF1, xmin, xmax, 6);
  fPL->SetParameters(A, mu, sigma, gamma, nuc, q);
  fPL->SetLineColor(kRed+1);
  fPL->SetLineWidth(3);

  // ---- Erf cutoff ----
  TF1 *fErf = new TF1("fErf", VoigtCutoffErf_TF1, xmin, xmax, 6);
  fErf->SetParameters(A, mu, sigma, gamma, nuc, sigc);
  fErf->SetLineColor(kBlue+1);
  fErf->SetLineWidth(3);

  fV->SetNpx(3000);
  fPL->SetNpx(3000);
  fErf->SetNpx(3000);

  // ---- Draw ----
  fPL->Draw();
  fPL->GetYaxis()->SetRangeUser(1.e-6,1.e-2);
  fV->Draw("same");
  fErf->Draw("same");

  gPad->SetLogy();
  gPad->SetGrid();

  // ---- Legend ----
  TLegend *leg = new TLegend(0.60, 0.70, 0.88, 0.88);
  leg->AddEntry(fV,   "Pure Voigt",        "l");
  leg->AddEntry(fPL,  "Voigt + power-law cutoff", "l");
  leg->AddEntry(fErf, "Voigt + erf cutoff","l");
  leg->Draw();

  c->Update();
}
