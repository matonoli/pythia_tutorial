#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace Pythia8;

int main() {

  // -----------------------
  // Pythia configuration
  // -----------------------
  Pythia pythia;

  // Setup the beams.
  pythia.readString("Beams:idA = 1000080160");
  pythia.readString("Beams:idB = 1000080160"); // Oxygen ions.
  pythia.readString("Beams:eCM = 200.0");  
  pythia.readString("Beams:frameType = 1");

  // Initialize the Angantyr model to fit the total and semi-inclusive cross sections.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typically suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = 2.15,17.24,0.33");
  // A simple genetic algorithm is run for 20 generations to fit the parameters.
  pythia.readString("HeavyIon:SigFitNGen = 20");

  // Centrality bin lower limits in summed forward transverse energy (assumed).
  // At 200 GeV, activity is ~2-3x lower than at 2.76 TeV, so scale down accordingly.
  // Original O-O @ 2.76 TeV eyeball: {230, 185, 122, 79, 52, 31, 17, 9, 4}
  // Scaling by ~1/2.5 for lower energy:
  double genlim[] = {
      28.0, 22.0, 15.0, 10.0, 7.0,
      5.0, 3.5, 2.5, 1.5};

  // The upper edge of the corresponding percentiles (for labeling/reporting).
  double pclim[]  = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  if (!pythia.init()) return 1;

  // -----------------------
  // ROOT output
  // -----------------------
  TFile outFile("angantyr.root", "RECREATE");

  // Book histograms.
  // Centrality observable (forward sum ET).
  // At 200 GeV, reduce upper range proportionally.
  TH1D hSumEt("SumETfwd", "SumETfwd;#Sigma E_{T} (3.2<|#eta|<4.9) [GeV];Events (weighted)",
              200, 0.0, 80.0);  // Changed from 500 to 200

  // Number of wounded nucleons (same for O-O regardless of energy).
  TH1D hWounded("Nwounded", "Nwounded;N_{wounded};Events (weighted)",
                33, -0.5, 32.5);

  // *** NEW: Centrality-integrated mid-rapidity multiplicity histogram ***
  TH1D hMultAll("MultAll", "Centrality-integrated;N_{ch}(|#eta|<0.5);Events (weighted)",
                100, -0.5, 99.5);  // Reduced upper range for 200 GeV
  hMultAll.Sumw2();

  // Centrality-binned histograms: eta distributions + multiplicity distributions.
  const int nCent = 9;
  std::vector<TH1D*> hEtaDist(nCent, nullptr);
  std::vector<TH1D*> hMultL (nCent, nullptr);
  std::vector<TH1D*> hMultH (nCent, nullptr);

  for (int i = 0; i < nCent; ++i) {
    const std::string idx = std::to_string(i);
    hEtaDist[i] = new TH1D(("EtadistC" + idx).c_str(),
                           ("EtadistC" + idx + ";#eta;dN_{ch}/d#eta (weighted)").c_str(),
                           54, -2.7, 2.7);
    hMultH[i]   = new TH1D(("MultCH" + idx).c_str(),
                           ("MultCH" + idx + ";N_{ch}(|#eta|<0.5);dN/dN_{ch} (weighted)").c_str(),
                           75, -0.5, 74.5);  // Reduced from 299.5 for lower energy
    hMultL[i]   = new TH1D(("MultCL" + idx).c_str(),
                           ("MultCL" + idx + ";N_{ch}(|#eta|<0.5);dN/dN_{ch} (weighted)").c_str(),
                           75, -0.5, 74.5);  // Reduced from 299.5
    // Store sum of squares for proper ROOT errors if scaled later.
    hEtaDist[i]->Sumw2();
    hMultH[i]->Sumw2();
    hMultL[i]->Sumw2();
  }
  hSumEt.Sumw2();
  hWounded.Sumw2();

  // Profiles stored as binned ROOT histograms (bin content = mean, bin error = error on mean).
  TH1D hNch0("Nch0", "Average central N_{ch};Centrality upper edge (%);<#it{N}_{ch}>(|#eta|<0.5)",
             nCent, 0.0, double(nCent));
  TH1D hNw0("Nwc", "Average wounded nucleons;Centrality upper edge (%);<#it{N}_{wounded}>",
            nCent, 0.0, double(nCent));
  hNch0.Sumw2();
  hNw0.Sumw2();
  for (int i = 0; i < nCent; ++i) {
    const int pct = int(pclim[i] * 100.0 + 0.5);
    hNch0.GetXaxis()->SetBinLabel(i + 1, (std::to_string(pct) + "%").c_str());
    hNw0.GetXaxis()->SetBinLabel(i + 1, (std::to_string(pct) + "%").c_str());
  }

  // -----------------------
  // Accumulators (match original logic)
  // -----------------------
  std::vector<double> gensumw(nCent, 0.0), gensumn(nCent, 0.0);
  std::vector<double> cmult(nCent, 0.0), cmult2(nCent, 0.0);
  std::vector<double> wound(nCent, 0.0), wound2(nCent, 0.0);

  // Track all (ET, weight) pairs to recompute actual class limits.
  std::multimap<double, double> gencent;

  // Total sum of event weights
  double sumw = 0.0;

  auto start = std::chrono::high_resolution_clock::now();

  // -----------------------
  // Event loop
  // -----------------------
  const int nEvents = 20000;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // Centrality measure: sum forward transverse energy in 3.2<|eta|<4.9
    // plus a trigger requiring at least one charged forward AND backward.
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    int nc = 0;

    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (!p.isFinal()) continue;

      const double eta = p.eta();

      if (p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84) trigfwd = true;
      if (p.isCharged() && p.pT() > 0.1 && eta >  2.09 && eta <  3.84) trigbwd = true;

      if (p.pT() > 0.1 && std::abs(eta) > 3.2 && std::abs(eta) < 4.9) etfwd += p.eT();
      if (p.isCharged() && p.pT() > 0.1 && std::abs(eta) < 0.5) ++nc;
    }

    // Skip if not triggered
    if (!(trigfwd && trigbwd)) continue;

    // Weight
    const double weight = pythia.info.weight();
    sumw += weight;

    // *** Fill centrality-integrated multiplicity histogram ***
    hMultAll.Fill(nc, weight);

    // Fill centrality observable + keep full distribution for limit finding
    hSumEt.Fill(etfwd, weight);
    gencent.insert(std::make_pair(etfwd, weight));

    // Wounded nucleons (absorptive + diffractive, target + projectile)
    const int nw = pythia.info.hiInfo->nAbsTarg() + pythia.info.hiInfo->nDiffTarg()
                 + pythia.info.hiInfo->nAbsProj() + pythia.info.hiInfo->nDiffProj();
    hWounded.Fill(nw, weight);

    // Determine centrality class index from assumed limits (higher ET => more central).
    int genidx = -1;
    for (int ic = 0; ic < nCent; ++ic) {
      if (etfwd >= genlim[ic]) { genidx = ic; break; }
    }
    if (genidx < 0) continue;

    // Fill class histograms and sums
    gensumw[genidx] += weight;
    gensumn[genidx] += 1.0;

    hMultH[genidx]->Fill(nc, weight);
    hMultL[genidx]->Fill(nc, weight);

    cmult [genidx] += nc * weight;
    cmult2[genidx] += double(nc) * double(nc) * weight;

    wound [genidx] += nw * weight;
    wound2[genidx] += double(nw) * double(nw) * weight;

    // Fill eta distribution for this class
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (p.isFinal() && p.isCharged() && std::abs(p.eta()) < 2.7 && p.pT() > 0.1) {
        hEtaDist[genidx]->Fill(p.eta(), weight);
      }
    }
  }

  auto stop = std::chrono::high_resolution_clock::now();
  const double elapsed = std::chrono::duration<double>(stop - start).count();

  // -----------------------
  // Normalization (mimic original factors)
  // -----------------------
  if (sumw > 0.0) {
    // Original code: sumet /= sumw*2.0; wounded /= sumw*2.0;
    hSumEt.Scale(1.0 / (sumw * 2.0));
    hWounded.Scale(1.0 / (sumw * 2.0));
  }

  // Original code divides by (gensumw * binWidth) for these:
  //   hmultH: 40.0, hmultL: 4.0, etadist: 0.1
  for (int ic = 0; ic < nCent; ++ic) {
    if (gensumw[ic] <= 0.0) continue;

    hMultH[ic]->Scale(1.0 / (gensumw[ic] * 40.0));
    hMultL[ic]->Scale(1.0 / (gensumw[ic] *  4.0));
    hEtaDist[ic]->Scale(1.0 / (gensumw[ic] * 0.1));

    // Fill mean profiles (bin content = mean, bin error = error on mean)
    const double Nch = cmult[ic] / gensumw[ic];
    double varNch = (cmult2[ic] / gensumw[ic]) - (Nch * Nch);
    if (gensumn[ic] > 0.0) varNch /= gensumn[ic];

    const double Nw = wound[ic] / gensumw[ic];
    double varNw = (wound2[ic] / gensumw[ic]) - (Nw * Nw);
    if (gensumn[ic] > 0.0) varNw /= gensumn[ic];

    hNch0.SetBinContent(ic + 1, Nch);
    hNch0.SetBinError  (ic + 1, (varNch > 0.0) ? std::sqrt(varNch) : 0.0);

    hNw0.SetBinContent(ic + 1, Nw);
    hNw0.SetBinError  (ic + 1, (varNw > 0.0) ? std::sqrt(varNw) : 0.0);
  }

  // -----------------------
  // Recompute "actual" ET limits between centrality classes (as in original)
  // -----------------------
  std::vector<double> newlim(nCent, 0.0);
  double curr = 0.0, prev = 0.0, acc = 0.0;
  int idxa = nCent - 1;
  double lim = sumw * (1.0 - pclim[idxa]);

  for (auto it = gencent.begin(); it != gencent.end(); ++it) {
    prev = curr;
    curr = it->first;
    const double w = it->second;

    if (acc < lim && acc + w >= lim) {
      newlim[idxa--] = prev + (curr - prev) * (lim - acc) / w;
      if (idxa < 0) break;
      lim = sumw * (1.0 - pclim[idxa]);
    }
    acc += w;
  }

  // -----------------------
  // Write and report
  // -----------------------
  outFile.cd();
  hSumEt.Write();
  hWounded.Write();
  hMultAll.Write();  
  hNch0.Write();
  hNw0.Write();
  for (int ic = 0; ic < nCent; ++ic) {
    if (hMultH[ic])   hMultH[ic]->Write();
    if (hMultL[ic])   hMultL[ic]->Write();
    if (hEtaDist[ic]) hEtaDist[ic]->Write();
  }
  outFile.Close();

  pythia.stat();

  std::cout << "Wrote ROOT file: angantyr.root\n";
  std::cout << "Time for " << nEvents << " attempted events: " << elapsed << " s\n\n";

  std::cout << "Generated limits between centrality classes in this run:\n"
            << "   %   assumed    actual\n";
  for (int idx = 0; idx < nCent; ++idx) {
    std::cout << std::setw(4) << int(pclim[idx] * 100.0 + 0.5)
              << std::setw(10) << std::fixed << std::setprecision(1) << genlim[idx]
              << std::setw(10) << std::fixed << std::setprecision(1) << newlim[idx]
              << "\n";
  }

  return 0;
}