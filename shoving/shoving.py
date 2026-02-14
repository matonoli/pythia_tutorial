import time
import math

import ROOT
import pythia8


def wrap_dphi(dphi: float, phi_min: float, phi_max: float) -> float:
    """Wrap dphi into [phi_min, phi_max)."""
    twopi = 2.0 * math.pi
    while dphi < phi_min:
        dphi += twopi
    while dphi >= phi_max:
        dphi -= twopi
    return dphi


def main() -> int:
    # Initialize Pythia
    pythia = pythia8.Pythia()

    # Generator. Process selection (match main441.cc)
    pythia.readString("Beams:idA = 2212")
    pythia.readString("Beams:idB = 2212")
    pythia.readString("Beams:eCM = 7000.")
    pythia.readString("SoftQCD:nonDiffractive = on")
    pythia.readString("Next:numberShowEvent = 0")

    # Enabling string shoving, setting model parameters.
    pythia.readString("Ropewalk:RopeHadronization = on")
    pythia.readString("Ropewalk:doShoving = on")
    pythia.readString("Ropewalk:doFlavour = off")
    pythia.readString("Ropewalk:rCutOff = 10.0")
    pythia.readString("Ropewalk:limitMom = on")
    pythia.readString("Ropewalk:pTcut = 2.0")
    pythia.readString("Ropewalk:r0 = 0.41")
    pythia.readString("Ropewalk:m0 = 0.2")
    pythia.readString("Ropewalk:gAmplitude = 10.0")
    pythia.readString("Ropewalk:gExponent = 1.0")
    pythia.readString("Ropewalk:deltat = 0.1")
    pythia.readString("Ropewalk:tShove = 1.")
    pythia.readString("Ropewalk:deltay = 0.1")
    pythia.readString("Ropewalk:tInit = 1.5")

    # Enabling setting of vertex information.
    pythia.readString("PartonVertex:setVertex = on")
    pythia.readString("PartonVertex:protonRadius = 0.7")
    pythia.readString("PartonVertex:emissionWidth = 0.1")

    # Initialize Pythia
    if not pythia.init():
        return 1

    # ROOT output + histograms
    out = ROOT.TFile("shoving.root", "RECREATE")

    # Histogram range must match dPhi wrapping range: [-pi/2, 3pi/2)
    phi_min = -0.5 * math.pi
    phi_max = 1.5 * math.pi

    nbins = 50
    h1 = ROOT.TH1D("deltaPhi1", "dPhi, 0 < Nch < 20;#Delta#varphi;Pairs", nbins, phi_min, phi_max)
    h2 = ROOT.TH1D("deltaPhi2", "dPhi, 20 < Nch < 40;#Delta#varphi;Pairs", nbins, phi_min, phi_max)
    h3 = ROOT.TH1D("deltaPhi3", "dPhi, 40 < Nch < 60;#Delta#varphi;Pairs", nbins, phi_min, phi_max)
    h4 = ROOT.TH1D("deltaPhi4", "dPhi, 60 < Nch < 120;#Delta#varphi;Pairs", nbins, phi_min, phi_max)

    # Event loop
    n_events = 40000
    t0 = time.time()

    for _ in range(n_events):
        if not pythia.next():
            continue

        event = pythia.event

        # Particle selection + multiplicity
        mult = 0
        parts = []  # store indices into the event

        for i in range(event.size()):
            p = event[i]
            if (
                p.isFinal()
                and p.isCharged()
                and abs(p.eta()) < 2.5
                and p.pT() > 0.5
            ):
                mult += 1
                if 1.0 < p.pT() < 3.0:
                    parts.append(i)

        if mult < 2:
            continue

        # Pair loop
        nparts = len(parts)
        for a in range(nparts):
            ia = parts[a]
            pa = event[ia]
            eta_a = pa.eta()
            phi_a = pa.phi()

            for b in range(nparts):
                ib = parts[b]
                if ia == ib:
                    continue
                pb = event[ib]

                d_eta = abs(eta_a - pb.eta())
                if 2.0 < d_eta < 4.0:
                    d_phi = wrap_dphi(phi_a - pb.phi(), phi_min, phi_max)

                    if mult <= 20:
                        h1.Fill(d_phi)
                    elif mult <= 40:
                        h2.Fill(d_phi)
                    elif mult <= 60:
                        h3.Fill(d_phi)
                    elif mult <= 120:
                        h4.Fill(d_phi)
                    # else: outside the labelled range, ignore

    t1 = time.time()

    # Save and report
    pythia.stat()

    out.cd()
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    out.Close()

    print(f"Wrote shoving.root (TH1D histograms).")
    print(f"Generated {n_events} events in {t1 - t0:.2f} s.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())