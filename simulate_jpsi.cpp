#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <iostream>
#include <chrono>

using namespace Pythia8;
using namespace std;

int main()
{
    // Initialize Pythia
    Pythia8::Pythia pythia;

    // Set up LHC proton-proton collisions at 7 TeV
    pythia.readString("Beams:idA = 2212");  // Proton
    pythia.readString("Beams:idB = 2212");  // Proton
    pythia.readString("Beams:eCM = 7000."); // Collision energy in GeV

    // Enable J/psi production
    pythia.readString("Charmonium:all = on");

    // Ensure J/Psi is decaying only into dimuons
    pythia.readString("443:onMode = off");       // Turn off all J/Psi decays
    pythia.readString("443:onIfMatch = 13 -13"); // Enable only J/Psi -> mu+ mu-

    // Enable MPI, ISR, FSR
    pythia.readString("PartonLevel:MPI = on");
    pythia.readString("PartonLevel:ISR = on");
    pythia.readString("PartonLevel:FSR = on");

    // Enable hadron decays
    pythia.readString("HadronLevel:Decay = on");

    // Initialize Pythia
    pythia.init();

    // Create a ROOT file and tree
    TFile output_file("jpsi_production.root", "RECREATE");
    TTree tree("T", "J/psi Production Data");

    // Define event-level variables
    int nMPI;
    int nCh;

    // Define J/psi mother variables
    vector<float> jpsi_pT;
    vector<float> jpsi_y;
    vector<float> jpsi_phi;
    vector<int> jpsi_id;

    // Define muon daughter variables
    vector<float> mu_pT;
    vector<float> mu_eta;
    vector<float> mu_phi;
    vector<int> mu_charge;
    vector<int> mu_motherID;
    vector<int> mu_motherPDG;

    // Attach branches to the tree
    tree.Branch("nMPI", &nMPI);
    tree.Branch("nCh", &nCh);
    tree.Branch("jpsi_pT", &jpsi_pT);
    tree.Branch("jpsi_y", &jpsi_y);
    tree.Branch("jpsi_phi", &jpsi_phi);
    tree.Branch("jpsi_id", &jpsi_id);
    tree.Branch("mu_pT", &mu_pT);
    tree.Branch("mu_eta", &mu_eta);
    tree.Branch("mu_phi", &mu_phi);
    tree.Branch("mu_charge", &mu_charge);
    tree.Branch("mu_motherID", &mu_motherID);
    tree.Branch("mu_motherPDG", &mu_motherPDG);

    // Start the clock
    auto start_time = chrono::high_resolution_clock::now();

    // Event loop
    int nEvents = 15000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        if (!pythia.next())
            continue;

        // Clear vectors for new event
        jpsi_pT.clear();
        jpsi_y.clear();
        jpsi_phi.clear();
        jpsi_id.clear();
        mu_pT.clear();
        mu_eta.clear();
        mu_phi.clear();
        mu_charge.clear();
        mu_motherID.clear();
        mu_motherPDG.clear();

        // Get event-level quantities
        nMPI = pythia.info.nMPI();

        int nCharged = 0;
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            Particle &p = pythia.event[i];
            if (p.isFinal() && p.isCharged() && abs(p.eta()) < 2.5 && p.pT() > 0.15)
            {
                nCharged++;
            }
        }
        nCh = nCharged;

        // Loop over particles to find J/psi and its decay products
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            Particle &p = pythia.event[i];

            if (abs(p.id()) == 443)
            { // J/psi identification
                jpsi_pT.push_back(p.pT());
                jpsi_y.push_back(p.y());
                jpsi_phi.push_back(p.phi());
                jpsi_id.push_back(i);
            }

            // Find muon daughters of J/psi
            if (abs(p.id()) == 13 && p.isFinal())
            { // Muon or anti-muon
                mu_pT.push_back(p.pT());
                mu_eta.push_back(p.eta());
                mu_phi.push_back(p.phi());
                mu_charge.push_back(p.id() / abs(p.id())); // +1 for mu+, -1 for mu-
                mu_motherID.push_back(p.mother1());
                mu_motherPDG.push_back(pythia.event[p.mother1()].id());
            }
        }
        tree.Fill();
    }

    // Stop the clock
    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration<double>(end_time - start_time).count();

    // Save tree and close file
    output_file.Write();
    output_file.Close();

    // Print final statistics
    pythia.stat();

    cout << "\nSimulation complete. ROOT file 'jpsi_production.root' saved with J/psi data." << endl;
    cout << "Time taken to generate " << nEvents << " events: " << elapsed_time << " seconds." << endl;

    return 0;
}
