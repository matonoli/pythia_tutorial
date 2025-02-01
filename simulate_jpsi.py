import pythia8
import ROOT
import time

# Initialize Pythia
pythia = pythia8.Pythia()

# Set up LHC proton-proton collisions at 7 TeV
pythia.readString("Beams:idA = 2212")  # Proton
pythia.readString("Beams:idB = 2212")  # Proton
pythia.readString("Beams:eCM = 7000.")  # Collision energy in GeV

# Enable J/psi production
pythia.readString("Charmonium:all = on")

# Ensure J/Psi is decaying only into dimuons
pythia.readString("443:onMode = off")       # Turn off all J/Psi decays
pythia.readString("443:onIfMatch = 13 -13") # Enable only J/Psi -> mu+ mu-

# Enable MPI, ISR, FSR
pythia.readString("PartonLevel:MPI = on")
pythia.readString("PartonLevel:ISR = on")
pythia.readString("PartonLevel:FSR = on")

# Enable hadron decays
pythia.readString("HadronLevel:Decay = on")

# Initialize Pythia
pythia.init()

# Create a ROOT file and tree
output_file = ROOT.TFile("jpsi_production.root", "RECREATE")
tree = ROOT.TTree("T", "J/psi Production Data")

# Define event-level variables
nMPI = ROOT.std.vector('int')(1)
nCh = ROOT.std.vector('int')(1)

# Define J/psi mother variables
jpsi_pT = ROOT.std.vector('float')()
jpsi_y = ROOT.std.vector('float')()
jpsi_phi = ROOT.std.vector('float')()
jpsi_id = ROOT.std.vector('int')()

# Define muon daughter variables
mu_pT = ROOT.std.vector('float')()
mu_eta = ROOT.std.vector('float')()
mu_phi = ROOT.std.vector('float')()
mu_charge = ROOT.std.vector('int')()
mu_motherID = ROOT.std.vector('int')()
mu_motherPDG = ROOT.std.vector('int')()

# Attach branches to the tree
tree.Branch("nMPI", nMPI)
tree.Branch("nCh", nCh)
tree.Branch("jpsi_pT", jpsi_pT)
tree.Branch("jpsi_y", jpsi_y)
tree.Branch("jpsi_phi", jpsi_phi)
tree.Branch("jpsi_id", jpsi_id)
tree.Branch("mu_pT", mu_pT)
tree.Branch("mu_eta", mu_eta)
tree.Branch("mu_phi", mu_phi)
tree.Branch("mu_charge", mu_charge)
tree.Branch("mu_motherID", mu_motherID)
tree.Branch("mu_motherPDG", mu_motherPDG)

# Start the clock
start_time = time.time()

# Event loop
nEvents = 4000
for iEvent in range(nEvents):
    if not pythia.next():
        continue
    
    # Clear vectors for new event
    nMPI[0] = pythia.infoPython().nMPI()
    
    nCharged = 0
    for i in range(pythia.event.size()):
        p = pythia.event[i]
        if p.isFinal() and p.isCharged() and abs(p.eta()) < 2.5 and p.pT() > 0.15:
            nCharged += 1
    nCh[0] = nCharged
    
    jpsi_pT.clear()
    jpsi_y.clear()
    jpsi_phi.clear()
    jpsi_id.clear()
    mu_pT.clear()
    mu_eta.clear()
    mu_phi.clear()
    mu_charge.clear()
    mu_motherID.clear()
    mu_motherPDG.clear()
    
    # Loop over particles to find J/psi and its decay products
    for i in range(pythia.event.size()):
        p = pythia.event[i]
        
        if abs(p.id()) == 443:  # J/psi identification
            jpsi_pT.push_back(p.pT())
            jpsi_y.push_back(p.y())
            jpsi_phi.push_back(p.phi())
            jpsi_id.push_back(i)
        
        # Find muon daughters of J/psi
        if abs(p.id()) == 13 and p.isFinal():  # Muon or anti-muon
            mu_pT.push_back(p.pT())
            mu_eta.push_back(p.eta())
            mu_phi.push_back(p.phi())
            mu_charge.push_back(int(p.id() / abs(p.id())))  # +1 for mu+, -1 for mu-
            mu_motherID.push_back(p.mother1())
            mu_motherPDG.push_back(pythia.event[p.mother1()].id())
    
    tree.Fill()

# Stop the clock
end_time = time.time()
elapsed_time = end_time - start_time

# Save tree and close file
output_file.Write()
output_file.Close()

# Print final statistics
pythia.stat()

print(f"\nSimulation complete. ROOT file 'jpsi_production.root' saved with J/psi data.")
print(f"Time taken to generate {nEvents} events: {elapsed_time:.2f} seconds.")
