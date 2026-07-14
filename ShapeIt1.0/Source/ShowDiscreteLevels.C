/* ***********************************************************************
* Copyright (C) 2019-2025, Dennis Muecher.                               *
* All rights reserved.                                                   *
*                                                                        *
* This program is free software: you can redistribute it and/or modify   *
* it under the terms of the GNU General Public License as published by   *
* the Free Software Foundation, either version 3 of the License, or      *
* (at your option) any later version.                                    *
* You should have received a copy of the GNU General Public License      *
* along with this program. If not, see  http://www.gnu.org/licenses/.    *
*************************************************************************/

// Function to display discrete levels histogram
// To be called from WebShapeIt.cxx when SHOW_DISCRETE_LEVELS message is received

void ShowDiscreteLevels(TCanvas* canvas, ShapeSetting* sett) {
    
    // Check if discrete level file is set
    if (sett->discreteLevelFile.empty()) {
        std::cout << "No discrete level file specified in settings" << std::endl;
        return;
    }
    
    // Read the discrete levels data file (2 column: Ex [keV], number_of_states)
    std::ifstream infile(sett->discreteLevelFile);
    if (!infile.is_open()) {
        std::cout << "ERROR: Could not open discrete level file: " << sett->discreteLevelFile << std::endl;
        return;
    }
    
    std::vector<double> energies;
    std::vector<double> states;
    
    double ex, nstates;
    while (infile >> ex >> nstates) {
        energies.push_back(ex);
        states.push_back(nstates);
    }
    infile.close();
    
    if (energies.empty()) {
        std::cout << "ERROR: No data found in discrete level file" << std::endl;
        return;
    }
    
    // Determine bin width from data spacing
    double binWidth = sett->discreteBins;  // Default from settings
    if (energies.size() > 1) {
        // Use spacing between first two points as bin width
        binWidth = energies[1] - energies[0];
    }
    
    // Create histogram
    double eMin = energies.front() - binWidth/2.0;
    double eMax = energies.back() + binWidth/2.0;
    int nBins = (int)((eMax - eMin) / binWidth + 0.5);
    
    TH1D* hDiscrete = new TH1D("hDiscrete", "Discrete Level Density;E_{x} [keV];#rho [keV^{-1}]", 
                                nBins, eMin, eMax);
    
    // Fill histogram
    for (size_t i = 0; i < energies.size(); i++) {
        // Convert from number_of_states to level density (divide by bin width in keV)
        double rho = states[i] / binWidth;
        hDiscrete->Fill(energies[i], rho);
    }
    
    // Draw on canvas
    canvas->cd();
    hDiscrete->SetLineColor(kBlue);
    hDiscrete->SetLineWidth(2);
    hDiscrete->Draw("HIST");
    canvas->SetLogy();
    
    std::cout << "Displayed discrete levels from: " << sett->discreteLevelFile << std::endl;
    std::cout << "Energy range: " << eMin << " - " << eMax << " keV" << std::endl;
    std::cout << "Bin width: " << binWidth << " keV" << std::endl;
    std::cout << "Number of bins: " << nBins << std::endl;
}
