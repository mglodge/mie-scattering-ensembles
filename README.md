# Mie Scattering (Ensembles)

## Overview

This code uses Mie theory to calculate scattering phase functions for ensembles (size-distributions) of particles. It was designed to calculate the optical properties of KCl particles as part of the collaborative paper:

*"The Effects of Cuboid Particle Scattering on Reflected Light Phase Curves: Insights from Laboratory Data and Theory"* (2025) <br>
Hamill, C.D., Johnson, A.V., Lodge, M., Gao, P., Nag, R., Batalha, N., Christie, D.A. and Wakeford, H.R., The Astrophysical Journal, 987(2), p.176.

<p align="center">
  <img src="https://github.com/user-attachments/assets/59b17472-5c13-4cb9-8265-b97b74536fec" alt="drawing" width="800" />
</p>

In this paper, laser light was experimentally scattered through aerosolised solid KCl particles (representing exoplanet cloud analogues), and the reflected light was measured at each scattering angle. By measuring the size distributions of the KCl particles, theoretical calculations could be compared to the experimental results to determine the impact of particle shape on the measurements (e.g. whether the particles were spheres, regular cuboids, irregular cuboids). The experimental work was completed by C. Hamil and A. Johnson, and the Mie scattering code was developed by M. Lodge and is attached here, including a sample size distribution file. Mie theory was used for spherical particles, and the discrete dipole approximation was used to determine the optical properties of non-spherical particles. The equations used for Mie scattering calculations are in Chapter 2 of Matt Lodge's thesis (also attached), where more information about the project is also explained in Chapter 5.

The figure below (from attached thesis) gives examples of the experimentally-measured size distributions, the number of particles $n_i$ measured in each radius bin $r$. The example file attached here ("size_distribution.txt") represents the 'small' distribution.

<p align="center">
  <img src="https://github.com/user-attachments/assets/83084efe-852e-4afb-913b-fd1ceb4b97fe" alt="drawing" width="500" />
</p>

The goal of the code is to produce model phase functions that determine the unpolarised intensity of light $|S_0|^2$ that has been scattered at each angle ($\theta$), which can then be compared to the experimental data, as shown here:

<p align="center">
  <img src="https://github.com/user-attachments/assets/37c851ba-2023-4ac6-85c3-86c2597efc82" alt="drawing" width="800" />
</p>

# Compilation/Running the code

To compile the code, download the project into a directory, navigate to that directory in the terminal, and then compile using:

    gcc mie_scattering_ensembles.c -o mie_scattering_ensembles

This should produce an executable version of the file. Run it using:

    ./mie_scattering_ensembles

## Code Inputs/Outputs

INPUTS:

    - Variables that can be changed directly within "mie_scattering_ensembles.c"
    
        - refractive index of material (as a complex number)

        - wavelength of light (units: um)

    - "size_distribution.txt" file

        -  This is the size distribution of particles, using the direct output file from ExCESS (Exoplanet Cloud Ensemble Scattering System)
        -  An example file is provided here
        -  Lines 1-16: Header information
        -  Line 17: Number of bins
        -  Line 18: onwards: bin_start_diameter, bin_end_diameter, number of particles in bin

OUTPUTS:

    - "S1_and_S2_vs_angle.txt" file

        -  .csv file containing: 
            - angle (degrees)
            - average ensemble scattering amplitude S1^2
            - average ensemble scattering amplitude S2^2
            - total unpolarised average scattering magnitude
            - degree of linear polarisation

    - Tables of data for each particle size bin, including:

        - Q_sca (scattering efficiency)
        - Q_ext (extinction efficiency)
        - Q_abs (absorption efficiency)
        - S1 (cumulative scattering amplitude at an angle of 0 degrees)
        - S2 (cumulative scattering amplitude at an angle of 0 degrees)

    - Phase function graph: Total unpolarised scattering intensity of ensemble vs angle

    - Phase function graph: Polarised components S1 and S2, as well as the unpolarised average 

    - Graph of 'Degree of linear polarisation vs angle'
