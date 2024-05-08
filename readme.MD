# DPE_module v1.0
An open-source GNSS Direct Position Estimation Plug-in Module for Two-Step Positioning SDRs, integrated into SoftGNSS v3.0 by Borre et al. (2007)

[Intelligent Positioning and Navigation Laboratory (IPNL)](https://www.polyu.edu.hk/aae/ipn-lab/us/index.html) / PNT Signal Processing Laboratory

The Hong Kong Polytechnic University

## Related Publication
DPE_module v1.0 is to be formally introduced at the upcoming Institute of Navigation GNSS+ 2024 conference.
<https://www.ion.org/gnss/abstracts.cfm?paperID=13473>

## Introduction
DPE_module v1.0 is a Direct Position Estimation (DPE) plug-in module that can be integrated into existing two-step positioning (2SP) MATLAB SDRs. 2SP information, namely tracking code phase, signal transmission time, receiver local time, satellite position from Least Squares, satellite clock bias, and Least Squares position solution, are used as input for the plug-in module. 

Programmed in a user-friendly language, MATLAB, DPE_module v1.0 is aimed for better understanding and familiarity of a practical implementation of DPE. In this repository, DPE_module v1.0 is integrated into the SoftGNSS MATLAB GPS L1 C/A 2SP SDR from Borre et al. (2007), which is based on Scalar Tracking Loop (STL). But integration of DPE_module v1.0 is not restricted to an STL-based receiver and can be integrated with other 2SP architectures, such as Vector Tracking Loop (VTL). Additionally, the GNSS constellation is not restricted to GPS L1 C/A, but is also integrable with SDRs working with other BPSK-modulated GNSS signals. 

DPE_module v1.0 would use the grid-based method, by which of establishing a set of candidate position, velocity, and time (PVT) and obtaining the correlations from each satellite before finally non-coherently summing them up to obtain the candidate position with the highest signal correlation, which would be the PVT estimate of DPE. Position as well as receiver clock bias estimate from 2SP are also used as initialization for the grid of candidate position. 

![DPE correlogram illustration](/main/Fig1.tif)

Parts of SoftGNSS v3.0 have been adapted to allow the DPE_module v1.0 to be integrated into it. 

## Running the software
The software presented in this repository is a SoftGNSS v3.0 that has been integrated with DPE_Module v1.0, and using it follows the same steps as running a regular SoftGNSS v3.0. Further information on using the software can be found in the original SoftGNSS v3.0 readme `readme - SoftGNSS.txt` and ppt `GPS_L1_CA_SDR.pdf`.

## DPE_module configuration
All DPE_module v1.0 parameters can be edited from `initSettings.m`, which are listed below.

1. `settings.candPVT_spacing` = Grid spacing for the latitude-longitude-height estimation (Default = 1 meter)

2. `settings.DPE_latlong_span` = Span of latitude and longitude search space (Default = ±30 meters)

3. `settings.DPE_height_span` = Span of height search space (Default = ±50 meters)

4. `settings.DPE_clkBias_span` = Span of clock bias search space (Default = ±20 meters)

5. `settings.DPE_nonCohInt` = DPE non-coherent integration time (Default = 1 ms) to improve the performance of DPE as the satellite correlations would be more filtered

6. `settings.DPE_plotCorrelogram` = Output the correlograms, plotted at the estimated DPE height

7. `settings.gt_llh` = Ground truth coordinates in geodetic coordinates to output the positioning errors from both DPE and 2SP

## Dependencies
DPE_module v1.0 was developed with MATLAB2022a and it is recommended to run the program with the same version or later. Running the program with other MATLAB versions has yet to be tested.

Dependencies of the underlying SoftGNSS v3.0 can be found in its readme file `readme - SoftGNSS.txt`

## Author

DPE_module v1.0 was written by Sergio Vicenzo.

Sergio Vicenzo is currently a PhD student at the Department of Aeronautical and Aviation Engineering, The Hong Kong Polytechnic University. He received first-class honours in Bachelor of Engineering in Aviation Engineering from the same university in 2022. He is supervised by Assistant Professor Bing Xu, and co-supervised by Associate Professor Li-Ta Hsu. His research interests include GNSS urban navigation and positioning with direct position estimation.

Email: <seergio.vicenzo@connect.polyu.hk>

## Disclaimer
I, Sergio Vicenzo, hereby do not claim ownership nor any responsibility for the SoftGNSS v3.0 software by Prof. Borre and Prof. Akos and the `llh2xyz.m` script by Prof. Todd Walter. All changes made to the SoftGNSS v3.0 software have been labelled appropriately and its GNU license kept intact in each code file, in accordance to GNU General Public License. No changes were made to the `llh2xyz.m` script and its copyright statement is kept intact.

No copyright infringement is intended from any of the MATLAB scripts provided in this GitHub repository.

## References
Borre K, Akos DM, Bertelsen N, Rinder P, Jensen SH (2007) [A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach](https://link.springer.com/book/10.1007/978-0-8176-4540-3). Birkhäuser, Boston, Massachusetts.

## Note
This repository is still under development. I apologise if some information is incomplete or if you encounter any problems...

Last Updated: 7 May 2024

	   
