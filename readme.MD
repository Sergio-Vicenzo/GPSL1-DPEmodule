## Direct Position Estimation Plug-in Modules

This repository has three algorithms:

1. Conventional Direct Position Estimation (DPE) in `DPE_GPSL1`

In DPE, the objective function aims to find the PVT with the code delay and Doppler frequency that fits the best with that of the received signal. This way, the receiver PVT is directly the parameters to be estimated. Hence, PVT is estimated in a “single step” in DPE, and no estimation of code delays and Doppler frequencies is required.

2. Corr-DPE, a low-complexity alternative to DPE in `CorrDPE_GPSL1`

Corr-DPE aims to act as a viable low complexity alternative to DPE by using the multicorrelator values from two-step positioning (2SP) to save the computational cost of the correlation operation in DPE. Instead of computing the correlations for every candidate position, velocity, and time (PVT) candidates, the values are interpolated from the tracking multicorrelators. Since the maximum code delay tracking error from MP in 2SP tracking is essentially half of the Early-Late spacing, the multicorrelators for Corr-DPE need only to be computed between the Early and Late correlators. To further accelerate Corr-DPE computation, Accelerated Random Search (ARS) was also used to obtain the PVT based on the interpolated multicorrelator values.

3. MMT-integrated Corr-DPE, an upgraded Corr-DPE which uses Multipath Mitigation Technology (MMT) to compensate for Multipath (MP) reception in Corr-DPE in `CorrDPE_MMT_GPSL1`

Improving the performance of Corr-DPE, the multicorrelators for Corr-DPE are compensated for MP reception by using the MP parameters estimated by MMT. MMT was introduced as an efficient maximum likelihood estimator for accurate estimation of the code delays and carrier phase of the LOS and reflected signal (Weill, 2002). An MMT-DPE would incorporate the cost function of MMT into Corr-DPE / DPE. In essence, both Corr-DPE / DPE and MMT are the same maximum likelihood estimators that share similar cost functions, with the difference that in MMT, a reflected path is assumed to be received in addition to the LOS path. Replacing Corr-DPE / DPE cost function with that of MMT would virtually mean that the MP component is incorporated into Corr-DPE / DPE signal model. To save computational load, the multicorrelators are computed on the positioning side instead of tracking side like in Corr-DPE.

All the above algorithms are provided in a plug-in module form that can be easily integrated into existing two-step positioning (2SP) MATLAB SDRs. Following 2SP, ephemerides are decoded through conventional 2SP tracking to obtain satellite position, velocity, clock bias, and clock drift.

Programmed in a user-friendly language, MATLAB, this repository is aimed for better understanding and familiarity of a practical implementation of DPE, as well as introducing alternatives to DPE. All the algorithms here are built upon the SoftGNSS MATLAB GPS L1 C/A 2SP SDR from Borre et al. (2007), which is based on Scalar Tracking Loop (STL). But application of the algorithms here are not restricted to an STL-based receiver and can be integrated with other 2SP architectures, such as Vector Tracking Loop (VTL). Additionally, the GNSS constellation is not restricted to GPS L1 C/A, but is also integrable with SDRs working with other BPSK-modulated GNSS signals. 

Parts of SoftGNSS v3.0 have been adapted to allow the DPE_module to be integrated into it. 

## Related Publications

For DPE_GPSL1 and CorrDPE_GPSL1

[GNSS Positioning with Multicorrelators as an Alternative to Direct Position Estimation: How Good Can It Get?](https://www.ion.org/pnt/abstracts.cfm?paperID=16451)

For CorrDPE_MMT_GPSL1

[Open-Source Multipath Mitigation Technology-Integrated GNSS Direct Position Estimation](https://www.ion.org/publications/abstract.cfm?articleID=20441)

# Running the software
The algorithms presented in this repository are built upon SoftGNSS v3.0, and using them follows the same steps as running a regular SoftGNSS v3.0. Modifications were made to SoftGNSS v3.0 to also run on 16-bit samples (`int16` or `short` data types). 

When running with 16-bit data samples, use `int16` in `settings.dataType` (instead of `short`). On the other hand, use `schar` in `settings.dataType` (instead of `int8`) when running with 8-bit data samples.

Further information on using the software can be found in the original SoftGNSS v3.0 readme `readme - SoftGNSS.txt` and ppt `GPS_L1_CA_SDR.pdf`.

# Configurations
All the parameters can be edited from `initSettings.m`, which are listed below.

1. `settings.candPVT_spacing` = Grid spacing for the latitude-longitude-height estimation (Default = 1 meter)

2. `settings.DPE_latlong_span` = Span of latitude and longitude search space (Default = ±30 meters)

3. `settings.DPE_height_span` = Span of height search space (Default = ±50 meters)

4. `settings.DPE_clkBias_span` = Span of clock bias search space (Default = ±20 meters)

5. `settings.DPE_nonCohInt` = DPE non-coherent integration time (Default = 1) to improve the performance of DPE as the satellite correlations would be more filtered

6. `settings.DPE_CohInt` = DPE coherent integration time (Default = 20 ms) to improve the performance of DPE and 2SP against MP. It is highly recommended to use 20 ms for MMT-DPE.

7. `settings.DPE_plotCorrelogram` = Output the correlograms, plotted at the estimated DPE height

8. `settings.gt_llh` = Ground truth coordinates in geodetic coordinates to output the positioning errors from both DPE and 2SP

9. `settings.chipspacing_dpe_precalc` = Chip spacings between the pre-calculated correlations (Default = chips/sample)

10. `settings.MMT` = 1 (Activate MMT-DPE and MMT-2SP) or 0 (Perform conventional DPE and 2SP) - Note that for DPE_GPSL1 (Conventional DPE), MMT tracking is the default setting.

11. `settings.MMT_const ` = Relative amplitude constraint for MMT (Default = 0.8)

# Dependencies
All the algorithms were developed with MATLAB2022a and it is recommended to run the program with the same version or later. Running the program with other MATLAB versions has yet to be tested. No additional MATLAB toolbox is required i.e., it can run with just the basic MATLAB package and no additional toolboxes.

Dependencies of the underlying SoftGNSS v3.0 can be found in its readme file `readme - SoftGNSS.txt`

# Author

The algorithms introduced here were written by Sergio Vicenzo.

Sergio Vicenzo is currently a PhD candidate at the Department of Aeronautical and Aviation Engineering, The Hong Kong Polytechnic University. He received first-class honours in Bachelor of Engineering in Aviation Engineering from the same university in 2022. He is supervised by Assistant Professor Bing Xu, and co-supervised by Associate Professor Li-Ta Hsu. His research interests include GNSS urban navigation and positioning with direct position estimation.

Email: <seergio.vicenzo@connect.polyu.hk>

# Disclaimer
I, Sergio Vicenzo, hereby do not claim ownership nor any responsibility for the SoftGNSS v3.0 software by Prof. Kai Borre and Prof. Dennis M. Akos and the `llh2xyz.m` script by Prof. Todd Walter. All changes made to the SoftGNSS v3.0 software have been labelled appropriately and its GNU license has been kept intact in each code file, in accordance to GNU General Public License. No changes were made to the `llh2xyz.m` script and its copyright statement is kept intact.

No copyright infringement is intended from any of the MATLAB scripts provided in this GitHub repository.

# References
Borre K, Akos DM, Bertelsen N, Rinder P, Jensen SH (2007) [A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach](https://link.springer.com/book/10.1007/978-0-8176-4540-3). Birkhäuser, Boston, Massachusetts.

Weill LR (2002) Multipath Mitigation using Modernized GPS Signals: How Good Can it Get?. In: Proc. ION GPS 2002, Institute of Navigation, Portland, Oregon, USA, September 24 - 27, 493 - 505.

## Note
DPE_GPSL1 (Conventional DPE) is only configured to run on zero (0) Hertz intermediate frequency (IF) datasets. Modifications are required to run non-zero IF data.

This repository is still under development. I apologise if some information is incomplete or if you encounter any problems...

If you have any questions, drop an email at <seergio.vicenzo@connect.polyu.hk>!

Last Updated: 08 December 2025

	   








