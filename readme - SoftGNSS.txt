An Open Source Archive of GNSS Software Defined Radio receivers based on MATLAB
===============================================================================



Authors
-------------------------------------------------------------------------------
Dennis Akos  
E-Mail: <dma@colorado.edu>
HP: <http://www.colorado.edu/aerospace/dennis-akos>

Nagaraj Channarayapatna Shivaramaiah
E-Mail: <nagarajcs@colorado.edu>

Yafeng Li
E-Mail: <lyf8118@126.com>

Jakob Almqvist
E-Mail: <>

Daehee Won 
E-Mail: <daehee.won@colorado.edu>

Hinckley
E-Mail: <>



Features
-------------------------------------------------------------------------------
* GNSS signal processing functions written in MATLAB
    * Code generations
    * Signal acquisition / tracking (data + pilot)
    * Decoding navigation messages 
    * Pseudo-range mesurement generation
    * Calculation of position
* Support following signals
    * GPS L1 C/A
    * GPS L2C (data + pilot)
    * GPS L5C (data + pilot)
    * Galileo E1 (data + pilot)
    * Galileo E5a (data + pilot)
    * Galileo E5b (data + pilot)
    * GLONASS L1 C/A
    * GLONASS L2 C/A
    * Beidou Phase II B1
    * Beidou Phase II B2
* Support RF binary file for post processing
    * All the SDRs have been tested with IF signals collected by the NUT4NT 
	sampler of the Ou Amungo company 

	

Directory and Files
-------------------------------------------------------------------------------
./Doc                     Summary PowerPoint documents for each SDR receiver
./Common                  Common functions between differnt SDR receivers
./IF_Data_Set             Folder containing the IF data sets to be processed and 
                          corresponding Metadata files
./GPS_L1_CA               GPS L1 C/A SDR receiver
    ./include             Functions related to NAV data decoding and PVT computation  
    ./init.m              Strating function of the receiver  
    ./initSettings.m      Parameter configurations of the receiver    
    ./postProcessing.m    Top level processing function of the receiver  
    ./acquisition.m       Signal acquisition function  
    ./tracking.m          Signal tracking function 
    ./postNavigation.m    PVT computation

Software Dependencies
-------------------------------------------------------------------------------
* All the SDRs have been tested on MATLAB R2016b and R2013b. Some functions 
  of the MATLAB Communications System Toolbox are needed:
  -- Create Galois field array: gf()
  -- BCH decoder: bchdec()
  -- Detect errors in input data using CRC: comm.CRCDetector()
  -- Convert convolutional code polynomials to trellis description: poly2trellis()
  -- Convolutionally decode binary data using Viterbi algorithm: vitdec()
  -- Detect errors in input data using CRC: step() 
   
 
How to use
-------------------------------------------------------------------------------
* Step 1: Copy the IF data file into folder "IF_Data_Set";
* Step 2: Configue parameters related to IF data file in function "initSettings.m";
* Step 2: Start processing by runing function "init.m".


Implementation details
-------------------------------------------------------------------------------
See the summary PowerPoint documents for each SDRs:
    * GPS L1 C/A
	   -- GPS_L1_CA_SDR.pptx
	   
Test signal and parameter configurations
-------------------------------------------------------------------------------
    * GPS_L1_CA_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 14.58e6 Hz
	   -- sampling Frequency: 53e6 Hz
    * GPS_L2C_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 7.4e6 Hz
	   -- sampling Frequency: 53e6 Hz	
    * GPS_L5C_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 13.55e6 Hz
	   -- sampling Frequency: 99.375e6 Hz
    * Galileo_E1_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 14.58e6 Hz
	   -- sampling Frequency: 53e6 Hz
    * Galileo_E5a_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 13.55e6 Hz
	   -- sampling Frequency: 99.375e6 Hz
    * Galileo_E5b_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 17.14e6 Hz
	   -- sampling Frequency: 99.375e6 Hz
    * Glonass_L1_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- Frequency spacing: 562.5e3 Hz 
	   -- IF: 11e6 Hz
	   -- sampling Frequency: 53e6 Hz
    * Glonass_L2_IF_signal.bin
       -- dataType: schar
	   -- Frequency spacing: 437.5e3 Hz
	   -- fileType: 8 bit real samples
	   -- IF: 12e6 Hz
	   -- sampling Frequency: 53e6 Hz
    * Beidou_B1_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 28.902e6 Hz
	   -- sampling Frequency: 99.375e6 Hz
    * Beidou_B2_IF_signal.bin
       -- dataType: schar
	   -- fileType: 8 bit real samples
	   -- IF: 17.14e6 Hz
	   -- sampling Frequency: 99.375e6 Hz

	   