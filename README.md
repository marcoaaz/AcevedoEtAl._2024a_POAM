# AcevedoEtAl._2024a_POAM

This repository contains an open-source image analysis pipeline to convert a traditional petrographic (polarising) microscope into a micro-fabric analyser using internal and external software modules. The software is widely documented in Supplementary Material 1 of the submitted manuscript:

Acevedo Zamora, M. A.; Schrank, C. E.; Kamber, B. S. Using the traditional microscope for mineral grain orientation determination: a prototype image analysis pipeline for optic axis mapping (POAM). Journal of Microscopy in revision.

The workflow is organised as follows (the main MatLab script 'stack_spectra_leica_v18.m' corresponds to blue panels):

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024a_POAM/assets/61703106/74f43cd2-57df-469f-ac9d-58fe6823a42f" width=60% height=60%>

The POAM dependencies are (updated up to *_v14.m):

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024a_POAM/assets/61703106/5f52711b-aeb7-4e63-82cf-cc703caad5ff" width=60% height=60%>

The scripts and functions documentation are listed below (preliminary descriptions, updated 13-Jan-2024):

1. EBSD validation
  +	readEBSD_h5oiana_v5.m = Uses MTEX to plot EBSD map of harzburgite showing enstatite grains with spectral transmission
  +	enstatite_model.m = Uses MTEX to plot the default enstatite 3D crystal model with crystallographic and optic axes
  +	validateBiaxialAzimuth.m = compare biaxial mineral slow-axis azimuth between POAM and MTEX Birefringence package (Sorensen review)
  +	validateUniaxialAzimuth.m = compare uniaxial mineral optic-axis azimuth between POAM and MTEX Birefringence package (first paper 2 submission)
2. fix_algorithm1_v2.m = Following Axer et al. 2011 paper for obtaining coefficients
3. imageFourierS_optim_ver2.m = matricial implementation of algorithm 1
4. sinDescriptor_ver3.m = function with boosted performance using algorithm 1 Function handles as input to algorithm 2
5. ROI Tool (multipol)
  + ROI_modulation_data.m = function with updated to work with pixelFourierS_ver2.m
  + pixelFourierS_ver2.m = update of Algorithm 1 for faster computation using Fourier series coefficients (without fitting a Fourier2)
6. highestOrderPeak_ver2.m = 10x improved performance using function handles
7. stack_spectra_leica_v17.m = update incorporating all other updates and supporting algorithms 3 and 4 as a switch/case.
  - stack_spectra_leica_v15_algorithm3.m = update to allow pixelFourierS_ver2.m vectorised
  - stack_spectra_leica_v16_algorithm3.m = TRIAL only. update to avoid running algorithm 1
  - Improvements:
    -	sinDescriptorPlot_ver2.m = function with updated input list
    -	stackImportLoop.m = function to use imread() import loop of selected variable range from the image stack. The sRGB images are linearised. Two output options are available: greyscale stack (3D) or RGB stack (4D). The output is rescaled following a scalar value (downscaling speeds up POAM). 
    -	imgCropAndOverlay.m = function to obtain POAM maps common bounding box and provide an aesthetic image overlay to be plotted below the orientation map and objects.
    -	rearrangePeakImages.m = uses 'img_closest*.tif to rearrange the peak images (algorithm 4)


In future updates, we expect to increase the efficiency of Algorithm 2 and provide explanatory videos for new users who would like to learn to acquire their images (microscope stage rotation and capturing) and manage the files for effectively use POAM.

To be able to run the scripts, you require installing:

The used MatLab version was:
  +	MATLAB Version: 9.14.0.2239454 (R2023a) Update 1
  +	MATLAB License Number: 729365
  +	Operating System: Microsoft Windows 10 Enterprise Version 10.0 (Build 19045)
  +	Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
, with toolboxes:
  +	Computer Vision Toolbox, Version 10.4
  +	Curve Fitting Toolbox, Version 3.9
  +	Deep Learning Toolbox, Version 14.6
  +	Fixed-Point Designer, Version 7.6
  +	Image Processing Toolbox, Version 11.7
  +	MATLAB Compiler, Version 8.6
  +	Mapping Toolbox, Version 5.5
  +	Parallel Computing Toolbox, Version 7.8
  +	Signal Processing Toolbox, Version 9.2
  +	Statistics and Machine Learning Toolbox, Version 12.5
  + Symbolic Math Toolbox, Version 9.3
  + Wavelet Toolbox, Version 6.3


Thanks.
:) 
