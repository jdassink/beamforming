# beamforming

Seismo-acoustic array processing routines written in Fortran 90. All programs can be compiled with GNU Fortran compiler by invoking the typical 'make' command. Individual programs can be compiled by invoking 'make <program name>', e.g. 'make timefisher'. The program uses the Fastest Fourier Transform in the West (FFTW) which can be downloaded on http://www.fftw.org. Make sure to change the INCLUDE and LIB_INC environmental variables to link to the FFTW3 headers (fftw3.f) and library files.
  
A short explanation of the routines is included below. More information on the usage can be found by executing the program. 
A detailed decription of the algorithms can be found in Evers, 2008. Recent examples of the use of the codes can be found in Assink et al., 2016 and Evers et al., 2018. An effort is underway to port these algorithms to Python.

---

**timefisher**
Time-Domain Fisher Detector (Melton and Bailey, 1957) - command line tool. Estimates the correlation as a function of frequency over a grid of horizontal slowness values.

**freqfisher**
Frequency-Domain Fisher Detector (Smart and Flinn, 1971) - command line tool. Estimates the coherency as a function of frequency over a grid of horizontal slowness values.

**tdoa**
Time-delay of arrival (e.g., Szuberla and Olson, 2004), estimates back azimuth and apparent velocity from cross-correlation functions of sensor pairs 

**ccts**
Cross-correlation Trace Stacking algorithm (e.g., Gibbons, 2015), beamforms pair-wise cross-correlation functions over a horizontal grid of slowness values

--

**bibliography**

Assink, J. D., G. Averbuch, P. S. M. Smets, and L. G. Evers (2016), On the infrasound detected from the 2013 and 2016 DPRK's underground nuclear tests, Geophys. Res. Lett., 43, 3526–3533, https://doi.org/10.1002/2016GL068497 

Evers, L.G. (2008), "The inaudible symphony: On the detection and source identification of atmospheric infrasound", Ph.D Thesis, Delft University of Technology, https://repository.tudelft.nl/islandora/object/uuid%3A4de38d6f-8f68-4706-bf34-4003d3dff0ce

Evers, L.G., J D Assink, P SM Smets (2018); Infrasound from the 2009 and 2017 DPRK rocket launches, Geophysical Journal International, Volume 213, Issue 3, 1 June 2018, Pages 1785–1791, https://doi.org/10.1093/gji/ggy092

Gibbons, S.J., Vladimir Asming, Lars Eliasson, Andrei Fedorov, Jan Fyen, Johan Kero, Elena Kozlovskaya, Tormod Kværna, Ludwik Liszka, Sven Peter Näsholm, Tero Raita, Michael Roth, Timo Tiira, Yuri Vinogradov; The European Arctic: A Laboratory for Seismoacoustic Studies. Seismological Research Letters ; 86 (3): 917–928. doi: https://doi.org/10.1785/0220140230

Melton, B.S. and L. F. Bailey (1957). ”MULTIPLE SIGNAL CORRELATORS.” GEOPHYSICS, 22(3), 565-588.
https://doi.org/10.1190/1.1438390

Smart, E. and Flinn, E. A. (1971), Fast Frequency‐Wavenumber Analysis and Fisher Signal Detection in Real‐Time Infrasonic Array Data Processing. Geophysical Journal of the Royal Astronomical Society, 26: 279-284. http://dx.doi.org/10.1111/j.1365-246X.1971.tb03401.x

Szuberla, C. and J. Olson (2004), "Uncertainties associated with parameter estimation in atmospheric infrasound arrays", The Journal of the Acoustical Society of America 115:1, 253-258, https://doi.org/10.1121/1.1635407
