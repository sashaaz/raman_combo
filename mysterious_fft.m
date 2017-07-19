clear all
c = 300;
amplitude_band=[1.5,1,1,1,.1];


% Center of each mode in nm, mode_wavelength=[ lambda1, lambda2, ... ]
mode_wavelength = [915, 840, 780, 730, 680];

% Number of frequency points
freq_num = 400;

% Bandwidth (FWHM) of each mode in nm, mode_bandwidth=[band1, band2, ... ]
mode_bandwidth = [30, 25, 20, 10, 10];

% beta = 2*GDD of each mode in fs^2
beta = [-10000, 0, 0, 0, 0];

% Convert to wavelength to frequency
mode_frequency = 2*pi*c./mode_wavelength;

% Convert FWHM bandwidth in nm to FWHM in 2*pi/fs
mode_bandwidth_freq = 2*pi*c.*mode_bandwidth./mode_wavelength.^2;


% Find max, min of mode frequencies
freq_max = max(mode_frequency);

band_max = max(mode_bandwidth_freq);

% Redefine mode_bandwidth to equal 2*sigma^2
mode_bandwidth_freq = 2*mode_bandwidth_freq/2.35482;

% Define frequency axis
% ********************* Note ************************
% The frequency axis must include zero for the FFT algorithm to 
% recognize the correct distance of the modes from zero. Not including
% zero will result in error, as the FFT will associate smallest frequency
% component with zero.
frequency = linspace(0, freq_max+2*band_max, freq_num);

%Define phase axis
% phase_range=10*pi;
% phase_num=200;
phase_range=100*pi;
phase_num=2000;

phase=linspace(0,phase_range,phase_num);

[phi,ww]=meshgrid(phase,frequency);

E1=amplitude_band(1)*exp(-((ww-mode_frequency(1))/mode_bandwidth_freq(1)).^2).*exp(1i*beta(1)*(ww-mode_frequency(1)).^2);

E1_gaussian = abs(E1(:, 10));
E1_gaussian_fft = fftshift(fft(E1_gaussian));

E1_time_wrong=ifftshift(ifft(E1,10000,1));
E1_time =ifftshift(ifft(E1,10000,2));

