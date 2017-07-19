% Center of each mode in nm, mode_wavelength=[ lambda1, lambda2, ... ]
mode_wavelength = [915, 840, 780, 730, 680];

% Bandwidth (FWHM) of each mode in nm, mode_bandwidth=[band1, band2, ... ]
mode_bandwidth = [30, 25, 20, 10, 10];

% beta = 2*GDD of each mode in fs^2
beta = [-10000, 0, 0, 0, 0];

% Convert to wavelength to frequency
mode_frequency = 2*pi*c./mode_wavelength;

% Convert FWHM bandwidth in nm to FWHM in 2*pi/fs
mode_bandwidth_freq = 2*pi*c.*mode_bandwidth./mode_wavelength.^2;

% Redefine mode_bandwidth to equal 2*sigma^2
mode_bandwidth_freq = 2*mode_bandwidth_freq/2.35482;

E1=amplitude_band(1)*exp(-((ww-mode_frequency(1))/mode_bandwidth_freq(1)).^2).*exp(1i*beta(1)*(ww-mode_frequency(1)).^2);

E1_gaussian = E1_gaussian = abs(E1(:, 10));
E1_gaussian_fft = fftshift(fft(E1_guassian));

E1_time_wrong=ifftshift(ifft(E1,10000,1));
E1_time =ifftshift(ifft(E1,10000,2));

