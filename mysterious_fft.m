clear all

% 0) DEFINE INITIAL CONDITIONS

c = 300;
amplitude_band=[1.5,1,1,1,.1];

% Number of zeros in zero padding for FFT algorithm
% Note, the more zeros are appended after the defined frequency range,
% the better the results will be, but the longer the calculation will take.
% The current value of 5000 seems to yield satisfactory results, although
% this value may be optimized by the user. Each FFT must be optimized. As
% the current algorithm contains two FFT's, zeros are padded twice.
num_zeros = 10000;

% Number of frequency modes
num_modes = 5;

% Center of each mode in nm, mode_wavelength=[ lambda1, lambda2, ... ]
mode_wavelength = [915, 840, 780, 730, 680];
delay_mode = 2; % I.e. delay_mode 2 gives picks mode_wavelength(2) - 840 nm.
% Define wavelength range of interest, in nm, for FROG trace

lambda_low = 300;
lambda_high = 500;
% Convert wavelenth range for FROG trace to frequency in 2*pi/fs
frequency_low = 2*pi*c/lambda_high;
frequency_high = 2*pi*c/lambda_low;


% Number of frequency points
freq_num = 6000;

% Bandwidth (FWHM) of each mode in nm, mode_bandwidth=[band1, band2, ... ]
mode_bandwidth = [30, 25, 20, 10, 10];

% beta = 2*GDD of each mode in fs^2
%beta = [-10000, 0, 0, 0, 0];
beta = [0, 0, 0, 0, 0];

% Convert to wavelength to frequency
mode_frequency = 2*pi*c./mode_wavelength;

% Convert FWHM bandwidth in nm to FWHM in 2*pi/fs
%mode_bandwidth_freq = 2*pi*c.*mode_bandwidth./mode_wavelength.^2;
mode_wavelength_low = mode_wavelength - mode_bandwidth/2;
mode_wavelength_high = mode_wavelength + mode_bandwidth/2;

mode_bandwidth_freq = 2*pi*c.*(1./mode_wavelength_low - 1./mode_wavelength_high);

% Find max, min of mode frequencies
freq_max = max(mode_frequency);

band_max = max(mode_bandwidth_freq);

% Redefine mode_bandwidth to equal sqrt(2)*sigma
mode_bandwidth_freq = sqrt(2)*mode_bandwidth_freq/2.35482;
% Not sure I agree with redefining it instead of making a new variable to
% correspond to the original FWHM, but we will leave it as is.

% Define frequency axis
% ********************* Note ************************
% The frequency axis must include zero for the FFT algorithm to 
% recognize the correct distance of the modes from zero. Not including
% zero will result in error, as the FFT will associate smallest frequency
% component with zero.

% This frequency axis works to plot E1! THIS IS ANGULAR FREQUENCY.
frequency = linspace(0, freq_max+2*band_max, freq_num); 

%Define phase (i.e. delay) axis.
% phase_range=10*pi;
% phase_num=200;
phase_range=100*pi;
phase_num=2000;

phase=linspace(0,phase_range,phase_num);

[phi,ww]=meshgrid(phase,frequency);
dw=mean(diff(frequency));
Fs_pulse = dw/(2*pi)*freq_num;
phase_axis = phi.*frequency';
dphase = phase_axis(2, 2) - phase_axis(2, 3);
dt = -dphase/(dw);


% 1) Build pulse in frequency domain
E1=amplitude_band(1)*exp(-((ww-mode_frequency(1))/mode_bandwidth_freq(1)).^2).*exp(1i*beta(1)*(ww-mode_frequency(1)).^2);
% In Weiner's notation, this corresponds to A(omega)^2
% So t_p = sqrt(2)/mode_bandwidth_freq 

if (num_modes>1)
    new_field = zeros(freq_num, phase_num, num_modes-1);
    for k=2:num_modes
    new_field(:, :, k-1)=amplitude_band(k)*exp(-((ww-mode_frequency(k))./mode_bandwidth_freq(k)).^2).*exp(1i.*beta(k).*(ww-mode_frequency(k)).^2);
    end
end

% 2) Add delay to pulse in frequency domain
% E1_delayed = E1 .* exp(-1i*phi.*frequency');

% Placeholders to check delay
% PLACEHOLDER
new_field_originalff = new_field(:, :, delay_mode-1);
%
% REALITY
new_field(:, :, delay_mode-1) = new_field(:, :, delay_mode-1) .* exp(-1i*phi.*frequency');
interfere_field = E1 + sum(new_field, 3);
%

% PLACEHOLDER
new_field_delayedff = new_field(:, :, delay_mode+1);
%

% 3) Fourier Transform from pulse freq to pulse time

delay_time = 0:dt:(phase_num-1)*dt;

% delay_time = -phase_num/2:phase_num/2 - 1;
E1_time = ifftshift(ifft(E1,num_zeros));
E1_interfered =ifftshift(ifft(interfere_field,num_zeros));
E1_interfered = ifftshift(E1_interfered, 2);

Esig = (E1_interfered).^2;

pulse_time = (-num_zeros/2:num_zeros/2-1)*freq_num/(2*num_zeros)/(Fs_pulse);

% PLACEHOLDERS
new_field_originalt = ifftshift(ifft(new_field_originalff,num_zeros));
new_field_delayedt =ifftshift(ifft(new_field_delayedff,num_zeros));
new_field_delayedt = ifftshift(new_field_delayedt, 2);


% 4) Fourier Transform back to pulse freq domain
% Defining the axes for this step goes really poorly for some reason.

E1_freq_delayed = fftshift(fft(Esig, num_zeros));
E1_freq_delayed = ifftshift(E1_freq_delayed, 2);

% DEFINING WAVELENGTH AXIS OF SPECTROGRAM

frequency_final = 0:dw:(num_zeros-1)*dw;

wavelength_final = 2*pi*c./frequency_final;
wavelength_final = wavelength_final(1001:num_zeros); 
% The first 1001 points are junk (basically DC and noice anyway).
% Unfortuantely, as a result you still get 1/x, which is a nonlinear axis.
% We need a liner one for imagesc.

wavelength_equal = linspace(2*pi*c/frequency_final(1001), ...
    2*pi*c/frequency_final(end), length(wavelength_final));

E1_freq_delayed_final = E1_freq_delayed(1001:num_zeros, :);

% INTERPOLATE FOR LINEAR AXIS.
E1_field_transfinal=interp1(wavelength_final, E1_freq_delayed_final, ...
    wavelength_equal, 'pchip');

imagesc(delay_time, wavelength_equal, abs(E1_field_transfinal));




% PLACEHOLDERS
new_field_transformed = fft(new_field_originalt, num_zeros);
new_field_dtransformed = fft(new_field_delayedt, num_zeros);

frequency_final = 0:dw:(num_zeros-1)*dw;

new_field_transhort = new_field_transformed(1001:num_zeros, :);
new_field_dtranshort = new_field_dtransformed(1001:num_zeros, :);

wavelength_final = 2*pi*c./frequency_final;
wavelength_final = wavelength_final(1001:num_zeros);

wavelength_equal = linspace(2*pi*c/frequency_final(1001), ...
    2*pi*c/frequency_final(end), length(wavelength_final));


new_field_transfinal=interp1(wavelength_final, new_field_transhort, ...
    wavelength_equal, 'pchip');


new_field_dtransfinal=interp1(wavelength_final, new_field_dtranshort, ...
    wavelength_equal, 'pchip');


