%function [] = SHG_FROG_multimode()
clear
% Author: Benjamin Strycker, Texas A&M University
% contact: ben.strycker@gmail.com
% Copied and modified codes from Jake Cohen and Vikrant Chauhan
%     from Georgia Tech University.

% This program calculates the SHG FROG trace of a multimodal pulse

% Specify unit of time as fs, unit of frequency as 1/fs
% Speed of light c is then
c = 300;

%*********************************************************************
%*********************************************************************

% Define wavelegnth range of interest, in nm, for FROG trace
lambda_low = 300;
lambda_high = 500;

% Number of frequency points
freq_num = 400;

% Number of temporal translation delay points
delay_num = 200;

% Range of temporal translation delay in fs
delay_range = 300;

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

% Bandwidth (FWHM) of each mode in nm, mode_bandwidth=[band1, band2, ... ]
mode_bandwidth = [30, 25, 20, 10, 10];

% Temporal delay of each mode in fs. One mode delay must be zero.
mode_delay = [0, 0, 0, 0, 0];

% Absolute temporal phase of each mode
mode_phase = [0*pi, 0*pi, 0.0*pi, 0*pi, 0.*pi];

% beta = 2*GDD of each mode in fs^2
beta = [-10000, 0, 0, 0, 0];

amplitude_band=[1.5,1,1,1,1];


%*********************************************************************
%*********************************************************************
% Convert wavelenth range for FROG trace to frequency in 2*pi/fs
frequency_low = 2*pi*c/lambda_high;
frequency_high = 2*pi*c/lambda_low;

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
% Find frequency step
dw=mean(diff(frequency));


% Define delay axis
%delay = linspace(-delay_range/2, delay_range/2, delay_num);
%[tau, ww]=meshgrid(delay, frequency);

%Define phase axis
phase_range=10*pi;
phase_num=200;
phase=linspace(0,phase_range,phase_num);

[phi,ww]=meshgrid(phase,frequency);
    
% Define electric field of pulses in frequency space
E1=amplitude_band(1)*exp(-((ww-mode_frequency(1))/mode_bandwidth_freq(1)).^2).*exp(i*beta(1)*(ww-mode_frequency(1)).^2)...
    .*exp(i.*(ww-mode_frequency(1)).*mode_delay(1)).*exp(i*mode_phase(1));
if (num_modes>1)
    for k=2:num_modes
    new_field=amplitude_band(k)*exp(-((ww-mode_frequency(k))./mode_bandwidth_freq(k)).^2).*exp(i.*beta(k).*(ww-mode_frequency(k)).^2)...
    .*exp(i.*(ww-mode_frequency(k)).*mode_delay(k)).*exp(i*mode_phase(k));
     if (k==3)%%%%%%%%%%%%%%%%%%%%%%%%%%scan asi (k=)
        %E1=E1+new_field;
        E1=E1+new_field.*exp(i.*phi);
     else
        E1=E1+new_field;
        %E2=E2+new_field.*exp(i*(ww-mode_frequency(k)).*tau);
        end
     end
end

E1_time=ifftshift(ifft(E1,10000,1));
   

dt = 2*pi/(num_zeros*dw);
for k=1:num_zeros
        temporal_axis(k) = 0+(k-1)*dt;
end
    
    
% Multiply the complex electric fields in time domain for Esig
Esig = (E1_time).^2;
    
        
% Transform signal back into the frequency domain

Isig0 = abs(fftshift((fft(Esig,(2*num_zeros),1)),2)).^2;
    
dw2 = 2*pi/(2*num_zeros*dt);
frequency2=zeros(1,(2*num_zeros));
for k=1:(2*num_zeros)
    frequency2(k) = 0+(k-1)*dw2;
end
    
    
% The current Isig contains much extraneous data than we want. We want
% just the frequency range of interest. 

[test1, min_index] = min(abs(frequency2-frequency_low));
[test2, max_index] = min(abs(frequency2-frequency_high));
frequency3=zeros(1,(max_index-min_index));
for k=1:(max_index-min_index)
    frequency3(k) = frequency2(min_index+(k-1));
end    
lambda2=2*pi*c./frequency3;
lambda_equal=linspace(2*pi*c/frequency2(min_index), ...
    2*pi*c/frequency2(max_index), length(frequency3));

  
% Interpolate to get equally spaced wavelength data
    
    
Isig1=interp1(frequency2, Isig0, frequency3, 'cubic');
Isig2=interp1(lambda2, Isig1, lambda_equal, 'cubic').*10e6;
phase=phase./pi;

figure
imagesc(phase, lambda_equal, Isig2),xlabel('phase'),ylabel('wavelength (nm)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%frequency%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% length(frequency3) is the pixel number of wavelength %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Save the data, the path could be added before the file name%%%%%%%%%%%%%%%%%%
save('Isig.dat','-ascii','Isig2');
%%%%%%Save the phase%%%%%%%%%%%%
save('phase.dat','-ascii','phase');
%%%%%%%%%%Save the lambda_equal%%%%%%%%
save('lambda_equal.dat','-ascii','lambda_equal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
down=1;
up=length(frequency3);
b=zeros(1,phase_num);
for m=1:phase_num
    for ij=down:up
    b(m)=b(m)+Isig2(ij,m);
    end
    
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%subplot(7,1,1),plot(phase,b),ylabel('Total')
 
 
%%%%%%%%%Pick up the right wave length%%%%%%
k1=175;
wavelength_shg=num2str(lambda_equal(k1));
a=Isig2(k1,:);
subplot(7,1,1), plot(phase,a),ylabel('AS 1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SHG_As1.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2=238;
wavelength_shg=num2str(lambda_equal(k2));
a=Isig2(k2,:);
subplot(7,1,2), plot(phase,a),ylabel('AS12')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SFG_As1AS2.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k3=294;
wavelength_shg=num2str(lambda_equal(k3));
a=Isig2(k3,:);
subplot(7,1,3),plot(phase,a),ylabel('AS2');%title(wavelength_shg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SHG_As2.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=346;
wavelength_shg=num2str(lambda_equal(k));
a=Isig2(k,:);
subplot(7,1,4), plot(phase,a),ylabel('AS23');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SFG_As2As3.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=395;
wavelength_shg=num2str(lambda_equal(k));
a=Isig2(k,:);
subplot(7,1,5), plot(phase,a),ylabel('AS3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SHG_As3.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=440;
wavelength_shg=num2str(lambda_equal(k));
a=Isig2(k,:);
subplot(7,1,6), plot(phase,a),xlabel('AS2 phase'),ylabel('AS34');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SFG_As3As4.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=482;
wavelength_shg=num2str(lambda_equal(k));
a=Isig2(k,:);
subplot(7,1,7),plot(phase,a),xlabel('AS4'),ylabel('AS4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SHG_As4.dat','-ascii','a');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=521;
wavelength_shg=num2str(lambda_equal(k));
a=Isig2(k,:);
%subplot(4,2,8),plot(phase,a),xlabel('AS45'),ylabel('Arb.unit'),title(wavelength_shg);

%k=560;
%wavelength_shg=num2str(lambda_equal(k));
%a=Isig2(k,:);
%figure
%plot(phase,a),xlabel('AS5'),ylabel('Arb.unit'),title(wavelength_shg);



%end