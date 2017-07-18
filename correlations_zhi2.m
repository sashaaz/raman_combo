% Here we now apply what we've learned to Miaochan's OL (screenshot attached)
% Basically, instead of applying delay via t-tau, we apply delay via the
% phase delay (which is actually presumably what happened in the paper).



c = 2.99792e8;

% Setting up all beams included in the paper.
wavelengths = [930e-9, 780e-9, 720e-9, 680e-9]; % AS 1, 3, 4, 5,
% wavelengths = [1064e-9, 807e-9, 650e-9, 468e-9]; % AS 1, 3, 4, 5,

% converting to frequency space
omega = 2*pi*c*(wavelengths).^-1;

% Only scanning one wavelength - this defines that one wavelength.
wavelength_scan =  840e-9;
omega_scan = 2*pi*c*(wavelength_scan).^-1;

% Original pulse length given in the paper
t_p = 100e-15;

% Added the possibility to include chirp - don't think this will be used so
% it is set to an absurdly high number.
t_c = 250e-9;

% Time grid of the pulse (presumably)
delta_t = 0.5e-16; 
t = -100e-15:delta_t:(100e-15-delta_t);

% Time grid of the delay time.
tau = -100e-15:delta_t:(100e-15-delta_t);
% [T, Tau] = meshgrid(t, tau);

Fs = 1/delta_t;

time_d = size(t);
time_d = time_d(1, 2);

phase = exp(1i*tau' .* 1*omega_scan);
% NON MESHGRID
a_t = real(exp(-t.^2/t_p^2).*exp(1i*t.^2/t_c^2));

% a_tau = exp(-(t-tau).^2/t_p^2).*exp(1i*(t-tau).^2/t_c^2);

% MESHGRID
% a_t = real(exp(-T(2000, :).^2/t_p^2).*exp(1i*T(2000, :).^2/t_c^2));
% a_tau = exp(-(T-Tau).^2/t_p^2).*exp(1i*(T-Tau).^2/t_c^2);
% y_array = real(a_t .* exp(1i*(omega' .* T(2000, :))));

y_array = real(a_t .* exp(1i*(omega' .* t)));
y_array = y_array .* [1, 1, 1, .1]';

pulse_1 = sum(y_array);
pulse_2 = 1*real(a_t .* exp(1i*(omega_scan .* t))).*phase;

% pulse_im = abs((pulse_1 + pulse_2).^2).^2;
pulse_im = abs(pulse_1).^2 .* abs(pulse_2).^2;
%intensity_1 = abs(pulse_1).^2;
%intensity_2 = abs(pulse_2).^2;

%pulse_cc = mean(intensity_1 .* intensity_2);
pulse_cc_sasha = mean(pulse_im);
pulse_cc_peter = mean(pulse_im, 2);

pulse_envelope = envelope(pulse_cc_sasha,30, 'peak');
% tau*omega_scan/(12*pi)
% plot(tau/2, pulse_envelope);
% hold on
% plot(tau/2, pulse_cc_peter);
% legend('Sasha, mean 1, more or less gets result in paper', 'Peter, mean 2');
% plot(tau/2, pulse_envelope);
%plot(t, pulse_1);

figure
plot(tau/16 * omega_scan, pulse_envelope);
hold on
plot(tau/16 * omega_scan, pulse_cc_peter);
