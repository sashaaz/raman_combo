% This was my first attempt to recreate the correlation plots found in
% earlier papers by people other than Sokolov (so that we can kind of trust
% them more).

% This program recreates the results in Shverdin's 2005 PRL, included in
% the shared directory. Even here, we run into the "mean" problem that
% plagues our own correlation stuff.

% Did not initially use a meshgrid for this - integrating meshgrid is
% proving to be a bit tricky. Leaving it as is and fixing later.

c = 3e8;
% Setting up all beams included in the paper.
wavelengths = [1.56e-6, 1.064e-6, 807e-9, 650e-9, 544e-9, 468e-9, 410e-9];

% converting to frequency space
omega = 2*pi*c*(wavelengths).^-1;

% Original pulse length given in the paper
t_p = 15e-9;

% Added the possibility to include chirp - don't think this will be used so
% it is set to an absurdly high number.
t_c = 60e-15;

% Wavelength that is being scanned is 1064 nm.
omega_0 = 2*pi*c/1.064e-6;
% Raman shift - seperates all beams. Included for testing purposes.
omega_a = 2994 *10^2 * c;

% Time grid of the pulse (presumably)
t = -30e-15:0.5e-16:(30e-15-0.1e-16);
% Time grid of the delay time.
tau = -30e-15:0.5e-16:(30e-15-0.1e-16);
% assumed same parameters as pulse time.
% [T, Tau] = meshgrid(t, tau);

% Fourier stuff - included but not needed.
% Fs = 1/5e-16;
% 
% Time grid size, needed to accurately recombine the beams later.
time_d = size(t);
time_d = time_d(1, 2);

tau_d = size(tau);
tau_d = tau_d(1, 2);

% Pulse envelope
a_t = real(exp(-t.^2/t_p^2).*exp(1i*t.^2/t_c^2));
% a_t = real(exp(-T(:, 600).^2/t_p^2).*exp(1i*T(:, 600).^2/t_c^2));

% Giving each beam a pulse enevelope and osillations.
y_array = real(a_t .* exp(1i*(omega' .* t)));
% y_array = real(a_t .* exp(1i*(omega .* T(:, 600))));

% In the paper, Shverdin et al delayed all the even beams with respect to
% the odd beams. This array sets that up.
y_even = zeros(4, time_d, tau_d); % 3D array necessary. 
% One axis for each beam, one for the pulse time, and the last for the
% delay time.

y_even2 = zeros(4, tau_d); % 2D array where pulse-time averaged beam will go.

y_odd = zeros(3, tau_d);

phase = exp(1i*tau' .* 1*omega); % Converting delay time to added phase.
% phase = exp(1i*Tau(600, :)' .* 1*omega);
% Splitting up the beams into even and odd arrays, giving each even beam a
% delay defined by phase above.

for i = 1:4
    y_even(i, 1:time_d, 1:tau_d) = y_array(2*i-1, 1:time_d).*phase(1:time_d, 2*i-1);
    % y_even2(i, 1:time_d) = y_array(2*i-1, 1:time_d);
%     y_even(i, 1:time_d, 1:time_d) = y_array(1:time_d, 2*i-1).*phase(1:time_d, 2*i-1)';
%     y_even2(i, 1:time_d) = y_array(1:time_d, 2*i-1);
end
for i = 1:3
    y_odd(i, 1:time_d) = y_array(2*i, 1:time_d);
%     y_odd(i, 1:time_d) = y_array(1:time_d, 2*i);
end

% Averaging all beams into a single array for delayed even beams
% (pulse_1), undelayed even beams (pulse_12), and odd beams (pulse_2).
pulse_1 = sum(y_even);
pulse_12 = sum(y_even2);
pulse_2 = sum(y_odd);

% Add delayed even pulse and odd pulse together.

pulse_1clean = zeros(time_d, time_d);
pulse_1clean(1:time_d, 1:time_d) = pulse_1(1, 1:time_d, 1:time_d);

intensity_1 = abs(pulse_1clean).^2;
intensity_2 = abs(pulse_2).^2;

% Ok, as far as I understand, this is the core point where me and you
% disagree.
% If you plot it along time axis 1, you get the result of Shverdin's 2005
% included in the shared directory. If you do
% it along axis 2, then you get just gaussian peaks.

pulse_cc_sasha = mean(intensity_1 .* intensity_2);
hold on
pulse_cc_peter = mean(intensity_1 .* intensity_2, 2);

plot(tau, pulse_cc_sasha);

hold on
plot(tau, pulse_cc_peter);
legend('Mean along axis 1 - correctly reproduces paper', 'Mean along axis 2')


