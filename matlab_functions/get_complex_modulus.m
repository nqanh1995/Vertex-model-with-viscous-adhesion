function [Gprime, Gprimeprime, omega1, omega2] = get_complex_modulus(stress, strain, omega, dt)
% GET_COMPLEX_MODULUS  Extract the storage and loss moduli (G', G'') from
% stress and strain time series under oscillatory (harmonic) shear.
%
%   [Gprime, Gprimeprime, omega1, omega2] = GET_COMPLEX_MODULUS(stress,
%   strain, omega, dt) discards the first half of the signals (to remove
%   transient/start-up effects), takes the FFT of both stress and strain,
%   locates the drive frequency omega in the spectrum, and computes the
%   in-phase (storage, G') and out-of-phase (loss, G'') moduli from the
%   phase lag between stress and strain at that frequency.
%
% Inputs:
%   stress, strain - time series of measured stress and applied strain
%   omega          - nominal (angular) driving frequency of the shear
%   dt             - simulation time step (sampling interval)
%
% Outputs:
%   Gprime      - storage modulus (elastic response, in-phase with strain)
%   Gprimeprime - loss modulus (viscous response, out-of-phase)
%   omega1      - angular frequency identified from the strain spectrum's
%                 dominant peak nearest to "omega"
%   omega2      - angular frequency bin nearest to "omega" directly

% Drop the first half of the trajectory to discard the transient response
% and only analyze the steady-state oscillatory regime.
stress = stress(round(length(stress)/2):end);
strain = strain(round(length(strain)/2):end);

sampling_freq = 1 / dt;
frequency = sampling_freq * (0:length(strain)-1)' / length(strain);

sigma = fft(stress)';
gamma = fft(strain)';

% Identify the true drive frequency by finding the strain spectrum's
% three largest peaks and picking whichever is closest to the nominal
% "omega" (guards against picking a numerical/DC artifact instead of the
% actual drive frequency).
[~, k] = maxk(abs(gamma(1:end)), 3);
[~, id] = min(abs(frequency(k)*2*pi - omega));
k1 = k(id);
omega1 = frequency(k1) * 2 * pi;

% Alternative: simply the frequency bin nearest to the nominal omega.
[~, k1] = min(abs(frequency*2*pi - omega));
omega2 = frequency(k1) * 2 * pi;

% Phase lag between stress and strain at the identified frequency bin.
phi_g = angle(gamma(k1));
phi_s = angle(sigma(k1));
delta = phi_g - phi_s;

% In-phase (storage) and out-of-phase (loss) components of the complex
% modulus, from the stress/strain amplitude ratio and phase lag.
Gprime = abs(sigma(k1)) / abs(gamma(k1)) * cos(delta);
Gprimeprime = abs(sigma(k1)) / abs(gamma(k1)) * sin(delta);
end
