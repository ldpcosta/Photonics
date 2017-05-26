function out = quantNoise(time,in,fc,Anoise)
    
c0 = 2.99792458e8;  % Speed of light in vacuum
constPlanck   = 6.62606947e-34; % Planck constant 
constBoltsman = 1.380662e-23;   % Boltzman constant

lambdac = c0/fc;
dTime = time(2) - time(1);
nSamp = length(time);

Eph = constPlanck*c0/lambdac;	% Energy of a photon        
Pph = Eph/dTime;                % Power of a photon
% Anoise = 1;

rng('shuffle');
noise_t = Anoise*sqrt(Pph).*(randn(size(in))+1j.*randn(size(in)))./2; % 1/2 photon per temporal mode is considered 
% noise_f = dTime*fftshift(fft(ifftshift(noise_t))); % Quantum noise in spectral domain 
% in_f = dTime*fftshift(fft(ifftshift(in)));
% out_f = in_f + noise_f;
% out = 1/dTime*fftshift(ifft(ifftshift(out_f))); 

out = in + noise_t;
