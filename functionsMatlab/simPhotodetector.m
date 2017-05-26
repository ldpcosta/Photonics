function [ PtracePD ] = simPhotodetector( trace11, f_tr, fc, BW_PD, Nsg)
%Simulates a photodetector Photodetector(signal, trace_freq_vector, trace center frequency, bandwidth, supergaussian filter order/2[OPTIONAL])

if nargin < 5
    Nsg = 4;
end

B0 = BW_PD/(2*(log(2))^(1/(2*Nsg)));
spect_PD = exp(-0.5*((f_tr-fc)/(2*B0)).^(2*Nsg));

Ptrace    = abs(trace11).^2;                    
Ptracef   = fftshift(fft(ifftshift(Ptrace)));
PtracefPD = Ptracef .* spect_PD;
PtracePD  = fftshift(ifft(ifftshift(PtracefPD)));

end

