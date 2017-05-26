function Pmeas_f = zeroPadAndFFT(Pmeas, nSamp)
    Pmeas = [Pmeas zeros(size(Pmeas,1),nSamp-size(Pmeas,2))]; 
    Pmeas_f = fftshift(fft(ifftshift(Pmeas.'))).';
end


