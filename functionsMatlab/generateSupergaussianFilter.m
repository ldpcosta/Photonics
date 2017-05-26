function BPF = generateSupergaussianFilter(cFreq, fShift, orderSupergaussian, nr)
    BW_filt = 3e9.*(nr==2)+1.5e9.*(nr==1); 
    centralFrequency = cFreq(nr);
    B0 = BW_filt/(2*(log(2))^(1/(2*orderSupergaussian)));  % De onde vem estas contas?
    BPF = exp(-0.5*((fShift-centralFrequency)/(2*B0)).^(2*orderSupergaussian)); %e estas
end

