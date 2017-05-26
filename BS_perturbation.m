% RAYLEIGH BACKSCATTERING SIMULATION AND SSFM
%
% =========================================================================

clc, clear all, close all, 
format long,

%% Settings
% --------
DimEcran = get(0,'ScreenSize'); 
V = 0.15; H = 0.15;
PosFig  = [V*DimEcran(3),(H)*DimEcran(4),(1-V)*DimEcran(3),(0.9-H)*DimEcran(4)];
set(0,'DefaultFigurePosition',PosFig),
AFS = 24;   LFS = 30;   TFS = 30;       % Axis, label and title font size
set(0,'DefaultLineLineWidth',2)

%% Functions
% ---------

% Description of the cases
cases = ['100ns'];      % LUIS: Nao e usado.

% Physical variables
c0 = 2.99792458e8; 
lambda0 = 1550e-9;  % 1560e-9;
fc = c0/lambda0;    % LUIS: Converter de lambda em f
omegac = 2*pi*fc;   % LUIS: Converter de Hz em Rad
n0 = 1.46;          % Indice fibra
ng = 1.01*n0;       % group index
energy = 3e-8;    % Pulses energy (J)


%% Time and Frequency Vectors
% --------------------------

% Define frequency window.
windowF = 200e9;    nSamp  = 2^16;        % LUIS: 200e9 ~ 1.6nm around 1550nm
fMax   = fc + windowF/2; 
fMin = fc - windowF/2;
freq   = linspace(fMin, fMax, nSamp+1);     freq   = freq(1:end-1);  % Optical frequency vector
interv = freq(end) - freq(1);                   % Teoricamente = a windowF
dFreq  = freq(2)-freq(1);
         
lambda = c0./freq;                           
lambda = linspace(lambda(end),lambda(1),nSamp); % Optical wavelength vector 
fShift = freq - fc;                             % Shifted frequency vector

dTime  = 1/interv;
time   = fix(0:nSamp-1)*dTime;                  % Time vector
timec  = fix(-nSamp/2:nSamp/2-1)*dTime;         % Time vector center at origin
omeg   = 2*pi*freq;


%% Characteristics of the sensing optical fiber
% --------------------------------------------
% SMF-28
Lspool = 2000;           % Medium length [m]
Gamma  = 1.3e-3;        % Nonlinear coefficient of the medium [1/W/km] (1e-3)
TR     = 3e-15;         % Intra-pulse Raman coefficient 
alpha_dBKm = 0.2;       Alpha = alpha_dBKm*1e-3/4.343; % 0.2
Beta2 = -2.28e-26;      D = -Beta2*(2*pi*c0/lambda0^2);         % Beta2_smf = -2.17e-26 s.s/m;   D_smf = 17e-6 s/m.m = 17 ps/nm.km
% -2.17e-26                       % Beta2 = -D/(2*pi*c0/lambda0^2);
Beta3 = 10e-41;         S = Beta3*(2*pi*c0/lambda0^2)^2;        % S_hnlf = 32 s/m^2.m = 0.032 ps/nm^2.km % Me sale S_smf = 54.7 - 58.5 s/m^2.m 
 % 8.9e-41;                        % Beta3 = S/(2*pi*c0/lambda0^2)^2;
                     

%% Definition input pulse
% ----------------------
FWHM = 20e-9;	% FWHM     
P0 = 0.2;    
Nsg = 10;       % Order of the super-gaussian pulse envelope

FW = FWHM*(log(100)/log(2))^(1/(2*Nsg));   % Full width at 1% of the peak value
T0 = FWHM/(2*(log(2))^(1/(2*Nsg)));        % FWHM (1/e)
          
chirp = 0;
freqInst = chirp*timec/FWHM; 
phaseIn = cumtrapz(timec,freqInst); %cumtrapz -> integral

eaux_t = sqrt(P0).*exp(-0.5*(timec/T0).^(2*Nsg)).*exp(2j*pi*phaseIn);
Ein_t = quantNoise(time,eaux_t,fc,10);          % Input amplitude          quantNoise introduz ruido em amp e fase
Pin_t = abs(Ein_t).^2;                          % Imput power
Ein_f = fftshift(fft(ifftshift(Ein_t)));        % Input spectrum
Pin_f = (abs(Ein_f)).^2;                        % Input power spectral density (in frequency)

Energy = trapz(time, Pin_t);

% figure,plot(timec/1e-9,abs(Ein_t)/max(abs(Ein_t)),'-','Color',[0 0.7 0]),
% xlabel('Time [ns]','fontsize',LFS), ylabel('Input amplitude [a.u.]','fontsize',LFS)
% figure,plot(freq/1e12,abs(Ein_f)/max(abs(Ein_f)),'-','Color',[0 0 0.7]),
% xlabel('Frequency [THz]','fontsize',LFS), ylabel('Input spectrum [a.u.]','fontsize',LFS)




%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCHRODINGER EQUATION RESOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalError = 1e-4;  
dZ_start = 1; % [m] 200!!!
Lcons = 0;   counter = 1;
step_inc = [];
 
% Dispersion operator
Dop	= -(1i/2)*Beta2*(omeg-omegac).^2 -(1i/6)*Beta3*(omeg-omegac).^3 - 0.5*Alpha;    

% Pulse propagation
At_dep = Ein_t;      % Set row 1 of E(z,t) = E(0,t)
Af_dep = Ein_f;      % Set row 1 of E(z,w) = FT of Ein
    
At_inc = []; Af_inc = [];
while Lcons <= Lspool
          
    % Coarse solution
    dz = 2*dZ_start;
    Position = Lcons + dZ_start;
                 
    DispExp = exp(Dop*dz/2);    % OJO: dz/2
    Af = Af_dep.*DispExp;  
    At = fftshift(ifft(ifftshift(Af)));

    PInst = abs(At).^2;
    dPInst = diff(PInst)./diff(time);       dPInst = [dPInst dPInst(end)];
    dPinAt = diff(PInst.*At)./diff(time);   dPinAt = [dPinAt dPinAt(end)];
    Nop = -1i*Gamma*(PInst + 1j./(omegac.*At).*dPinAt - TR.*dPInst);    % Non-linear operator
    At = At.*exp(dz*Nop);  At(isnan(At)) = eps;    % OJO: dz
    Af = fftshift(fft(ifftshift(At)));
    
    Af = Af.*DispExp;     % OJO: dz/2
    coarseSol_f = Af;

    % Fine solution
    dz = dZ_start;
    
    DispExp = exp(Dop*dz/2);    % OJO: dz/2
    Af_fine = Af_dep;

    for k = 1:2  
        Af = Af_fine.*DispExp;  
        At = fftshift(ifft(ifftshift(Af)));

        PInst = abs(At).^2;
        dPInst = diff(PInst)./diff(time);       dPInst = [dPInst dPInst(end)];
        dPinAt = diff(PInst.*At)./diff(time);   dPinAt = [dPinAt dPinAt(end)];
        Nop = -1i*Gamma*(PInst + 1j./(omegac.*At).*dPinAt - TR.*dPInst);    % Non-linear operator
        At = At.*exp(dz*Nop);  At(isnan(At)) = eps;
        Af = fftshift(fft(ifftshift(At)));
        
        Af = Af.*DispExp;     % OJO: dz/2 
        Af_fine = Af;
    end
    fineSol_f = Af_fine;

    % Calculation of the local error
    fineSol_t = fftshift(ifft(ifftshift((fineSol_f))));
    coarseSol_t = fftshift(ifft(ifftshift((coarseSol_f))));
    localError = sqrt(sum(abs(fineSol_t-coarseSol_t).^2)./sum(abs(fineSol_t).^2));
    localErrorReg(counter) = localError;

    if localError > globalError*2 % We discard the soltion and divide dZ_start / 2
        Af_dep = Af_dep;
        dZ_start = dZ_start/2;
        Lcons = Lcons; 
    elseif localError <= globalError*2 
        if localError > globalError*0.5 % We save the soltion and divide h / 2^1/3
            Af_dep = fineSol_f;
            Lcons = Lcons + 2*dZ_start;
            step_inc(counter) = 2*dZ_start;
            dZ_start = dZ_start/2^(1/3);
        elseif localError <= globalError*0.5 % We save the soltion and multiply h * 2^1/3
            Af_dep = fineSol_f;
            Lcons = Lcons + 2*dZ_start;
            step_inc(counter) = 2*dZ_start;
            dZ_start = dZ_start*2^(1/3);
        end
        At_inc(counter,:) = fineSol_t;  
        Af_inc(counter,:) = fineSol_f;  
        counter = counter+1;
    end
     
%         figure(99),plot(freq/1e9,20*log10(abs(Af_dep))),axis([(fc-windowF/2)/1e9,(fc+windowF/2)/1e9 -260 -110])
% %         figure(99),plot(10*log(Pout_f)),% axis([3000,7000 -inf inf])
%         drawnow,
%         
    clc; PercentSSFM = fix(100*Lcons/Lspool)
end

dZinc = [0 cumsum(step_inc)];

Eout_f = Af_dep;
Pout_f = dFreq*(abs(Eout_f)).^2;
Eout_t = fftshift(ifft(ifftshift(Eout_f))); 
Pout_t = abs(Eout_t).^2;

phaseOut_t = unwrap(angle(Eout_t));    phaseOut_t = phaseOut_t - phaseOut_t(1);
chirpOut = diff(phaseOut_t)/(dTime);   chirpOut(nSamp)=0;

% % figure,plot(step_inc)
% 
% figure,set(gca,'FontSize',AFS)
% subplot(2,1,1),plot(timec/1e-12,Pin_t,'b-+',timec/1e-12,Pout_t,'g-')
% % axis([-2e-10 2e-10 0 1.4])
% legend('In','Out')
% xlabel('Time [ps]','fontsize',LFS), ylabel('Power [W]','fontsize',LFS)
% subplot(2,1,2),plot(freq/1e12,Pin_f,'b-*',freq/1e12,Pout_f,'g+-')
% % axis([fc-1e12 fc+1e12 0 1.4e-13])
% xlabel('Frequency [THz]','fontsize',LFS), ylabel('PSD [W/Hz]','fontsize',LFS)
% 
% % figure,imagesc(20*log10(abs(Af_inc)))
% figure,mesh(freq(3700:end-3700)/1e12,dZinc(1:end-1).'/1e3,10*log10(abs(Af_inc(:,3700:end-3700))))



%% Scattering properties
% ---------------------
Lsensor = dZinc(end);
densScatt = round(n0/(dTime*c0));   % Scattering centers density (points/m)   (1/vf * dt) = 1/(dist percorrida pela luz em dt)

scatPulse = floor(FW/dTime);        % Number of scattering points whithin the pulse width   NAO IMPORTA PARA NADA ACHO EU
        
S   = 1/dTime;                      % Sampling ratio (Hz)
Msc = fix(densScatt*Lsensor);       % Number of scattering point along the sensor fiber (Lsensor/dist percorrida pela luz a cada dT = um ponto a cada dT).
% de  = 25;                   % Strain induced (um/m)
% dn  = de*1e-6*n0;           % Index change induced    
% reflectorFlag = 0;          % Generate a reflector at the end of the fiber







%% Calculation of the backscattered signal
% ---------------------------------------
pos = 0:Lsensor/Msc:Lsensor-Lsensor/Msc; % Position of scattering centers 
dPos = Lsensor/Msc;
ma = sqrt(dTime/FW);                     % Scatter's amplitude (no absolute signal levels)

da = ma/4;                              % LUIS: porque sim  -> VARIANCE OF SCATTERING POINTS REFLECTION COEFFICIENT

r = abs(ma+da*randn(1,Msc));             % Scattering centers reflectivity (Juan: 2.5e-9^rand)
phaRand = 2*pi*rand(1,Msc);              % Random phase associated to scattering centers
%     r(end) = 1e-7*(reflectorFlag==1);        % Reflector at the end of the fiber




% Time and frequency vectors for the trace
dt_tr = 1/S;                                                        % S = 1/dTime <=> dt_tr == dTime
t_tr  = (0:(nSamp+2*(Msc-1)-1))*dt_tr;                              % Trace time vector (s)   ->   LUIS: 2 * M e para upsampling para ter em conta a reflexao tambem. 
        % LUIS: Porque len=nSamp +2*(Msc-1)-1??

df_tr = S/(nSamp+2*(Msc-1));                                        % Frequency sampling rate (Hz)            
f_tr  = (-(nSamp+2*(Msc-1))/2:(nSamp+2*(Msc-1))/2-1)*df_tr + fc;	% Trace frequency vector (Hz)
zPos = c0*t_tr/n0/2;

% generate perturbation
startPerturbation = 1200000; % what index does the perturbation start
endPerturbation = 1200051;   % what index does the perturbation end
measurement_time = 0:20;     % how long to simulate
perturbation_dn = 1e-6*n0;
perturbation_period = 10;
perturbation = n0 + perturbation_dn * sin(measurement_time * (2*pi) / perturbation_period);

traceMatrix = zeros(size(measurement_time,2),nSamp+(Msc-1)*2);

for measurement_dT = measurement_time+1

    nVect = n0*ones(1,Msc);
    nVect(startPerturbation:endPerturbation) = perturbation(measurement_dT);
    %     nVectS(round(Msc/2):round(Msc/2+1*FW*S)) = nVectS(round(Msc/2):round(Msc/2+1*FW*S))+dn;  

    ph  = phaRand+2*pi.*cumsum(nVect)/lambda0.*dPos;     % Total phase of scattering centers   
    %     phS = phaRand+2*pi.*cumsum(nVectS)/lambda0.*dPos;    % Total phase of scattering centers (with strain)
    prodRP  = r.*exp(1j*ph).* exp(-0.5*Alpha*pos);   % Complex r(z)=r*exp(i*ph);
    %     prodRPS = r.*exp(1j*phS).* exp(-0.5*Alpha*pos);  % Complex r(z)=r*exp(i*ph);

    changePulse = Lsensor/10;  % Propagation length at which the pulse changes its shape in simulation (m)
    nChanges = fix(Lsensor/changePulse);
    lengthChange = fix(Msc/nChanges);


    %% Code:conv %
    % ========= %
    prodRP2 = upsample(prodRP,2);           % prodRPS2 = upsample(prodRPS,2);      prodRP is the reflection function for each point in the fiber

    trace11  = zeros(1,nSamp+(Msc-1)*2);   %  trace2  = zeros(1,nSamp+(Msc-1)*2);
    trace11f = zeros(1,nSamp+(Msc-1)*2);   %  trace2f = zeros(1,nSamp+(Msc-1)*2);
    NconvAux = 2*lengthChange;

    % convolve beam with fiber


    for pp = 1:nChanges
        [void,indC] = isclose2(dZinc(1:end-1),pos(pp*lengthChange)); %encontrar o ponto da fibra mais perto do ponto definido para actualizar o pulso ??

        %     in = fliplr(At_inc(indC,:));
        in_f = Af_inc(indC,:);
        in_fD = in_f.*exp(-2j*(pi*fShift).^2*Beta2*pos(pp*lengthChange));   % aplicar a dispersao ao pulso
        in = fliplr(fftshift(ifft(ifftshift(in_fD))));  % convol pulso com disperso

        trace11((1+NconvAux*(pp-1)):(NconvAux*pp+nSamp-1)) = trace11((1+NconvAux*(pp-1)):(NconvAux*pp+nSamp-1)) + ...
                                                             conv(prodRP2((1+NconvAux*(pp-1)):NconvAux*pp),in);


        clc; PercentBS = fix(100*pp/nChanges)
    end
    traceMatrix(measurement_dT,:) = trace11;
end
%figure,plot(zPos,abs(trace11))


%% Photodetection
% --------------
% seno1 = [];
% seno2 = [];

BW_PD = 500e6;
LO = 0.2


for measurement_dT = measurement_time+1
    clc; disp(measurement_dT);
    trace11 = traceMatrix(measurement_dT,:);
    BW_PD = 500e6; 

    LO = 0.2*ones(1,length(trace11));

    PtracePDI1 = simPhotodetector(trace11 - LO, f_tr, fc, BW_PD);
    PtracePDI2 = simPhotodetector(1j*trace11 + 1j*LO, f_tr, fc, BW_PD);

    PtracePDQ1 = simPhotodetector(1j*trace11 - LO, f_tr, fc, BW_PD);
    PtracePDQ2 = simPhotodetector(-trace11 + 1j*LO, f_tr, fc, BW_PD);

    I_ = PtracePDI2 - PtracePDI1; 
    Q_ = PtracePDQ1 - PtracePDQ2;

    % Generate Noise
    Inoise = (4.62 * 1e-3 * randn(1,length(I_)));
    Qnoise = (4.62 * 1e-3 * randn(1,length(Q_)));

    I = I_ + Inoise;
    Q = Q_ + Qnoise;
    if measurement_dT == 1

        Magnitude_ref = sqrt(I.^2 + Q.^2); 
        %Phase_ref = unwrap(atan2(Q,I));% - median(unwrap(atan2(Q,I)));
        Phase_ref = atan2(Q,I);
        
        %angleTheo_ref = unwrap(angle(trace11));% - median(unwrap(angle(trace11)));
        angleTheo_ref = unwrap(angle(trace11));
        magTheo_ref = abs(trace11);
 
    else
        Magnitude = sqrt(I.^2 + Q.^2); 
        %Phase = unwrap(atan2(Q,I));% - median(unwrap(atan2(Q,I)));
        Phase = atan2(Q,I);
        
        %angleTheo = unwrap(angle(trace11));% - median(unwrap(angle(trace11)));
        angleTheo = unwrap(angle(trace11));
        magTheo = abs(trace11);
        
        %plot(zPos/1e3,Phase-Phase_ref, zPos/1e3, angleTheo-angleTheo_ref);
         figure(1), plot(zPos/1e3,Phase-Phase_ref), hold on;
         figure(2), plot(zPos/1e3, angleTheo-angleTheo_ref), hold on;
        drawnow;

        seno1 = [seno1 (angleTheo(1400051)-angleTheo_ref(1400051))];
        seno2 = [seno2 (Phase(1400051)-Phase_ref(1400051))];
    end
    %figure(1), plot(zPos/1e3, Magnitude-Magnitude_ref, zPos/1e3, magTheo-magTheo_ref);

end
% 
% figure, plot(seno1);
% figure, plot(seno2);
% %
% figure,plot(zPos/1e3,abs(PtracePD))




