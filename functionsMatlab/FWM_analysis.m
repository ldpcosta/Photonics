clc, clear all, close all, warning off
   
%% Define Physical constants
c0 = 2.99792458e8;
n0 = 1.452;


%% Define global variables.
fc = 193e12;        % ! Perguntar porque
FWHMin = 100e-9;    % ! Perguntar o que e
lambda0 = fc/c0;
omegac = 2*pi*fc;


%% Parse files
path = 'C:\Users\W10\Desktop\09052017';
fileExtension = 'txt';

cd(path)
[filenames_, type] = parseFilenames(fileExtension, path);


%% File data
types = {'Het';'FWM'};
nSegments = 512;

powerIn_dB = 17:-2:5;
nSim = length(powerIn_dB);
typesOfMeasurements = 3;

filenames = reshape(filenames_, nSim, typesOfMeasurements);


%% Define values and scales for time and frequency.
nSamp = nSamplesPow2(filenames{1});
windowF = 40e9;                       % sampling frequency
fMax   = fc + windowF/2;                % frequencia maxima do espaco de frequencia
fMin = fc - windowF/2;                  % frequencia minima do espaco de frequencia
freq   = linspace(fMin, fMax, nSamp+1); % cria o eixo dos x (basicamente valores de fmin ate fmax equidistantes, de modo a termos nSamp+1 valores)do espaco de frequencias
interv = freq(end) - freq(1);           % da-nos o verdadeiro span de frequencia
dFreq  = freq(2)-freq(1);               % intervalo entre frequencias
freq   = freq(1:end-1);                 % Optical frequency vector ( x axis frequency)           
lambda = c0./freq;                      % wavelength scale
lambda = linspace(lambda(end),lambda(1),nSamp); % Optical wavelength vector (x axis, wavelength)
fShift = freq - fc;                             % Shifted frequency vector (shifted so frequency is from -f to +f)
dTime  = 1/interv;
time   = fix(0:nSamp-1)*dTime;                  % Time vector
timec  = fix(-nSamp/2:nSamp/2-1)*dTime;         % Time vector center at origin
omeg   = 2*pi*freq;



%% Define (or find) center frequencies and signal time-zones
cFreq = [5.7e9 2*5.7e9];

nIn  = [2100 2100; 1800 1800]; 
nInS = nIn;

nFin = [5700 5700; 5400 5400]; 
nFinS = nFin;

for i = 1:typesOfMeasurements
    filelist = filenames(:,i);

    for filename = 1:length(filelist)
        if filelist{filename}(1) == 'H', nBeats = 1;
        else, nBeats = 2;
        end
        
        for nr = 1:nBeats
            pumpPower_dBm = 20-getPowerFromFilename(filelist{filename});
            
          %open file, alloc to dataArray
            file = filelist{filename};
            file_path = strcat(path,file);
            delimiter = '\t'; formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            fileID = fopen(file,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            fclose(fileID);

            Pmeas = divideInSegments(dataArray, 512);
            clearvars delimiter formatSpec fileID dataArray ans;

            orderSuperGaussian = 4;
            spect_BPF1 = generateSupergaussianFilter(cFreq, fShift, orderSuperGaussian, nr);
            
            title(strcat(filelist{filename}(1:3),' | P_{pump}:', num2str(pumpPower_dBm), ' | \omega:', num2str(nr)));

            Pmeas_f = zeroPadAndFFT(Pmeas, nSamp);
            Eout_filt_f = Pmeas_f .* repmat(spect_BPF1,nSegments,1); %applies the filter to all segments in the spectral domain
            Eout_filt = fftshift( ifft( ifftshift( Eout_filt_f.'))).';
            
            
            fcFWM = cFreq(nr);
            Eout_dc = Eout_filt .* repmat(exp(-2j*pi*fcFWM*time),nSegments,1);  % ! PORQUE ?

            Eout_dc_f = fftshift(fft(ifftshift(Eout_dc.'))).';
            plot(Eout_dc_f(7,:));
            
            
        end
    end
end

% for nType = 1:numberOfMeasurements
%     nHarm = find(cellfun('length',regexp(types,filenames{1,nType}(1:3)) == 1)); % 1 for FWM, 2 for Het
%     
%     for nsim = 1:nSim
%         

%         
%         power = getPowerFromFilename(file)+20;
%    
%         Pmeas = divideInSegments(dataArray, nSegments);
%         
%         clearvars delimiter formatSpec fileID dataArray ans;
%         
%         % Generate time and frequency vectors
%         nSamp  = 2^(nextpow2(size(Pmeas,2))); % diz o numero de dados tirado por amostra, arredondado (para cima) a potencia de 2 mais proxima
%         windowF = 40e9;                       % sampling frequency
%         fMax   = fc + windowF/2;                % frequencia maxima do espaco de frequencia
%         fMin = fc - windowF/2;                  % frequencia minima do espaco de frequencia
%         freq   = linspace(fMin, fMax, nSamp+1); % cria o eixo dos x (basicamente valores de fmin ate fmax equidistantes, de modo a termos nSamp+1 valores)do espaco de frequencias
%         interv = freq(end) - freq(1);           % da-nos o verdadeiro span de frequencia
%         dFreq  = freq(2)-freq(1);               % intervalo entre frequencias
%         freq   = freq(1:end-1);                 % Optical frequency vector ( x axis frequency)           
%         lambda = c0./freq;                      % wavelength scale
%         lambda = linspace(lambda(end),lambda(1),nSamp); % Optical wavelength vector (x axis, wavelength)
%         fShift = freq - fc;                             % Shifted frequency vector (shifted so frequency is from -f to +f)
%         dTime  = 1/interv;
%         time   = fix(0:nSamp-1)*dTime;                  % Time vector
%         timec  = fix(-nSamp/2:nSamp/2-1)*dTime;         % Time vector center at origin
%         omeg   = 2*pi*freq;
%         %
%         
%         Pmeas_f = zeroPadAndFFT(Pmeas, nSamp);
%         figure,plot(fShift,abs(Pmeas_f)) %% PLOT
%         
%         orderSuperGaussian = 4;
%         spect_BPF1 = generateFilter(cFreq, fShift, orderSuperGaussian);
%         hold on,plot(fShift,10*abs(spect_BPF1)) %% PLOT
%         
%         
%         Eout_filt_f = Pmeas_f.*repmat(spect_BPF1,nSegments,1); %applies the filter to all segments in the spectral domain
%         Eout_filt = fftshift(ifft(ifftshift(Eout_filt_f.'))).';
%         Eout_dc = Eout_filt.*repmat(exp(-2j*pi*fcFWM*time),nSegments,1);
%         Eout_dc_f = fftshift(fft(ifftshift(Eout_dc.'))).';
%         figure,plot(abs(Eout_dc(1,:)))
%         
%         vkh = nHarm;    
%         auxSNR = [];
%         for rr = 1:nSegments
%         sig2w = abs(Eout_dc(rr,:))/50; 
%             levS2w = mean(sig2w(nIn(vkh,nr):nFin(vkh,nr))).^2;  
%             levN2w = std(sig2w(nInS(vkh,nr):nFinS(vkh,nr))).^2;    
%             auxSNR(rr) = 10*log10((levS2w-0*levN2w)/levN2w);
%             if vkh == 1 || vkh == 2
%                 auxComp(vkh,:,rr,ns,nr) = sig2w;
%             end
%         end
%         SNR_fin(ns,vkh,nr) = mean(auxSNR);
%     end
% end
