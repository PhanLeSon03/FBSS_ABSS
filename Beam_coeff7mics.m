clear
clc
close all
c = 340.0;
fs = 48000;
fl = 1200;                    % lower cutoff frequency
fu = 8000;                    % higher cutoff frequency
fTest = 4000;
N = 512;                      % number of samples in one frame
N_RIR = 4800;
fStep = fs/N;

klow = round(fl/fStep);       % low index
kup = round(fu/fStep);        % high index
f = fStep*(0:N/2);            % frequencies used to compute W
nf = length(f);

lam = c/fu;
dmics = 0.02; 

mu = 0.001;                     % parameter for beam-forming                     
null = 1;                     % parameter for beam-forming

M = 7;                        % number of microphones
L = 256;
window =hanning(L);
window(L/2) = 0.9999;
window(L/2+1) = 0.9999;
alpha1 = 0.01;
alpha2 = 0.01;
xMics = (-(M-1)/2:(M-1)/2)*dmics;

% Room Impulse Response
xRoom = 1.5 + xMics;
r = [xRoom' 2*ones(M,1) 1*ones(M,1)];                     % Receiver position [x y z] (m)


s1 = [1.5+0.5*cosd(70)  2+0.5*sind(70)   1];              % Source position [x y z] (m)
s2 = [1.5+0.5*cosd(150) 2+0.5*sind(150)  1];              % Source position [x y z] (m)
s3 = [1.5+0.5*cosd(0)   2+0.5*sind(0)    1];              % Source position [x y z] (m), SOI

incidentAngle1 = [20 ;0];  %10° azimuth and 90° elevation, 20 = 90 -70:70
incidentAngle2 = [-60 ;0]; % -60 = 90 -150:150
incidentAngle3 = [90 ;0];  % source of interset 90 = 90 -0: 0

R = [3.5 6 3];              % Room dimensions [x y z] (m)
beta = 0.3;                 % Reverberation time (s)
mtype = 'omnidirectional';  % Type of microphone
order = 2;                  % Reflection order
dim = 3;                    % Room dimension
orientation = [0 0];        % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

% free-field
free_field=0;
% diffuse noise
diffuse_noise =0;



[H1,beta1_hat] = rir_generator(c, fs, r, s1, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s1);%
[H2,beta2_hat] = rir_generator(c, fs, r, s2, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s2);% 
[H3,beta3_hat] = rir_generator(c, fs, r, s3, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s3);%


transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);

t = 0:1/fs:1;


% Simulate a chirp signal with a 500 Hz bandwidth.
%Source1 = chirp(t,5000,0.5,5000);%0.5*randn(1,8001);%

Nx=(N/2)-2; % the number of tones
mindex=0:Nx;
phi=pi*mindex.*mindex/Nx;
k0 = 5;
omega=2*pi*(k0*ones(Nx+1,1)+mindex')*fs/N;
Source1=sum(sin(omega*t+phi'*ones(1,length(t))),1)/Nx;
Source1 = 0.5*Source1/max(Source1);


Source2 =  0.5*chirp(t,2000,1,2000);%rand(1,16001);%
Source3 = zeros(size(Source1));
%Source3(6001:6500) = sin(2*pi*3000/fs*(0:500-1));


Source3(6001) = 1;Source3(6500) = 1;


Source1 = bandpass(Source1,[fl fu],fs);
Source2 = bandpass(Source2,[fl fu],fs);
Source3 = bandpass(Source3,[fl fu],fs);

% Create an incident wave arriving at the array. Add gaussian noise to the wave.
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',N);


Source_C = 10*chirp(t,fl,1,fu,'linear');%randn(1,16001);%
pspectrum(Source_C,fs,'spectrogram','TimeResolution',0.1, ...
        'OverlapPercent',99,'Leakage',0.85)
if free_field==1
    signal_C = collector(Source_C.' ,incidentAngle3);
        %collector(Source_C.' ,incidentAngle1)+...
        %collector(Source_C.' ,incidentAngle2); 
else
    signal_C = zeros(length(Source_C),M);
    for iMic=1:M
        signal_C(:,iMic) = filter(H3(iMic,:),1,Source_C);
         %                  filter(H1(iMic,:),1,Source_C)+...
         %                  filter(H2(iMic,:),1,Source_C);
    end
end

nFrame = floor(length(Source_C)/N);

Rx = zeros(M,M,N);
for iFrame = 1:nFrame
   signal_F = fft(signal_C((iFrame-1)*N+1:iFrame*N,:),[],1);
   
   for k = 1:N
       normVal = abs(signal_F(k,:)')*abs(signal_F(k,:));
       Rx(:,:,k) = Rx(:,:,k) +  signal_F(k,:)'*signal_F(k,:);
   end
end
Rx =  Rx/nFrame;

% free-field
if free_field==1
    signal1 = collector(Source1.' ,incidentAngle1);
    signal2 = collector(Source2.' ,incidentAngle2);
    signal3 = collector(Source3.' ,incidentAngle3); % SOI
else
    signal1 = zeros(length(Source1),M);
    signal2 = zeros(length(Source2),M);
    signal3 = zeros(length(Source3),M);

    % Room Impulse Response
    for iMic=1:M
        signal1(:,iMic) = filter (H1(iMic,:),1,Source1);
        signal2(:,iMic) = filter (H2(iMic,:),1,Source2);
        signal3(:,iMic) = filter (H3(iMic,:),1,Source3); % SOI
    end
    %signal1 = collector(Source1.' ,incidentAngle1);
    %signal2 = collector(Source2.' ,incidentAngle2);
end
signal = signal3  +  signal2 + signal1;

if diffuse_noise==0
   SNR = 10;
   noise = randn(size(signal));
   noise = noise/max(max(abs(noise)));
   noise = (10^(-SNR/20))*noise*max(max(signal3));
   recsignal = signal+noise;

else
    params.fs  = fs;
    params.c  = c;
    params.N_phi = 360;
    signal_diff = sinf_1D(xMics,length(Source1),params)';
    signal_diff = bandpass(signal_diff,[fl fu],fs)/max(max(abs(signal_diff)));
    SNR1 = 10;
    recsignal = signal + (10^(-SNR1/20))*max(max(signal3))*signal_diff;
end


% pspectrum(recsignal(:,4),fs,'spectrogram','TimeResolution',0.1, ...
%       'OverlapPercent',99,'Leakage',0.85)


% MATLAB toolbox for GSC---------------------------------------------------
gscbeamformer = phased.GSCBeamformer('SensorArray',array, ...
    'PropagationSpeed',c,'SampleRate',fs,'DirectionSource','Input port', ...
    'FilterLength',L,'LMSStepSize',alpha1);

ygsc = gscbeamformer(recsignal,[90;0]); %% 90,0


% plotting
pos = [0.0 0.0 0.45 0.45];
figure('numbertitle','on','name','Full wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,signal3(:,4),'g','Linewidth',2)  % Source3
plot(t*1000,ygsc,'r')
xlabel('Time (ms)')
%ylabel('Amplitude')
if free_field==1
    legend('Received signal','SOI','GSC','Location','Best')
else
    legend('Received signal','SOI+Reverberation','GSC','Location','Best')
end
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);

pos(2) = pos(2)+0.4;
figure('numbertitle','on','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 5800:7200;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)
plot(t(idx)*1000,ygsc(idx),'r')
xlabel('Time (ms)')
if free_field==1
    legend('Received signal','SOI','GSC','Location','Best')
else
    legend('Received signal','SOI+Reverberation','GSC','Location','Best')
end
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
% % %--------------------------------------------------------------------------

% Synchronize signal toward look direction in FFT domain


%% Propose method

phi = pi/180*(0:1:360);
Nphi = length(phi);
WNG_SD = zeros(nf,1);
WNG_SDSS = zeros(nf,1);
WNG_SDSS_C = zeros(nf,1);

DF_SD = zeros(nf,1);
DF_SDSS = zeros(nf,1);

bpdB_SD = zeros(nf,Nphi);
bpdB_SDSS = zeros(nf,Nphi);

W_SDSS = zeros(M,nf);


% for main-lobe beam
phi_desired = 0;
phi_zero = [70 150];

fd = [1 zeros(size(phi_zero))];  % resonse in desired directions

phi3 = [phi_desired phi_zero];
theta3 = 90*ones(1,length(fd));

% find optimum frequency-domain weight vector      
Gamma = zeros(length(xMics),length(xMics));
for i=1:length(xMics)
   Gamma(i,:) =  abs(xMics(i)-xMics);   
end


Wup = zeros(M,N/2+1)/M;

% Delay and Sum as remaining frequencies
for k=1:N/2+1
    beta = 2*pi*f(k)/c;          % wave number
    d = exp(-1j*beta*xMics);      % steering vector
    Wup(:,k) = d'/M; 
end

Wup(:,klow+1:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow:kup),mu,null);
%Wup(:,klow+1:N/2+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow:N/2),mu,null);


kx1 = round(300/fStep);%-klow+1; % index threshold for spatial aliasing
kx2 = kup+2; % index threshold for spatial aliasing
for k = klow+1:kup+1
    
    [R,theta,p] = array_pattern_fft(xMics',Wup,f(k),k); 

    % for side-lobe beam
    if k < kx2
        phi_desired_s = 180;
        phi_zero_s = [0 phi_zero];
    else
      if k < 65 % 64 is 6 kHz
        [~,locs]=findpeaks(1./R(1:length(R)/2));
      end 
      phi_desired_s = 180;
      phi_zero_s = [0 (locs-1)];  
    end
    fd_s = [1 zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_lo = bf_coefs(xMics',theta3_s,phi3_s,fd_s,f,mu,null);
    
    % scale factor
    beta = 2*pi*f(k)/c;
    df = exp(1j*beta*xMics'*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
    scale = Wup(:,k)'*df;

    D = exp(1j*beta*xMics'*cos(p));                    % steering matrix at a frequency
    bp = Wup(:,k)'*D;
    bpdB_SD(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    W_SDSS(:,k) =  1*W_lo(:,k);
    bp = W_SDSS(:,k)'*D;
%     W_SDSS(:,k) =  W_s(:,k)/max(abs(bp));
%     bp = W_SDSS(:,k)'*D;
    bpdB_SDSS(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    ds = exp(1j*beta*xMics'*cosd(phi_desired(1)));         % steering vector at look direction
    WNG_SD(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Wup(:,k)));
    WNG_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*W_SDSS(:,k)));
    
    Rdif = spatio_spect_corr(beta,xMics');
    DF_SD(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Rdif*Wup(:,k)));
    DF_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*Rdif*W_SDSS(:,k)));
end



Wlow = W_SDSS;


Hup = fftshift(irfft(Wup,[],2),2);
Hlow = fftshift(irfft(Wlow,[],2),2);

% fvtool((Hup(4,:)))
% fvtool((Hlow(4,:)))

% myFig = figure('numbertitle','off','name','FIR coefficients','Units','normal',...
%    'Position',pos);
% subplot(2,1,1);
% imagesc(Hup)
% subplot(2,1,2);
% imagesc(Hlow)
%%

out_ABSS = zeros(length(recsignal(:,1)),1);
out_FBSS = zeros(length(recsignal(:,1)),1);
yFill_up = zeros(L,1);
yFill_low = zeros(L,1);
h_u = zeros(L,1);
h_u(L/2) = 1;
h_l = zeros(L,1);

coef = 0;
for iLoop = 1:length(recsignal(:,1))- N + 1
   y_beam1 = sum(sum(Hup.*recsignal(iLoop:iLoop+N-1,:)',2),1);
   yFill_up = circshift(yFill_up,1);
   yFill_up(1) = y_beam1;
   y_up = sum(h_u.*yFill_up);
   
   if 0 % adaptive for null beamforming
       sig_F = fft(recsignal(iLoop:iLoop+N-1,:),N,1);
       
       for k = klow:kup+1
            Rx1 = sig_F(k,:)'*sig_F(k,:);
            fd_s = [1 zeros(size(phi_zero))];  % resonse in desired directions
            phi3_s = [phi_desired_s, phi_zero];
            theta3_s =  90*ones(1,length(phi3_s));
            W_lo_C = bf_coefs_R(xMics',theta3_s,phi3_s,fd_s,f,mu,0,Rx1);

            % scale factor
            beta = 2*pi*f(k)/c;
            df = exp(1j*beta*xMics'*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
            scale = Wup(:,k)'*df;

            W_SDSS_C(:,k) =  scale*W_lo_C(:,k);

            D = exp(1j*beta*xMics'*cos(p));                    % steering matrix at a frequency


            bp = W_SDSS_C(:,k)'*D;
            bpdB_SDSS_C(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);

            ds = exp(1j*beta*xMics'*cosd(phi_desired(1)));         % steering vector at look direction
            WNG_SDSS_C(k) = 10*log10(1/(W_SDSS_C(:,k)'*W_SDSS_C(:,k)));

            Rdif = spatio_spect_corr(beta,xMics');
            DF_SDSS_C(k) = 10*log10(1/(W_SDSS_C(:,k)'*Rdif*W_SDSS_C(:,k)));
       end
        Hlow_C = fftshift(irfft(W_SDSS_C,[],2),2);
   end
   y_beam2 = sum(sum(Hlow.*recsignal(iLoop:iLoop+N-1,:)',2),1);
   yFill_low = circshift(yFill_low,1);
   yFill_low(1) = y_beam2;
   y_low = sum(h_l.*yFill_low);
   
   y_out = y_up - y_low;
    
   out_ABSS(N/2-L/2+iLoop) =  y_out;
   out_FBSS(N/2-L/2+iLoop) =   yFill_up(L/2+1) - yFill_low(L/2+1);
   
   % update filter coefficient
   Pk0 = sum(yFill_up.*yFill_up);
   Pk = sum(yFill_low.*yFill_low);
   mu = alpha2/Pk0;
   delta = mu*y_out*yFill_low;
%    if delta > 0.1
%        delta = 0.1;
%    elseif delta < -0.1
%        delta = -0.1;
%    end
   h_l = h_l + delta;
   coef = [coef, h_l(10)];
   
end
figure(1)
plot(coef)

pos = [0.5 0.0 0.45 0.45];

figure('numbertitle','off','name','Wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,signal3(:,4),'g','Linewidth',2) % Source3 
plot(t*1000,out_ABSS,'r')
%plot(t*1000,noise(:,4))
xlabel('Time (ms)')
if free_field==1
    legend('Received signal','SOI','ABSS','Location','Best')
else
    legend('Received signal','SOI+Reverberation','ABSS','Location','Best')
end
grid on
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
    
    
pos(2) = pos(2) +0.3;
figure('numbertitle','on','name','Wave signals','Units','normal',...
       'Position',pos);
idx = 5800:7200;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)/20
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
if free_field==1
    legend('Received signal','SOI','ABSS','Location','Best')
else
    legend('Received signal','SOI+Reverberation','ABSS','Location','Best')
end
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);


pos(2) = pos(2) - 0.3;
figure('numbertitle','on','name','Wave signales','Units','normal',...
       'Position',pos);
if free_field==1
    idx = 5800:7200;
else
    idx = 5000:8000;
end
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)/20 
hold on
plot(t(idx)*1000,out_FBSS(idx),'b') % ygsc(idx)
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
if free_field==1
    legend('SOI','FBSS','ABSS','Location','Best')
else
    legend('SOI+Reverberation','FBSS','ABSS','Location','Best')
end
    set(gcf,'color','w'); 
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);

grid on
%%
e1 = mean((ygsc(idx) - signal3(idx,4)).^2);
e2 = mean((out_ABSS(idx) - signal3(idx,4)).^2);
%%
% pos = [0.055 0.04 0.4 0.38];
% figure('numbertitle','off','name','WNG',...
%        'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
% plot(fStep*(klow:kup),WNG_SD(klow+1:kup+1),'-b')
% hold on
% plot(fStep*(klow:kup),WNG_SDSS(klow+1:kup+1),'-r')
% legend('Up','AGSC1','Location','SouthEast');
% grid on;
% set(gcf,'color','w');
% xlabel('f (Hz)');
% ylabel('WNG (dB)');
% 
% pos(1) = pos(1) + 0.3;
% figure('numbertitle','off','name','DF',...
%        'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
% plot(fStep*(klow:kup),DF_SD(klow+1:kup+1),'-b')
% hold on
% plot(fStep*(klow:kup),DF_SDSS(klow+1:kup+1),'-r')
% legend('Up','AGSC1','Location','SouthEast');
% grid on;
% set(gcf,'color','w');
% xlabel('f (Hz)');
% ylabel('DF (dB)');

%%
pos = [0.045 0.45 0.3 0.35];
figure('numbertitle','on','name','Upper path''s beam pattern',...
       'Units','normal','Position',pos);
ph = 180/pi*p; 
  
imagesc(ph(1:180),fStep*(klow:kup),bpdB_SD((klow+1:kup+1),1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('f in Hz');
zlabel('Magnitude');
title('Array response at the upper path (preserving the main-lobe)');
set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'normal');
%zlim([-dBmax 0])
set(gcf,'color','w');
grid on
colorbar 
%colormap jet

pos(1) = pos(1)+0.32;
figure('numbertitle','on','name','Lower path''s beam pattern',...
       'Units','normal','Position',pos);
  
imagesc(ph(1:180),fStep*(klow:kup), bpdB_SDSS((klow+1:kup+1),1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('f in Hz');
zlabel('Magnitude');
title('Array response at the lower path (suppressing the main-lobe)');
set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'normal');
%zlim([-dBmax 0])
set(gcf,'color','w');
grid on
colorbar 
%colormap jet

%%


fvtool((Hup(4,:)),1) 
fvtool((Hlow(1,:)),1)

hMAXup = max(Hup(:))

 filename = './beam_weights/FIR.txt';
 fileID = fopen(filename,'w');
for iMic=1:M
    % convert float to 24bit 
    FIR24b = zeros(3*N,1,'uint8');
    for iEle=1:N
        intVal = int32(Hup(iMic,iEle)*(2^23));
        FIR24b((iEle-1)*3+1) = uint8(bitand(bitshift(intVal,-16),255));
        FIR24b((iEle-1)*3+2) = uint8(bitand(bitshift(intVal,-8),255));
        FIR24b((iEle-1)*3+3) = uint8(bitand(intVal,255));
    end
   
    formatSpec = '0x%x';
    fprintf(fileID,'const byte hFIRUp%d[] = {',iMic); 
    for iEle = 1:3*N
        fprintf(fileID,formatSpec,FIR24b(iEle));
        if (iEle <3*N)
            fprintf(fileID,',');
        end
    end
    fprintf(fileID,'};\n');
end


hMAXlow = max(Hlow(:))

for iMic=1:M
    % convert float to 24bit 
    FIR24b = zeros(3*N,1,'uint8');
    for iEle=1:N
        intVal = int32(Hlow(iMic,iEle)*(2^23));
        FIR24b((iEle-1)*3+1) = uint8(bitand(bitshift(intVal,-16),255));
        FIR24b((iEle-1)*3+2) = uint8(bitand(bitshift(intVal,-8),255));
        FIR24b((iEle-1)*3+3) = uint8(bitand(intVal,255));
    end
    formatSpec = '0x%x';
    fprintf(fileID,'const byte hFIRLow%d[] = {',iMic); 
    for iEle = 1:3*N
        fprintf(fileID,formatSpec,FIR24b(iEle));
        if (iEle <3*N)
            fprintf(fileID,',');
        end
    end
    fprintf(fileID,'};\n');
end


fprintf(fileID,'const byte hAFIR[] = {');
for iEle = 1:3*L
    fprintf(fileID,'0x00');
    if (iEle <3*L)
        fprintf(fileID,',');
    end
end
fprintf(fileID,'};\n');

formatSpec = '0x%x';
fprintf(fileID,'const byte window[] = {');
% convert float to 24bit 
window24b = zeros(3*L,1,'uint8');
for iEle=1:L
    intVal = int32(window(iEle)*(2^23));
    window24b((iEle-1)*3+1) = uint8(bitand(bitshift(intVal,-16),255));
    window24b((iEle-1)*3+2) = uint8(bitand(bitshift(intVal,-8),255));
    window24b((iEle-1)*3+3) = uint8(bitand(intVal,255));
end
for iEle = 1:3*L
    fprintf(fileID,formatSpec,window24b(iEle));
    if (iEle <3*N)
        fprintf(fileID,',');
    end
end
fprintf(fileID,'};\n');






