clear
clc
close all
c = 340.0;
fs = 16000;
fl = 1000;                    % lower cutoff frequency
fu = 6000;                   % higher cutoff frequency
fTest = 4000;
lam = c/fu;
dmics = 0.02; 
N = 512;                    % number of samples in one frame
fStep = fs/N;
M = 7;                      % number of microphones
L = 101;
xMics = (-(M-1)/2:(M-1)/2)*dmics;

% Room Impulse Response
xRoom = 1 + (0:M-1)'*dmics;
r = [xRoom ones(M,1) ones(M,1)];              % Receiver position [x y z] (m)

s1 = [1+(M-1)*dmics/2+2*cosd(60)  1+2*sind(60)   1];              % Source position [x y z] (m)
s2 = [1+(M-1)*dmics/2+2*cosd(150) 1+2*sind(150)  1];              % Source position [x y z] (m)
s3 = [1+(M-1)*dmics/2+2*cosd(0)   1+2*sind(0)    1];              % Source position [x y z] (m)

R = [5 4 4];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reverberation time (s)

h1_1 = rir_generator(c, fs, r(1,:), s1, R, beta, N);
h2_1 = rir_generator(c, fs, r(2,:), s1, R, beta, N);
h3_1 = rir_generator(c, fs, r(3,:), s1, R, beta, N);
h4_1 = rir_generator(c, fs, r(4,:), s1, R, beta, N);
h5_1 = rir_generator(c, fs, r(5,:), s1, R, beta, N);
h6_1 = rir_generator(c, fs, r(6,:), s1, R, beta, N);
h7_1 = rir_generator(c, fs, r(7,:), s1, R, beta, N);

h1_2 = rir_generator(c, fs, r(1,:), s2, R, beta, N);
h2_2 = rir_generator(c, fs, r(2,:), s2, R, beta, N);
h3_2 = rir_generator(c, fs, r(3,:), s2, R, beta, N);
h4_2 = rir_generator(c, fs, r(4,:), s2, R, beta, N);
h5_2 = rir_generator(c, fs, r(5,:), s2, R, beta, N);
h6_2 = rir_generator(c, fs, r(6,:), s2, R, beta, N);
h7_2 = rir_generator(c, fs, r(7,:), s2, R, beta, N);

h1_3 = rir_generator(c, fs, r(1,:), s3, R, beta, N);
h2_3 = rir_generator(c, fs, r(2,:), s3, R, beta, N);
h3_3 = rir_generator(c, fs, r(3,:), s3, R, beta, N);
h4_3 = rir_generator(c, fs, r(4,:), s3, R, beta, N);
h5_3 = rir_generator(c, fs, r(5,:), s3, R, beta, N);
h6_3 = rir_generator(c, fs, r(6,:), s3, R, beta, N);
h7_3 = rir_generator(c, fs, r(7,:), s3, R, beta, N);

transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);
alpha1 = 0.01;
alpha2 = 0.01;
t = 0:1/fs:0.5;

incidentAngle1 = [20 ;0]; %10° azimuth and 90° elevation, 20 = 90 -70:70
incidentAngle2 = [-60 ;0]; % -60 = 90 -150:150
incidentAngle3 = [90 ;0]; % source of interset 90 = 90 -0: 0

% Simulate a chirp signal with a 500 Hz bandwidth.
Source1 = chirp(t,5000,0.5,5000);%0.5*randn(1,8001);%
% 
% Nx=(N/2)-10; % the number of tones
% mindex=0:Nx;
% phi=pi*mindex.*mindex/Nx;
% k0 = 5;
% omega=2*pi*(k0*ones(Nx+1,1)+mindex')*fs/N;
% Source1=sum(sin(omega*t+phi'*ones(1,length(t))),1)/Nx;
% pspectrum(Source1,fs,'spectrogram','TimeResolution',0.1, ...
%         'OverlapPercent',99,'Leakage',0.85)

Source2 =  randn(1,8001);%0.5*chirp(t,2000,0.5,2000);%
Source3 = zeros(size(Source1));
%Source3(6001:6500) = 5*ones(500,1);
Source3(6001) = 5;Source3(6500) = 5;
Source1 = bandpass(Source1,[fl fu],fs);
Source2 = bandpass(Source2,[fl fu],fs);
Source3 = bandpass(Source3,[fl fu],fs);


% Create an incident wave arriving at the array. Add gaussian noise to the wave.
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',N);


% free-field
%signal1 = collector(Source1.' ,incidentAngle1);
%signal2 = collector(Source2.' ,incidentAngle2);
%signal3 = collector(Source3.' ,incidentAngle3);

signal1 = zeros(length(Source1),M);
signal2 = zeros(length(Source2),M);
signal3 = zeros(length(Source3),M);

% Room Impulse Response
signal1(:,1) = filtfilt(h1_1,1,Source1);
signal1(:,2) = filtfilt(h2_1,1,Source1);
signal1(:,3) = filtfilt(h3_1,1,Source1);
signal1(:,4) = filtfilt(h4_1,1,Source1);
signal1(:,5) = filtfilt(h5_1,1,Source1);
signal1(:,6) = filtfilt(h6_1,1,Source1);
signal1(:,7) = filtfilt(h7_1,1,Source1);

signal2(:,1) = filtfilt(h1_2,1,Source2);
signal2(:,2) = filtfilt(h2_2,1,Source2);
signal2(:,3) = filtfilt(h3_2,1,Source2);
signal2(:,4) = filtfilt(h4_2,1,Source2);
signal2(:,5) = filtfilt(h5_2,1,Source2);
signal2(:,6) = filtfilt(h6_2,1,Source2);
signal2(:,7) = filtfilt(h7_2,1,Source2);

signal3(:,1) = filtfilt(h1_3,1,Source3);
signal3(:,2) = filtfilt(h2_3,1,Source3);
signal3(:,3) = filtfilt(h3_3,1,Source3);
signal3(:,4) = filtfilt(h4_3,1,Source3);
signal3(:,5) = filtfilt(h5_3,1,Source3);
signal3(:,6) = filtfilt(h6_3,1,Source3);
signal3(:,7) = filtfilt(h7_3,1,Source3);

signal = signal1 + signal2 + signal3;

SNR = 40;
noise = (10^(-SNR/20))*randn(size(signal))*mean(mean(signal));
recsignal = signal + noise;

% pspectrum(recsignal(:,4),fs,'spectrogram','TimeResolution',0.1, ...
%       'OverlapPercent',99,'Leakage',0.85)


% MATLAB toolbox for GSC---------------------------------------------------
gscbeamformer = phased.GSCBeamformer('SensorArray',array, ...
    'PropagationSpeed',c,'SampleRate',fs,'DirectionSource','Input port', ...
    'FilterLength',L,'LMSStepSize',alpha1);

ygsc = gscbeamformer(recsignal,[90;0]);


% plotting
pos = [0.0 0.0 0.45 0.45];
figure('numbertitle','off','name','Full wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,signal3(:,4),'g','Linewidth',2)  % Source3
plot(t*1000,ygsc,'r')
xlabel('Time (ms)')
ylabel('Amplitude')
legend('Received signal','SOI','GSC','Location','Best')
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);

pos(2) = pos(2)+0.4;
figure('numbertitle','off','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 5800:7200;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)
plot(t(idx)*1000,ygsc(idx),'r')
xlabel('Time (ms)')
legend('Received signal','SOI','GSC','Location','Best')
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
% % %--------------------------------------------------------------------------

% Synchronize signal toward look direction in FFT domain


%% Propose method

klow = round(fl/fStep);     % low index
kup = round(fu/fStep);      % high index
f = fStep*(0:N/2);       % frequencies used to compute W
nf = length(f);
phi = pi/180*(0:1:360);
Nphi = length(phi);
WNG_SD = zeros(nf,1);
WNG_SDSS = zeros(nf,1);
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
mu = 0.1;                    % parameter for beam-forming                     
null = 1;                    % parameter for beam-forming
% find optimum frequency-domain weight vector      
Gamma = zeros(length(xMics),length(xMics));
for i=1:length(xMics)
   Gamma(i,:) =  abs(xMics(i)-xMics);   
end
nguy = 0.1;

Wup = ones(M,N/2+1)/M;
for k = 2:klow
    beta = 2*pi*(k-1)*fStep/c;
    d0 = exp(1j*beta*xMics'*cosd(phi_desired));        % steering vector at peak of side-lobe
    Shi = (sin(beta*Gamma)./(beta*Gamma));
    Shi(logical(eye(size(Shi)))) = 1;       
    Wup(:,k) =(Shi*(1-nguy) + nguy*eye(length(xMics)))^-1*d0 / (d0'*(Shi*(1-nguy) + nguy*eye(length(xMics)))^-1*d0);
    Wup(:,k) = d0/M;
end
Wup(:,klow+1:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow:kup),mu,null);
for k = kup+2:N/2+1
    beta = 2*pi*(k-1)*fStep/c;
    df = exp(1j*beta*xMics'*cosd(phi_desired));        % steering vector at peak of side-lobe
    Wup(:,k) = df/M;
end

kx1 = round(300/fStep);%-klow+1; % index threshold for spatial aliasing
kx2 = klow+1; % index threshold for spatial aliasing
for k = klow:kup+1
    [R,theta,p] = array_pattern_fft(xMics',Wup,f(k),k); 

    % for side-lobe beam
    if k < kx2
        phi_desired_s = 180;
        phi_zero_s = [0 phi_zero];
    else
      [~,locs]=findpeaks(1./R(1:length(R)/2));
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
    
    W_SDSS(:,k) =  scale*W_lo(:,k);
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



%Wup(:,klow+1:kup+1) =  W;

Wlow = W_SDSS;
Wlow(:,1:klow) = 0; 
Wlow(:,kup+2:end) = 0; 

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
h_u((L-1)/2+1) = 1;
h_l = zeros(L,1);
for iLoop = 1:length(recsignal(:,1))- N + 1
   y_beam1 = sum(sum(Hup.*recsignal(iLoop:iLoop+N-1,:)',2),1);
   yFill_up = circshift(yFill_up,1);
   yFill_up(1) = y_beam1;
   y_up = sum(h_u.*yFill_up);
   
   y_beam2 = sum(sum(Hlow.*recsignal(iLoop:iLoop+N-1,:)',2),1);
   yFill_low = circshift(yFill_low,1);
   yFill_low(1) = y_beam2;
   y_low = sum(h_l.*yFill_low);
   
   y_out = y_up - y_low;
    
   out_ABSS(N/2-(L-1)/2+iLoop) =  y_out;
   out_FBSS(N/2-(L-1)/2+iLoop) =   yFill_up((L-1)/2+1) - yFill_low((L-1)/2+1);
   
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
   
end

pos = [0.5 0.0 0.45 0.45];

figure('numbertitle','off','name','Wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,signal3(:,4),'g','Linewidth',2) % Source3
plot(t*1000,out_ABSS,'r')
xlabel('Time (ms)')
legend('Received signal','SOI','ABSS','Location','Best')
grid on
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
    
    
pos(2) = pos(2) +0.3;
figure('numbertitle','off','name','Wave signals','Units','normal',...
       'Position',pos);
idx = 5800:7200;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
legend('Received signal','SOI','ABSS','Location','Best')
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);


pos(2) = pos(2) - 0.3;
figure('numbertitle','off','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 5800:7000;
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx) 
hold on
plot(t(idx)*1000,ygsc(idx),'b')
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
legend('SOI','GSC','ABSS','Location','Best')
    set(gcf,'color','w'); 
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);

grid on
%%
e1 = mean((ygsc(idx) - signal3(idx,4)).^2)
e2 = mean((out_ABSS(idx) - signal3(idx,4)).^2)
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
% pos = [0.045 0.45 0.3 0.35];
% figure('numbertitle','off','name','Upper path''s beam pattern',...
%        'Units','normal','Position',pos);
% ph = 180/pi*p; 
%   
% imagesc(ph(1:180),fStep*(klow:kup),bpdB_SD((klow+1:kup+1),1:180));
% axis tight
% %set(gca,'XTick',[0 45 90 135 180]);
% %view([25,50]);
% xlabel('Incident angle \phi in °');
% ylabel('f in Hz');
% zlabel('Magnitude');
% title('Array response at the upper path (preserving the main-lobe)');
% set(gca, 'xdir', 'reverse');
% set(gca, 'ydir', 'normal');
% %zlim([-dBmax 0])
% set(gcf,'color','w');
% grid on
% colorbar 
% %colormap jet
% 
% pos(1) = pos(1)+0.32;
% figure('numbertitle','off','name','Lower path''s beam pattern',...
%        'Units','normal','Position',pos);
%   
% imagesc(ph(1:180),fStep*(klow:kup), bpdB_SDSS((klow+1:kup+1),1:180));
% axis tight
% %set(gca,'XTick',[0 45 90 135 180]);
% %view([25,50]);
% xlabel('Incident angle \phi in °');
% ylabel('f in Hz');
% zlabel('Magnitude');
% title('Array response at the lower path (suppressing the main-lobe)');
% set(gca, 'xdir', 'reverse');
% set(gca, 'ydir', 'normal');
% %zlim([-dBmax 0])
% set(gcf,'color','w');
% grid on
% colorbar 
% %colormap jet























































