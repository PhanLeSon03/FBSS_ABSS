clear
clc
close all
c = 340.0;
fs = 16000;
fl = 500;                    % lower cutoff frequency
fu = 6000;                   % higher cutoff frequency
fTest = 4000;
N = 512;                    % number of samples in one frame
N_RIR = 4800;
fStep = fs/N;

klow = round(fl/fStep);     % low index
kup = round(fu/fStep);      % high index
f = fStep*(0:N/2);       % frequencies used to compute W
nf = length(f);

lam = c/fu;
dmics = 0.02; 

mu = 0.1;                    % parameter for beam-forming                     
null = 1;                    % parameter for beam-forming

M = 7;                      % number of microphones
L = 511;
alpha1 = 0.01;
alpha2 = 0.01;
xMics = (-(M-1)/2:(M-1)/2)*dmics;

% Room Impulse Response
xRoom = 1.5 + xMics;
r = [xRoom' 2*ones(M,1) 1*ones(M,1)];              % Receiver position [x y z] (m)


s1 = [1.5+0.5*cosd(90)  2+0.5*sind(90)   1];              % Source position [x y z] (m)
s2 = [1.5+0.5*cosd(180) 2+0.5*sind(180)  1];              % Source position [x y z] (m)
s3 = [1.5+0.5*cosd(0)   2+0.5*sind(0)    1];              % Source position [x y z] (m), SOI


R = [3.5 6 3];                % Room dimensions [x y z] (m)
beta = 0.3;                 % Reverberation time (s)
mtype = 'omnidirectional';  % Type of microphone
order = 2;                  % Reflection order
dim = 3;                    % Room dimension
orientation = [0 0];            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

% diffuse noise
diffuse_noise =1;



[H1,beta1_hat] = rir_generator(c, fs, r, s1, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s1);%
[H2,beta2_hat] = rir_generator(c, fs, r, s2, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s2);% 
[H3,beta3_hat] = rir_generator(c, fs, r, s3, R, beta, N_RIR, mtype, order, dim, orientation, hp_filter); %rir(fs, r(iMic,:), 12, beta, R, s3);%


transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);

t = 0:1/fs:0.5;


Source1=chirp(t,fl,0.5,fu);

Source2 =  randn(1,8001);%0.5*chirp(t,2000,0.5,2000);%


Source3 = chirp(t,fTest-2000,0.5,fTest-2000);
Source3(1:end)=0;
% Source3(6501:end)=0;
% Source3(6001:6500) = 5*ones(500,1);
Source3(6001) = 100;Source3(6500) = 100;

Source1 = bandpass(Source1,[fl fu],fs);
Source2 = bandpass(Source2,[fl fu],fs);
Source3 = bandpass(Source3,[fl fu],fs);

idx = 6000:6501;
sigma_s = std(Source3(idx));

Source1 = sigma_s*Source1/std(Source1);
Source2 = sigma_s*Source2/std(Source2);
% Source3 = Source3/std(Source3);

% Create an incident wave arriving at the array. Add gaussian noise to the wave.
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',N);

signal1 = zeros(length(Source1),M);
signal2 = zeros(length(Source2),M);
signal3 = zeros(length(Source3),M);

% Room Impulse Response
for iMic=1:M
    signal1(:,iMic) = filter (H1(iMic,:),1,Source1);
    signal2(:,iMic) = filter (H2(iMic,:),1,Source2);
    signal3(:,iMic) = filter (H3(iMic,:),1,Source3); % SOI
end

signal = signal3  +  signal2 + signal1;

if diffuse_noise==0
   SNR = -10;
   noise = randn(size(signal));
   noise = noise/std(noise);
   noise = (10^(-SNR/20))*noise;
   recsignal = signal+noise;

else
    params.fs  = fs;
    params.c  = c;
    params.N_phi = 360;
    signal_diff = sinf_1D(xMics,length(Source1),params)';
    signal_diff = bandpass(signal_diff,[fl fu],fs);
    
    SNR1 = -10;
    recsignal = signal + (10^(-SNR1/20))*signal_diff/std(signal_diff);
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
legend('Received signal','SOI+Reverberation','GSC','Location','Best')
axis tight
grid on
set(gcf,'color','w');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);

pos(2) = pos(2)+0.4;
figure('numbertitle','on','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 6000:6800;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)
plot(t(idx)*1000,ygsc(idx),'r')
xlabel('Time (ms)')
axis tight
legend('Received signal','SOI+Reverberation','GSC','Location','Best')

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
DF_SDSS_C = zeros(nf,1);
bpdB_SD = zeros(nf,Nphi);
bpdB_SDSS = zeros(nf,Nphi);
bpdB_SDSS_C = zeros(nf,Nphi);
W_SDSS = zeros(M,nf);
W_SDSS_C =  zeros(M,nf);

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
nguy = 0.1;

Wup = zeros(M,N/2+1)/M;

Wup(:,klow+1:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow:kup),mu,null);


kx1 = round(300/fStep);%-klow+1; % index threshold for spatial aliasing
kx2 = klow+1; % index threshold for spatial aliasing
for k = klow+1:kup+1
    
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

% if Calibration step enables
for k = klow:kup+1
    [R,theta,p] = array_pattern_fft(xMics',Wup,f(k),k); 

    % for side-lobe beam
    if k < kx2
        phi_desired_s = 180;
        phi_zero_s = phi_zero;
    else
      [~,locs]=findpeaks(1./R(1:length(R)/2));
       phi_desired_s = 180;
       phi_zero_s = phi_zero;   %(locs-1)
    end
    fd_s = [1 zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_lo_C = bf_coefs(xMics',theta3_s,phi3_s,fd_s,f,mu,null);
    
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

Wlow = W_SDSS;


Hup = fftshift(irfft(Wup,[],2),2);
Hlow = fftshift(irfft(Wlow,[],2),2);
Hlow_C = fftshift(irfft(W_SDSS_C,[],2),2);

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
Pk0Old = 0.0;
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
    
   out_ABSS(N/2-(L-1)/2+iLoop) =  y_out;
   out_FBSS(N/2-(L-1)/2+iLoop) =   yFill_up((L-1)/2+1) - yFill_low((L-1)/2+1);
   
   % update filter coefficient
   Pk0 = sum(yFill_up.*yFill_up);
   Pk = sum(yFill_low.*yFill_low);
   Pk0 = 0.9*Pk0 + 0.1*Pk0Old;
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
%plot(t*1000,noise(:,4))
xlabel('Time (ms)')

legend('Received signal','SOI+Reverberation','ABSS','Location','Best')

grid on
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
    
    
pos(2) = pos(2) +0.3;
figure('numbertitle','on','name','Wave signals','Units','normal',...
       'Position',pos);
idx = 6000:6800;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)/20
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')

legend('Received signal','SOI+Reverberation','ABSS','Location','Best')
axis tight
grid on
set(gcf,'color','w');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);


pos(2) = pos(2) - 0.3;
figure('numbertitle','on','name','Wave signales','Units','normal',...
       'Position',pos);

idx = 6000:6800;

plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2) % Source3(idx)/20 
hold on
plot(t(idx)*1000,ygsc(idx),'b') % ygsc(idx)
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')

legend('SOI+Reverberation','GSC','ABSS','Location','Best')

set(gcf,'color','w'); 
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);
axis tight
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

% %%
% pos = [0.045 0.45 0.3 0.35];
% figure('numbertitle','on','name','Upper path''s beam pattern',...
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
% figure('numbertitle','on','name','Lower path''s beam pattern',...
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
% 
% pos(1) = pos(1)+0.32;
% figure('numbertitle','on','name','Lower path''s beam pattern',...
%        'Units','normal','Position',pos);
% 
% imagesc(ph(1:180),fStep*(klow:kup), bpdB_SDSS_C((klow+1:kup+1),1:180));
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























































