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
M = 7;                      % number of microphones
L = 511;
xMics = (-(M-1)/2:(M-1)/2)*dmics;

transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);
alpha1 = 0.01;
alpha2 = (M-1)*0.01;
t = 0:1/fs:0.5;

incidentAngle1 = [20 ;0]; %10° azimuth and 90° elevation, 20 = 90 -70
incidentAngle2 = [-60 ;0]; % -60 = 90 -150
incidentAngle3 = [90 ;0]; % source of interset

% Simulate a chirp signal with a 500 Hz bandwidth.
Source1 = chirp(t,fl,0.5,fu);%0.5*randn(1,8001);%
% 
% Nx=(N/2)-10; % the number of tones
% mindex=0:Nx;
% phi=pi*mindex.*mindex/Nx;
% k0 = 5;
% omega=2*pi*(k0*ones(Nx+1,1)+mindex')*fs/N;
% Source1=sum(sin(omega*t+phi'*ones(1,length(t))),1)/Nx;
% pspectrum(Source1,fs,'spectrogram','TimeResolution',0.1, ...
%         'OverlapPercent',99,'Leakage',0.85)


Source2 =randn(1,8001);% 0.5*chirp(t,2000,0.5,2000);%


Source3 = chirp(t,fTest-2000,0.5,fTest-2000);
Source3(1:end)=0;
% Source3(6501:end)=0;
% Source3(6001:6500) = 5*ones(500,1);
Source3(6001) = 10;Source3(6500) = 10;

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


signal1 = collector(Source1.' ,incidentAngle1);
signal2 = collector(Source2.' ,incidentAngle2);
signal3 = collector(Source3.' ,incidentAngle3);


signal = signal1 + signal2 + signal3;

SNR = -10;
noise = randn(size(signal));
noise = sigma_s*sqrt(10^(-SNR/10))*noise/std(noise);
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
figure('numbertitle','off','name','Wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,Source3,'g','Linewidth',2)
plot(t*1000,ygsc,'r')
xlabel('Time (ms)')
ylabel('Amplitude')
legend('Received signal','SOI','GSC','Location','Best')
grid on
    set(gcf,'color','w');
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
axis tight

pos(2) = pos(2)+0.4;
figure('numbertitle','off','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 5800:7200;
plot(t(idx)*1000,recsignal(idx,4))
hold on
plot(t(idx)*1000,Source3(idx),'g','Linewidth',2)
plot(t(idx)*1000,ygsc(idx),'r')
xlabel('Time (ms)')
legend('Received signal','SOI','GSC','Location','Best')
grid on
set(gcf,'color','w');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);
axis tight
% % %--------------------------------------------------------------------------

% % Synchronize signal toward look direction in FFT domain
% fStep = fs/N;
% D = zeros(M,N/2+1);
% for iBin = 2:N/2+1
%     fBin = fStep*(iBin-1);
%     beta = 2*pi*fBin/c;        % wavenumber
%     df = exp(-1j*beta*xMics*sind(incidentAngle3(1))); 
%     D(:,iBin) = df;
% end
% 
% firD = real(fftshift(irfft(D,N,2)));
% %  fvtool((firD(4,:)),1)
% signal_D = recsignal;
% for iMic=1:M
%    signal_D(:,iMic) = conv(recsignal(:,iMic),firD(iMic,:),'same');
% end
% 
% 
% % figure()
% % plot(signal(:,1))
% % hold on
% % plot(signal_D(:,1),'r')
% 
% 
% y_B = zeros(M-1,L);
% h_up = zeros(L,1);
% h_up((L-1)/2+1) = 1;
% h_low = zeros(M-1,L);
% 
% 
% out = signal_D(:,4);
% % upper path
% for iLoop = 1:length(signal_D(:,1))- L + 1
%    y_c = sum(signal_D(iLoop:iLoop+L-1,:),2)/M;
%    y_up = sum(h_up.*y_c);
%    % blocking matrix
%    for iB = 1:M-1
%       y_B(iB,:) =  signal_D(iLoop:iLoop+L-1,iB) - signal_D(iLoop:iLoop+L-1,iB+1);  
%    end
% 
%    y_low = sum(sum(y_B.*h_low,2),1);
% 
%    y_out = y_up - y_low;
% 
%    Pk = sum(sum(y_B.*y_B,2),1);
%    mu = alpha1/Pk;
% 
%    delta = mu*y_out*y_B;
% %    if delta > 0.1
% %        delta = 0.1;
% %    elseif delta < -0.1
% %        delta = -0.1;
% %    end
%    h_low = h_low + delta; 
% 
%    out((L-1)/2+iLoop) = -y_out;
% 
% end

% pos = [0.045 0.045 0.45 0.45];
% 
% figure('numbertitle','off','name','Wave signales','Units','normal',...
%        'Position',pos);
% plot(t*1000,recsignal(:,4))
% hold on
% plot(t*1000,Source3,'g')
% plot(t*1000,out,'r')
% xlabel('Time (ms)')
% legend('Received signal','interest signal','GSC beamformed signal','Location','Best')
% 
% pos(2) = pos(2) +0.2;
% figure('numbertitle','off','name','Wave signales','Units','normal',...
%        'Position',pos);
% idx = 5800:7200;
% plot(t(idx)*1000,recsignal(idx,4))
% hold on
% plot(t(idx)*1000,Source3(idx),'g')
% plot(t(idx)*1000,out(idx),'r')
% xlabel('Time (ms)')
% legend('Received signal','interest signal','GSC beamformed signal','Location','Best')

%% Propose method

% Fixed beamformings design
fStep = fs/N;
klow = round(fl/fStep);     % low index
kup = round(fu/fStep);      % high index
f = fStep*(0:N/2);       % frequencies used to compute W
nf = length(f);
phi = pi/180*(0:1:360);
Nphi = length(phi);
WNG_up = zeros(nf,1);
WNG_low = zeros(nf,1);
DF_SD = zeros(nf,1);
DF_SDSS = zeros(nf,1);
bpdB_Up = zeros(nf,Nphi);
bpdB_low = zeros(nf,Nphi);
Wlow = zeros(M,nf);
Wup = zeros(M,nf);

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


% for k = 2:klow
%     beta = 2*pi*(k-1)*fStep/c;
%     d0 = exp(1j*beta*xMics'*cosd(phi_desired));        % steering vector at peak of side-lobe
%     Shi = (sin(beta*Gamma)./(beta*Gamma));
%     Shi(logical(eye(size(Shi)))) = 1;       
%     Wup(:,k) =(Shi + nguy*eye(length(xMics)))^-1*d0 / (d0'*(Shi + nguy*eye(length(xMics)))^-1*d0);
%     % Wup(:,k) = d0/M;
% end
Wup(:,klow:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow-1:kup),mu,null);


% for k = kup+2:N/2+1
%     beta = 2*pi*(k-1)*fStep/c;
%     df = exp(1j*beta*xMics'*cosd(phi_desired));        % steering vector at peak of side-lobe
%     Wup(:,k) = df/M;
% end


for k = klow:kup+1
    fTest = fStep*k
    [R,theta,p] = array_pattern_fft(xMics',Wup,f(k),k); 

    % for look direction finding for side-lobe beam
   [~, locs] = findpeaks(1 ./ R(1:floor(length(R)/2)));
   
   [~, locs_max] = max(R(phi_zero(1) + 1:phi_zero(2) + 1));

   phi_desired_s = locs_max + phi_zero(1) - 1;

   
   phi_zero_s = [0 (locs-1)];    
   % phi_zero_s = [0 phi_zero];
    
   fd_s = [1 zeros(size(phi_zero_s))];  % resonse in desired directions
   phi3_s = [phi_desired_s, phi_zero_s];
   theta3_s =  90*ones(1,length(phi3_s));
   W_s = bf_coefs(xMics',theta3_s,phi3_s,fd_s,f,mu,null);
    
   % scale factor
   beta = 2*pi*f(k)/c;
   df = exp(1j*beta*xMics'*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
   scale = Wup(:,k)'*df;

    D = exp(1j*beta*xMics'*cos(p));                    % steering matrix at a frequency
    bp = Wup(:,k)'*D;
    bpdB_Up(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    Wlow(:,k) = scale*W_s(:,k);
    bp = Wlow(:,k)'*D;
    bpdB_low(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    ds = exp(1j*beta*xMics'*cosd(phi_desired(1)));         % steering vector at look direction
    WNG_up(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Wup(:,k)));
    WNG_low(k) = 10*log10(1/(Wlow(:,k)'*Wlow(:,k)));
    
    Rdif = spatio_spect_corr(beta,xMics');
    DF_SD(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Rdif*Wup(:,k)));
    DF_SDSS(k) = 10*log10(1/(Wlow(:,k)'*Rdif*Wlow(:,k)));

    
end


Hup = fftshift(irfft(Wup,[],2),2);
Hlow = fftshift(irfft(Wlow,[],2),2);


pos = [0.045 0.45 0.3 0.35];
myFig = figure('numbertitle','off','name','Upper-path beam pattern',...
       'Units','normal','Position',pos);
ph = 180/pi*p; 

imagesc(ph(1:180),f,bpdB_Up(:,1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('f in Hz');
zlabel('Magnitude');
title('Directivity response of the array');
set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'normal');
%zlim([-dBmax 0])
set(gcf,'color','w');
grid on
colorbar 
%colormap jet
set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);


pos(1) = pos(1)+0.32;
myFig = figure('numbertitle','off','name','Lower-path beam pattern',...
       'Units','normal','Position',pos);

imagesc(ph(1:180),f, bpdB_low(:,1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('f in Hz');
zlabel('Magnitude');
title('Directivity response of the array');
set(gca, 'xdir', 'reverse');
set(gca, 'ydir', 'normal');
%zlim([-dBmax 0])
set(gcf,'color','w');
grid on
colorbar 
%colormap jet
set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);


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
   Pk0Old = Pk0;
   mu = alpha2/Pk0;
   delta = mu*y_out*yFill_low;
%    if delta > 0.1
%        delta = 0.1;
%    elseif delta < -0.1
%        delta = -0.1;
%    end
   h_l = h_l+ delta;
   
end

pos = [0.5 0.0 0.45 0.45];

figure('numbertitle','off','name','Wave signals','Units','normal',...
       'Position',pos);
plot(t*1000,recsignal(:,4))
hold on
plot(t*1000,Source3,'g','Linewidth',2)
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
plot(t(idx)*1000,Source3(idx),'g','Linewidth',2)
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
legend('Received signal','SOI','ABSS','Location','Best')
grid on
set(gcf,'color','w');
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);
axis tight

pos(2) = pos(2) - 0.3;
figure('numbertitle','off','name','Wave signales','Units','normal',...
       'Position',pos);
idx = 5800:7000;
plot(t(idx)*1000,signal3(idx,4),'g','Linewidth',2)
hold on
plot(t(idx)*1000,ygsc(idx),'b')
plot(t(idx)*1000,out_ABSS(idx),'r')
xlabel('Time (ms)')
legend('SOI','GSC','ABSS','Location','Best')
    set(gcf,'color','w'); 
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
axis tight
grid on
%%
e1 = mean((ygsc(idx) - signal3(idx,4)).^2)
e2 = mean((out_ABSS(idx) - signal3(idx,4)).^2)

 SOI_std = std(signal3(idx,4));
 oSINR_GSC =  20*log10(SOI_std /(eps +std(ygsc(idx) - signal3(idx,4) )))
 oSINR_ABSS =  20*log10(SOI_std /(eps +std(out_ABSS(idx) - signal3(idx,4))))
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



















































% comupte GSC in Frequency-domain
% % numFrame = floor(length(recsignal(:,1))/N);
% % y_B_FFT = zeros(M-1,1);
% % h_low_FFT = zeros(M-1,N);
% % out_FFT = zeros(N,1);
% % oGSB = zeros(length(recsignal(:,1)),1);
% % for iFrame = 1:numFrame
% %     signal_FFT = fft(recsignal((iFrame-1)*N+1:iFrame*N,:),N,1);
% %  
% %     for iBin = 1:N
% %         fBin = fStep*(iBin-1);
% %         beta = 2*pi*fBin/c;        % wavenumber
% %         ds = exp(-1j*beta*xMics'*sind(incidentAngle3(1))); % steering vector 
% %         up_FFT = signal_FFT(iBin,:)*ds;
% %         
% %         sync_FFT = signal_FFT(iBin,:).*ds;
% %         
% %         for iB = 1:M-1
% %             y_B_FFT(iB) =  sync_FFT(iB) - sync_FFT(iB+1);  
% %         end
% %         
% %         low_FFT = y_B_FFT'*h_low_FFT(:,iBin);
% %         
% %         out_FFT(iBin) = up_FFT - low_FFT;
% %         
% %         Pk = y_B_FFT'*y_B_FFT;
% %         
% %         mu = alpha/Pk;
% %         
% %         delta = mu*out_FFT(iBin)*y_B_FFT;
% %         
% %         h_low_FFT(:,iBin) = h_low_FFT(:,iBin) + delta;
% %         
% %     end
% %     oGSB((iFrame-1)*N+1:iFrame*N) = real(ifft(out_FFT));
% %     
% % end





