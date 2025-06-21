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
% xMics = [-4,-3,-1,0,2,4,6]*dmics;

transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);
alpha1 = 0.01;
alpha2 = 0.01;
t = 0:1/fs:0.5;

incidentAngle1 = [0 ;0]; %10° azimuth and 90° elevation, 0 = 90 -90
incidentAngle2 = [-90 ;0]; % -90 = 90 -180
incidentAngle3 = [90 ;0]; % source of interset


% Create an incident wave arriving at the array. Add gaussian noise to the wave.
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',N);


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
%%
Alg = 1;
Monte=25;
SIR = (-30:5:10);

E1 = zeros(length(SIR),5);
E2 = zeros(length(SIR),5);


K =200;
yStack = zeros(L,K);
y_up_stack = zeros(K,1);


for iRate = 1:3
    idxE = 1;
    alpha2 = 0.05*iRate;
    for iSIR = SIR
        sigma_i = sqrt(0.5*10^(-iSIR/10));
        for iMonte = 1:Monte

            Source1=chirp(t,fl,0.5,fu);

            Source2 =  randn(1,8001);
            

            Nx=(N/2)-10; % the number of tones
            mindex=0:Nx;
            phi=pi*mindex.*mindex/Nx;
            k0 = 5;
            omega=2*pi*(k0*ones(Nx+1,1)+mindex')*fs/N;
            Source3 = sum(sin(omega*t+phi'*ones(1,length(t))),1)/Nx; % multi-tone signal
            Source3(6001) = Source3(6001) + 10;
            Source3(6500) = Source3(6500)+ 10;

            Source1 = bandpass(Source1,[fl fu],fs);
            Source2 = bandpass(Source2,[fl fu],fs);
            Source3 = bandpass(Source3,[fl fu],fs);


            Source1 = sigma_i*Source1/std(Source1);
            Source2 = sigma_i*Source2/std(Source2);
            Source3 = Source3/std(Source3);

            signal1 = collector(Source1.' ,incidentAngle1);
            signal2 = collector(Source2.' ,incidentAngle2);
            signal3 = collector(Source3.' ,incidentAngle3);


            signal = signal1 + signal2 + signal3;

            SNR = -10;
            noise = randn(size(signal));
            noise = sqrt(10^(-SNR/10))*noise/std(noise);

            recsignal = signal + noise;

            out_ABSS = zeros(length(recsignal(:,1)),1);
            out_ABSS_base = zeros(length(recsignal(:,1)),1);
            yFill_up = zeros(L,1);
            yFill_low = zeros(L,1);
            h_u = zeros(L,1);
            h_u((L-1)/2+1) = 1;
            h_l = zeros(L,1);
            h_l((L-1)/2+1) = 1;

            h_l_base = zeros(L,1);
            h_l_base((L-1)/2+1) = 1;

            gOld = 0;
            PkOld = 0.1;
            Pk0Old = 0.1;
            Q = 0.01*eye(L,L);
            b = yFill_low;
            for iLoop = 1:length(recsignal(:,1))- N + 1
               
               y_beam1 = sum(sum(Hup.*recsignal(iLoop:iLoop+N-1,:)',2),1);
               yFill_up = circshift(yFill_up,1);
               yFill_up(1) = y_beam1;
               y_up = sum(h_u.*yFill_up);
               
               y_beam2 = sum(sum(Hlow.*recsignal(iLoop:iLoop+N-1,:)',2),1);
               yFill_low = circshift(yFill_low,1);
               yFill_low(1) = y_beam2;

               y_low = sum(h_l.*yFill_low);
               y_low_base = sum(h_l_base.*yFill_low);

               y_out = y_up - y_low;
               y_out_base = y_up - y_low_base;

               out_ABSS(N/2-(L-1)/2+iLoop) =  y_out;
               out_ABSS_base(N/2-(L-1)/2+iLoop) =  y_out_base;
        
               % update filter coefficient
               Pk0 = sum(yFill_up.*yFill_up);
               Pk = sum(yFill_low.*yFill_low);
               Pk = 0.9*Pk + 0.1*PkOld;
               PkOld = Pk;
               Pk0 = 0.9*Pk0 + 0.1*Pk0Old;
               Pk0Old = Pk0;

               g = y_out*yFill_low; 
               g_base = y_out_base*yFill_low; 


               mu = alpha2/Pk0;
               delta = mu*g;

               mu = alpha2/Pk;
               delta_base = mu*g_base;
               
               h_l = h_l + delta;
               h_l_base = h_l_base + delta_base;
                          
            end


            idx = 800:6800;
            SOI_std = std(signal3(idx,4));
            E1(idxE,iRate) = E1(idxE,iRate) + 20*log10(SOI_std /(eps +std(out_ABSS_base(idx) - signal3(idx,4) )))
            E2(idxE,iRate) = E2(idxE,iRate) + 20*log10(SOI_std /(eps +std(out_ABSS(idx) - signal3(idx,4))))
            
        end
        E1(idxE,iRate) = E1(idxE,iRate)/Monte;
        E2(idxE,iRate) = E2(idxE,iRate)/Monte;
        idxE = idxE +1;
    end
end
%%
pos = [0.5 0.0 0.45 0.45];

figure('numbertitle','off','name','oSIR: learning rate','Units','normal',...
       'Position',pos);
hold on;
% === Main Plot ===
% Plot for base method (e.g., GSC)
plot(SIR, E1(:,1), '-ob', 'LineWidth', 1.5); % beta = 0.05
plot(SIR, E1(:,2), '-^b', 'LineWidth', 1.5); % beta = 0.1
plot(SIR, E1(:,3), '-sb', 'LineWidth', 1.5); % beta = 0.15


% Plot for proposed method (e.g., ABSS adaptive)
plot(SIR, E2(:,1), '-og', 'LineWidth', 1.5); % beta = 0.05
plot(SIR, E2(:,2), '-^g', 'LineWidth', 1.5); % beta = 0.1
plot(SIR, E2(:,3), '-sg', 'LineWidth', 1.5); % beta = 0.15


xlabel('SIR (dB)')
ylabel('oSINR (dB)')
legend(...
    'Base, \beta = 0.05', 'Base, \beta = 0.1', 'Base, \beta = 0.15', ...
    'Proposed, \beta = 0.05', 'Proposed, \beta = 0.1', 'Proposed, \beta = 0.15', ...
    'Location','Best')

set(gca, 'FontSize', 15)
set(gcf, 'color', 'w')
grid on;


% === Define zoomed region ===
x1 = 4.8; x2 = 5.2;
y1 = 16.5; y2 = 19;

% === Draw circle to indicate zoom region ===
theta = linspace(0, 2*pi, 100);
xc = (x1 + x2)/2;
yc = (y1 + y2)/2;
rx = (x2 - x1)/2;
ry = (y2 - y1)/2;
plot(xc + rx*cos(theta), yc + ry*sin(theta), 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% === Inset axes at safe position (avoid overlap) ===
inset_axes = axes('Position',[0.2 0.7 0.2 0.22]); % Change as needed
box on;
hold on;
plot(SIR, E1(:,1), '-ob', 'LineWidth', 1.5); % beta = 0.05
plot(SIR, E1(:,2), '-^b', 'LineWidth', 1.5); % beta = 0.1
plot(SIR, E1(:,3), '-sb', 'LineWidth', 1.5); % beta = 0.15


% Plot for proposed method (e.g., ABSS adaptive)
plot(SIR, E2(:,1), '-og', 'LineWidth', 1.5); % beta = 0.5
plot(SIR, E2(:,2), '-^g', 'LineWidth', 1.5); % beta = 0.1
plot(SIR, E2(:,3), '-sg', 'LineWidth', 1.5); % beta = 0.15
x1 = 4.8; x2 = 5.2;
y1 = 17.7; y2 = 17.9;
xlim([x1 x2])
ylim([y1 y2])
set(gca,'FontSize',10)
grid on;

























































