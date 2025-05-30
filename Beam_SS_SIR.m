clear
clc
close all
c = 340.0;
fs = 16000;
fl = 500;                    % lower cutoff frequency
fu = 6000;                   % higher cutoff frequency
fTest = 4000;
lam = c/fu;
dmics = 0.02; 
N = 512;                    % number of samples in one frame
M = 7;                      % number of microphones
L = 511;
xMics = (-(M-1)/2:(M-1)/2)*dmics;
%xMics = (0:M-1)*dmics;
transducer = phased.OmnidirectionalMicrophoneElement; %('FrequencyRange',[20 20000])
array = phased.ULA('Element',transducer,'NumElements',M,'ElementSpacing',dmics);
alpha1 = 0.01;
alpha2 = 0.01;
t = 0:1/fs:0.5;

incidentAngle1 = [0 ;0]; %10° azimuth and 90° elevation, 20 = 90 -70
incidentAngle2 = [-90 ;0]; % -60 = 90 -150
incidentAngle3 = [90 ;0]; % source of interset


% Create an incident wave arriving at the array. Add gaussian noise to the wave.
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',N);



%% Fixed beamformings design
fStep = fs/N;
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
    Wup(:,k) =(Shi + nguy*eye(length(xMics)))^-1*d0 / (d0'*(Shi + nguy*eye(length(xMics)))^-1*d0);
    Wup(:,k) = d0/M;
end
Wup(:,klow:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow-1:kup),mu,null);
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
    W_s = bf_coefs(xMics',theta3_s,phi3_s,fd_s,f,mu,null);
    
    % scale factor
    beta = 2*pi*f(k)/c;
    df = exp(1j*beta*xMics'*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
    scale = Wup(:,k)'*df;

    D = exp(1j*beta*xMics'*cos(p));                    % steering matrix at a frequency
    bp = Wup(:,k)'*D;
    bpdB_SD(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    W_SDSS(:,k) =  scale*W_s(:,k);
    bp = W_SDSS(:,k)'*D;

    bpdB_SDSS(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    ds = exp(1j*beta*xMics'*cosd(phi_desired(1)));         % steering vector at look direction
    WNG_SD(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Wup(:,k)));
    WNG_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*W_SDSS(:,k)));
    
    Rdif = spatio_spect_corr(beta,xMics');
    DF_SD(k) = 10*log10(abs(Wup(:,k)'*ds)/(Wup(:,k)'*Rdif*Wup(:,k)));
    DF_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*Rdif*W_SDSS(:,k)));

    
end

Wlow = W_SDSS;
Wlow(:,1:klow) = 0; 
Wlow(:,kup+2:end) = 0; 

Hup = fftshift(irfft(Wup,[],2),2);
Hlow = fftshift(irfft(Wlow,[],2),2);


pos = [0.045 0.45 0.3 0.35];
myFig = figure('numbertitle','off','name','LCMV beam pattern',...
       'Units','normal','Position',pos);
ph = 180/pi*p; 
  
imagesc(ph(1:180),f,bpdB_SD(:,1:180));
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
myFig = figure('numbertitle','off','name','MLMV beam pattern',...
       'Units','normal','Position',pos);
  
imagesc(ph(1:180),f, bpdB_SDSS(:,1:180));
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
SIR = (-30:2.5:10);

E1 = zeros(length(SIR),2);
E2 = zeros(length(SIR),2);
E3 = zeros(length(SIR),2);

K =200;
yStack = zeros(L,K);
y_up_stack = zeros(K,1);

for iSNR = 1:2
    idxE = 1;
    for iSIR = SIR
        sigma_i = sqrt(0.5*10^(-iSIR/10));
        for iMonte = 1:Monte
            % Simulate a chirp signal with a 500 Hz bandwidth.
            %Source1 = sigma_i*chirp(t,5000,0.5,5000);%0.5*randn(1,8001);%
            
            Nx=(N/2)-10; % the number of tones
            mindex=0:Nx;
            phi=pi*mindex.*mindex/Nx;
            k0 = 5;
            omega=2*pi*(k0*ones(Nx+1,1)+mindex')*fs/N;
            Source1=sum(sin(omega*t+phi'*ones(1,length(t))),1)/Nx;

            Source2 =  sigma_i*randn(1,8001);%0.5*chirp(t,2000,0.5,2000);%
            
            Source3 = zeros(size(Source1));
            Source3(6001) = 10;
            Source3(6500) = 10;
            Source3 = Source3 + chirp(t,fTest-2000,0.5,fTest-2000);

            Source1 = bandpass(Source1,[fl fu],fs);
            Source2 = bandpass(Source2,[fl fu],fs);
            Source3 = bandpass(Source3,[fl fu],fs);


            signal1 = collector(Source1.' ,incidentAngle1);
            signal2 = collector(Source2.' ,incidentAngle2);
            signal3 = collector(Source3.' ,incidentAngle3);

            signal = signal1 + signal2 + signal3;

            SNR = iSNR*10;
            noise = sqrt(10^(-SNR/10))*randn(size(signal));
            recsignal = signal + noise;

            % pspectrum(recsignal(:,4),fs,'spectrogram','TimeResolution',0.1, ...
            %       'OverlapPercent',99,'Leakage',0.85)


            % MATLAB toolbox for GSC---------------------------------------------------
            gscbeamformer = phased.GSCBeamformer('SensorArray',array, ...
                'PropagationSpeed',c,'SampleRate',fs,'DirectionSource','Input port', ...
                'FilterLength',L,'LMSStepSize',alpha1);

            ygsc = gscbeamformer(recsignal,[90;0]);


            out_ABSS = zeros(length(recsignal(:,1)),1);
            out_FBSS = zeros(length(recsignal(:,1)),1);
            yFill_up = zeros(L,1);
            yFill_low = zeros(L,1);
            h_u = zeros(L,1);
            h_u((L-1)/2+1) = 1;
            h_l = zeros(L,1);
            h_l((L-1)/2+1) = 1;
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

               y_out = y_up - y_low;

               out_ABSS(N/2-(L-1)/2+iLoop) =  y_out;
               out_FBSS(N/2-(L-1)/2+iLoop) =   yFill_up((L-1)/2+1) - yFill_low((L-1)/2+1);
                   % either NLMS algorithm or steepest decent 
                   if Alg == 1
                               % update filter coefficient
                               Pk0 = sum(yFill_up.*yFill_up)+0.00001;
                               Pk = sum(yFill_low.*yFill_low)+0.00001;
                               %Pk = 0.95*Pk + 0.05*PkOld;
                               %PkOld = Pk;
                               Pk0 = 0.9*Pk0 + 0.1*Pk0Old;
                               Pk0Old = Pk0;

                               g = y_out*yFill_low; 
                               %g = 0.95*g + 0.05*gOld;
                               %gOld = g;

                               mu = alpha2/Pk0;
                               delta = mu*g;
                               
%                                for iD = i:L
%                                    if delta(iD) > 0.5
%                                        delta(iD) = 0.5;
%                                    elseif delta(iD) < -0.5
%                                        delta(iD) = -0.5;
%                                    end
%                                end
                               h_l = h_l*(1-0.2*mu) + delta;
                   else
                               yStack(:,mod(iLoop,K)+1) = yFill_low; 
                               y_up_stack(mod(iLoop,K)+1) = y_up;

                               Q = (yStack*yStack')/K ; %  covariance matrix 
                               b = (yStack*y_up_stack)/K;       %  correlation between upper path and lower path
                               %h_l = Q^-1*b;
                               g =Q*h_l -b;
                               alpha_n = (g'*g/(g'*Q*g+0.01));
                               delta = - alpha2*alpha_n*g;

                %                for iD = i:L
                %                    if delta(iD) > 0.5
                %                        delta(iD) = 0.5;
                %                    elseif delta(iD) < -0.5
                %                        delta(iD) = -0.5;
                %                    end
                %                end


                               h_l = h_l + delta;
                   end          
            end


    
            idx = 5800:7000;
            E1(idxE,iSNR) = E1(idxE,iSNR) + mean((ygsc(idx) - signal3(idx,4)).^2)
            E2(idxE,iSNR) = E2(idxE,iSNR) + mean((out_ABSS(idx) - signal3(idx,4)).^2)
            E3(idxE,iSNR) = E3(idxE,iSNR) + mean((out_FBSS(idx) - signal3(idx,4)).^2);
        end
        E1(idxE,iSNR) = E1(idxE,iSNR)/Monte;
        E2(idxE,iSNR) = E2(idxE,iSNR)/Monte;
        E3(idxE,iSNR) = E3(idxE,iSNR)/Monte;
        idxE = idxE +1;
    end
end
pos = [0.5 0.0 0.45 0.45];

figure('numbertitle','off','name','Signal Error','Units','normal',...
       'Position',pos);

plot(SIR,E1(:,1),'-sb','Linewidth',1.5)
hold on
plot(SIR,E2(:,1),'-sg','Linewidth',1.5)
plot(SIR,E3(:,1),'-sr','Linewidth',1.5)

plot(SIR,E1(:,2),'-db','Linewidth',1.5)
plot(SIR,E2(:,2),'-dg','Linewidth',1.5)
plot(SIR,E3(:,2),'-dr','Linewidth',1.5)

xlabel('SIR (dB)')
ylabel('SE')

legend('GSC(SNR=10 dB)','ABSS(SNR=10 dB)','FBSS(SNR=10 dB)', 'GSC(SNR=20 dB)','ABSS(SNR=20 dB)','FBSS(SNR=20 dB)','Location','Best')
set(gcf,'color','w');
grid on
set(gcf,'defaultAxesFontSize',15)
set(gca,'FontSize', 15);























































