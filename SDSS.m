clc; % clearing the screen  
clear all; % clear all the variables 
close all; % close all the figures 
echo off

vs = 340;                    % sound speed
Nf = 512;                    % FFT length
Fs = 16000;                  % sampling frequency
Nh = Nf/2+1;                 % number of frequency points in [0,Fs/2]
Fsh = Fs/2;                  % half of sampling frequency
fl = 300;                    % lower cutoff frequency
fu = 5000;                   % higher cutoff frequency
klow = round(fl/Fsh*Nh);     % low index
kup = round(fu/Fsh*Nh);      % high index
f = Fsh/Nh*(klow:kup);       % frequencies used to compute W
mu = 0.01;                    % parameter for beam-forming                     
null = 1;                    % parameter for beam-forming

M = 7;                  % number of sensors
x_array = ones(1,M);    % uniform array geometry
dmics = 0.02;          % inter-dsitance  

% for main-lobe beam
phi_desired = 0;
phi_zero = 180;

theta_desired = 90;
theta_zero = 90;

fd = [1 zeros(size(phi_zero))];  % resonse in desired directions

phi3 = [phi_desired phi_zero];
theta3 = [theta_desired theta_desired];


i1 = find(x_array == 1);     % positions of active sensors in array
i1 = i1(:);
M = length(i1);              % number of active sensors
r = dmics*(i1(:)-i1(1));     % active sensor locations in m
mics = r - r(M)/2;           % x = 0 in middle of sensor array

% find optimum frequency-domain weight vector      
W = bf_coefs(mics,theta3,phi3,fd,f,mu,0);

% plot microphones positions
pos = [0.045 0.045 0.45 0.15];

tstr = ['M = ',int2str(M),', Nf = ',int2str(Nf)];

figure('numbertitle','off','name','Array geometry (cm)','Units','normal',...
       'Position',pos,'MenuBar','none');

plot(mics*100,zeros(M,1),'o','MarkerEdgeColor','k','MarkerFaceColor','r',...
     'MarkerSize',6);
set(gca, 'YDir'); 
grid on
set(gcf,'color','w');
grid minor
%set(gca,'xtick',[d_mics*100/2:5*100*d_mics:25*100*d_mics + 100*d_mics/2])
%set(gca,'ytick',[d_mics*100/2:5*100*d_mics:25*100*d_mics + 100*d_mics/2])

daspect([1 1 1]);
hold off

% plot beam response
fp = [500 1000 2000 4000 5000];
kf = round(fp/Fsh*Nh);
fp = Fsh/Nh*kf;              % frequencies rounded to FFT grid
kf = kf-klow+1;              % index used in W matrix corresponding to fp          

dBmax = 40;
t1 = pi/180*(-90:2.5:90);
pos = [0.045 0.045 0.3 0.6];

for k = 1:length(fp)

    figure('numbertitle','off','name',...
      'Fixed Array radiation pattern (dB)','Units','normal','Position',pos);
    
    %Plotting beam-pattern for 1st looking direction 
    [R,t,p] = array_pattern_fft(mics,W,fp(k),kf(k)); 
    [~,locs]=findpeaks(1./R(1:length(R)/2));
    [pks,locss] = findpeaks(R)
    idxPeaks = find(locss<=181);
    [~,idxMax] = max(pks);
    

    
    % for side-lobe beam
    phi_desired_s = 180;
    phi_zero_s = [0 locs];
    fd_s = [ones(size(phi_desired_s)) zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_s = bf_coefs(mics,theta3_s,phi3_s,fd_s,f,mu,null);
    % scale factor
    %scale = R(181) %max(pks)
    d = exp(1j*2*pi*fp(k)*mics*cosd(phi_desired_s(1))/vs);        % steering vector at peak of side-lobe
    scale = W(:,kf(k))'*d

    % side-lobe beam 
    [R_s,t,p] = array_pattern_fft(mics,W_s,fp(k),kf(k)); 
    RdB_s = max(0,10*log10(R_s+eps)+dBmax);
    RB_90_s = RdB_s(end,:);
    
    [R_ss,t,p] = array_pattern_fft(mics,W+(scale)*W_s,fp(k),kf(k)); 
    
    Noise_Power_Reduce =sum(sum(R))/(size(R,1)*size(R,2));
    SNR_Improve = -10*log10(Noise_Power_Reduce);
    
    RdB = max(0,10*log10(R+eps)+dBmax);
    
    % plotting ------------------------------------------------------------
    subplot(2,1,1); 
    
    polar(p,RdB); 
    ght = title(sprintf('f = %3.0f Hz', fp(k)));
    set(ght,'color','black');
    set(gcf,'color','w');
      
    hold on
    RdB_ss = max(0,10*log10(R_ss+eps)+dBmax);
    polar(p,RdB_ss,'--r');
    
    legend('SD','SDSS','Location','Best')
    hold off
    
    subplot(2,1,2); 
    polar(p,RdB_s);
    ght = title(sprintf('Side-lobe beam pattern, f = %3.0f Hz', fp(k)));
    
    pos(1) = pos(1) + 0.15;
end

%% WNG and DF
nf = length(f);
Nphi = length(p);
WNG_SD = zeros(nf,1);
WNG_SDSS = zeros(nf,1);
DF_SD = zeros(nf,1);
DF_SDSS = zeros(nf,1);
bpdB_SD = zeros(nf,Nphi);
bpdB_SDSS = zeros(nf,Nphi);
W_SDSS = zeros(M,nf)

for k = 1:nf
    [R,t,p] = array_pattern_fft(mics,W,f(k),k); 
    [~,locs]=findpeaks(1./R(1:length(R)/2));
    [pks,locss] = findpeaks(R);
    idxPeaks = find(locss<=181);
    [~,idxMax] = max(pks);
    
    % for side-lobe beam
    phi_desired_s = 180;
    phi_zero_s = [0 locs];
    fd_s = [ones(size(phi_desired_s)) zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_s = bf_coefs(mics,theta3_s,phi3_s,fd_s,f,mu,null);
    
    % scale factor
    beta = 2*pi*f(k)/vs;
    df = exp(1j*beta*mics*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
    scale = W(:,k)'*df;
    
    W_SDSS(:,k) = W(:,k)+(scale)*W_s(:,k);
    
    
    WNG_SD(k) = 10*log10(1/(W(:,k)'*W(:,k)));
    WNG_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*W_SDSS(:,k)));
    
    Rdif = spatio_spect_corr(beta,mics);
    DF_SD(k) = 10*log10(1/(W(:,k)'*Rdif*W(:,k)));
    DF_SDSS(k) = 10*log10(1/(W_SDSS(:,k)'*Rdif*W_SDSS(:,k)));
    
    D = exp(1j*beta*mics*cos(p));                    % steering matrix at a frequency
    bp = W(:,k)'*D;
    bpdB_SD(k,:) = min(20*log10(abs(bp)+eps),dBmax);
    bp = W_SDSS(:,k)'*D;
    bpdB_SDSS(k,:) = min(20*log10(abs(bp)+eps),dBmax);
end


pos = [0.055 0.04 0.4 0.38];
figure('numbertitle','off','name','WNG',...
       'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
plot(f,WNG_SD,'-b')
hold on
plot(f,WNG_SDSS,'-r')
legend('SD','SDSS','Location','SouthEast');
grid on;
set(gcf,'color','w');
xlabel('f (Hz)');
ylabel('WNG (dB)');

pos(1) = pos(1) + 0.3;
figure('numbertitle','off','name','DF',...
       'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
plot(f,DF_SD,'-b')
hold on
plot(f,DF_SDSS,'-r')
legend('SD','SDSS','Location','SouthEast');
grid on;
set(gcf,'color','w');
xlabel('f (Hz)');
ylabel('DF (dB)');
