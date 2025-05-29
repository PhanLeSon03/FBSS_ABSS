clc; % clearing the screen  
clear all; % clear all the variables 
close all; % close all the figures 
echo off

vs = 340;                    % sound speed
Nf = 512;                    % FFT length
Fs = 16000;                  % sampling frequency
Nh = Nf/2+1;                 % number of frequency points in [0,Fs/2]
Fsh = Fs/2;                  % half of sampling frequency
fl = 500;                    % lower cutoff frequency
fu = 6000;                   % higher cutoff frequency
klow = round(fl/Fsh*Nh);     % low index
kup = round(fu/Fsh*Nh);      % high index
f = Fsh/Nh*(klow:kup);       % frequencies used to compute W
mu = 0.1;                    % parameter for beam-forming                     
null = 1;                    % parameter for beam-forming

M = 7;                  % number of sensors
x_array = ones(1,M);    % uniform array geometry
dmics = 0.02;          % inter-dsitance  

% for main-lobe beam
phi_desired = 0;
phi_zero = [70 150];

theta_desired = 90;
theta_zero = [90 90];

fd = [1 zeros(size(phi_zero))];  % resonse in desired directions

phi3 = [phi_desired phi_zero];
theta3 = [theta_desired theta_zero];


i1 = find(x_array == 1);     % positions of active sensors in array
i1 = i1(:);
M = length(i1);              % number of active sensors
r = dmics*(i1(:)-i1(1));     % active sensor locations in m
mics = r- r(M)/2;           % x = 0 in middle of sensor array

%  x_sparse = [-0.1775 ;  -0.0925 ;  -0.0575  ; -0.0425 ;  -0.0225 ;  ...
%      0   ; 0.0225  ;  0.0425 ;   0.0575 ;   0.0925 ;   0.1775]; % array 0
%  mics = x_sparse;
%  M = length(mics);
% find optimum frequency-domain weight vector      
W = bf_coefs(mics,theta3,phi3,fd,f,mu,null);

% comparison: Superdirective beamforming with multiple constraints  
phi_desired_SDM = 0;
phi_zero_SDM = [70 150 180];

theta_desired_SDM = 90;
theta_zero_SDM = [90 90 90];

fd_SDM = [1 zeros(size(phi_zero_SDM))];  % resonse in desired directions

phi3_SDM = [phi_desired_SDM phi_zero_SDM];
theta3_SDM = [theta_desired_SDM theta_zero_SDM];

W_SDM = bf_coefs(mics,theta3_SDM,phi3_SDM,fd_SDM,f,mu,null);

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
fp = [1000 2000 4000 5500];
kf = round(fp/Fsh*Nh);
fp = Fsh/Nh*kf;              % frequencies rounded to FFT grid
kf = kf-klow+1;              % index used in W matrix corresponding to fp          

dBmax = 50;
t1 = pi/180*(-90:2.5:90);
pos = [0.045 0.045 0.45 0.45];

for k = 1:length(fp)

    myFig = figure('numbertitle','off','name',...
      'Fixed Array radiation pattern (dB)','Units','normal','Position',pos);
    
    % upper path
    [R,t,p] = array_pattern_fft(mics,W,fp(k),kf(k));     
    [Phase,~,~] = array_phase_fft(mics,W,fp(k),kf(k));  
    
    % lower beamforming design
    [~,locs]=findpeaks(1./R(1:length(R)/2));
    phi_desired_s = 180;
    phi_zero_s = [0 (locs-1)];  
       
    % lower path weight compute
    fd_s = [ones(size(phi_desired_s)) zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_s = bf_coefs(mics,theta3_s,phi3_s,fd_s,f,mu,null);
    % scale factor
    d = exp(j*2*pi*fp(k)*mics*cosd(phi_desired_s(1))/vs);        % steering vector at peak of side-lobe
    scale = W(:,kf(k))'*d


    % lowerpath beam 
    [R_s,t,p] = array_pattern_fft(mics,W_s,fp(k),kf(k)); 
    [Phase_s,~,~] = array_phase_fft(mics,scale*W_s,fp(k),kf(k)); 
    
    RdB_s = max(0,10*log10(R_s+eps)+dBmax);
    RB_90_s = RdB_s(end,:);
    
    % final beam
    [R_ss,t,p] = array_pattern_fft(mics,W-(scale)*W_s,fp(k),kf(k)); 
    [Phase_ss,~,~] = array_phase_fft(mics,W-(scale)*W_s,fp(k),kf(k)); 
    
    RdB = max(0,10*log10(R+eps)+dBmax);
    
    %comparison
    [R_SDM,t,p] = array_pattern_fft(mics,W_SDM,fp(k),kf(k));  
    
    RdB_SDM = max(0,10*log10(R_SDM+eps)+dBmax);
    
    % plotting ------------------------------------------------------------
    subplot(1,2,1); 
    
    polarplot(p,RdB,'b','LineWidth',1.5); 
    ght = title(sprintf('f = %3.0f Hz', fp(k)),'FontSize',15);
    set(ght,'color','black');
    set(gcf,'color','w');  
    set(gca,'FontSize', 15);
    hold on

    RdB_ss = max(0,10*log10(R_ss+eps)+dBmax);
    polarplot(p,RdB_ss,'-r','LineWidth',1.5);
    
    polarplot(p,RdB_SDM,'--g','LineWidth',1.5);
    
    legend('Upper path','FBSS','SDM','Location','southoutside')
    hold off
    
    subplot(1,2,2); 
    polarplot(p,RdB_s,'b','LineWidth',1.5);
    %ght = title(sprintf('Lower path''s beam pattern, f = %3.0f Hz', fp(k)),'FontSize',10);
    ght = title(sprintf('f = %3.0f Hz', fp(k)),'FontSize',15);
    legend('Lower path','FBSS','Location','southoutside')
    pos(1) = pos(1) + 0.15;
    set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal');
    set(gca,'FontSize', 15);
end

%% WNG and DF
nf = length(f);
Nphi = length(p);
WNG_SD = zeros(nf,1);
WNG_SDSS = zeros(nf,1);
WNG_SDM = zeros(nf,1);

DF_SD = zeros(nf,1);
DF_SDSS = zeros(nf,1);
DF_SDM = zeros(nf,1);

bpdB_SD = zeros(nf,Nphi);
bpdB_SD_low = zeros(nf,Nphi);
bpdB_SDSS = zeros(nf,Nphi);


W_SDSS = zeros(M,nf);

kx = round(1000/Fsh*Nh)-klow+1; % index threshold for spatial aliasing
for k = 1:nf
    [R,t,p] = array_pattern_fft(mics,W,f(k),k); 

    
    % for side-lobe beam
    if k < kx
        phi_desired_s = 180;
        phi_zero_s = [0 phi_zero];
    else
      [~,locs]=findpeaks(1./R(1:length(R)/2));
       phi_desired_s = 180;
       phi_zero_s = [0 (locs-1)];  
    end
    fd_s = [ones(size(phi_desired_s)) zeros(size(phi_zero_s))];  % resonse in desired directions
    phi3_s = [phi_desired_s, phi_zero_s];
    theta3_s =  90*ones(1,length(phi3_s));
    W_s = bf_coefs(mics,theta3_s,phi3_s,fd_s,f,mu,null);
    
    % scale factor
    beta = 2*pi*f(k)/vs;
    df = exp(j*beta*mics*cosd(phi_desired_s(1)));        % steering vector at peak of side-lobe
    scale = W(:,k)'*df;
    
    W_SDSS(:,k) = W(:,k)-(scale)*W_s(:,k);
    
    ds = exp(1j*beta*mics*cosd(phi_desired(1)));         % steering vector at look direction
    WNG_SD(k) = 10*log10(abs(W(:,k)'*ds)/(W(:,k)'*W(:,k)));
    WNG_SDSS(k) = 10*log10(abs(W_SDSS(:,k)'*ds)/(W_SDSS(:,k)'*W_SDSS(:,k)));
    WNG_SDM(k) = 10*log10(abs(W_SDM(:,k)'*ds)/(W_SDM(:,k)'*W_SDM(:,k)));
    
    Rdif = spatio_spect_corr(beta,mics);
    
    DF_SD(k) = 10*log10(abs(W(:,k)'*ds)/(W(:,k)'*Rdif*W(:,k)));
    DF_SDSS(k) = 10*log10(abs(W_SDSS(:,k)'*ds)/(W_SDSS(:,k)'*Rdif*W_SDSS(:,k)));
    DF_SDM(k) = 10*log10(abs(W_SDM(:,k)'*ds)/(W_SDM(:,k)'*Rdif*W_SDM(:,k)));
    
    D = exp(1j*beta*mics*cos(p));                    % steering matrix at a frequency
    bp = W(:,k)'*D;
    bpdB_SD(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    bp = W_SDSS(:,k)'*D;
    bpdB_SDSS(k,:) = abs(bp);%min(20*log10(abs(bp)+eps),dBmax);
    
    bp = (scale)*W_s(:,k)'*D;
    bpdB_SD_low(k,:) = abs(bp)
end

%% plot phase of 2 beamforming
pos = [0.055 0.04 0.4 0.38];
figure('numbertitle','off','name','Phases',...
       'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
plot(0:1:180,Phase*180/pi,'-b')
hold on
plot(0:1:180,Phase_s*180/pi,'-r')
legend('Upper path','Lower path','Location','SouthEast');
grid on;
set(gcf,'color','w');
xlabel('Arrival direction');
ylabel('Phase');


%%
pos = [0.055 0.04 0.4 0.38];
figure('numbertitle','off','name','WNG',...
       'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
plot(f,WNG_SD,'-sb','LineWidth',1.0)
hold on
plot(f,WNG_SDSS,'-dr','LineWidth',1.0)
plot(f,WNG_SDM,'-^g','LineWidth',1.0)
legend('Upper path','FBSS','SDM','Location','SouthEast');
grid on;
set(gcf,'color','w');
xlabel('Frequency (Hz)');
ylabel('WNG (dB)');
axis tight
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);


pos(1) = pos(1) + 0.3;
figure('numbertitle','off','name','DF',...
       'Units','normal','Position',pos); %%sop1hc ,'Menubar','none'
plot(f,DF_SD,'-sb','LineWidth',1.0)
hold on
plot(f,DF_SDSS,'-dr','LineWidth',1.0)
plot(f,DF_SDM,'-^g','LineWidth',1.0)
legend('Upper path','FBSS','SDM','Location','SouthEast');
grid on;
set(gcf,'color','w');
xlabel('Frequency (Hz)');
ylabel('DF (dB)');
axis tight
    set(gcf,'defaultAxesFontSize',15)
    set(gca,'FontSize', 15);
%%
pos = [0.045 0.45 0.3 0.35];
myFig = figure('numbertitle','off','name','LCMV beam pattern',...
       'Units','normal','Position',pos);
ph = 180/pi*p; 
  
imagesc(ph(1:180),f,bpdB_SD(:,1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('Frequency (Hz)');
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


pos(1) = pos(1)+0.16;
myFig = figure('numbertitle','off','name','MLMV beam pattern',...
       'Units','normal','Position',pos);
  
imagesc(ph(1:180),f, bpdB_SD_low(:,1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('Frequency (Hz)');
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


pos(1) = pos(1)+0.32;
myFig = figure('numbertitle','off','name','MLMV beam pattern',...
       'Units','normal','Position',pos);
  
imagesc(ph(1:180),f, bpdB_SDSS(:,1:180));
axis tight
%set(gca,'XTick',[0 45 90 135 180]);
%view([25,50]);
xlabel('Incident angle \phi in °');
ylabel('Frequency (Hz)');
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


