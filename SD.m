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
fu = 5500;                   % higher cutoff frequency
klow = round(fl/Fsh*Nh);     % low index
kup = round(fu/Fsh*Nh);      % high index
f = Fsh/Nh*(klow:kup);       % frequencies used to compute W
mu = 0.01;                    % parameter for beam-forming                     
null = 0;                    % parameter for beam-forming

M = 9;                  % number of sensors
x_array = ones(1,M);   % uniform array geometry
dmics = 0.025; 

phi_desired = 0;
phi_zero = 180;

theta_desired = 90;
theta_zero = 90;

phi1 = phi_desired;
phi2 = phi_zero(:).';

fd = [1 zeros(size(phi2))];  % resonse in desired directions

theta1 = theta_desired;
theta2 = theta_zero(:).';

phi3 = [phi1 phi2];
theta3 = [theta1 theta2];



i1 = find(x_array == 1);  % positions of active sensors in thinned array
i1 = i1(:);
M = length(i1);           % number of active sensors
r = dmics*(i1(:)-i1(1));  % active sensor locations in m
mics = r - r(M)/2;           % x = 0 in middle of sensor array

% find optimum frequency-domain weight vector      
W = bf_coefs(mics,theta3,phi3,fd,f,mu,null);

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
fp = [500 2000 5000];
kf = round(fp/Fsh*Nh);
fp = Fsh/Nh*kf;              % frequencies rounded to FFT grid
kf = kf-klow+1;              % index used in W matrix corresponding to fp          

dBmax = 40;
t1 = pi/180*(-90:2.5:90);

Plot_Legend = {sprintf('%d Hz',round(fp(1))), sprintf('%d Hz',round(fp(2))), sprintf('%d Hz',round(fp(3)))};
Plot_Color = {'r', 'g', 'b', 'k'};
pos = [0.045 0.045 0.45 0.45];
figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
                  'Units','normal','Position',pos);
for k = 1:length(fp)

    %Plotting beam-pattern for 1st looking direction 
    [R,t,p] = array_pattern_fft(mics,W,fp(k),kf(k)); 
    R = R/max(R(:));

    Noise_Power_Reduce =sum(sum(R))/(size(R,1)*size(R,2));
    SNR_Improve = -10*log10(Noise_Power_Reduce);
    
    RdB = max(0,10*log10(R+eps)+dBmax);
    RB_90 = RdB(end,:);
    BW  = find(RB_90 > dBmax -3 );
    BW = fftshift(BW);
    %subplot(2,2,1); 
    polarplot(p,RdB,Plot_Color{k},'LineWidth',1.5); 
    grid on
    hold on
    

  
%     polar(p(BW(1))*ones(1,length(p)),linspace(0,dBmax,length(p)),'--g');
%     polar(p(BW(end))*ones(1,length(p)),linspace(0,dBmax,length(p)),'--g');
%     hold off
    %legend('25 mics','625 mics');
    
    % 3D plot
%     Xc = RdB .* (sin(t')*cos(p));
%     Yc = RdB .* (sin(t')*sin(p));
%     Zc = RdB .* (cos(t')*ones(1,length(p)));
%     mesh(Xc,Yc,Zc,RdB);
%     axis([-dBmax dBmax -dBmax dBmax 0 dBmax]);
%     view([25,50]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');

    % draw radial grid lines

%     hold on;
%     plot3(dBmax*cos(p),dBmax*sin(p),zeros(length(p),1),'b--');
%     plot3(dBmax*sin(t1),zeros(length(t1),1),dBmax*cos(t1),'b--');
%     plot3(zeros(length(t1),1),dBmax*sin(t1),dBmax*cos(t1),'b--');
%     hold off;
     
      %ght = title(sprintf('%s, f = %3.0f Hz, beam-width=%.1f, SNR improve=%.1f',tstr, fp(k),length(BW)*2.5,SNR_Improve));
%       ght = title(sprintf('f = %3.0f Hz, beam-width=%.1f °', fp(k),length(BW)));
%       set(ght,'color','black');
%       pos(1) = pos(1) + 0.305;
end
set(gcf,'color','w');
legend(Plot_Legend)
