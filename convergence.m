
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


Wup(:,klow:kup+1) =  bf_coefs(xMics',theta3,phi3,fd,fStep*(klow-1:kup),mu,null);


for k = klow:kup+1
   fTest = fStep*k;
   [R,~,p] = array_pattern_fft(xMics',Wup,f(k),k); 

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

%%
Monte   = 5;               % at least 3–10 recommended
ITER    = 50;              % number of full adaptation passes over data
DUR     = 1024;            % number of output samples per pass (~0.15 s)

alpha1  = 0.001;            % step-size for upper path (if adapted – usually fixed)
alpha2s = [0.001, 0.01, 0.05];   % candidates for lower path

SNR     = -10;             % dB
SIR_vec = -5:10:-5;         % multiple SIR conditions to compare

% Storage for learning curves
EB_all     = zeros(length(alpha2s), ITER, length(SIR_vec), Monte);
EP_all     = zeros(length(alpha2s), ITER, length(SIR_vec), Monte);
oSINR_all  = zeros(length(alpha2s), ITER, length(SIR_vec), Monte);
oSINR_all_base = zeros(length(alpha2s), ITER, length(SIR_vec), Monte);
WH_dist_all= zeros(length(alpha2s), ITER,  length(SIR_vec), Monte);   
WH_dist_all_base= zeros(length(alpha2s), ITER,  length(SIR_vec), Monte); 

rng(123);  % for reproducibility

% ────────────────────────────────────────────────
%               Main Monte-Carlo loop
% ────────────────────────────────────────────────

for i_alpha = 1:length(alpha2s)
    alpha2 = alpha2s(i_alpha);
    
    for i_sir = 1:length(SIR_vec)
        SIR_dB = SIR_vec(i_sir);
        sigma_i = sqrt(0.5 * 10^(-SIR_dB/10));

        
        for mc = 1:Monte
            
            % ─── Generate signals (same as your code) ───
            t = (0 : DUR+N-1)/fs;
            
            s1 = chirp(t,1000,0.5,1000); 
            s1 = bandpass(s1,[fl fu],fs);
            s1 = sigma_i * s1 / std(s1);
            
            s2 = randn(size(t));
            s2 = bandpass(s2,[fl fu],fs);
            s2 = sigma_i * s2 / std(s2);
            
            % s3 = your desired source (tone cloud or speech ...)
            Nx = (N/2)-10;
            mindex = 0:Nx;
            phi    = pi*mindex.*mindex/Nx;
            k0     = 5;
            omega  = 2*pi*(k0 + mindex')*fs/N;
            s3     = sum(sin(omega.*t + phi'*ones(1,length(t))),1) / Nx;
            s3     = bandpass(s3,[fl fu],fs);
            s3     = s3 / std(s3);
            
            % Collect array signals
            sig1 = collector(s1(:), [  0; 0]);
            sig2 = collector(s2(:), [-90; 0]);
            sig3 = collector(s3(:), [ 90; 0]);   % look direction
            
            rec  = sig1 + sig2 + sig3;
            rec  = rec + sqrt(10^(-SNR/10)) * randn(size(rec));
            % ─── Fixed beamformer outputs ───
            y_up   = zeros(DUR,1);
            y_low  = zeros(DUR,1);
            
            yFill_up = zeros(L,1);
            yFill_low= zeros(L,1);
            
            for n = 1:DUR
                frame = rec(n:n+N-1,:)';
                y_beam_up  = sum(sum(Hup .* frame, 2));
                y_beam_low = sum(sum(Hlow.* frame, 2));
                
                yFill_up  = [y_beam_up;  yFill_up(1:end-1)];
                yFill_low = [y_beam_low; yFill_low(1:end-1)];
                
                y_up(n)  = yFill_up'  * ones(L,1);   
                y_low(n) = yFill_low' * ones(L,1);
            end
            
            % ─── Adaptive part ─────────────────────────────────────
            h_l       = zeros(L,1);  h_l((L+1)/2) = 1;
            h_l_base  = zeros(L,1);  h_l_base((L+1)/2) = 1;
            
            PkOld     = 0.1;
            Pk0Old    = 0.1;
            
            e2_base   = zeros(ITER,1);
            e2_prop   = zeros(ITER,1);
            oSINR_buf = zeros(ITER,1);
            oSINR_buf_base = zeros(ITER,1);
            
            
            % ─── Multiple adaptation passes ───
            for iter = 1:ITER
                y_out      = zeros(DUR,1);
                y_out_base = zeros(DUR,1);
                
                yFill_up(:)  = 0;
                yFill_low(:) = 0;
                
                e_buf      = zeros(DUR,1);   % instantaneous error (prop)
                e_buf_base = zeros(DUR,1);
                
                for n = 1:DUR
                    frame = rec(n:n+N-1,:)';
                    y_beam_up  = sum(sum(Hup .* frame,2));
                    y_beam_low = sum(sum(Hlow.* frame,2));
                    
                    yFill_up  = [y_beam_up;  yFill_up(1:end-1)];
                    yFill_low = [y_beam_low; yFill_low(1:end-1)];
                    
                    y_u   = yFill_up'  * ones(L,1);     
                    y_low = yFill_low' * h_l;
                    y_low_base = yFill_low' * h_l_base;
                    
                    outIdx = N/2 - (L-1)/2 + n;
                    y_out(outIdx)      = y_u - y_low;
                    y_out_base(outIdx) = y_u - y_low_base;
                    
                    % ─── NLMS-like updates ───
                    Pk0 = yFill_up' * yFill_up;   Pk0 = 0.9*Pk0Old + 0.1*Pk0;  Pk0Old = Pk0;
                    Pk  = yFill_low'* yFill_low;  Pk  = 0.9*PkOld  + 0.1*Pk;   PkOld  = Pk;
                    
                    g      = y_out(outIdx) * yFill_low;
                    g_base = y_out_base(outIdx) * yFill_low;
                    
                    h_l      = h_l      + (alpha2/Pk0) * g;
                    h_l_base = h_l_base + (alpha2/Pk)  * g_base;
                    
                    % error for convergence plot
                    e_buf(n)      = y_out(outIdx)      - s3(outIdx);  
                    e_buf_base(n) = y_out_base(outIdx) - s3(outIdx);
                end
                
                % Aggregate metrics per iteration
                e2_prop(iter) = mean(e_buf.^2);
                e2_base(iter) = mean(e_buf_base.^2);
                
                sig_var = var(s3(N/2- (L-1)/2+1 : N/2- (L-1)/2+DUR));
                oSINR_buf(iter) = 10*log10( sig_var / (mean(e_buf.^2)+1e-14) );
                oSINR_buf_base(iter) = 10*log10( sig_var / (mean(e_buf_base.^2)+1e-14) );

                
            end   
            
            % Save results
            EB_all(    i_alpha, :, i_sir, mc) = e2_base;
            EP_all(    i_alpha, :, i_sir, mc) = e2_prop;
            oSINR_all( i_alpha, :, i_sir, mc) = oSINR_buf;
            oSINR_all_base( i_alpha, :, i_sir, mc) = oSINR_buf_base;

        end   % Monte
    end   % SIR
end   % alpha

%% ===== Convergence over iterations (all α₂, different colors & markers) =====
pos = [0.5 0.0 0.45 0.45];

colors  = lines(length(alpha2s));              
markers = {'o','s','^','d','v','>','<','p','h','x','+'};

for i_sir = 1:length(SIR_vec)

    fig = figure('numbertitle','off', ...
        'Name', sprintf('Error convergence | SIR = %d dB, SNR = %d dB', SIR_vec(i_sir), SNR), ...
        'Units','normal','Position',pos);
    clf(fig); set(fig,'color','w');
    hold on;

    it = 1:ITER;

    for i_alpha = 1:length(alpha2s)
        avg_e2_base = mean(EB_all(i_alpha,:,i_sir,:),4);
        avg_e2_prop = mean(EP_all(i_alpha,:,i_sir,:),4);

        % --- BASELINE (dashed) ---
        plot(it, avg_e2_base, '--', ...
            'Color', 'b', ...
            'Marker', markers{mod(i_alpha-1,length(markers))+1}, ...
            'LineWidth',1.5, ...
            'DisplayName', sprintf('Base, \\beta = %.3f', alpha2s(i_alpha)));

        % --- PROPOSED (solid + marker) ---
        plot(it, avg_e2_prop, '-', ...
            'Color', 'g', ...
            'Marker', markers{mod(i_alpha-1,length(markers))+1}, ...
            'LineWidth',1.8, ...
            'MarkerSize',6, ...
            'DisplayName', sprintf('Proposed, \\beta = %.3f', alpha2s(i_alpha)));
    end

    grid on;
    xlabel('Iteration');
    ylabel('Mean squared error');
    title(sprintf('Error energy convergence | SIR = %d dB, SNR = %d dB', SIR_vec(i_sir), SNR));
    legend('Location','Best');
    set(gca,'FontSize',15);

end









