%% =========================================================
%  Full Script: Proposed ABSS vs GSC
%  Sections:
%   1. Parameters & array setup
%   2. Fixed beamformer design (Wup, Wlow -> Hup, Hlow)
%   3. GSC blocking matrix
%   4. Single trial: run both adaptive filters
%   5. Figure: 2-panel weight comparison
%      Panel 1: converged h (tap domain)
%      Panel 2: learning curve — proposed h_l + GSC h_1..h_6 individually
%   6. Figure: single zoomed plot around the pulse
% =========================================================
clear; clc; close all;

%% =========================================================
%  1. PARAMETERS
% =========================================================
c     = 340.0;
fs    = 16000;
fl    = 1000;
fu    = 6000;
N     = 512;
M     = 7;
L     = 511;
dmics = 0.02;
xMics = (-(M-1)/2:(M-1)/2) * dmics;

t = 0:1/fs:0.5;

incidentAngle1 = [0;  0];
incidentAngle2 = [-90; 0];
incidentAngle3 = [90;  0];   % SOI

transducer = phased.OmnidirectionalMicrophoneElement;
array      = phased.ULA('Element',transducer,'NumElements',M,...
                        'ElementSpacing',dmics);
collector  = phased.WidebandCollector('Sensor',array,...
                'PropagationSpeed',c,'SampleRate',fs,...
                'ModulatedInput',false,'NumSubbands',N);

%% =========================================================
%  2. FIXED BEAMFORMER DESIGN
% =========================================================
fStep = fs/N;
klow  = round(fl/fStep);
kup   = round(fu/fStep);
f     = fStep*(0:N/2);

phi_desired = 0;
phi_zero    = [70 150];
fd          = [1 zeros(size(phi_zero))];
phi3        = [phi_desired phi_zero];
theta3      = 90*ones(1,length(fd));
mu_bf       = 0.1;
null_bf     = 1;

Wup  = zeros(M, length(f));
Wlow = zeros(M, length(f));
Wdsb = zeros(M, length(f));

Wup(:,klow:kup+1)  = bf_coefs(xMics',theta3,phi3,fd,...
                               fStep*(klow-1:kup),mu_bf,null_bf);
Wdsb(:,klow:kup+1) = dsb_coefs(xMics',theta3,phi3,fStep*(klow-1:kup));

for k = klow:kup+1
    [R,~,p] = array_pattern_fft(xMics', Wup, f(k), k);

    [~, locs]     = findpeaks(1 ./ R(1:floor(length(R)/2)));
    [~, locs_max] = max(R(phi_zero(1)+1:phi_zero(2)+1));
    phi_desired_s = locs_max + phi_zero(1) - 1;
    phi_zero_s    = [0 (locs-1)];
    fd_s          = [1 zeros(size(phi_zero_s))];
    phi3_s        = [phi_desired_s, phi_zero_s];
    theta3_s      = 90*ones(1,length(phi3_s));
    W_s           = bf_coefs(xMics',theta3_s,phi3_s,fd_s,f,mu_bf,null_bf);

    beta_bf = 2*pi*f(k)/c;
    df      = exp(1j*beta_bf*xMics'*cosd(phi_desired_s(1)));
    scale   = Wup(:,k)' * df;
    Wlow(:,k) = scale * W_s(:,k);
end

Hup  = fftshift(irfft(Wup, [], 2), 2);
Hlow = fftshift(irfft(Wlow,[], 2), 2);
Hdsb = fftshift(irfft(Wdsb,[], 2), 2);

fprintf('Filter banks ready: Hup and Hlow [%dx%d]\n', size(Hup,1), size(Hup,2));

%% =========================================================
%  3. GSC BLOCKING MATRIX
% =========================================================
B_gsc = zeros(M-1, M);
for i = 1:M-1
    B_gsc(i,i)   =  1;
    B_gsc(i,i+1) = -1;
end
fprintf('B_gsc [%dx%d], rank=%d, max|row-sum|=%.1e\n',...
    size(B_gsc,1),size(B_gsc,2),rank(B_gsc),max(abs(sum(B_gsc,2))));

phi_look_deg = 0;
tau = round(xMics * cosd(phi_look_deg) / c * fs);
fprintf('Time-alignment delays (samples): ');
fprintf('%d ', tau); fprintf('\n');

%% =========================================================
%  4. PARAMETERS FOR SINGLE TRIAL
% =========================================================
iSIR_val  = -5;
SNR       = -10;
beta      = 0.01;
beta_ABSS = (M-1)*beta;
sigma_i   = sqrt(0.5 * 10^(-iSIR_val/10));

%% --- Signal generation ---
Source1 = chirp(t, fl, 0.5, fu);
Source2 = randn(1, 8001);

Nx     = (N/2) - 10;
midx   = 0:Nx;
phi_mt = pi*midx.*midx/Nx;
k0     = 5;
omega  = 2*pi*(k0*ones(Nx+1,1)+midx')*fs/N;
Source3 = sum(sin(omega*t + phi_mt'*ones(1,length(t))),1)/Nx;
Source3(6001) = Source3(6001) + 10;
Source3(6500) = Source3(6500) + 10;

Source1 = bandpass(Source1,[fl fu],fs);
Source2 = bandpass(Source2,[fl fu],fs);
Source3 = bandpass(Source3,[fl fu],fs);
Source1 = sigma_i * Source1 / std(Source1);
Source2 = sigma_i * Source2 / std(Source2);
Source3 = Source3  / std(Source3);

sig1 = collector(Source1.',incidentAngle1);
sig2 = collector(Source2.',incidentAngle2);
sig3 = collector(Source3.',incidentAngle3);

signal    = sig1 + sig2 + sig3;
noise     = randn(size(signal));
noise     = sqrt(10^(-SNR/10)) * noise / std(noise);
recsignal = signal + noise;

Nframes  = length(recsignal(:,1)) - N + 1;
offset   = N/2 - (L-1)/2;
mic4_idx = 4;
soi_ref  = sig3(:, mic4_idx);

%% =========================================================
%  4a. PROPOSED ABSS
% =========================================================
yFill_up  = zeros(L,1);
yFill_low = zeros(L,1);
h_u       = zeros(L,1);  h_u((L-1)/2+1) = 1;
h_l       = zeros(L,1);  h_l((L-1)/2+1) = 1;
Pk0Old    = 0.1;

norm_prop = zeros(Nframes,1);
out_prop  = zeros(length(recsignal(:,1)),1);

for iLoop = 1:Nframes
    frame = recsignal(iLoop:iLoop+N-1,:)';

    y_up_beam  = sum(sum(Hup  .* frame, 2), 1);
    yFill_up   = circshift(yFill_up,  1);
    yFill_up(1) = y_up_beam;

    y_low_beam = sum(sum(Hlow .* frame, 2), 1);
    yFill_low  = circshift(yFill_low, 1);
    yFill_low(1) = y_low_beam;

    y_out = sum(h_u.*yFill_up) - sum(h_l.*yFill_low);
    out_prop(offset+iLoop) = y_out;

    Pk0    = 0.9*sum(yFill_up.*yFill_up) + 0.1*Pk0Old;
    Pk0Old = Pk0;
    h_l    = h_l + (beta_ABSS/(Pk0+eps)) * y_out * yFill_low;

    norm_prop(iLoop) = norm(h_l);
end
h_prop_final = h_l;
fprintf('Proposed converged: ||h_l|| = %.4f\n', norm(h_prop_final));

%% =========================================================
%  4b. GSC — store per-channel norms
% =========================================================
g_f      = zeros(L,1);  g_f((L-1)/2+1) = 1;
yDSB_buf = zeros(L,1);
H_gsc    = zeros(M-1,L);  H_gsc(:,(L-1)/2+1) = 1;
xB_buf   = zeros(M-1,L);

norm_gsc_ch = zeros(Nframes, M-1);   % ||h_j(n)|| for each channel j
out_gsc     = zeros(length(recsignal(:,1)),1);

for iLoop = 1:Nframes
    frame = recsignal(iLoop:iLoop+N-1,:)';

    x_align  = sum(Hdsb .* frame, 2);
    y_dsb    = sum(x_align, 1);
    yDSB_buf = circshift(yDSB_buf, 1);
    yDSB_buf(1) = y_dsb;
    y_upper  = sum(yDSB_buf .* g_f);

    x_B_n  = M * B_gsc * x_align;
    xB_buf = circshift(xB_buf, 1, 2);
    xB_buf(:,1) = x_B_n;

    y_lower = sum(sum(H_gsc .* xB_buf));
    y_gsc   = y_upper - y_lower;
    out_gsc(offset+iLoop) = y_gsc;

    denom  = sum(sum(xB_buf.*xB_buf, 2), 1);
    lambda = (M-1)*beta / (denom + eps);
    H_gsc  = H_gsc + lambda * y_gsc * xB_buf;

    % store norm of each channel separately
    for j = 1:M-1
        norm_gsc_ch(iLoop,j) = norm(H_gsc(j,:));
    end
end
H_gsc_final = H_gsc;
fprintf('GSC converged:     ||H||_F = %.4f\n', norm(H_gsc_final,'fro'));

%% =========================================================
%  5. FIGURE: Weight comparison (2 panels)
% =========================================================
c0       = (L-1)/2 + 1;
half_win = 50;
idx_show = c0-half_win : c0+half_win;
tap_ax   = -half_win : half_win;
cmap     = lines(M-1);

figure('Name','Adaptive Weights: Proposed vs GSC',...
    'Units','normalized','Position',[0.04 0.08 0.90 0.55]);

%-- Panel 1: Converged weight vectors in tap domain --------
subplot(1,2,1);
plot(tap_ax, real(h_prop_final(idx_show)), 'g-', 'LineWidth',2.2);
hold on;
for j = 1:M-1
    plot(tap_ax, real(H_gsc_final(j,idx_show)), '--',...
        'Color', cmap(j,:), 'LineWidth',1.0);
end

xlabel('Tap index');
ylabel('Amplitude');
title('\bf Converged weights h  (centre \pm50 taps)');
leg1 = {['Proposed h_l  (1\times' num2str(L) ')']};
for j = 1:M-1
    leg1{end+1} = sprintf('GSC h_%d', j);
end
leg1{end+1} = 'GSC mean h';
legend(leg1, 'Location','best','FontSize',8);
grid on;  set(gca,'FontSize',11);

%-- Panel 2: Learning curve — h_l + each GSC h_j -----------
subplot(1,2,2);

% Proposed: single green line (absolute norm)
plot(norm_prop, 'g-', 'LineWidth',2.2);
hold on;

% GSC: one coloured dashed line per channel h_j
leg2 = {sprintf('Proposed h_l  (L=%d)', L)};
for j = 1:M-1
    plot(norm_gsc_ch(:,j), '--', 'Color', cmap(j,:), 'LineWidth',1.3);
    leg2{end+1} = sprintf('GSC h_%d', j);
end

xlabel('Frame index n');
ylabel('||\mathbf{h}_j(n)||');
title('\bf Learning Curve: Proposed h_l  vs  GSC h_1 \ldots h_6');
legend(leg2, 'Location','best','FontSize',8);
grid on;  set(gca,'FontSize',11);

sgtitle(sprintf(['Adaptive Weight Comparison  |  '...
    'Proposed: 1\\times%d DOF   vs   GSC: %d\\times%d = %d DOF  '...
    '|  SIR=%ddB,  \\beta=%.2f'], ...
    L, M-1, L, (M-1)*L, iSIR_val, beta), ...
    'FontSize',11, 'FontWeight','bold');
set(gcf,'color','w');

%% =========================================================
%  6. FIGURE: Single zoomed plot around the pulse
% =========================================================
Nsig         = length(recsignal(:,1));
t_sig        = (0:Nsig-1) / fs * 1000;   % ms
mic4_sig     = recsignal(:, mic4_idx);
pulse_sample = 6001;
half_win_sig = 200;
idx_pulse    = max(1, pulse_sample-half_win_sig) : ...
               min(Nsig, pulse_sample+half_win_sig);
t_pulse_ms   = t_sig(idx_pulse);

figure('Name','Signal around pulse: Mic4 vs ABSS vs GSC',...
    'Units','normalized','Position',[0.10 0.20 0.75 0.45]);

plot(t_pulse_ms, mic4_sig(idx_pulse), 'Color',[0.6 0.6 0.6], 'LineWidth',1.2); hold on;
plot(t_pulse_ms, soi_ref(idx_pulse),  'k-',  'LineWidth',1.5);
plot(t_pulse_ms, out_prop(idx_pulse), 'g-',  'LineWidth',1.8);
plot(t_pulse_ms, out_gsc(idx_pulse),  'r--', 'LineWidth',1.8);
xline(pulse_sample/fs*1000, 'b:', 'LineWidth',1.2, 'HandleVisibility','off');

xlabel('Time (ms)');
ylabel('Amplitude');
legend(sprintf('Mic %d (received)', mic4_idx), ...
       'SOI ref (clean)', 'Proposed ABSS', 'GSC', ...
       'Location','best','FontSize',11);
title(sprintf('Signal around pulse (sample %d)  |  SIR=%ddB, SNR=%ddB, \\beta_{GSC}=%.2f, \\beta_{ABSS}=%.2f',...
    pulse_sample, iSIR_val, SNR, beta, beta_ABSS));
grid on;
set(gca,'FontSize',12);
set(gcf,'color','w');

%% ---- oSINR + weight summary ----------------------------
idx_eval   = 800:min(6800,Nsig);
SOI_std    = std(soi_ref(idx_eval));
oSINR_abss = 20*log10(SOI_std/(eps+std(out_prop(idx_eval)-soi_ref(idx_eval))));
oSINR_gsc  = 20*log10(SOI_std/(eps+std(out_gsc(idx_eval) -soi_ref(idx_eval))));

fprintf('\n======= Signal Quality =======\n');
fprintf('oSINR  Proposed ABSS : %.2f dB\n', oSINR_abss);
fprintf('oSINR  GSC           : %.2f dB\n', oSINR_gsc);
fprintf('\n======= Adaptive Weight DOF =======\n');
fprintf('                Proposed        GSC\n');
fprintf('Vectors         1               %d\n',     M-1);
fprintf('Taps/vector     %-7d         %-7d\n',      L, L);
fprintf('Total DOF       %-7d         %-7d\n',      L, (M-1)*L);
fprintf('||h|| final     %-10.4f      %-10.4f\n',...
    norm(h_prop_final), norm(H_gsc_final,'fro'));
fprintf('====================================\n');