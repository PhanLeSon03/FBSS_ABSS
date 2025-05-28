function [R,theta,phi] = array_pattern_fft(mics,W,f,k)
% compute radiation pattern of 1 or 2 dim. microphone array
%
% mics      matrix of mic x-, y-coordinates, (z-coord. = 0)
% W         weight matrix of fixed beamformer
%           (rows are FFT coeff. for each channel)
% f         frequency in Hz (must correspond to frequency index kf of W)
% k         frequency index corresponding to f
% R         squared magnitude of rad. pattern
%           1 dim. array: R(phi)   
%           2 dim. array: R(theta,phi)   
% phi       azimuth vector in radians
% theta     elevation vector in radians (pi/2, if 1 dim. array)

% Db, 24.1.02, 18.9.06

   vs = 340;
   beta = 2*pi*f/vs;                  % wave number

   phi = pi/180*(0:1:360);
   
   [N,K] = size(mics);
   
   if (K == 1)             % 1-dim Array?
      theta = pi/2;
      D = exp(j*beta*mics*cos(phi));        % matrix of steering vectors
      R = abs(W(:,k)'*D).^2;
   else
      theta = pi/180*(0:2.5:90);
      Nphi = length(phi);
      Ntheta = length(theta);
      U = j*beta*mics;
      V =[cos(phi) ; sin(phi)];  
      R = zeros(Ntheta,Nphi);
      for m = 1:Ntheta
         er = sin(theta(m))*V;
         D = exp(U*er);               % matrix of steering vectors
         R(m,:) = abs(W(:,k)'*D).^2;
      end
   end
