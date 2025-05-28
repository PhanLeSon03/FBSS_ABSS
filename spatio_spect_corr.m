function Sxx = spatio_spect_corr(beta,x)
%function Sxx = spatio_spect_corr(beta,x)
%
% compute spatio-spectral correlation matrix of 1d (linear) array 
% with mic. positions x, beta = 2*pi*f/vs (wave number)
% (isotropic noise field, spherical mic radiation pattern)

% compute distance matrix of all microphones

x = x(:);
N = length(x);
x = x(:,ones(N,1));
dx = x - x.';
dR = sqrt(dx.^2);
%Sxx = sinc(beta/pi*dR);      % spatio-spectral correlation matrix
 Sxx = (sin(beta*dR)./(beta*dR));
 Sxx(logical(eye(size(Sxx)))) = 1;
end 