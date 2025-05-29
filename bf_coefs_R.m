
function W = bf_coefs_R(mics,theta_d,phi_d,resp,f,mue, null,R)
%
% compute fixed beamformer weight in frequency domain 
%
% mics     (x,y) coordinates of array (see mics_config.m)
% theta_d  elevations in deg. of desired direction followed by direction of nulls
% phi_d    azimuths in deg., of desired direction followed by direction of nulls 
% resp     beamformer response in desired direction followed by direction of nulls
% f        frequency vector in Hz at which W is computed
% mue      regularization parameter
% null     using spatial null or not
% W        N x Nhigh-Nlow+1 array of beamformer weights


    vs = 340;

    % compute matrices used for optimzation

    theta_d = theta_d.* pi / 180;
    phi_d = phi_d.* pi / 180;

    [N,K] = size(mics);
    if (K < 1) || (N < 1)
       error('bad microphone positions');
    end

    if K == 2                      % 2 dim. array 
       mics(:,1) = mics(:,1);  
       mics(:,2) = mics(:,2); 
       rn = [mics zeros(N,1)];
    else
       rn = [mics zeros(N,1) zeros(N,1)];
    end

    Ntheta = length(theta_d);
    if (Ntheta ~= length(phi_d))
       error('angle vectors must have same sizes');
    end

    er = [sin(theta_d).*cos(phi_d) ; sin(theta_d).*sin(phi_d) ; cos(theta_d)];  % steering vector
    Rc = rn*er;            % used to compute matrix C = exp(j*beta*rn*er)
                           % size of NxNtheta


% perform optimization for each frequency to find reference values
nf = length(f);
W = zeros(N,nf);

for l = 1:nf

   if null==1
       beta = 2*pi*f(l)/vs;     % wave number
       C = exp(1j*beta*Rc);      % constraints matrix, size of NxNtheta

       Amue = R + mue*eye(size(R));   % regularization, size of NxN
       B = Amue\C;                % size of NxNtheta
       Lambda = (C'*B)\resp(:);   % Lagrange multiplicator
       W(:,l) = B*Lambda;         % optimum coefficient vector at given frequency
   else
       beta = 2*pi*f(l)/vs;     % wave number
       C = exp(1j*beta*Rc);      % constraints matrix

       Amue = R + mue*eye(size(R));   % regularization
       B = Amue\C(:,1); % B = Rx^-1*a_d
       W(:,l) = B/(C(:,1)'*B);  

   end
end 


