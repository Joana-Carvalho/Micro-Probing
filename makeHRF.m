function hrf = makeHRF(t,p)
% makeHRF - make HRF as a difference of two gamma functions
%
% hrf = makeHRF(time_vector,params_vector)
%
% Input:
%    time_vector: a range of latencies
%    params_vector: five parameters defining the two gamma functions
%                   default: [5.4 5.2 10.8 7.35 0.35]
%       params(1): peak gamma 1
%       params(2): fwhm gamma 1
%       params(3): peak gamma 2
%       params(4): fwhm gamma 2
%       params(5): dip
%                  
% Final hrf is:   gamma1/max(gamma1)-dip*gamma2/max(gamma2)
%
% from Glover, NeuroImage, 9:416-429
%
% 2009/03 SD: wrote it.

if ~exist('p','var') || isempty(p),
    p = [5.4 5.2 10.8 7.35 0.35];
end;

% params
peak1 = p(1);
fwhm1 = p(2);
peak2 = p(3);
fwhm2 = p(4);
dip   = p(5);

% sanity check
if(peak1 == 0 || fwhm1 ==0),
    fprintf('[%s]: zero params',mfilename);
    return;
end;

% first gamma function:
alpha1=peak1^2/fwhm1^2*8*log(2);
beta1=fwhm1^2/peak1/8/log(2);
gamma1=(t/peak1).^alpha1.*exp(-(t-peak1)./beta1);

if peak2>0 && fwhm2>0
    % second gamma function:
    alpha2=peak2^2/fwhm2^2*8*log(2);
    beta2=fwhm2^2/peak2/8/log(2);
    gamma2=(t/peak2).^alpha2.*exp(-(t-peak2)./beta2);
else
    gamma2=min(abs(t-peak2))==abs(t-peak2);
end

% final hrf
hrf = gamma1-dip*gamma2;

return;