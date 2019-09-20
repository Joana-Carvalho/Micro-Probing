function [ p,xi,varE, fitted_tseries] = computing_mcmc_tiny_ica (Y,x,y,radius,hrf,stimulus2d,l_rho,l_theta,sigma_init,r_min,l_beta,how_beta)
%
% p will be the likelyhood of the fit 
% xi will be the parameters as calculated from the latent variables. it
% will be a column vector containing: [rho;theta;sigma;beta]
%
%-------------------------------------------------------------------------
% Functions definition


rho=@(l_rho,radius) radius*normcdf(l_rho,0,1);

theta=@(l_theta) 2*pi.*normcdf(l_theta,0,1)-pi; 

beta=@(l_beta) exp(l_beta);

%sigma=@(l_sigma,radius_sigma, r_min)(radius_sigma-r_min).*normcdf(l_sigma,0,1)+r_min;

miu_x0=@(l_rho,radius,l_theta) rho(l_rho,radius).*cos(theta(l_theta));
       
miu_y0=@(l_rho,radius,l_theta) rho(l_rho,radius).*sin(theta(l_theta));

pRF =@(l_rho,radius,l_theta,r_min,sigma_init) exp(((x-miu_x0(l_rho,radius,l_theta)).^2+(y-miu_y0(l_rho,radius,l_theta)).^2)/(-2*sigma_init.^2));
    
pTime_series =@(pRF) (stimulus2d*pRF(:))';

%-------------------------------------------------------------------------
pRF_1=pRF(l_rho,radius,l_theta,r_min,sigma_init);

if how_beta==0;
    
    pTime_series_1 =pTime_series(pRF_1);
    
    % Convolve with hrf
    pTime_series_1 = conv(pTime_series_1,hrf); % Convolve
    pTime_series_1 = pTime_series_1(1:size(stimulus2d,1)); % Clip data to correct length
    
    % Estimate beta using classical GLM and OLS
    varBase = ones(size(Y,1),1);
    X = [pTime_series_1' varBase];
    B_hat = pinv(X)*Y;
    E= Y-(X*B_hat);
    fitted_tseries=X*B_hat;
    varE=var(E);
    
    
    xi=[rho(l_rho,radius);theta(l_theta);sigma_init;B_hat(1)];
               
elseif how_beta==1;
    
    pTime_series_1 =l_beta*pTime_series(pRF_1);
    pTime_series_1 = conv(pTime_series_1,hrf);  % Convolve
    pTime_series_1 = pTime_series_1(1:size(stimulus2d,1)); % Clip data to correct length
    pTime_series_1=detrend(pTime_series_1,'constant');
    Y_demean=detrend(Y,'constant');
    E=transpose(Y_demean)-pTime_series_1;
    fitted_tseries=pTime_series_1;
    varE=var(E);
   
    xi=[rho(l_rho,radius);theta(l_theta);sigma_init; l_beta];
    
end

% Find likelihood
[muhat,sigmahat] = normfit(E);
loglikelihood=log(normpdf(-abs(E),muhat,sigmahat)); 
prior_r=normpdf(l_rho, 0,1);
prior_t=normpdf(l_theta, 0,1);
%prior_s=normpdf(l_sigma, 0,1);

if how_beta==0
    %a=sum(loglikelihood);
    %p=sum(loglikelihood)+log(prior_r)+log(prior_t)+log(prior_s);
    
    p=sum(loglikelihood)+log(prior_r)+log(prior_t);
    
elseif how_beta==1
    prior_b=normpdf(l_beta,-2,5);
    %a=sum(loglikelihood);
    %p=sum(loglikelihood)+log(prior_r)+log(prior_t)+log(prior_s)+log(prior_b);
    
    p=sum(loglikelihood)+log(prior_r)+log(prior_t)+log(prior_b);
end

end
