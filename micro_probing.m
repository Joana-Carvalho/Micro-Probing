function [data_bayes_pRF] = micro_probing(stimulus, tSeries_data, k , j, i)
% runs the MP pipeline using the following input arguments:
%stimulus: binary representation of the stimulus (1 frame per TR)
%tSeries_data: BOLD tseries per voxel, hemisphere, condition and ROI
%k: hemisphere
%j:condition
%i: ROI
variables_definition

display(['k : ' num2str(k)]); 
display(['j : ' num2str(j)]);
display(['i : ' num2str(i)]);

for u=1:size(tSeries_data{k,j}{1,i},2)
    variables_definition_v2
    n_start=1;
    Y=tSeries_data{k,j}{1,i}(:,u);
    size(Y)
    
    stimulus_f=stimulus(:,:,1:size(Y,1));

    Y=Y./max(Y);
    varY = var(Y);% compute variance of Y
    err = varY;% set error matrix (U) to variance of Y
    varBase = ones(size(Y,1),1);% make vector of all ones for variable baseline
    
    
    [a b c]=size(stimulus_f);
    
    stimulus2d = reshape(stimulus_f,[a*b c])';

    stim_time = 1:1.5:size(stimulus2d,1)*1.5;
    
    
    if exist('posterior','var');
        [nparam,col]=size(posterior); % allow to run the script twice to just double the number of samples.
        posterior=cat(2,posterior,zeros(nparam,n_samples));
    else
        posterior=zeros(nparam,n_samples);%each row will contain one parameter each collumn one (accepted) simulation result
        col=1;
    end
    
    if exist('posterior_latent','var');
        [nparam,col]=size(posterior_latent); % allow to run the script twice to just double the number of samples.
        posterior_latent=cat(2,posterior_latent,zeros(nparam,n_samples));
    else
        posterior_latent=zeros(nparam,n_samples);%each row will contain one parameter each collumn one (accepted) simulation result
        col=1;
    end
    
    
    
    
    for kk=1:n_samples-1
        %
         if rem(kk,n_s_sp)==0
            n_start=n_start+1;
            l_rho=rand(1);
            l_theta=randn(1);
        end
        
        proposal_width=abs(normrnd(2, 0.5));
        
        
        if how_beta==0
            
            l_rho_proposal=normrnd(l_rho, proposal_width);
            l_theta_proposal=normrnd(l_theta, proposal_width);
            
            [p_current,xi_current,var_current, fitted_tseries_current] = computing_mcmc_tiny_ica (Y,x,y,radius,hrf,stimulus2d, l_rho,l_theta, sigma_init,r_min,l_beta,how_beta);
            [p_proposal,xi_proposal,var_proposal, fitted_tseries_proposal] = computing_mcmc_tiny_ica (Y,x,y,radius,hrf,stimulus2d,l_rho_proposal,l_theta_proposal,sigma_init,r_min,l_beta,how_beta);
            
            
            pvalue_current(kk)=p_current;
            pvalue_proposal(kk)=p_proposal;
            var_curr(kk)=var_current;
            var_pro(kk)=var_proposal;
            
            posterior_allvalues(:,kk)=xi_proposal;
            
            p_accept(kk)=exp(p_proposal-p_current);
            accept = normrnd(0,1);
            accepted(kk)=abs(accept)<p_accept(kk);
            
            if accepted(kk)==1
                
                l_rho = l_rho_proposal;
                l_theta=l_theta_proposal;
                %l_sigma=l_sigma_proposal;
                
                pstore(kk)=p_proposal;
                var_u(kk)=var_proposal;
                posterior_latent(:,col+kk-1)=[l_rho;l_theta; sigma_init;l_beta];
                posterior(:,col+kk-1)=xi_proposal;
                
                fitted_tseries_all(:,kk)= fitted_tseries_proposal;
                
            else
                pstore(kk)=p_current;
                var_u(kk)=var_current;
                posterior_latent(:,col+kk-1)=[l_rho;l_theta; sigma_init;l_beta];
                posterior(:,col+kk-1)=xi_current;
                
                fitted_tseries_all(:,kk)= fitted_tseries_current;
                
                
            end
            
        elseif how_beta==1
            
            l_rho_proposal=normrnd(l_rho, proposal_width);
            l_theta_proposal=normrnd(l_theta, proposal_width);
            l_beta_proposal=normrnd(l_beta, proposal_width);
            
            [p_current,xi_current,var_current, fitted_tseries_current] = computing_mcmc_tiny_ica (Y,x,y,radius,hrf,stimulus2d, l_rho,l_theta, sigma_init,r_min,l_beta,how_beta);
            [p_proposal,xi_proposal,var_proposal, fitted_tseries_proposal] = computing_mcmc_tiny_ica (Y,x,y,radius,hrf,stimulus2d,l_rho_proposal,l_theta_proposal,sigma_init,r_min,l_beta,how_beta);
            
            
            pvalue_current(kk)=p_current;
            pvalue_proposal(kk)=p_proposal;
            var_curr(kk)=var_current;
            var_pro(kk)=var_proposal;
            
            p_accept(kk)=exp(p_proposal-p_current);
            accept = normrnd(0,1);
            accepted(kk)=abs(accept)<p_accept(kk);
            
            posterior_allvalues(:,kk)=xi_proposal;
            
            if accepted(kk)==1
                
                l_rho = l_rho_proposal;
                l_theta=l_theta_proposal;
                l_beta=l_beta_proposal;
                
                pstore(kk)=p_proposal;
                var_u(kk)=var_proposal;
                posterior_latent(:,col+kk-1)=[l_rho;l_theta; sigma_init;l_beta];
                posterior(:,col+kk-1)=xi_proposal;
                
            else
                pstore(kk)=p_current;
                var_u(kk)=var_current;
                posterior_latent(:,col+kk-1)=[l_rho;l_theta; sigma_init;l_beta];
                posterior(:,col+kk-1)=xi_current;
                
            end
        end
        
    end
    
    [p_max{u},p_avg{u},varExpl{u},pstore_b{u},posterior_latent_b{u},posterior_b{u}, varExpl_all_subs{u}, rho_all{u}, theta_all{u}, fitted_tseries_all_bi{u}] = compute_burn_in_tiny_block_ica(Y,var_u,posterior_latent,posterior,pstore,radius,r_min,n,how_beta,burn_in,how_pRF, sigma_init, n_s_sp, n_start, fitted_tseries_all);
    
    
    
    data_bayes_pRF = struct('p_max', p_max, 'p_avg', p_avg, 'varExpl', varExpl, 'pstore_b', pstore_b, 'posterior_latent_b', posterior_latent_b, 'varExpl_all_subs', varExpl_all_subs, 'rho_all', rho_all, 'theta_all', theta_all);
    %  data_fitted_tseries=struct('fitted_tseries_all', fitted_tseries_all_bi);
    clear var_u; clear posterior_latent; clear posterior; clear pstore; clear fitted_tseries_all;
    %     clear p_max; clear varExpl; clear pstore_b; clear varExpl_all_subs;
    %     clear sigma_all; clear x0_all; clear y0_all; clear posterior_latent_b;
    
end

