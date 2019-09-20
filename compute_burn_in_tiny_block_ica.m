function [p_max,p_avg,varExpl,pstore,posterior_latent,posterior, varExpl_all, rho_all, theta_all, fitted_tseries_all] = compute_burn_in_tiny_block_ica (Y,var_u,posterior_latent,posterior,pstore,radius,r_min,n,how_beta,burn_in,how_pRF, sigma_init, n_s_sp,start_points,fitted_tseries_all)
%-------------------------------------------------------------------------
clear x0_all; clear y0_all;

% Variables definition
p_max = zeros(1,size(posterior_latent,1));
p_avg = zeros(1,size(posterior_latent,1));

% Functions definition
rho=@(l_rho,radius) radius*normcdf(l_rho,0,1);

theta=@(l_theta) 2*pi.*normcdf(l_theta,0,1)-pi;

beta=@(l_beta) exp(l_beta);

alpha=@(l_alpha)  2*pi.*normcdf(l_alpha,0,1)-pi;

%sigma=@(l_sigma,radius_sigma, r_min)(radius_sigma-r_min).*normcdf(l_sigma,0,1)+r_min;

mu_x0=@(l_rho,radius,l_theta) rho(l_rho,radius).*cos(theta(l_theta));

mu_y0=@(l_rho,radius,l_theta) rho(l_rho,radius).*sin(theta(l_theta));

%-------------------------------------------------------------------------
if burn_in==0
    % Burn-in - remove the first 1000 of elements on 10000 runs
    if size(posterior_latent,2)<n
        burn_in=1;
    else
%         
        n_burn=round(n_s_sp/n);
%         pstore_mat = reshape(pstore,[n_s_sp,start_points]);
%         pstore_mat=pstore_mat(n_burn+1:end,:);
%         pstore=reshape(pstore_mat,[1,size(pstore_mat,1)*size(pstore_mat,2)]);
          
          rem_index=linspace(1,n_burn,n_burn);
        for i=2:start_points-1
          rem_index=[rem_index linspace(n_s_sp*i+1,n_s_sp*i+n_burn,n_burn)];
        end
        
        
         pstore(rem_index)=[];
         posterior_latent(:,rem_index)=[];
         var_u(:,rem_index)=[];
        posterior(:,rem_index)=[];
        
        fitted_tseries_all(:,rem_index)=[];
        
        varExpl_all= 1 - var_u./var(Y);
        
        A=find(pstore==max(pstore));
        max_A=posterior_latent(:,A(1));
        b=mean(posterior_latent(:,A),2);
        max_position=A(1);
        varExpl = 1 - var_u(max_position)./var(Y);
        
        p_max(1)=rho(max_A(1),radius);
        p_max(2)=theta(max_A(2));
        p_max(3)=sigma_init;
        p_avg(1)=rho(b(1),radius);
        p_avg(2)=theta(b(2));
        p_avg(3)=sigma_init;
        
        
        
        rho_all=rho(posterior_latent(1,:),radius);
        theta_all=theta(posterior_latent(2,:));
        %sigma_all=sigma(posterior_latent(3,:),radius_sigma,r_min);
        
        
        if how_pRF==1
            p_max(4)=sigma(max_A(4),radius_sigma,r_min);
            p_avg(4)=sigma(b(4),radius_sigma,r_min);
            p_max(5)=alpha(max_A(5));
            p_avg(5)=alpha(b(5));
        end
        if how_pRF==2
            p_max(4)=sigma_init;
            p_avg(4)=sigma_init;
        end
        
        if how_beta==0
            p_max(end)=max_A(end);
            p_avg(end)=b(end);
        elseif how_beta==1
            p_max(end)=beta(max_A(end));
            p_avg(end)=beta(b(end));
        end
        
    end
    
elseif burn_in==1
    
    
    varExpl_all= 1 - var_u./var(Y);
    
    A=find(pstore==max(pstore));
    max_A=posterior_latent(:,A(1));
    max_position=A(1);
    b=mean(posterior_latent(:,A),2);
    max_position=A(1);
    
    
    
    
    varExpl = 1 - var_u(max_position)./var(Y);
    
    p_max(1)=mu_x0(max_A(1),radius,max_A(2));
    p_max(2)=mu_y0(max_A(1),radius,max_A(2));
    p_max(3)=sigma_init;
    p_avg(1)=mu_x0(b(1),radius,b(2));
    p_avg(2)=mu_y0(b(1),radius,b(2));
    p_avg(3)=sigma_init;
    
    if how_pRF==1
        p_max(4)=sigma(max_A(4),radius_sigma,r_min);
        p_avg(4)=sigma(b(4),radius_sigma,r_min);
        p_max(5)=alpha(max_A(5));
        p_avg(5)=alpha(b(5));
    end
    if how_pRF==2
        p_max(4)=sigma_init;
        p_avg(4)=sigma_init;
    end
    if how_beta==0
        p_max(end)=max_A(end);
        p_avg(end)=b(end);
    elseif how_beta==1
        p_max(end)=beta(max_A(end));
        p_max(end)=beta(b(end));
    end
    
end

end
