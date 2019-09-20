% initialization of  variables 

f=rng;
r_min=0.1; % minimum position to be estimated
radius=7;  % stimulus radius
sigma_init=0.05;

nparam=4;
n_samples=10000;
how_beta=0; % 0=estimate beta using classical glm or 1=get a posterior distribution of beta
burn_in=0; % 0= apply burn-in, 1=no burn-in
how_pRF=0; % 0=classical pRF; 1=elliptical; 2=DoG.
n=10; %10
X=linspace(-radius,radius,1000);
[x, y] = meshgrid(linspace(-radius,radius,101)); % define visual space grid
accepted= false(1,n_samples);
% init latent parameters

l_sigma=1;
l_beta=1;

hrf = makeHRF(0:1.5:24);

l_rho=-0.5;
l_theta=0;
n_s_sp=50;
n_start=1;



