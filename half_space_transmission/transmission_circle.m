clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

%% kernel parameters
omega = 1+sqrt(2);
ep0 = 1; mu0 = 3; zk0 = omega*sqrt(mu0*ep0); % exterior
ep1 = 2; mu1 = 4; zk1 = omega*sqrt(mu1*ep1); % interior

%% geometry parameters and construction

cparams = [];
cparams.nover = 4;
cparams.ifclosed = 1;
cparams.eps = 1.0e-9;

pref = []; 
pref.k = 16;

circle = @(t,ctr,rad) chnk.curves.bymode(t,rad*ones(1,1),ctr);

chnkr = chunkerfunc(@(t) circle(t,[0;0],1),cparams,pref); 
chnkr.n = chnk.normal2d(chnkr);

%% boundary data generation
% 1) generate the scattered field on the interior and exterior.
% 2) calculate the data of scattered field on the boundary. 
% 3) calculate the incident field on the boundary. 
% 4) solve the BIE to get the densities for the boundary integral representation
%       for the scattered fields
% 5) varify the solution. 

% initialize sources

ns = 1;
ts0 = 2*pi*rand(ns,1);
ts1 = 2*pi*rand(ns,1);
sources0 = circle(ts0,[2;2],1);
sources1 = circle(ts0,[0;0],0.2);
strengths1 = randn(ns,1);
strengths0 = randn(ns,1);

% plot sources

% scatter(sources0(1,:),sources0(2,:),'o')
% hold on
% scatter(sources1(1,:),sources1(2,:),'x')
% axis equal 


% eval incident field on the boundary

% exterior

srcinfo0 = []; srcinfo0.r = sources0;

kerns0 = @(s,t) chnk.helm2d.kern(zk0,s,t,'s');
kernsprime0 = @(s,t) chnk.helm2d.kern(zk0,s,t,'sprime');
f_u0 = @(targs) kerns0(srcinfo0,targs)*strengths0;
f_u0n = @(targs) kernsprime0(srcinfo0,targs)*strengths0;

u0bdr = f_u0(chnkr);
u0nbdr = f_u0n(chnkr);

% interior

srcinfo1 = []; srcinfo1.r = sources1;

kerns1 = @(s,t) chnk.helm2d.kern(zk1,s,t,'s');
kernsprime1 = @(s,t) chnk.helm2d.kern(zk1,s,t,'sprime');
f_u1 = @(targs) kerns1(srcinfo1,targs)*strengths1;
f_u1n = @(targs) kernsprime1(srcinfo1,targs)*strengths1;

u1bdr = f_u1(chnkr);
u1nbdr = f_u1n(chnkr);

% the incident field at the boundary,

% u_inc = u1bdr - u0bdr;
u_inc = u1bdr;

% un_inc = ep0/ep1 * u1nbdr - u0nbdr;
un_inc = ep0/ep1 * u1nbdr;
rhs = [-u_inc; -un_inc/ep0];

%% generating and solving the BIE

% The Boundary Integral Equation is 
% (ep1+ep0)/2 * rho + (ep0 D0 - ep1D1) rho + (ep0^2 S0 - ep1^2 S1) sigma 
%   = -u_inc
% (ep1+ep0)/2 * sigma - (D'0 - D'1) rho - (ep0 S'0 - ep1 S'1) sigma
%   = -un_inc / ep0

[~, npts] = size(chnkr.r);

A = zeros(2*npts,2*npts,'like', 1j);

fkern_d_diff = @(s,t) ep0*chnk.helm2d.kern(zk0,s,t,'D')...
    - ep1*chnk.helm2d.kern(zk1,s,t,'D');
fkern_s_diff = @(s,t) ep0^2 * chnk.helm2d.kern(zk0,s,t,'S')...
    - ep1^2 * chnk.helm2d.kern(zk1,s,t,'S');
fkern_dprime_diff = @(s,t) chnk.helm2d.kern(zk0,s,t,'dprime')...
    - chnk.helm2d.kern(zk1,s,t,'dprime');
kern_sprime_diff = @(s,t) ep0*chnk.helm2d.kern(zk0,s,t,'sprime') ...
    - ep1*chnk.helm2d.kern(zk1,s,t,'sprime');

A(1:npts,1:npts) = chunkermat(chnkr,fkern_d_diff) ...
    + 0.5*(ep0+ep1)*eye(npts,'like',1j);
A(1:npts,npts+1:2*npts) = chunkermat(chnkr,fkern_s_diff);
A(npts+1:2*npts,1:npts) = -chunkermat(chnkr,fkern_dprime_diff);
A(npts+1:2*npts,npts+1:2*npts) = -chunkermat(chnkr,kern_sprime_diff) ...
    + 0.5*(ep0+ep1)*eye(npts,'like',1j);


[x,flag,relres,iter,resvec] = gmres(A,rhs,1,1e-6,50);

rho = x(1:npts);
sigma = x(npts+1:2*npts);


%% verifying the correctness of the solution

%% maybe first I should try to evaluate the solutions on the boundary

% evaluating the solutions
fd0 = @(s,t) ep0 * chnk.helm2d.kern(zk0,s,t,'D');
fs0 = @(s,t) ep0^2 * chnk.helm2d.kern(zk0,s,t,'S');

fd1 = @(s,t) ep1 * chnk.helm2d.kern(zk1,s,t,'D'); 
fs1 = @(s,t) ep1^2 * chnk.helm2d.kern(zk1,s,t,'S');


h = 0.001;
uscatbdr_sol = chunkerkerneval(chnkr,fd0,rho,chnkr.r + h*chnkr.n)...
    + chunkerkerneval(chnkr,fs0,sigma,chnkr.r + h*chnkr.n);



%%

% generating the targets in the interior and exterior

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 100;
xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot); 
ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);
targetsinfo = [];
targetsinfo.r = targets;

in = chunkerinterior(chnkr,targets);
out = ~in;

% the exact evaluation of scattered field 
uscat_exact = f_u0(targetsinfo); uscat_exact = uscat_exact(out);
uin_exact = f_u1(targetsinfo);   uin_exact = uin_exact(in);

% evaluating the solutions
fd0 = @(s,t) ep0 * chnk.helm2d.kern(zk0,s,t,'D');
fs0 = @(s,t) ep0^2 * chnk.helm2d.kern(zk0,s,t,'S');

fd1 = @(s,t) ep1 * chnk.helm2d.kern(zk1,s,t,'D'); 
fs1 = @(s,t) ep1^2 * chnk.helm2d.kern(zk1,s,t,'S');

uscat_sol1 = chunkerkerneval(chnkr,fd0,rho,targets(:,out));
uscat_sol2 = chunkerkerneval(chnkr,fs0,sigma,targets(:,out));
uscat_sol = uscat_sol1 + uscat_sol2;

uin_sol1 = chunkerkerneval(chnkr,fd1,rho,targets(:,in));
uin_sol2 = chunkerkerneval(chnkr,fs1,sigma,targets(:,in));
uin_sol = uin_sol1 + uin_sol2;



%% plotting

maxin = max(abs(uscat_exact(:)));
maxsc = max(abs(uin_exact(:)));
maxtot = max(abs(uin_exact(:)));
maxu = max(max(maxin,maxsc),maxtot);

figure(1)
clf

subplot(1,2,1)
zztarg = nan(size(xxtarg));
zztarg(out) = uscat_exact;
zztarg(in) = uin_exact;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'LineWidth',2)
axis equal
axis tight
colormap(redblue)
clim([-maxu,maxu])
title('$u_{scat}$ exact','Interpreter','latex','FontSize',24)

colorbar


subplot(1,2,2)
zztarg = nan(size(xxtarg));
zztarg(out) = uscat_sol;
zztarg(in) = uin_sol;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'LineWidth',2)
axis equal
axis tight
colormap(redblue)
clim([-maxu,maxu])
title('$u_{scat}$ soln','Interpreter','latex','FontSize',24)

colorbar
