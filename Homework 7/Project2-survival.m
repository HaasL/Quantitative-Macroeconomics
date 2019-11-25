% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% project 2
% solution of household problem for T = \infty


clear all
close all

global betta tetta r g gridx vpfun epsi probepsi ne nx min_cons
tic;
% -------------------------------------------------------------------------
% SETTINGS
maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;
nj=80;
jr=45;

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 5;
betta = 1/(1+rho);

% grid
nx=nj;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';
sr1=readfile([],'MR.txt',3);
sr=sr1(21:end,2)

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
%[epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
%mat=[[1:ne]',epsi,probepsi];
%save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% initial guess
cons = gridx;               % consumption
vpfun = margutil(cons);     % derivative of value function

% iteration
for it=1:maxit,
    disp(['iteration # ', num2str(it)]);
    consm1=cons;
    
    for xc=nj-1:-1:1,
        % check binding constraint:
        mincons=gridx(xc);
        mu = foc(mincons,gridx(xc));
        if (mu>=0.0),
            cons(xc)=mincons;
        else,
            [cons(xc),fval] = fzero('foc',cons(xc),[],gridx(xc));
            cons(xc)=max(cons(xc),min_cons);
        end;
%          if (cons(xc)>gridx(xc))
%              cons(xc)=gridx(xc);
%          end;
    end;    

    % update vpfun
    vpfun = margutil(cons);
    
    % check convergence
    dist=cons-consm1;
    maxdist=max(abs(dist));
    if (maxdist<tol),
        disp(['I am so happy! The thing has converged in iteration ', num2str(it)]);
        break;
    else,
        disp(['current maximum norm is ', num2str(maxdist)]);
    end;
end;
if (it==maxit),
    warning('increase # of iters for specified tolerance');
end;

% -------------------------------------------------------------------------

figure;
plot(gridx,cons,'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy');

%[gridx,cons]



% -------------------------------------------------------------------------
% SIMULATION
ct = zeros(nt,1);
xt = zeros(nt,1);
at = zeros(nt,1);
yt = zeros(nt,1);
et = zeros(nt,1);
eulert = zeros(nt,1);

% random numbers:
rand('seed',0);
probet=rand(nt,1);
indet=ceil(probet*ne);

outx = 0;
for tc=1:nt,
    yt(tc) = epsi(indet(tc));
    xt(tc) = at(tc)+yt(tc);
    ct(tc) = func_intp(gridx,cons,xt(tc));
    chkoutx=false;
    if (ct(tc)==max(cons)),
        outx=outx+1;
        chkoutx=true;
    end;
    if (tc<nt),
        at(tc+1)=(xt(tc)-ct(tc))*(1+r)/(1+g);
    end;
    
    % error evaluation
    if (ct(tc)<xt(tc) && chkoutx==false),
        margu=margutil(ct(tc));
        et(tc)=abs(foc(ct(tc),xt(tc))/margu);
    end;
end;
% -------------------------------------------------------------------------

fracoutx=outx/nt;
%if ( routx>0.01 ),
%    beep; beep; beep;
%    warning('grid too small, enlarge your grid');
    disp(['fraction of points outside grid is ', num2str(fracoutx)]);
%end;

nnt=nt-dt+1;
disp(['maximum Euler equation error is : ', num2str(max(et(dt+1:nt)))]);
disp(['mean Euler equation error is : ', num2str(mean(et(dt+1:nt)))]);
I = (et==0);
rbc = sum(I(dt+1:nt))/nnt;
disp(['borrowing constraint binds in ', num2str(rbc), ' cases']);

figure;
plot([dt+1:nt],ct(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('c_t');
title('consumption over time');

figure;
plot([dt+1:nt],xt(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('x_t');
title('cash-on-hand over time');

figure;
plot([dt+1:nt],et(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('e_t');
title('euler error over time');
time=toc
disp(['Time is : ', num2str(time)]);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
