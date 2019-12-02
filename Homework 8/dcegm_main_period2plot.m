function dcegm_main
% Illustrates solution to a simple 3-period model with discrete labor
% decision in the last period (see write-up)
% Timing: in the code initial period is t=1 (refers to t=0 in write-up)

clc
clear all
close all

% options
opt.plt      = 1;               % option for plotting graphs
opt.pr       = 0;               % option for printing of graphs
opt.meth     = 3;               % 1 - global solution, 2 - DC-EGM, 3 - exogm (rootfinder)
opt.comp     = 1;               % compare solution methods (global vs. dc-egm or exogm)
opt.det      = 0;               % deterministic problem (this option is relevant only for global solution => it solves then a static problem; not relevant for dc-egm)

% deep parameters
sig2e        = sqrt(eps);             % variance of wage shocks (when increasing the variance, make sure to set grids large enough )
param.sigma  = sqrt(eps);             % variance of taste shocks (closed form solution available only for sigma=1) 
param.phi    = 1.0;             % preference for leisure   - DO NOT CHANGE (if want comparison with analytical solution)
param.Gam    = 2.0;             % time endowment           - DO NOT CHANGE (if want comparison with analytical solution)

% computational parameters
param.nm     = 3;               % number of solution methods available
param.na     = 41;              % gridsize assets/coh
param.ne     = 20;              % gridsize shocks
param.nc     = 20;              % gridsize consumption (for global solution)
param.nd     = 2;               % number of discrete states/choices
param.nt     = 2;               % number of periods minus one (i.e. number of periods with continuous decisions; discrete decision ALWAYS ONLY in the last period)
param.tol    = 1.0e-08;         % tolerance criterion
param.curv   = 3.0;             % curvature parameter for grid
param.np1aux = 6;               % number of gridpoints+1 in bc region for dcegm
max_ass      = 10.0;             % maximium assets
min_ass      = 0.0; %1.0e-02; 

% -------------------------------------------------------------------------

% correction of options
if (opt.det==1) % if deterministic option (static problem) is chosen, variance of shocks should be zero)
    opt.meth=1;
    sig2e=0.0;
    param.sigma=0.0;
end;

if (sig2e>param.tol)  
    if (param.sigma<param.tol),
        warning('if income uncertainty is present, then solution with DC-EGM should be done WITH TASTE SHOCK');
        warning('reason: may get spurious discountinuities in policies, cf. Iskhakov et al 2017 for further details');
        ans = input('PRESS 1 TO CONTINUE ');
        if (ans~=1)
            return
        end;
    end;
end;
    
% draw shocks
if (sig2e<=sqrt(eps))
    param.eta=zeros(param.ne,1);
    param.wght=zeros(param.ne,1);
    param.eta(:)=1.0;
    param.wght(:)=param.ne^(-1);
else
    [eta,param.wght] = qnwnorm(param.ne,-sig2e/2.0,sig2e);
    param.eta = exp(eta);
    mu1_eta = param.wght'*param.eta;
    dist = abs(mu1_eta-1.0);
    if (dist>param.tol)
        error('something is wrong with the shocks');
    end
end

% Preallocate grids and policies (dimensions are: na - assets/coh, nt - time, nm - solution method)
% grids
grid.coh=zeros(param.na,param.nt,param.nm);
grid.coh_anal=zeros(param.na,1); % for illustrating analytical solution
grid.ass=zeros(param.na,1);
grid.sav=zeros(param.na,param.nt); % age-dependent savings grid (egm)
% cons, sav policy and value fun
fun.cons=zeros(param.na,param.nt,param.nm);
fun.val=zeros(param.na,param.nt,param.nm);
fun.sav=zeros(param.na,param.nt,param.nm);
fun.assets=zeros(param.na,param.ne,param.nt+1); % current period assets

% set up exogenous grids
grid.ass(:) = makegrid(min_ass,max_ass,param.na,param.curv); % exogenous assets 
param.min_coh = min_ass + min(param.eta);
param.max_coh = max_ass + max(param.eta);

% coh grid for illustrating analytical solution
grid.coh_anal(:) = makegrid(param.min_coh,param.max_coh,param.na,param.curv);

% exogenous coh for global solution 
for tc=1:param.nt
    % for global method, keep coh grid constant
    grid.coh(:,tc,1) = makegrid(param.min_coh,param.max_coh,param.na,param.curv);
    % same grid for exogm with a rootfinder
    grid.coh(:,tc,3) = grid.coh(:,tc,1);
end

% Numerical solution of dynamic program (with a chosen method)
tic;
[fun,grid] = fun_solve(fun,grid,param,opt);
disp(['time elapsed for method ', num2str(opt.meth), ' is: ', num2str(toc)]);

% Closed form solution of dynamic program (fully deterministic + var of taste shock = 1)
if ( sig2e<=sqrt(eps) && ( param.sigma<=sqrt(eps) || (param.sigma==1.0) ) ),
    
    fun = fun_closed(grid,fun,param);
    
    % one cannot do this comparison bec the objects do not live on the same grid,
    % first one would have to interpolate
        % one cannot do this comparison bec the objects do not live on the same grid,
    % first one would have to interpolate
    c0_intp=zeros(param.na,1);
    c1_intp=zeros(param.na,1);
    for ac=1:param.na
        c0_intp(ac) = interp1(grid.coh(:,1,opt.meth),fun.cons(:,1,opt.meth),grid.coh_anal(ac,1),'linear');
        c1_intp(ac) = interp1(grid.coh(:,2,opt.meth),fun.cons(:,2,opt.meth),grid.coh_anal(ac,1),'linear');
        
    end
    dist_c0=abs(fun.c0_closed./c0_intp-1.0);
    dist_c1=abs(fun.c1_closed./c1_intp-1.0); 
    max_dist_c0=max(dist_c0);
    max_dist_c1=max(dist_c1);
    disp(['maximum distance from closed form solution: ', num2str(max_dist_c0)]);
    disp(['maximum distance from closed form solution: ', num2str(max_dist_c1)]);
    
    if (opt.plt)
        xla='$x$';
        yla='$c$';
        leg={'$c_{0}$ - analytical','$c_{0}$ - numerical'};
        tit=['comparision of c(0): analytical and numerical solution with $\sigma$ = ',num2str(param.sigma)] ;
        outpath=[];
        grname='c0_anal_num';
        ax=[];
        onesideplt_marker([grid.coh_anal,grid.coh(:,1,opt.meth)],[fun.c0_closed,fun.cons(:,1,opt.meth)],xla,yla,leg,tit,outpath,grname,ax,opt.pr)
        
        
        xla='$x$';
        yla='$c$';
        leg={'$c_{1}$ - analytical','$c_{1}$ - numerical'};
        tit=['comparision of c(1): analytical and numerical solution with $\sigma$ = ',num2str(param.sigma)] ;
        outpath=[];
        grname='c1_anal_num';
        ax=[];
        onesideplt_marker([grid.coh_anal,grid.coh(:,2,opt.meth)],[fun.c1_closed,fun.cons(:,2,opt.meth)],xla,yla,leg,tit,outpath,grname,ax,opt.pr)
    end
    
end

% Compare solution methods
if (opt.comp==1 && opt.det==0)
    % solve numerical problem with an alternative method (if method 1 => compare with method 2; if method 2/3 => compare with method 1)
    if (opt.meth==1)
        opt.meth=2;
        opt.meth_comp=2;
    else
        opt.meth_comp=opt.meth;
        opt.meth=1;
    end
    
    tic;
    [fun,grid] = fun_solve(fun,grid,param,opt);
    disp(['time elapsed second method: ', num2str(toc)]);
    
    % plot consumption policy for period 0 (in the code period 1) obtained
    % from two solution methods
    if (opt.plt)
        xla='$x$';
        yla='$c$';
        leg={'$c_{0}$ - global ',['$c_{0}$ - method ', num2str(opt.meth_comp)]};
        tit=['comparision: solution from two methods with $\sigma$ = ',num2str(param.sigma)] ;
        outpath=[];
        grname=['c0_global_meth', num2str(opt.meth_comp)];
        ax=[];
        onesideplt_marker([grid.coh(:,1,opt.meth),grid.coh(:,1,opt.meth_comp)],[fun.cons(:,1,opt.meth),fun.cons(:,1,opt.meth_comp)],xla,yla,leg,tit,outpath,grname,ax,opt.pr)
    end
end    

% compute cross-sectional measure and aggregate variables
if (opt.det==0)
    [grid,agg]=fun_aggr(grid,fun,param,opt);
    disp(['aggregate coh, cons and assets: ', num2str([agg.coh,agg.cons,agg.ass])]);
end;

disp(' ');  
disp('I am so incredibly happy that I am done!');
disp(' ');

end     % end function snapshot_dynprg_dcegm_3per
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [fun] = fun_closed(grid,fun,param)

    disp('get closed form solutions for comparison');

    fun.c0_closed=zeros(param.na,1);
    fun.c1_closed=zeros(param.na,1);
        
    if (param.sigma<=sqrt(eps))  % option w/o taste shocks
        
        Gam_tilde=param.Gam/(param.Gam-1.0);
        temp=Gam_tilde^param.phi;
        param.temp=temp;
    
        func=@fun_thrs;
        w0=1.0;
        thres_w0=fzero(func,w0,[],param);
        thres_x0=thres_w0-1.0;

        for ac=1:param.na
            
            if (param.phi<eps)
                fun.c0_closed(ac)=2/3+grid.coh_anal(ac)/3;
                %fun.c1_closed(ac)=1/3+grid.coh_anal(ac)/3;
               
            else
            
                if (grid.coh_anal(ac)<thres_x0)
                    fun.c0_closed(ac)=2/3+grid.coh_anal(ac)/3;
                    %fun.c1_closed(ac)=1/3+grid.coh_anal(ac)/3;
                else
                    fun.c0_closed(ac)=1/3+grid.coh_anal(ac)/3;
                    %fun.c1_closed(ac)=2/3+grid.coh_anal(ac)/3;
                end
            end
            fun.a1_closed=grid.coh_anal-fun.c0_closed;
            fun.a2_closed=grid.coh_anal+1.0-2.0*fun.c0_closed;
        fun.c1_closed = 1.0+2.0*fun.c0_closed
        end
    
        
    elseif (param.sigma == 1.0) % option w/ taste shocks, assuming sigma=1
        
        fun.a1_closed = 1/3 * grid.coh_anal - 4/9;
        fun.a2_closed = 1/3 * grid.coh_anal + 1/9;
        fun.c0_closed = 1/3 * grid.coh_anal + 4/9;
        fun.c1_closed = 1/3 * grid.coh_anal + 4/9;
        fun.prob_closed = (10 + 3.0*grid.coh_anal)./(12 + 9.0*grid.coh_anal);
        
    end    

end   % end function fun_closed
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [fun,grid] = fun_solve(fun,grid,param,opt)
% solution of dynamic program

disp('solve household problem numerically');

if (opt.det==1) % solve a static problem
    % evaluate value function on grid:
    val_cons = zeros(param.nc,1);
    min_cons = param.tol;
    for ac = 1:param.na
        max_cons  = grid.coh(ac,1,1) + min(param.eta) - param.tol;     % maximum possible current period consumption 
        grid_cons = makegrid(min_cons,max_cons,param.nc,1.0); % equidistant grid
        for cc = 1:param.nc
            val_cons(cc) = fun_val_det(grid_cons(cc),grid.coh(ac,1,1),param);
        end
        % pick maximum on that grid
        [fun.val(ac,1,1),max_cc] = max(val_cons(:));
        % search for optimum in neighborhood of max using golden search;
        if ( max_cc == 1 )
            a = grid_cons(1);
            c = grid_cons(2);
        elseif (max_cc == param.nc)
            a = grid_cons(param.nc);
            c = grid_cons(param.nc-1);
        else
            a = grid_cons(max_cc-1);
            c = grid_cons(max_cc+1);
        end

        [fun.cons(ac,1,1),fun.val(ac,1,1)] = golden(@fun_val_det,a,c,grid.coh(ac,1,1),param);
    end
    % compute savings 
    fun.sav(:,1,1) = grid.coh(:,1,1)-fun.cons(:,1,1);
  
else % solve hh problem recursively
    
    % final period
    fun.val_final=zeros(param.na,param.ne,param.nd);
    fun.cons_final = zeros(param.na,param.ne,param.nd);
    fun.prob_final = zeros(param.na,param.ne);
    grid.coh_final = zeros(param.na,param.ne);
    
    % compute choice specific policy functions
    for ac=1:param.na
        for ec=1:param.ne
            for dc=1:param.nd
                if (dc==1)
                    lab=1;
                else 
                    lab=0;
                end;
                fun.cons_final(ac,ec,dc) = grid.ass(ac) + (lab)*param.eta(ec);
                fun.val_final(ac,ec,dc)=fun_util(fun.cons_final(ac,ec,dc),lab,param);                
            end;
            
            % compute choice probabilities:
            prob_temp=fun_chprob(param,fun.val_final(ac,ec,:));
            fun.prob_final(ac,ec)=prob_temp(1);
            
            grid.coh_final(ac,ec)= grid.ass(ac)+(param.Gam)*param.eta(ec);
            fun.ass(ac,ec,param.nt+1)=grid.ass(ac);
        end;
    end; 
    
    % other periods (iterate backwards)
    for tc=param.nt:-1:1
        
        if (tc==param.nt) % last period before final discrete decision
            opt_discr=1;
        else
            opt_discr=0;
        end;
    
        if (opt.meth==1) % global
            
            val_cons = zeros(param.nc,1); % consumption grid
            min_cons = param.tol;
            
            % evaluate value function on grid:
            for ac = 1:param.na
            
                max_cons  = grid.coh(ac,tc,opt.meth);     % maximum possible consumption today 
                grid_cons = makegrid(min_cons,max_cons,param.nc,1.0);
                
                for cc = 1:param.nc
                    val_cons(cc) = fun_val(grid_cons(cc),grid.coh(ac,tc,opt.meth),param,fun,grid,tc,opt_discr);
                end
                
                % pick maximum on that grid
                [fun.val(ac,tc,opt.meth),max_cc] = max(val_cons(:));
                
                % search for optimum in neighborhood of max using golden search;
                if ( max_cc == 1 )
                    a = grid_cons(1);
                    c = grid_cons(2);
                elseif (max_cc == param.nc)
                    a = grid_cons(param.nc);
                    c = grid_cons(param.nc-1);
                else
                    a = grid_cons(max_cc-1);
                    c = grid_cons(max_cc+1);
                end

                [fun.cons(ac,tc,opt.meth),fun.val(ac,tc,opt.meth)] = golden(@fun_val,a,c,grid.coh(ac,tc,opt.meth),param,fun,grid,tc,opt_discr);
            end

        elseif (opt.meth==2) % dc-egm
            
            % set up exogenous savings grid
            max_sav = min(fun.ass(param.na,:,tc+1)) ;
            min_sav = max(fun.ass(1,:,tc+1));
            grid.sav(param.np1aux:param.na,tc) = makegrid(min_sav,max_sav,param.na-param.np1aux+1,param.curv);
            grid.sav(1:param.np1aux-1,tc) = 0.0;
            
            % solution on savings grid (endogenous part of coh grid)
            for ac=param.na:-1:param.np1aux
                sav = grid.sav(ac,tc);
                [fun.cons(ac,tc,opt.meth),grid.coh(ac,tc,opt.meth),fun.val(ac,tc,opt.meth)]=fun_dcegm(sav,grid,fun,param,ac,tc,opt_discr,0);
            end;
            
            % construct exogenous part of coh grid (linear)
            min_coh = min(param.eta);
            if (min_coh>grid.coh(param.np1aux,tc,opt.meth) )
                min_coh = 0.99*grid.coh(param.np1aux,tc,opt.meth);
            end  
            coh_aux = makegrid(min_coh,grid.coh(param.np1aux,tc,opt.meth),param.np1aux,1.0);
            % solution on savings grid (exogenous part of coh grid)
            for ac=param.np1aux-1:-1:1
                sav = grid.sav(ac,tc);
                grid.coh(ac,tc,opt.meth) = coh_aux(ac);
                [fun.cons(ac,tc,opt.meth),grid.coh(ac,tc,opt.meth),fun.val(ac,tc,opt.meth)]=fun_dcegm(sav,grid,fun,param,ac,tc,opt_discr,1);
            end
            
            % check monotonicity of endogenous coh grid, and if necessary apply correction (see Iskhakov et al. 2017 paper for details)    
            coh_raw = grid.coh(:,tc,opt.meth); % store raw coh grid
            if sum(diff(coh_raw)<0)>0
                % determine common grid
                min_coh_common=param.min_coh; %min(grid.coh(:,tc,opt.meth));
                max_coh_common=max(grid.coh(:,tc,opt.meth));
                coh_common(1:param.na,1) = makegrid(min_coh_common,max_coh_common,param.na,param.curv);
                
                % take raw solution and kick out invalid ones
                cons_raw = fun.cons(:,tc,opt.meth);
                val_raw = fun.val(:,tc,opt.meth);
                % call upper envelope
                [coh_ref,cons_ref,val_ref] = fun_UpperEnvelope(coh_raw,cons_raw,val_raw,coh_common,param);
                
                % assign refined values to coh grid, value and policy functions
                grid.coh(:,tc,opt.meth)=coh_ref;
                fun.val(:,tc,opt.meth)=val_ref;
                fun.cons(:,tc,opt.meth)=cons_ref;
            end
            
        elseif (opt.meth==3) % exogm (rootfinding)
            
            for ac=1:param.na
                % disp(['state ac=',num2str(ac)]);
                % solve for cons with a univariate solver
                coh = grid.coh(ac,tc,opt.meth);
                x0=0.5*grid.coh(ac,tc,opt.meth);
              
                f = @(x)fun_foc(x,coh,grid,fun,param,tc);
                cons = fzero(f,x0);
                %cons = max(eps,cons);
               
                % check if the borrowing constraint is binding
                sav = grid.coh(ac,tc,opt.meth)-cons;
                
                 if (sav<0.0 || cons<1.0e-02) % bc binding
                     fun.cons(ac,tc,opt.meth)=grid.coh(ac,tc,opt.meth);
                 else
                     fun.cons(ac,tc,opt.meth)=cons;
                 end
               
                fun.val(ac,tc,opt.meth)=fun_val_exogm(fun.cons(ac,tc,opt.meth),coh,param,fun,grid,tc,opt_discr);
             
            end
          
        end;
        
        fun.sav(:,tc,opt.meth) = grid.coh(:,tc,opt.meth)-fun.cons(:,tc,opt.meth);
        for ec=1:param.ne
            fun.ass(:,ec,tc)=grid.coh(:,tc,opt.meth)-param.eta(ec);
        end
    end;
end 

% plot of policy functions (period 0)
if (opt.plt)
    xla='$x$';
    yla='$c(x),s(x)$';
    leg={'consumption','savings'};
    tit=['consumption and savings policies for t=0, solution method ',num2str(opt.meth)];
    outpath=[];
    grname=['polfun_conssav_meth', num2str(opt.meth)];
    ax=[];
    onesideplt_marker(squeeze(grid.coh(:,1,opt.meth)),[ fun.cons(:,1,opt.meth),fun.sav(:,1,opt.meth)],xla,yla,leg,tit,outpath,grname,ax,opt.pr)
    
    % plot of value function
    xla='$x$';
    yla='$v(x)$';
    leg=[];
    tit=['value function for t=0, solution method ',num2str(opt.meth)];
    grname=['valfun_meth', num2str(opt.meth)];
    onesideplt_marker(grid.coh(:,1,opt.meth),fun.val(:,1,opt.meth),xla,yla,leg,tit,outpath,grname,ax,opt.pr)
end

end    % end function fun_solve
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [res] = fun_foc(x,coh,grid,fun,param,tc)

cons=x;

muc = fun_margu(cons);
sav = coh-cons;
sav = max(0.0,sav);

vp_coh_tp1=zeros(param.nd,1);
v_tp1=zeros(param.nd,1);
evp1_coh=0.0;
evp1=0.0;
for ecc=1:param.ne % next period income
    if (tc+1==param.nt+1)
        
        coh_tp1 = sav + (param.Gam)*param.eta(ecc) ; % coh next period
      
        for dcc=1:param.nd
            temp = func_intp(grid.coh_final(:,ecc),fun.cons_final(:,ecc,dcc),coh_tp1,0);
            vp_coh_tp1(dcc) = 1.0/temp;
            v_tp1(dcc) = func_intp(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh_tp1,0);
        end
        
        if (param.sigma <= sqrt(eps))
            [val1,maxd] = max(v_tp1(:));
            prob(maxd,1)=1.0;
            evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1(maxd);

        else
            val1 = fun_logsum(param,v_tp1);
            prob = fun_chprob(param,v_tp1);
            
            for dcc=1:param.nd
                evp1_coh = evp1_coh + param.wght(ecc)*prob(dcc)*vp_coh_tp1(dcc);
            end;
        end
        
        evp1 = evp1 + param.wght(ecc)*val1;  
        
    else % other (intermediate) periods
        
        coh_tp1 = sav + param.eta(ecc); % next period coh
        
        temp = func_intp(grid.coh(:,tc+1,3),fun.cons(:,tc+1,3),coh_tp1,0);

        vp_coh_tp1 = 1.0/temp;

        v_tp1 = func_intp(grid.coh(:,tc+1,3),fun.val(:,tc+1,3),coh_tp1,0);   

        evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1;

        evp1 = evp1 + param.wght(ecc)*v_tp1;   
           
    end
end
   
lhs = muc;
rhs = evp1_coh;

res = rhs-lhs;
        
end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function res=fun_rhs(x,coh,grid,fun,param,tc)


cons=x;

muc = fun_margu(cons);
sav = coh-cons;

vp_coh_tp1=zeros(param.nd,1);
v_tp1=zeros(param.nd,1);
evp1_coh=0.0;
evp1=0.0;
for ecc=1:param.ne % next period income
    if (tc+1==param.nt+1)
        
        coh_tp1 = sav + (param.Gam)*param.eta(ecc) ; % coh next period
      
        for dcc=1:param.nd
            temp = func_intp(grid.coh_final(:,ecc),fun.cons_final(:,ecc,dcc),coh_tp1,0);
            vp_coh_tp1(dcc) = 1.0/temp;
            v_tp1(dcc) = func_intp(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh_tp1,0);
        end
        
        if (param.sigma <= sqrt(eps))
            [val1,maxd] = max(v_tp1(:));
            prob(maxd,1)=1.0;
            evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1(maxd);

        else
            val1 = fun_logsum(param,v_tp1);
            prob = fun_chprob(param,v_tp1);
            
            for dcc=1:param.nd
                evp1_coh = evp1_coh + param.wght(ecc)*prob(dcc)*vp_coh_tp1(dcc);
            end;
        end
        
        evp1 = evp1 + param.wght(ecc)*val1;  
        
    else % other (intermediate) periods
        
        coh_tp1 = sav + param.eta(ecc); % next period coh
        
        temp = func_intp(grid.coh(:,tc+1,3),fun.cons(:,tc+1,3),coh_tp1,0);

        vp_coh_tp1 = 1.0/temp;

        v_tp1 = func_intp(grid.coh(:,tc+1,3),fun.val(:,tc+1,3),coh_tp1,0);   

        evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1;

        evp1 = evp1 + param.wght(ecc)*v_tp1;   
           
    end
end
   
res = evp1_coh;

end %fun_rhs
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function val = fun_val_exogm(c,x,param,fun,grid,tc,opt_discr)
% compute value function (for given cons and coh), 
% opt_discr=1 => next period (which is final) contains discrete choice

u  = fun_util(c,0.0,param);    % utility for labor = 0;
val = u;
ap  = x - c;     % budget constraint today

if (opt_discr==1) % solution for last period before discrete choice
    
    v_tp1 = zeros(param.nd,1);
    for ecc = 1:param.ne    % shock tomorrow
        for dcc = 1:param.nd % discrete choice tomorrow
        
            coh_tp1 = ap + (param.Gam)*param.eta(ecc) ;
        
            % interpolate on tomorrow's choice-specific value function
            v_tp1(dcc) = func_intp(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh_tp1,0);
        end

        if (param.sigma <= sqrt(eps))
            val1 = max(v_tp1(:));
        else
            val1 = fun_logsum(param,v_tp1);
        end
        % today's value function
        val = val + param.wght(ecc)*val1;
    end
    
else % solution for other periods
    
    for ecc = 1:param.ne % shock tomorrow
        % next period coh
        coh_tp1 = ap+param.eta(ecc);
     
        % interpolate on next period value function
        v_tp1 = func_intp(grid.coh(:,tc+1,3),fun.val(:,tc+1,3),coh_tp1,0);
        % compute expectation
        val = val + param.wght(ecc)*v_tp1;
    end
end

end     % end function fun_val_exogm
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [cons,coh,val]=fun_dcegm(sav,grid,fun,param,ac,tc,opt_discr,indbc)
% compute solution of hh problem, no discrete decision in the current period
% opt_discr=1 => discrete decision in the next period (use dc-egm), otherwiese basic egm

evp1_coh=0;
evp1=0;

if (opt_discr==1) % use dc-egm
    
    v_tp1 = zeros(param.nd,1);
    prob = zeros(param.nd,1);
    vp_coh_tp1 = zeros(param.nd,1);

    for ecc=1:param.ne % next period income shock states
        for dcc=1:param.nd % next period discrete choice
          
            coh_tp1 = sav + (param.Gam)*param.eta(ecc) ; % coh next period
            
            if (coh_tp1<grid.coh_final(1,ecc))
                disp('x too small, how can this be?')
            end

            temp = func_intp(grid.coh_final(:,ecc),fun.cons_final(:,ecc,dcc),coh_tp1,0);

            vp_coh_tp1(dcc) = 1.0/temp;

            v_tp1(dcc) = func_intp(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh_tp1,0);
            if (v_tp1(dcc)==-inf)
                v_tp1(dcc)=-90000;
            end
            if (vp_coh_tp1(dcc)==inf)
                vp_coh_tp1(dcc)=90000;
            end
        end;    

        if (param.sigma <= sqrt(eps))
            [val1,maxd] = max(v_tp1(:));
            prob(maxd,1)=1.0;
            evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1(maxd);

        else
            val1 = fun_logsum(param,v_tp1);
            prob = fun_chprob(param,v_tp1);
            
            for dcc=1:param.nd
                evp1_coh = evp1_coh + param.wght(ecc)*prob(dcc)*vp_coh_tp1(dcc);
            end;
        end

        evp1 = evp1 + param.wght(ecc)*val1;    
    end;

else % standard egm
  
    for ecc=1:param.ne % next period income shock states

        coh_tp1 = sav + param.eta(ecc); % next period coh
        
        if (coh_tp1<grid.coh(1,tc+1,2))
             disp('x too small, how can this be?')
        end

        temp = func_intp(grid.coh(:,tc+1,2),fun.cons(:,tc+1,2),coh_tp1,0);

        vp_coh_tp1 = 1.0/temp;

        v_tp1 = func_intp(grid.coh(:,tc+1,2),fun.val(:,tc+1,2),coh_tp1,0);   

        evp1_coh = evp1_coh + param.wght(ecc)*vp_coh_tp1;

        evp1 = evp1 + param.wght(ecc)*v_tp1;    
    end;
end;

% consumption policy
if (indbc==0)
    cons = 1.0/evp1_coh;
    % endogenous coh
    coh = cons + sav;
else
    cons = grid.coh(ac,tc,2);
    coh= grid.coh(ac,tc,2);
end

% utility
util = fun_util(cons,0.0,param);

% value function
val = util + evp1;

end  % end function fun_dcegm
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function val = fun_val(c,x,param,fun,grid,tc,opt_discr)
% compute value function (for given cons and coh), 
% opt_discr=1 => next period (which is final) contains discrete choice

u  = fun_util(c,0.0,param);    % utility for labor = 0 (no labor disutility in all but last period);
val = u;
ap  = x - c;     % budget constraint today

if (opt_discr==1) % solution for last period before discrete choice
    
    v_tp1 = zeros(param.nd,1);
    for ecc = 1:param.ne    % shock tomorrow
        for dcc = 1:param.nd % discrete choice tomorrow
        
            coh_tp1 = ap + (param.Gam)*param.eta(ecc) ;
        
            % interpolate on tomorrow's choice-specific value function
            v_tp1(dcc) = func_intp(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh_tp1,0);
        end

        if (param.sigma <= sqrt(eps))
            val1 = max(v_tp1(:));
        else
            val1 = fun_logsum(param,v_tp1);
        end
        % today's value function
        val = val + param.wght(ecc)*val1;
    end
    
else % solution for other periods
    
    for ecc = 1:param.ne % shock tomorrow
        % next period coh
        coh_tp1 = ap+param.eta(ecc);
     
        % interpolate on next period value function
        v_tp1 = func_intp(grid.coh(:,tc+1,1),fun.val(:,tc+1,1),coh_tp1,0);
        % compute expectation
        val = val + param.wght(ecc)*v_tp1;
    end
end

end     % end function fun_val
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function val = fun_val_det(c0,x0,param)
% compute value function (for given cons and coh), next period discrete choice
epsi = 1.0e-6;

u0  = fun_util(c0,1.0,param);    % utility for labor = 1;
val = 2*u0;
a1  = x0 - c0;
c1  = c0;
a2  = a1 + 1.0 - c1;

c2_temp   = zeros(param.nd,1);
u2_temp   = zeros(param.nd,1);
lab2_temp = zeros(param.nd,1);
val_temp  = zeros(param.nd,1);

for dcc = 1:2

    if (dcc == 1) % employment
        lab2 = 1.0;
    else % unemployment
        lab2 = 0.0;
    end

    lab2_temp(dcc) = lab2;

    c2_temp(dcc) = a2 + lab2_temp(dcc);     % budget constraint tomorrow

    if (c2_temp(dcc)<epsi)
        u2_temp(dcc) = -90000.0;
    else    
        u2_temp(dcc) = fun_util(c2_temp(dcc),lab2_temp(dcc),param);
    end
    
    val_temp(dcc) = val+u2_temp(dcc);
end

if (param.sigma <= sqrt(eps))
    val = max(val_temp(:));
else
    val = fun_logsum(param,val_temp);
end

end     % end function fun_val_det
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [grid,agg]=fun_aggr(grid,fun,param,opt)

disp('compute aggregation and cross-sectional measure');

grid.Phi = zeros(param.na,param.nt+1);          % distribution of assets by age
grid.Phi_temp = zeros(param.na,param.ne,param.nt+1); % age and income shock
grid.Phi_temp2 = zeros(param.na,param.ne,param.nd); % final period - discr choice and income shock
grid.PhiAss = zeros(param.na,1);             % distribution of assets

% initial distribution over cash at hand
for ec=1:param.ne
    
    % income in current period/age:
    inc=param.eta(ec);
    
    % initial cash-on-hand:
    cohini=inc;
    
    frac=1.0*param.wght(ec);
    
    [vals,inds]=basefun(grid.coh(:,1,opt.meth),cohini,param.na);
    
    for wc=1:2
        grid.Phi_temp(inds(wc),ec,1)=grid.Phi_temp(inds(wc),ec,1)+vals(wc)*frac;
    end
    
end;

for ac=1:param.na
    grid.Phi(ac,1)=sum(grid.Phi_temp(ac,:,1));
end

% distribution for periods 2,3,..
for tc=2:param.nt+1 
    
    if (tc<param.nt+1) % intermediate periods w/o discrete decision (in default case only one such period)  
   
        for ac=1:param.na, 

            for ecc=1:param.ne, % income today

                % income in current period/age:
                inc=param.eta(ecc);

                % cash on hand: x=a+y ;
                coh =inc + fun.sav(ac,tc-1);

                frac = grid.Phi(ac,tc-1);
                if (frac==0) 
                    continue 
                end
                frac = grid.Phi(ac,tc-1)*param.wght(ecc);

                [vals,inds]=basefun(grid.coh(:,tc,opt.meth),coh,param.na);
                for wc=1:2
                    grid.Phi_temp(inds(wc),ecc,tc)=grid.Phi_temp(inds(wc),ecc,tc)+frac*vals(wc);
                end;

            end;
        end;
        
        for acc=1:param.na
            grid.Phi(acc,tc)=sum(grid.Phi_temp(acc,:,tc));
        end
    
    else % final period - discrete choice
       
        for ac=1:param.na,

            for ecc=1:param.ne, % income today

                prob=zeros(param.nd,1);
                v_temp=zeros(param.nd,1);
                vals_temp = zeros(2,param.nd);
                inds_temp = zeros(2,param.nd);

                frac = grid.Phi(ac,tc-1);

                if (frac==0)
                    continue
                end

                % coh doesn't depend on discrete choice: x=a+y;
                coh = fun.sav(ac,tc-1)+ (param.Gam)*param.eta(ecc);

                for dcc=1:param.nd, % discrete choice
                    [v_temp(dcc),vals_temp(:,dcc),inds_temp(:,dcc)]=func_intp_mod(grid.coh_final(:,ecc),fun.val_final(:,ecc,dcc),coh,0);
                end

                % obtain probabilities
                if (param.sigma <= sqrt(eps))
                    [val1,maxd] = max(v_temp(:));
                    prob(maxd,1)=1.0;
                else
                    prob = fun_chprob(param,v_temp);
                end

                for dcc=1:param.nd,    

                    frac = grid.Phi(ac,tc-1)*param.wght(ecc)*prob(dcc);

                    for wc=1:2
                        grid.Phi_temp2(inds_temp(wc,dcc),ecc,dcc)=grid.Phi_temp2(inds_temp(wc,dcc),ecc,dcc)+frac*vals_temp(wc,dcc);
                    end;

                end;

            end;    

        end;
        
        for acc=1:param.na
            for ecc=1:param.ne
                grid.Phi_temp(acc,ecc,tc)=sum(grid.Phi_temp2(acc,ecc,:));
            end
            grid.Phi(acc,tc)=sum(sum(grid.Phi_temp(acc,:,tc)));
        end

    end;
end;    % end for tc

% Check that for each period distribution sums to 1
for tc=1:param.nt+1
    sumprob=sum(sum(sum(grid.Phi(:,tc))));
    if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
        beep; beep; beep;
        warning(['distribution is bad in period ',num2str(tc)]);
    end;
end

% Check if Grid is Big enough
sumprob=sum(sum(grid.Phi(param.na,:)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning(['grid too small -- increase your grid']);
    pause
end;

agg.ass=0.0;
agg.cons=0.0;
agg.coh=0.0;

% aggregation
for tc=1:param.nt+1,
    
    for ac=1:param.na,
        
        grid.PhiAss(ac)=grid.PhiAss(ac)+grid.Phi(ac,tc);
        
        for ec=1:param.ne
            agg.ass=agg.ass+grid.Phi_temp(ac,ec,tc)*fun.ass(ac,ec,tc); 
        end
           
        if (tc==param.nt+1)
            for ec=1:param.ne
                agg.coh=agg.coh+grid.Phi_temp(ac,ec,tc)*grid.coh_final(ac,ec);
                for dc=1:param.nd
                    agg.cons=agg.cons+grid.Phi_temp2(ac,ec,dc)*fun.cons_final(ac,ec,dc);
                end                
            end
        else
            agg.coh=agg.coh+grid.Phi(ac,tc)*grid.coh(ac,tc,opt.meth);
            agg.cons=agg.cons+grid.Phi(ac,tc)*fun.cons(ac,tc,opt.meth); 
        end
    end;

end;

% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup(grid_x,x,0);
        
        if ( (i+1)>nx),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------


end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function util = fun_util(c,lab,param)
% computes the utility function

util = log(c) + param.phi*log(param.Gam-lab);

end     % end function fun_util
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function margu = fun_margu(c)
% computes marginal utility 

margu = 1.0./c; 

end     % end function fun_util
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)

% makes curved grid according to curvature parameter c
scale  = x2-x1;
grd    = zeros(n,1);
grd(1) = x1;
grd(n) = x2;
for i = 2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end

end          % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function val_thrs = fun_thrs(x,param)

val_thrs = (1.0-param.temp)*x^3 + 3*x^2 + 3*x + 1.0;

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function val=fun_logsum(param,x)

if (param.sigma <= sqrt(eps))
    val = max(x);
else    
    %logsum by columns
    mx  = max(x,[],1);
    mxx = x-repmat(mx,size(x,1),1);
    val = mx+param.sigma*log(sum(exp(mxx/param.sigma),1));
end

end  % end function fun_logsum
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function prob=fun_chprob(param,x)
% choice probabilities (multinomial logit formula)

% logsum by columns
mx  = max(x,[],1);
mxx = x-repmat(mx,size(x,1),1);
logsum = mx+param.sigma*log(sum(exp(mxx/param.sigma),1));

% choice prob
prob    = exp( (x - repmat(logsum,1))./param.sigma );

end % function fun_chprob
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [uM,uC,uVt,Increase, Fall,CommonMgrid,CommonVt,CommonC]= fun_UpperEnvelope(M,C,Vt,gridcom,param)
% This function is taken from Jorgensen's replication package of Iskhakov et al (2017) 
% only minor adjustments were undertaken
% input:
%   M: raw grid
%   C: raw control
%   Vt: raw value function
%   gridcom: common grid
% output: 
%   uM: corrected grid
%   uC: corrreced control
%   uVt: corrected value function

% First find ALL drops and increases
Fall = []; % initialize with empty and then add the last point below while-loop
Increase = 1; % Initialize such that the lowest point is the first grid point
i = 2; % Initialize
while (i<=param.na-1)
    if ((M(i+1)<M(i) && M(i)>M(i-1)) || (Vt(i)<Vt(i-1) && M(i)>M(i-1)) ) % The last point before resources fall: i   (or the point where the value function decreased) 
        Fall = [Fall;i]; % add the point to the vector of points
        % find the point where resources again is increasing:
        k = i;
        while (M(k+1)<M(k))
            k = k + 1;
        end
        Increase = [Increase;k];
        i = k; % Set the index to the point where resources again is increasing
    end
    i = i +1;
end
Fall = [Fall ;param.na]; % Add the last point to the vector. This makes the follwing easier since both end-points now are included in the vectors
NumKinks    = length(Fall);

% Use these segments to sequentially find upper envelopes
NumP        = param.na;                             % Number of points used to interpolate on common grid
CommonMgrid = gridcom;

CommonVt    = -(1.0e10)*ones(NumP,NumKinks)*NaN;   % Initialize interpolated value function to vey low number (or NaN)
CommonC     = NaN(NumP,NumKinks);
% TAKE THE FIRST ONE BY HAND: prevent all the NaN-stuff..
for j = 1:NumKinks
    below = (M( Increase(j) ) >= CommonMgrid);
    above = (M( Fall(j) ) <= CommonMgrid);
    InRange = (above+below)==0;
    CommonVt(InRange,j) = interp1(M(Increase(j):Fall(j)),Vt(Increase(j):Fall(j)), CommonMgrid(InRange) );
    CommonC(InRange,j)  = interp1(M(Increase(j):Fall(j)),C(Increase(j):Fall(j)),  CommonMgrid(InRange) ); % Interpolat econsumption also. May not be nesserary
end

CommonVt(param.na,NumKinks)=Vt(param.na);
CommonC(param.na,NumKinks)=C(param.na);

% Now take the max of all these functions. Since the CommonVt
% is either NaN or very low number outside the range of the actual line-segment this works "globally"
[uVt,id]    = max(CommonVt,[],2);
uM          = CommonMgrid;
% Ad the zero point in the bottom
if isnan(uVt(1)) 
   uVt(1) = 0; % Since M=0 here
   CommonC(1,1)  = uM(1);%
end
% Extrapolate if NaNs are introduced due to the common grid
% going outside all the sub-line segments
IsNaN = isnan(uVt);
uVt(IsNaN) = interp1(uM(IsNaN==0),uVt(IsNaN==0),uM(IsNaN));
LastBeforeNaN = [diff(IsNaN)>0 ; 0];
LastId = LastBeforeNaN'*id; % Find last id-number 
id(IsNaN) = LastId;

LinInd      = cumsum(ones(NumP,1)) + (id - 1)*NumP; % Linear index used to get optimal consumption based on "id"  from max
uC          = CommonC(LinInd);
uC(IsNaN)   = interp1(uM(IsNaN==0),uC(IsNaN==0),uM(IsNaN));

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

        