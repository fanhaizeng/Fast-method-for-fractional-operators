function fast_conv_fode_Y18m02d19
% fast convolution for system of FODEs with nonsmooth solutions
% Gauss quadrature is used to discretize the kernel function
close all;    

% nonlinear_eq_fast_conv_corrections;
nonlinear_eq_fast_conv_corrections_jsc_revision;
%END


function nonlinear_eq_fast_conv_corrections_jsc_revision
clr = {'r+-','kd-','mo-','b>-','k*-','c*-','b-'};
m = 4;   % number of correction terms
B = 5;  % basis
vps = 1e-10;  % accuracy of Gauss quadrature
solu = 2; % different analytical solutions
dT = 1;   % tau <= dT < T
[uv,fv,tau,mu,T,ut0,alf,sgm,Hv] = solu_system_nonlinear_system_1(solu);
% uv(t) : analytical solution, handle function
% ut0 : initial value
% fv(u,t) : nonlinar term
% tau : stepsiz
% alf : fractional order
% sgm : correction indices
% Hv : Jacobian of Hv = partial u of fv(u,t)
tau = 1/256;
sgm = [0.1 0.2 0.3 0.4];
% s = sgm(1:m); 
% tau = 0.01; dT = 1; T = 10;
len = length(sgm)+1; 
for k = 1:len    
    s = sgm(1:k-1);
% [uh,ee] = fast_L1_method(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT);
% [uh,ee] = fast_L1_method_trap(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT);
% [uh,ee] = fode_nonlinear_quadratic_p2_sisc(uv,fv,ut0,tau,mu,T,alf,Hv,s);
[uh,ee] = fast_p2_method(uv,fv,ut0,tau,mu,T,alf,Hv,s,B,dT,vps);    
% [uh,ee] = L1_method_corection_m(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
    eee(k,:) = ee;
end
t = 0:tau:T; p0 = 200;
for k = 1:len
    semilogy(t(6:p0:end),eee(k,6:p0:end),clr{k}); hold on;
    eee;
end
xlabel('t'); ylabel('Error'); 
legend('m = 0','m = 1','m = 2','m = 3','m = 4');
%END

function nonlinear_eq_fast_conv_corrections
m = 4;   % number of correction terms
B = 5;  % basis
vps = 1e-10;  % accuracy of Gauss quadrature
solu = 2; % different analytical solutions
dT = 1;   % tau <= dT < T
[uv,fv,tau,mu,T,ut0,alf,sgm,Hv] = solu_system_nonlinear_system_1(solu);
% uv(t) : analytical solution, handle function
% ut0 : initial value
% fv(u,t) : nonlinar term
% tau : stepsiz
% alf : fractional order
% sgm : correction indices
% Hv : Jacobian of Hv = partial u of fv(u,t)
sgm = [0.1 0.2 0.3 0.4 0.5];
s = sgm(1:m); 
% tau = 0.01; dT = 1; T = 10;
len = length(tau);   cputime = zeros(len,1); 
for k = 1:len
    tic       
% [uh,ee] = fast_L1_method(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT);
% [uh,ee] = fast_L1_method_trap(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT);
[uh,ee] = fode_nonlinear_quadratic_p2_sisc(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
% [uh,ee] = fast_p2_method(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT,vps);    
% [uh,ee] = L1_method_corection_m(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
    e(k,1) = norm(ee,inf);      e2(k,1) = ee(end);
    cputime(k) = toc;
end
fprintf('\nalpha = ');
fprintf('%-f   ',alf); fprintf('\n');
fprintf('T = %f\n',T);
fprintf('sigma(X) = %s\n',num2str(s));
odr = myorder(tau(1:len),tau(1:len),e);  odr2 = myorder(tau(1:len),tau(1:len),e2);  
s = {'T/tau','max-err (X)','order','error at T (X)','order','cputime'};
fmts = {'%-8s','%-16s','%-10s','%-16s','%-8s','%12s'};
 fmt = {'%-8d','%-16.4e','%-10.4f','%-16.4e','%-8.4f','%12.6f'};
output2(fmts,s,fmt,round(1./tau(1:len)),e,odr,e2,odr2,cputime);
p0 = 50;  k = len;
myplot0(uh,tau(k),T,p0);
myplot0(ee,tau(k),T,40);
% nT,cpu'
% fprintf('%-12.4e',ee5); fprintf('\n');
%END


function nonlinear_eq_fast_conv_p2_correction_guass_system_1_N
clr = {'r+-','kd-','mo-','b>-','k*-','c*-','b-'};
m = 3;   % number of correction terms
B = 1e8;  % basis
vps = 1e-10;

dT = 1;
solu = 2; % different analytical solutions
[uv,fv,tau,mu,T,ut0,alf,sgm,Hv] = solu_system_nonlinear_system_1(solu);
tau = 0.001*[1 1 1 1 1];   
N = [64 128 256 512 1024]'*8; 
% N = N(end:-1:1);
len = length(tau);
len = length(tau);   
cputime = zeros(len,1); 
s = sgm(1:m); 
m0 = 4; d0 = 2^m0;
T = 2000;  

% N = 1024; dT = 0.1;
% tmax = dT*(1/(1-vps^(1/2/N))-1);

for k = 1:len
    tic 
%     [uh,ee] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
    
%     [uh,ee] = fode_fast_conv_p2_gauss_correction_s1(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,NN0,dTT0);   
% dT = 1*tau(k);
% [uh,ee] = fast_L1_method(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT);
% [uh1,ee] = fode_nonlinear_quadratic_p2_sisc(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
[uh,ee] = fast_p2_method(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,dT,vps,N(k));

%     [uh,ee] = fast_L1_method_bak(uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B);
%     [uh,ee] = L1_method_corection_1(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
%     [uh2,ee] = L1_method_corection_m(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
%     eeee = norm(uh - uh1,inf)
%     eeee = uh - uh1;
%     [uh,ee] = fode_fast_conv_p2_gauss_correction_system_11....
%         (uv,fv,ut0,tau(k),mu,T,alf,Hv,s,B,NN0,TT0,dTT0);      
%     [uh,ee] = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);        
%     [uh,ee] = system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);  
%       [uh,ee] = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
%       [uh,ee] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau(k),mu,T,alf,Hv,s);
    if isempty(uv)
        ee  = abs(uh0(1:2^(len+m0-k):end) - uh);
    end
    e(k,1) = norm(ee,inf);      e2(k,1) = ee(end);
    eee{k} = ee;   u{k} = uh; 
    cputime(k) = toc;
end
fprintf('\nalpha = ');
fprintf('%-f   ',alf); fprintf('\n');
fprintf('T = %f\n',T);
fprintf('sigma(X) = %s\n',num2str(s));
odr = myorder(tau(1:len),tau(1:len),e);  odr2 = myorder(tau(1:len),tau(1:len),e2);  
s = {'T/tau','max-err (X)','order','error at T (X)','order','cputime'};
fmts = {'%-8s','%-16s','%-10s','%-16s','%-8s','%12s'};
 fmt = {'%-8d','%-16.4e','%-10.4f','%-16.4e','%-8.4f','%12.6f'};
output2(fmts,s,fmt,round(T./tau(1:len)),e,odr,e2,odr2,cputime);
p0 = 50;  k = len;
myplot1(u{k},u{k},tau(k),T,p0);
myplot2(eee,eee,tau,T,40);
xlabel('t'); ylabel('Error'); title('');
legend('N = 64','N = 128','N = 256','N = 512','N = 1028');
% legend('N = 512','N = 1028','N = 2048','N = 4096','N = 8192');
%END



function [uh,e,e2,ee] = fdm_multi_term_fode(uv,fv,ut0,tau,mu,T,alf0,nu,sgm)
% D_{0,t}^{alf}u + nv*D_{0,t}^{beta}u = mu*u + f
% the method is exact for u = t^smg, alf=<smg<2 is not an integer. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
% [w,B,C,w1,w2] = weght_2(ne,alf,sgm);
alf = alf0(1); alf1 = alf0(2);
[w0,B0,~,C0] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgm);
w = nu*tau^(1-alf)*w0;  B = nu*tau^(1-alf)*B0;  C = C0*nu*tau^(1-alf);
w1 = tau^(1-alf1)*w10;  B1 = tau^(1-alf1)*B10;  C1 = C10*tau^(1-alf1);
w = w+w1; B = B+B1; C = C+C1;
% C(:,round(ne/4):end) = 0;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = w(1) - mu*tau;
uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh0 = ut0;
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);  
no = length(sgm);
for k = 1:no
    uh(k+1) = uuh0(k) + uh0; 
end
uh(1) = uh0;  rhs1 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    rhs0 = B(n)*uh0 - uuh0*C(:,n+1);
    rhs = rhs0  + tau*fv(t(n+1));
    if n>1
        rhs1 =  uh(2:n)*w(n:-1:2); 
    end
    rhs = rhs - rhs1;
    u2 = rhs/arrL; 
    uh(n+1) = u2;    
end
%----------------------------------------------------- 
ue = uh;
for k = 1:nt+1
    ue(k) = uv(t(k));
end
ee = abs(ue-uh);
e2 = norm(ue-uh,inf);
e = abs(ue(nt)-uh(nt));
%END

function [uh,vh,eeu,eev] = fdm_system_fode_linear(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(t)
% No correction terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,~] = my_weght_2(ne,alf,[]);
[w10,B10,~,~] = my_weght_2(ne,alf1,[]);
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
taualf = tau^alf;  taualf1 = tau^alf1;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = 1:nt      
    rhs =   B(n+1)*uh0 + taualf*fv(t(n+1));
    rhs1 = B1(n+1)*vh0 + taualf1*gv(t(n+1));
%     if n>1
        rhs =   rhs  - uh(1:n)*w(n+1:-1:2);
        rhs1 = rhs1 -  vh(1:n)*w1(n+1:-1:2);
%     end
    u2 = arrL\[rhs;rhs1];
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END
function [uh,vh,eeu,eev] = fdm_system_fode_linear_c(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,sgm,sgmv)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(t)
% Correction terms are used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,C0] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgmv);
w = w0;     B = B0;   C = C0;
w1 = w10;   B1 = B10; C1 = C10;  
taualf = tau^alf;  taualf1 = tau^alf1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
uh00 = 0; vh00 = 0;
%--------------------------------------------
no2 = length(sgm); no1 = length(sgmv); 
no = max(no1,no2);
% One correction is applied with smaller stepsize to get starting values
if no>0 
    N = nt;
    [uh00,vh00] = ini_system_fode_linear_c(uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,sgm(1),sgmv(1));
    uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
    uh00 = uh(2:no2+1) - uh0;    vh00 = vh(2:no1+1) - vh0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    rhs =   B(n+1)*uh0 + taualf*fv(t(n+1))  - uh00*C(:,n+1);
    rhs1 = B1(n+1)*vh0 + taualf1*gv(t(n+1)) - vh00*C1(:,n+1);
    rhs =   rhs  - uh(1:n)*w(n+1:-1:2);
    rhs1 = rhs1 -  vh(1:n)*w1(n+1:-1:2);
    u2 = arrL\[rhs;rhs1];
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END
function [uh,vh,eeu,eev] = ini_system_fode_linear_c(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,sgm,sgmv)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(t)
% One term is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,C0] = my_weght_2(ne,alf,sgm(1));
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgmv(1));
w = w0;     B = B0;   C = C0;
w1 = w10;   B1 = B10; C1 = C10;  
taualf = tau^alf;  taualf1 = tau^alf1;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
%------------------ first step ----------------------
no = 1;  n = 1;
rhs =   B(n+1)*uh0 + taualf*fv(t(n+1))  + uh0*C(n+1);
rhs1 = B1(n+1)*vh0 + taualf1*gv(t(n+1)) + vh0*C1(n+1);
rhs =   rhs  - uh(1:n)*w(n+1:-1:2);
rhs1 = rhs1 -  vh(1:n)*w1(n+1:-1:2);
u2 = (arrL + diag([C(n+1),C1(n+1)]))\[rhs;rhs1];
uh(n+1) = u2(1);  vh(n+1) = u2(2);
uh00 = uh(2) - uh(1);   vh00 = vh(2) - vh(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    rhs =   B(n+1)*uh0 + taualf*fv(t(n+1))  - uh00*C(:,n+1);
    rhs1 = B1(n+1)*vh0 + taualf1*gv(t(n+1)) - vh00*C1(:,n+1);
    rhs =   rhs  - uh(1:n)*w(n+1:-1:2);
    rhs1 = rhs1 -  vh(1:n)*w1(n+1:-1:2);
    u2 = arrL\[rhs;rhs1];
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] = fdm_system_fode_nonlinear(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% No correction terms, fixed point method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,~] = my_weght_2(ne,alf,[]);
[w10,B10,~,~] = my_weght_2(ne,alf1,[]);
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
taualf = tau^alf;  taualf1 = tau^alf1;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = 1:nt      
    rhs =   B(n+1)*uh0 - uh(1:n)*w(n+1:-1:2) + taualf*fv(uh(n),vh(n),t(n+1));
    rhs1 = B1(n+1)*vh0 - vh(1:n)*w1(n+1:-1:2) + taualf1*gv(uh(n),vh(n),t(n+1));
    u2 = arrL\[rhs;rhs1];
    RHS0 = B(n+1)*uh0 - uh(1:n)*w(n+1:-1:2);
    RHS1 = B1(n+1)*vh0 - vh(1:n)*w1(n+1:-1:2);
    for kk=1:50
        rhs =  RHS0 + taualf*fv(u2(1),u2(2),t(n+1));
        rhs1 = RHS1 + taualf1*gv(u2(1),u2(2),t(n+1));
        uu2 = arrL\[rhs;rhs1];
        if norm(u2-uu2,inf)<1e-14
            break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] .....
    = fdm_system_fode_nonlinear_newton(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% No correction terms, Newton method is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,~] = my_weght_2(ne,alf,[]);
[w10,B10,~,~] = my_weght_2(ne,alf1,[]);
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
taualf = tau^alf;  taualf1 = tau^alf1;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = 1:nt      
    RHS0 = -B(n+1)*uh0 + uh(1:n)*w(n+1:-1:2);
    RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2);
    u2 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2(1),u2(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),t(n+1));
        rhs = -(arrL*u2 + [rhs0;rhs1]);
        J = arrL - E*Hv(u2(1),u2(2),t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] ....
    = fdm_system_fode_nonlinear_newton_c(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
  [w0,B0,~,C] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C1] = my_weght_2(ne,alf1,sgmv);
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
taualf = tau^alf;  taualf1 = tau^alf1;
% nt0 = round(1/tau);  C(:,nt0:end) = 0; C1(:,nt0:end) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
uh00 = 0; vh00 = 0;
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(no1,no2);
% One correction term with smaller stepsize is used to get the starting
% values
if no>0
    N = nt;
    [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
        (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
    uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
    uh00 = uh(2:no1+1) - uh0;   vh00 = vh(2:no2+1) - vh0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0 + uh(1:n)*w(n+1:-1:2) + uh00*C(:,n+1);
    RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) + vh00*C1(:,n+1);
    u2 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2(1),u2(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),t(n+1));
        rhs = -(arrL*u2 + [rhs0;rhs1]);
        J = arrL - E*Hv(u2(1),u2(2),t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] ....
    = system_fode_nonlinear_fast_conv(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% No correction terms, Newton method is applied, Talbot contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);      
taualf = tau^alf;        taualf1 = tau^alf1;
taualf_1 = tau^(alf-1);  taualf1_1 = tau^(alf1-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,lambda,lmd1,lmd2,lmd3,aphi1,aphi2] = ini_contour(tau,-alf,BB,nt);
[bw,lambda,lmd1,lmd2,lmd3,bphi1,bphi2] = ini_contour(tau,-alf1,BB,nt);
aphi1 = taualf*aphi1;   aphi2 = taualf_1*aphi2;
bphi1 = taualf1*bphi1;  bphi2 = taualf1_1*bphi2;
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);   u1 = uh(1) - uh0;  v1 = vh(1) - vh0;
for n = 1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi1*u1 - aphi2*(u1+uh0)....
        + taualf*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*aw{1}.*uy2{1}));
    RHS1 = bphi1*v1 - bphi2*(v1+vh0)....
        + taualf1*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*bw{1}.*vy2{1}));
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
    
    uy2{1} = lmd1{1}.*uy2{1} + lmd2{1}*u1 + lmd3{1}*u2;
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v1 + lmd3{1}*v2;
    for kk = 1:len
        uyy4{kk} = lmd1{kk}.*uyy4{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v1 + lmd3{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd1{1}.*uy4{1} + lmd2{1}*u1 + lmd3{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v1 + lmd3{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd1{2}.*uy4{2} + lmd2{2}*u1 + lmd3{2}*u2;
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v1 + lmd3{2}*v2;
    end
    u1 = u2;  v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi1*u1 - aphi2*(u1+uh0);
    RHS1 = bphi1*v1 - bphi2*(v1+vh0);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*aw{kk}.*uy2{kk}));
        RHS1 = RHS1 + taualf1*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*bw{kk}.*vy2{kk}));
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
% updation
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,u1,u2);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
        (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    u1 = u2;  v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] ....
    = fode_fast_conv_gauss(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% No correction terms, Newton method is applied, Gauss-Laguerre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);      
taualf = tau^alf;        taualf1 = tau^alf1;
taualf_1 = tau^(alf-1);  taualf1_1 = tau^(alf1-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% [awc,lambdac,lmd1c,lmd2c,lmd3c,aphi1c,aphi2c] = ini_contour(tau,-alf,BB,nt);
% [bwc,lambdac,lmd1c,lmd2c,lmd3c,bphi1c,bphi2c] = ini_contour(tau,-alf1,BB,nt);
[aw,lambda1,lmd10,lmd20,lmd30,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt);
[bw,lambda2,lmd1,lmd2,lmd3,bphi1,bphi2] = ini_contour_gauss(tau,-alf1,BB,nt);

aphi1 = taualf*aphi1;   aphi2 = taualf_1*aphi2;
bphi1 = taualf1*bphi1;  bphi2 = taualf1_1*bphi2;
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);   u1 = uh(1) - uh0;  v1 = vh(1) - vh0;
for n = 1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi1*u1 - aphi2*(u1+uh0)....
         + taualf*sum(exp(lambda1{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    RHS1 = bphi1*v1 - bphi2*(v1+vh0)....
        + taualf1*sum(exp(lambda2{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1});
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    uy2{1} = lmd10{1}.*uy2{1} + lmd20{1}*u1 + lmd30{1}*u2;
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v1 + lmd3{1}*v2;
    for kk = 1:len
        uyy4{kk} = lmd10{kk}.*uyy4{kk} + lmd20{kk}*u1 + lmd30{kk}*u2;
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v1 + lmd3{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd10{1}.*uy4{1} + lmd20{1}*u1 + lmd30{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v1 + lmd3{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd10{2}.*uy4{2} + lmd20{2}*u1 + lmd30{2}*u2;
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v1 + lmd3{2}*v2;
    end
    u1 = u2;  v1 = v2;
end
% ------------------------------------------------------------
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi1*u1 - aphi2*(u1+uh0);
    RHS1 = bphi1*v1 - bphi2*(v1+vh0);
    for kk = 1:L-1
        RHS0 = RHS0 +  taualf*sum(exp(lambda1{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(lambda2{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
% updation
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd10,lmd20,lmd30,u1,u2);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
        (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    u1 = u2;  v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev] ....
    = system_fode_nonlinear_fast_conv_p1_correction(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);      
taualf = tau^alf;        taualf1 = tau^alf1;
taualf_1 = tau^(alf-1);  taualf1_1 = tau^(alf1-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(no1,no2);
% Derive the starting values
uhm0 = 0; vhm0 = 0;
if no > 0
    if isempty(uv)
        N = nt;
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
    else
        uh(1:no+1) = uv(t(1:no+1));   vh(1:no+1) = vv(t(1:no+1));
    end
    uhm0 = uh(2:no1+1) - uh0;   vhm0 = vh(2:no2+1) - vh0;
end
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2] = ini_contour(tau,-alf,BB,nt);
[bw,blambda,blmd1,blmd2,blmd3,bphi1,bphi2] = ini_contour(tau,-alf1,BB,nt);
aphi1 = taualf*aphi1;   aphi2 = taualf_1*aphi2;
bphi1 = taualf1*bphi1;  bphi2 = taualf1_1*bphi2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};
end

T0 = 1; N0 = round(T0/tau)-1;
if N0 < no+1
    N0 = no+1;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32);
if no>0
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_p1(tau,T0,alambda{kk}(n),....
                almd2{kk}(n),almd3{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_p1(tau,T0,blambda{kk}(n),....
                blmd2{kk}(n),blmd3{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
        end
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
    end
end
clear wwww1 wwww2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);   u1 = uh(1) - uh0;  v1 = vh(1) - vh0;
for n = 1:no
    u2 = uh(n+1) - uh0;    v2 = vh(n+1) - vh0;
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u1 + almd3{kk}*u2 + acw1{kk}(n,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v1 + blmd3{kk}*v2 + bcw1{kk}(n,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u1 + almd3{2}*u2 + acw1{2}(n,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v1 + blmd3{2}*v2 + bcw1{2}(n,:);
    end
    u1 = u2;  v1 = v2;
end
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi1*u1 - aphi2*(u1+uh0)....
        + taualf*sum(imag(exp(alambda{1}*(t(n+1)-btau(1))).*aw{1}.*uy2{1}));
    RHS1 = bphi1*v1 - bphi2*(v1+vh0)....
        + taualf1*sum(imag(exp(blambda{1}*(t(n+1)-btau(1))).*bw{1}.*vy2{1}));
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
    
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u1 + almd3{kk}*u2 + acw1{kk}(n,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v1 + blmd3{kk}*v2 + bcw1{kk}(n,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u1 + almd3{2}*u2 + acw1{2}(n,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v1 + blmd3{2}*v2 + bcw1{2}(n,:);
    end
    u1 = u2;  v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi1*u1 - aphi2*(u1+uh0);
    RHS1 = bphi1*v1 - bphi2*(v1+vh0);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(imag(exp(alambda{kk}*(t(n+1)-btau(kk))).*aw{kk}.*uy2{kk}));
        RHS1 = RHS1 + taualf1*sum(imag(exp(blambda{kk}*(t(n+1)-btau(kk))).*bw{kk}.*vy2{kk}));
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;

%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,u1,u2);
%     [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
%         (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_2_p1......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_2_p1......
        (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,v1,v2,bcw1);
    u1 = u2;  v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh); 
%END

function [uh,vh,eeu,eev] ....
    = fode_fast_conv_p1_gauss_correction(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);      
taualf = tau^alf;        taualf1 = tau^alf1;
taualf_1 = tau^(alf-1);  taualf1_1 = tau^(alf1-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(no1,no2);
% Derive the starting values
uhm0 = 0; vhm0 = 0;
if no > 0
    if isempty(uv)
        N = nt;
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
    else
        uh(1:no+1) = uv(t(1:no+1));   vh(1:no+1) = vv(t(1:no+1));
    end
    uhm0 = uh(2:no1+1) - uh0;   vhm0 = vh(2:no2+1) - vh0;
end
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt);
[bw,blambda,blmd1,blmd2,blmd3,bphi1,bphi2] = ini_contour_gauss(tau,-alf1,BB,nt);
aphi1 = taualf*aphi1;   aphi2 = taualf_1*aphi2;
bphi1 = taualf1*bphi1;  bphi2 = taualf1_1*bphi2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};
end

T0 = 1; N0 = round(T0/tau)-1;
if N0 < no+1
    N0 = no+1;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32);
if no>0
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_p1(tau,T0,alambda{kk}(n),....
                almd2{kk}(n),almd3{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_p1(tau,T0,blambda{kk}(n),....
                blmd2{kk}(n),blmd3{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
        end
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
    end
end
clear wwww1 wwww2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);   u1 = uh(1) - uh0;  v1 = vh(1) - vh0;
for n = 1:no
    u2 = uh(n+1) - uh0;    v2 = vh(n+1) - vh0;
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u1 + almd3{kk}*u2 + acw1{kk}(n,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v1 + blmd3{kk}*v2 + bcw1{kk}(n,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u1 + almd3{2}*u2 + acw1{2}(n,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v1 + blmd3{2}*v2 + bcw1{2}(n,:);
    end
    u1 = u2;  v1 = v2;
end
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi1*u1 - aphi2*(u1+uh0)....
        + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    RHS1 = bphi1*v1 - bphi2*(v1+vh0)....
        + taualf1*sum(exp(blambda{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1});
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;
    
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u1 + almd3{kk}*u2 + acw1{kk}(n,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v1 + blmd3{kk}*v2 + bcw1{kk}(n,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u1 + almd3{1}*u2 + acw1{1}(n,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v1 + blmd3{1}*v2 + bcw1{1}(n,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u1 + almd3{2}*u2 + acw1{2}(n,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v1 + blmd3{2}*v2 + bcw1{2}(n,:);
    end
    u1 = u2;  v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi1*u1 - aphi2*(u1+uh0);
    RHS1 = bphi1*v1 - bphi2*(v1+vh0);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1) - uh0;   v2 = uu1(2) - vh0;

%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,u1,u2);
%     [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
%         (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_2_p1......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_2_p1......
        (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,v1,v2,bcw1);
    u1 = u2;  v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,eeu,eev]....
    = system_fode_nonlinear_fast_conv_p2_correction(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied, quatratic
% inerpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(max(no1,no2),2);
% use one correction term with smaller stepsize to get the starting values
N = nt;
if isempty(uv)
    if isempty(sgm)
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,[],[]);
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        uhm0 = 0; vhm0 = 0;
    else
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
    end
else
        uh(1:no+1) = uv(t(1:no+1));  vh(1:no+1) = vv(t(1:no+1));
        uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
end
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% uhm0 = 0;  vhm0 = 0;
% uh(1:3) = uv(t(1:3));  vh(1:3) = vv(t(1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% [aw,lambda,lmd1,lmd2,lmd3,aphi1,aphi2] = ini_contour(tau,-alf,BB,nt);
[aw,lambda,lmd1,lmd2,lmd3,lmd4,aphi0,aphi1,aphi2] = ini_contour_p2(tau,-alf,BB,nt);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,lmd11,lmd22,lmd33,lmd44] = ini_contour_p2_2(tau,-alf,BB,nt);
% aphi00 = tau^(-alf)*aphi00;  aphi11 = tau^(-alf)*aphi11; aphi22 = tau^(-alf)*aphi22;

% [bw,lambda,lmd1,lmd2,lmd3,bphi1,bphi2] = ini_contour(tau,-alf1,BB,nt);
[bw,lambda,lmd1,lmd2,lmd3,lmd4,bphi0,bphi1,bphi2] = ini_contour_p2(tau,-alf1,BB,nt);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
% [~,~,lmd11,lmd22,lmd33,lmd44,bphi00,bphi11,bphi22] = ini_contour_p2_2(tau,-alf1,BB,nt);
% bphi00 = tau^(-alf1)*bphi00;  bphi11 = tau^(-alf1)*bphi11;  bphi22 = tau^(-alf1)*bphi22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(lambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};
    acw2{kk} = acw1{kk};         bcw2{kk} = acw1{kk};
end

T0 = 0.4; 
N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(lambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,lambda{kk}(n),lmd2{kk}(n),lmd3{kk}(n),.....
                lmd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_11(tau,T0,lambda{kk}(n),lmd2{kk}(n),lmd3{kk}(n),....
                lmd4{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,lambda{kk}(n),lmd22{kk}(n),lmd33{kk}(n),.....
                lmd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,lambda{kk}(n),lmd22{kk}(n),lmd33{kk}(n),....
                lmd44{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww4(:,n) = vhm0*tmp;
        end
        acw2{kk} = wwww3;  bcw2{kk} = wwww4;
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
    end
end
clear wwww1 wwww2 wwww3 wwww4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  u1 = uh(2);  v1 = vh(2);
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);
uy2{1} = lmd11{1}.*uy2{1} + lmd22{1}*u0 + lmd33{1}*u1 + lmd44{1}*u2;
vy2{1} = lmd11{1}.*vy2{1} + lmd22{1}*v0 + lmd33{1}*v1 + lmd44{1}*v2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
vy2{1} = vy2{1} + bcw2{1}(1,:);
for kk = 1:len
    uyy4{kk} = lmd11{kk}.*uyy4{kk} + lmd22{kk}*u0 + lmd33{kk}*u1 + lmd44{kk}*u2;
    vyy4{kk} = lmd11{kk}.*vyy4{kk} + lmd22{kk}*v0 + lmd33{kk}*v1 + lmd44{kk}*v2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
    vyy4{kk} = vyy4{kk} + bcw2{kk}(1,:);
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = lmd11{1}.*uy4{1} + lmd22{1}*u0 + lmd33{1}*u1 + lmd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
    vy4{1} = lmd11{1}.*vy4{1} + lmd22{1}*v0 + lmd33{1}*v1 + lmd44{1}*v2 + bcw2{1}(1,:);
    vy44{1} = vy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = lmd11{2}.*uy4{2} + lmd22{2}*u0 + lmd33{2}*u1 + lmd44{2}*u2 + acw2{2}(1,:);
    vy4{2} = lmd11{2}.*vy4{2} + lmd22{2}*v0 + lmd33{2}*v1 + lmd44{2}*v2 + bcw2{2}(1,:);
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
clear lmd11 lmd22 lmd33 lmd44 acw2 bcw2; 
%
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);
    uy2{1} = lmd1{1}.*uy2{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2 + bcw1{1}(n-1,:);
    
    for kk = 1:len
        uyy4{kk} = lmd1{kk}.*uyy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v0 + lmd3{kk}*v1 + lmd4{kk}*v2 + bcw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd1{1}.*uy4{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd1{2}.*uy4{2} + lmd2{2}*u0 + lmd3{2}*u1 + lmd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v0 + lmd3{2}*v1 + lmd4{2}*v2 + bcw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
arrL = diag([aphi2,bphi2]) - diag([taualf,taualf1])*mu;
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*aw{1}.*uy2{1}));
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*bw{1}.*vy2{1}));
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1); v2 = uu1(2);
    
    uy2{1} = lmd1{1}.*uy2{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2 + bcw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = lmd1{kk}.*uyy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v0 + lmd3{kk}*v1 + lmd4{kk}*v2 + bcw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd1{1}.*uy4{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd1{2}.*uy4{2} + lmd2{2}*u0 + lmd3{2}*u1 + lmd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v0 + lmd3{2}*v1 + lmd4{2}*v2 + bcw1{2}(n-1,:);
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*aw{kk}.*uy2{kk}));
        RHS1 = RHS1 + taualf1*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*bw{kk}.*vy2{kk}));
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
%         if cond(J)>1e10
%             n,nt,uu1,J,return;
%         end
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1);  v2 = uu1(2);

    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,....
        uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_p2_correction(len,vy3,vy33,vy4,....
       vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,v0,v1,v2,bcw1);
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2......
%        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2);
%     [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2_p2......
%        (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,v0,v1,v2);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
eeu;
%END

function [uh,vh,wh,eeu,eev,eew] = fode_fast_conv_p2_gauss_correction_system_3......
    (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T,alf0,Hv,sgm,sgmv,sgmw,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + mu13*w + f(u,v,w,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + mu23*w + g(u,v,w,t)
% D_{0,t}^{gamma}w = mu31*u + mu32*v + mu33*w + g(u,v,w,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
% % % % % % % % % % % % % % % % % % % % % % % % 
alf3 = alf0(3);   taualf3 = tau^alf3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;      wh = uh;
uh0 = ut0;             vh0 = vt0;    wh0 = wt0;
uh(1) = uh0;           vh(1) = vh0;  wh(1) = wh0;
E = diag([taualf,taualf1,taualf3]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);  no3 = length(sgmw); 
no = max(max([no1,no2,no3]),2);
% use one correction term with smaller stepsize to get the starting values
if isempty(uv)
    N = nt;
    if isempty(sgm)
        [uh00,vh00,wh00] = ini_system_fode_nonlinear_newton_system_3.......
            (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau/N,mu,no*tau,alf0,Hv,[],[],[]);
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        wh(1:no+1) = wh00(1:N:end);
        uhm0 = 0; vhm0 = 0;  whm0 = 0;
    else
        [uh00,vh00,wh00] = ini_system_fode_nonlinear_newton_system_3.......
            (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1),sgmw(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        wh(1:no+1) = wh00(1:N:end);
        uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
        whm0 = wh(2:no3+1) - wh(1);
    end
else
    uh(1:no+1) = uv(t(1:no+1));  vh(1:no+1) = vv(t(1:no+1));
    wh(1:no+1) = wv(t(1:no+1));
    uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
    whm0 = wh(2:no3+1) - wh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isempty(uv) 
    uv = @(t) 0.*t;  vv = @(t) 0.*t;  wv = @(t) 0.*t; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
NN0 = 40;
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,tau,NN0);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,tau,NN0);

[bw,blambda,blmd1,blmd2,blmd3,blmd4,bphi0,bphi1,bphi2]....
    = ini_contour_p2_gauss_b(tau,-alf1,BB,nt,tau,NN0);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
[~,~,blmd11,blmd22,blmd33,blmd44] = ini_contour_p2_gauss_2_b(tau,-alf1,BB,nt,tau,NN0);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[cw,clambda,clmd1,clmd2,clmd3,clmd4,cphi0,cphi1,cphi2]....
    = ini_contour_p2_gauss_b(tau,-alf3,BB,nt,tau,NN0);
cphi0 = tau^alf3*cphi0;  cphi1 = tau^alf3*cphi1;  cphi2 = tau^alf3*cphi2;
[~,~,clmd11,clmd22,clmd33,clmd44] = ini_contour_p2_gauss_2_b(tau,-alf3,BB,nt,tau,NN0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
lenN1 = length(alambda{1}); lenN2 = length(blambda{1}); lenN3 = length(clambda{1});
lenN = max([lenN1,lenN2,lenN3]);
for kk = 1:len
    aw{kk}(lenN1+1:lenN) = 0; bw{kk}(lenN2+1:lenN) = 0; cw{kk}(lenN3+1:lenN) = 0;
    alambda{kk}(lenN1+1:lenN) = 0; blambda{kk}(lenN2+1:lenN) = 0; clambda{kk}(lenN3+1:lenN) = 0;
    almd1{kk}(lenN1+1:lenN) = 0;  blmd1{kk}(lenN2+1:lenN) = 0;  clmd1{kk}(lenN3+1:lenN) = 0;
    almd2{kk}(lenN1+1:lenN) = 0;  blmd2{kk}(lenN2+1:lenN) = 0;  clmd2{kk}(lenN3+1:lenN) = 0;
    almd3{kk}(lenN1+1:lenN) = 0;  blmd3{kk}(lenN2+1:lenN) = 0;  clmd3{kk}(lenN3+1:lenN) = 0;
    almd4{kk}(lenN1+1:lenN) = 0;  blmd4{kk}(lenN2+1:lenN) = 0;  clmd4{kk}(lenN3+1:lenN) = 0;
    almd11{kk}(lenN1+1:lenN) = 0; blmd11{kk}(lenN2+1:lenN) = 0; clmd11{kk}(lenN3+1:lenN) = 0;
    almd22{kk}(lenN1+1:lenN) = 0; blmd22{kk}(lenN2+1:lenN) = 0; clmd22{kk}(lenN3+1:lenN) = 0;
    almd33{kk}(lenN1+1:lenN) = 0; blmd33{kk}(lenN2+1:lenN) = 0; clmd33{kk}(lenN3+1:lenN) = 0;
    almd44{kk}(lenN1+1:lenN) = 0; blmd44{kk}(lenN2+1:lenN) = 0; clmd44{kk}(lenN3+1:lenN) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};  ccw1{kk} = acw1{kk};
    acw2{kk} = acw1{kk};         bcw2{kk} = acw1{kk};  ccw2{kk} = acw1{kk};
end

T0 = 1;   N0 = max(round(T0/tau)-1,2*BB-1);
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
for k = 1:no3
    s = sgmw(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0w{k} = x0; w0w{k} = w0;
end

[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_11(tau,T0,blambda{kk}(n),blmd2{kk}(n),blmd3{kk}(n),....
                blmd4{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
            tmp = correction_weights_11(tau,T0,clambda{kk}(n),clmd2{kk}(n),clmd3{kk}(n),....
                clmd4{kk}(n),x0w,w0w,sgmw,xx0,ww0);
            wwww2_c(:,n) = whm0*tmp;
            
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,blambda{kk}(n),blmd22{kk}(n),blmd33{kk}(n),....
                blmd44{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww4(:,n) = vhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,clambda{kk}(n),clmd22{kk}(n),clmd33{kk}(n),....
                clmd44{kk}(n),x0w,w0w,sgmw,xx0,ww0);
            wwww4_c(:,n) = whm0*tmp;
        end
        acw2{kk} = wwww3;   bcw2{kk} = wwww4;  ccw2{kk} = wwww4_c;
        acw1{kk} = zeros(N0,lenN);   bcw1{kk} = acw1{kk};  ccw1{kk} = acw1{kk};
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
        ccw1{kk}(1:N0-1,:) = wwww2_c(1:N0-1,:);
    end
end
clear wwww1 wwww2 wwww3 wwww4 wwww2_c wwww4_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 

wy2 = yz0; wy3 = yz0;  wy4 = yz0; wyy4 = yz0;
wy33 = yz0; wy44 = yz0;  wy66 = yz0; wy77 = yz0;  wy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  w0 = wh(1);
u1 = uh(2);  v1 = vh(2);  w1 = wh(2);
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);   w2 = wh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
vy2{1} = blmd11{1}.*vy2{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2;
wy2{1} = clmd11{1}.*wy2{1} + clmd22{1}*w0 + clmd33{1}*w1 + clmd44{1}*w2; 
uy2{1} = uy2{1} + acw2{1}(1,:);   
vy2{1} = vy2{1} + bcw2{1}(1,:);
wy2{1} = wy2{1} + ccw2{1}(1,:);
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    vyy4{kk} = blmd11{kk}.*vyy4{kk} + blmd22{kk}*v0 + blmd33{kk}*v1 + blmd44{kk}*v2;
    wyy4{kk} = clmd11{kk}.*wyy4{kk} + clmd22{kk}*w0 + clmd33{kk}*w1 + clmd44{kk}*w2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
    vyy4{kk} = vyy4{kk} + bcw2{kk}(1,:);
    wyy4{kk} = wyy4{kk} + ccw2{kk}(1,:);
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
    vy4{1} = blmd11{1}.*vy4{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2 + bcw2{1}(1,:);
    vy44{1} = vy4{1};
    wy4{1} = clmd11{1}.*wy4{1} + clmd22{1}*w0 + clmd33{1}*w1 + clmd44{1}*w2 + ccw2{1}(1,:);
    wy44{1} = wy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
    vy4{2} = blmd11{2}.*vy4{2} + blmd22{2}*v0 + blmd33{2}*v1 + blmd44{2}*v2 + bcw2{2}(1,:);
    wy4{2} = clmd11{2}.*wy4{2} + clmd22{2}*w0 + clmd33{2}*w1 + clmd44{2}*w2 + ccw2{2}(1,:);
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2; w0 = w1; w1 = w2;
clear almd11 almd22 almd33 almd44 blmd11 blmd22 blmd33 blmd44 acw2 bcw2; 
clear clmd11 clmd22 clmd33 clmd44 ccw2; 
%
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);  w2 = wh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    wy2{1} = clmd1{1}.*wy2{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
        wyy4{kk} = clmd1{kk}.*wyy4{kk} + clmd2{kk}*w0 + clmd3{kk}*w1 + clmd4{kk}*w2 + ccw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
        wy4{1} = clmd1{1}.*wy4{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
        wy44{1} = wy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
        wy4{2} = clmd1{2}.*wy4{2} + clmd2{2}*w0 + clmd3{2}*w1 + clmd4{2}*w2 + ccw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  v0 = v1;  v1 = v2;  w0 = w1;  w1 = w2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = diag([aphi2,bphi2,cphi2]) - diag([taualf,taualf1,taualf3])*mu;
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum(exp(blambda{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1});
    RHS2 = cphi0*w0 + cphi1*w1 - n^(-alf3)/gamma(1-alf3)*wh0.....
        + taualf3*sum(exp(clambda{1}*(t(n+1)-btau(1)-nT0(1))).*cw{1}.*wy2{1});
    uu1 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs2 = RHS2 - taualf3*hv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1;rhs2]);
        J = arrL - E*Hv(uu1(1),uu1(2),uu1(3),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);  wh(n+1) = uu1(3);
    u2 = uu1(1); v2 = uu1(2);  w2 = uu1(3);
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    wy2{1} = clmd1{1}.*wy2{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
        wyy4{kk} = clmd1{kk}.*wyy4{kk} + clmd2{kk}*w0 + clmd3{kk}*w1 + clmd4{kk}*w2 + ccw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
        wy4{1} = clmd1{1}.*wy4{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
        wy44{1} = wy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
        wy4{2} = clmd1{2}.*wy4{2} + clmd2{2}*w0 + clmd3{2}*w1 + clmd4{2}*w2 + ccw1{2}(n-1,:);
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;  w0 = w1;  w1 = w2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:N0      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    [wy2,wy3,wy4,b2,b3,d3,dd3] = update_1(wy2,wy3,wy4,wy6,wy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    RHS2 = cphi0*w0 + cphi1*w1 - n^(-alf3)/gamma(1-alf3)*wh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf *sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
        RHS2 = RHS2 + taualf3*sum(exp(clambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*cw{kk}.*wy2{kk});
    end
    uu1 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs2 = RHS2 - taualf3*hv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1;rhs2]);
        J = arrL - E*Hv(uu1(1),uu1(2),uu1(3),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);  wh(n+1) = uu1(3);
    u2 = uu1(1);  v2 = uu1(2);  w2 = uu1(3);

    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_p2_correction(len,vy3,vy33,vy4,....
       vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,blmd4,v0,v1,v2,bcw1);
    [wy3,wy33,wy4,wy44,wyy4,wy6,wy66,wy77] = update_p2_correction(len,wy3,wy33,wy4,....
       wy44,wyy4,wy6,wy66,wy77,b2,b3,d3,dd3,BB,n,clmd1,clmd2,clmd3,clmd4,w0,w1,w2,ccw1);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2; w0 = w1; w1 = w2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% t>=2*BB*tau, delete the correction terms
for kk = 1:len
    acw1{kk} = zeros(nt+1,1);    bcw1{kk} = zeros(nt+1,1);    
    ccw1{kk} = zeros(nt+1,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = N0+1:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    [wy2,wy3,wy4,b2,b3,d3,dd3] = update_1(wy2,wy3,wy4,wy6,wy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    RHS2 = cphi0*w0 + cphi1*w1 - n^(-alf3)/gamma(1-alf3)*wh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf *sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
        RHS2 = RHS2 + taualf3*sum(exp(clambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*cw{kk}.*wy2{kk});
    end
    uu1 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs2 = RHS2 - taualf3*hv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1;rhs2]);
        J = arrL - E*Hv(uu1(1),uu1(2),uu1(3),t(n+1));

        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);  wh(n+1) = uu1(3);
    u2 = uu1(1);  v2 = uu1(2);  w2 = uu1(3);

    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_p2_correction(len,vy3,vy33,vy4,....
       vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,blmd4,v0,v1,v2,bcw1);
    [wy3,wy33,wy4,wy44,wyy4,wy6,wy66,wy77] = update_p2_correction(len,wy3,wy33,wy4,....
       wy44,wyy4,wy6,wy66,wy77,b2,b3,d3,dd3,BB,n,clmd1,clmd2,clmd3,clmd4,w0,w1,w2,ccw1);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2; w0 = w1; w1 = w2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));  we = wv(t(1:nt+1)); 
eeu = abs(ue-uh);    eev = abs(ve-vh);    eew = abs(we-wh); 
eeu;
%END

function [uh,vh,wh,eeu,eev,eew] = fode_fast_conv_p2_gauss_correction_system_33......
    (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T,alf0,Hv,sgm,sgmv,sgmw,BB,T0)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + mu13*w + f(u,v,w,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + mu23*w + g(u,v,w,t)
% D_{0,t}^{gamma}w = mu31*u + mu32*v + mu33*w + g(u,v,w,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
% % % % % % % % % % % % % % % % % % % % % % % % 
alf3 = alf0(3);   taualf3 = tau^alf3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;       
uh = zeros(1,nt+1);    vh = uh;      wh = uh;
uh0 = ut0;             vh0 = vt0;    wh0 = wt0;
uh(1) = uh0;           vh(1) = vh0;  wh(1) = wh0;
E = diag([taualf,taualf1,taualf3]);
%--------------------------------------------
[uh5,vh5,wh5] = ini_system_fode_nonlinear_newton_system_3....
    (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T0,alf0,Hv,sgm,sgmv,sgmw);
nt0 = length(uh5);
% ------------------------------------------
no1 = length(sgm); no2 = length(sgmv);  no3 = length(sgmw); 
no = max(max([no1,no2,no3]),2);
% use one correction term with smaller stepsize to get the starting values
if isempty(uv)
    N = nt;
    if isempty(sgm)
        [uh00,vh00,wh00] = ini_system_fode_nonlinear_newton_system_3.......
            (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau/N,mu,no*tau,alf0,Hv,[],[],[]);
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        wh(1:no+1) = wh00(1:N:end);
        uhm0 = 0; vhm0 = 0;  whm0 = 0;
    else
        [uh00,vh00,wh00] = ini_system_fode_nonlinear_newton_system_3.......
            (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1),sgmw(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        wh(1:no+1) = wh00(1:N:end);
        uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
        whm0 = wh(2:no3+1) - wh(1);
    end
else
    uh(1:no+1) = uv(t(1:no+1));  vh(1:no+1) = vv(t(1:no+1));
    wh(1:no+1) = wv(t(1:no+1));
    uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
    whm0 = wh(2:no3+1) - wh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isempty(uv) 
    uv = @(t) 0.*t;  vv = @(t) 0.*t;  wv = @(t) 0.*t; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0] = ini_contour_p2_gauss(tau,-alf,BB,nt);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2(tau,-alf,BB,nt);

[bw,blambda,blmd1,blmd2,blmd3,blmd4,bphi0,bphi1,bphi2] = ini_contour_p2_gauss(tau,-alf1,BB,nt);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
[~,~,blmd11,blmd22,blmd33,blmd44] = ini_contour_p2_gauss_2(tau,-alf1,BB,nt);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[cw,clambda,clmd1,clmd2,clmd3,clmd4,cphi0,cphi1,cphi2] = ini_contour_p2_gauss(tau,-alf3,BB,nt);
cphi0 = tau^alf3*cphi0;  cphi1 = tau^alf3*cphi1;  cphi2 = tau^alf3*cphi2;
[~,~,clmd11,clmd22,clmd33,clmd44] = ini_contour_p2_gauss_2(tau,-alf3,BB,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};  ccw1{kk} = acw1{kk};
    acw2{kk} = acw1{kk};         bcw2{kk} = acw1{kk};  ccw2{kk} = acw1{kk};
end

T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
for k = 1:no3
    s = sgmw(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0w{k} = x0; w0w{k} = w0;
end

[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_11(tau,T0,blambda{kk}(n),blmd2{kk}(n),blmd3{kk}(n),....
                blmd4{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
            tmp = correction_weights_11(tau,T0,clambda{kk}(n),clmd2{kk}(n),clmd3{kk}(n),....
                clmd4{kk}(n),x0w,w0w,sgmw,xx0,ww0);
            wwww2_c(:,n) = whm0*tmp;
            
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,blambda{kk}(n),blmd22{kk}(n),blmd33{kk}(n),....
                blmd44{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww4(:,n) = vhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,clambda{kk}(n),clmd22{kk}(n),clmd33{kk}(n),....
                clmd44{kk}(n),x0w,w0w,sgmw,xx0,ww0);
            wwww4_c(:,n) = whm0*tmp;
        end
        acw2{kk} = wwww3;   bcw2{kk} = wwww4;  ccw2{kk} = wwww4_c;
        acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};  ccw1{kk} = acw1{kk};
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
        ccw1{kk}(1:N0-1,:) = wwww2_c(1:N0-1,:);
    end
end
clear wwww1 wwww2 wwww3 wwww4 wwww2_c wwww4_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 

wy2 = yz0; wy3 = yz0;  wy4 = yz0; wyy4 = yz0;
wy33 = yz0; wy44 = yz0;  wy66 = yz0; wy77 = yz0;  wy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  w0 = wh(1);
u1 = uh(2);  v1 = vh(2);  w1 = wh(2);
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);   w2 = wh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
vy2{1} = blmd11{1}.*vy2{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2;
wy2{1} = clmd11{1}.*wy2{1} + clmd22{1}*w0 + clmd33{1}*w1 + clmd44{1}*w2; 
uy2{1} = uy2{1} + acw2{1}(1,:);   
vy2{1} = vy2{1} + bcw2{1}(1,:);
wy2{1} = wy2{1} + ccw2{1}(1,:);
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    vyy4{kk} = blmd11{kk}.*vyy4{kk} + blmd22{kk}*v0 + blmd33{kk}*v1 + blmd44{kk}*v2;
    wyy4{kk} = clmd11{kk}.*wyy4{kk} + clmd22{kk}*w0 + clmd33{kk}*w1 + clmd44{kk}*w2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
    vyy4{kk} = vyy4{kk} + bcw2{kk}(1,:);
    wyy4{kk} = wyy4{kk} + ccw2{kk}(1,:);
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
    vy4{1} = blmd11{1}.*vy4{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2 + bcw2{1}(1,:);
    vy44{1} = vy4{1};
    wy4{1} = clmd11{1}.*wy4{1} + clmd22{1}*w0 + clmd33{1}*w1 + clmd44{1}*w2 + ccw2{1}(1,:);
    wy44{1} = wy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
    vy4{2} = blmd11{2}.*vy4{2} + blmd22{2}*v0 + blmd33{2}*v1 + blmd44{2}*v2 + bcw2{2}(1,:);
    wy4{2} = clmd11{2}.*wy4{2} + clmd22{2}*w0 + clmd33{2}*w1 + clmd44{2}*w2 + ccw2{2}(1,:);
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2; w0 = w1; w1 = w2;
clear almd11 almd22 almd33 almd44 blmd11 blmd22 blmd33 blmd44 acw2 bcw2; 
clear clmd11 clmd22 clmd33 clmd44 ccw2; 
%
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);  w2 = wh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    wy2{1} = clmd1{1}.*wy2{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
        wyy4{kk} = clmd1{kk}.*wyy4{kk} + clmd2{kk}*w0 + clmd3{kk}*w1 + clmd4{kk}*w2 + ccw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
        wy4{1} = clmd1{1}.*wy4{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
        wy44{1} = wy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
        wy4{2} = clmd1{2}.*wy4{2} + clmd2{2}*w0 + clmd3{2}*w1 + clmd4{2}*w2 + ccw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  v0 = v1;  v1 = v2;  w0 = w1;  w1 = w2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
% arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
arrL = diag([aphi2,bphi2,cphi2]) - diag([taualf,taualf1,taualf3])*mu;
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum(exp(blambda{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1});
    RHS2 = cphi0*w0 + cphi1*w1 - n^(-alf3)/gamma(1-alf3)*wh0.....
        + taualf3*sum(exp(clambda{1}*(t(n+1)-btau(1)-nT0(1))).*cw{1}.*wy2{1});
    uu1 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs2 = RHS2 - taualf3*hv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1;rhs2]);
        J = arrL - E*Hv(uu1(1),uu1(2),uu1(3),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);  wh(n+1) = uu1(3);
    u2 = uu1(1); v2 = uu1(2);  w2 = uu1(3);
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    wy2{1} = clmd1{1}.*wy2{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
        wyy4{kk} = clmd1{kk}.*wyy4{kk} + clmd2{kk}*w0 + clmd3{kk}*w1 + clmd4{kk}*w2 + ccw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
        wy4{1} = clmd1{1}.*wy4{1} + clmd2{1}*w0 + clmd3{1}*w1 + clmd4{1}*w2 + ccw1{1}(n-1,:);
        wy44{1} = wy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
        wy4{2} = clmd1{2}.*wy4{2} + clmd2{2}*w0 + clmd3{2}*w1 + clmd4{2}*w2 + ccw1{2}(n-1,:);
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;  w0 = w1;  w1 = w2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    [wy2,wy3,wy4,b2,b3,d3,dd3] = update_1(wy2,wy3,wy4,wy6,wy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    RHS2 = cphi0*w0 + cphi1*w1 - n^(-alf3)/gamma(1-alf3)*wh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf *sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
        RHS2 = RHS2 + taualf3*sum(exp(clambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*cw{kk}.*wy2{kk});
    end
    uu1 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs2 = RHS2 - taualf3*hv(uu1(1),uu1(2),uu1(3),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1;rhs2]);
        J = arrL - E*Hv(uu1(1),uu1(2),uu1(3),t(n+1));
%         if cond(J)>1e10
%             n,nt,uu1,J,return;
%         end
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);  wh(n+1) = uu1(3);
    u2 = uu1(1);  v2 = uu1(2);  w2 = uu1(3);

    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_p2_correction(len,vy3,vy33,vy4,....
       vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,blmd4,v0,v1,v2,bcw1);
    [wy3,wy33,wy4,wy44,wyy4,wy6,wy66,wy77] = update_p2_correction(len,wy3,wy33,wy4,....
       wy44,wyy4,wy6,wy66,wy77,b2,b3,d3,dd3,BB,n,clmd1,clmd2,clmd3,clmd4,w0,w1,w2,ccw1);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2; w0 = w1; w1 = w2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));  we = wv(t(1:nt+1)); 
eeu = abs(ue-uh);    eev = abs(ve-vh);    eew = abs(we-wh); 
eeu;
%END

function [uh,vh,eeu,eev] ....
    = fode_fast_conv_p2_gauss_correction(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end 
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(max(no1,no2),2);
% use one correction term with smaller stepsize to get the starting values
if isempty(uv)
    N = nt;
    if isempty(sgm)
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,[],[]);
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        uhm0 = 0; vhm0 = 0;
    else
        [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
            (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
        uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
        uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
    end
else
    uh(1:no+1) = uv(t(1:no+1));  vh(1:no+1) = vv(t(1:no+1));
    uhm0 = uh(2:no1+1) - uh(1);  vhm0 = vh(2:no2+1) - vh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0] = ini_contour_p2_gauss(tau,-alf,BB,nt);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2(tau,-alf,BB,nt);

[bw,blambda,blmd1,blmd2,blmd3,blmd4,bphi0,bphi1,bphi2] = ini_contour_p2_gauss(tau,-alf1,BB,nt);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
[~,~,blmd11,blmd22,blmd33,blmd44] = ini_contour_p2_gauss_2(tau,-alf1,BB,nt);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};
    acw2{kk} = acw1{kk};  bcw2{kk} = acw1{kk};
end

T0 = 1; 
N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
for k = 1:no2
    s = sgmv(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64); 
    x0v{k} = x0; w0v{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_11(tau,T0,blambda{kk}(n),blmd2{kk}(n),blmd3{kk}(n),....
                blmd4{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww2(:,n) = vhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,blambda{kk}(n),blmd22{kk}(n),blmd33{kk}(n),....
                blmd44{kk}(n),x0v,w0v,sgmv,xx0,ww0);
            wwww4(:,n) = vhm0*tmp;
        end
        acw2{kk} = wwww3;  bcw2{kk} = wwww4;
        acw1{kk} = zeros(nt,lenN);   bcw1{kk} = acw1{kk};
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  bcw1{kk}(1:N0-1,:) = wwww2(1:N0-1,:);
    end
end
clear wwww1 wwww2 wwww3 wwww4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  u1 = uh(2);  v1 = vh(2);
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
vy2{1} = blmd11{1}.*vy2{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
vy2{1} = vy2{1} + bcw2{1}(1,:);
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    vyy4{kk} = blmd11{kk}.*vyy4{kk} + blmd22{kk}*v0 + blmd33{kk}*v1 + blmd44{kk}*v2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
    vyy4{kk} = vyy4{kk} + bcw2{kk}(1,:);
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
    vy4{1} = blmd11{1}.*vy4{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2 + bcw2{1}(1,:);
    vy44{1} = vy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
    vy4{2} = blmd11{2}.*vy4{2} + blmd22{2}*v0 + blmd33{2}*v1 + blmd44{2}*v2 + bcw2{2}(1,:);
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
clear almd11 almd22 almd33 almd44 blmd11 blmd22 blmd33 blmd44 acw2 bcw2; 
%
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum(exp(blambda{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1});
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1); v2 = uu1(2);
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2 + bcw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2 + bcw1{1}(n-1,:);
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2 + bcw1{2}(n-1,:);
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
        RHS1 = RHS1 + taualf1*sum(exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk});
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
%         if cond(J)>1e10
%             n,nt,uu1,J,return;
%         end
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1);  v2 = uu1(2);

    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,....
        uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = update_p2_correction(len,vy3,vy33,vy4,....
       vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,blmd4,v0,v1,v2,bcw1);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
eeu;
%END

function [uh,vh,eeu,eev] ....
    = system_fode_nonlinear_fast_conv_p2(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied, quatratic
% inerpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end 
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(no1,no2) + 1;
% use one correction term with smaller stepsize to get the starting values
if no > 1
    N = nt;
    [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
        (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
    uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
else
     N = nt; no = 2;
    [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
        (uv,vv,fv,gv,ut0,vt0,tau/N,mu,2*tau,alf0,Hv,[],[]);
    uh(1:3) = uh00(1:N:end);  vh(1:3) = vh00(1:N:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% [aw,lambda,lmd1,lmd2,lmd3,aphi1,aphi2] = ini_contour(tau,-alf,BB,nt);
[aw,lambda,lmd1,lmd2,lmd3,lmd4,aphi0,aphi1,aphi2] = ini_contour_p2(tau,-alf,BB,nt);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,lmd11,lmd22,lmd33,lmd44] = ini_contour_p2_2(tau,-alf,BB,nt);
% aphi00 = tau^(-alf)*aphi00;  aphi11 = tau^(-alf)*aphi11; aphi22 = tau^(-alf)*aphi22;

% [bw,lambda,lmd1,lmd2,lmd3,bphi1,bphi2] = ini_contour(tau,-alf1,BB,nt);
[bw,lambda,lmd1,lmd2,lmd3,lmd4,bphi0,bphi1,bphi2] = ini_contour_p2(tau,-alf1,BB,nt);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
% [~,~,lmd11,lmd22,lmd33,lmd44,bphi00,bphi11,bphi22] = ini_contour_p2_2(tau,-alf1,BB,nt);
% bphi00 = tau^(-alf1)*bphi00;  bphi11 = tau^(-alf1)*bphi11;  bphi22 = tau^(-alf1)*bphi22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  u1 = uh(2);  v1 = vh(2);
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);
uy2{1} = lmd11{1}.*uy2{1} + lmd22{1}*u0 + lmd33{1}*u1 + lmd44{1}*u2;
vy2{1} = lmd11{1}.*vy2{1} + lmd22{1}*v0 + lmd33{1}*v1 + lmd44{1}*v2;
for kk = 1:len
    uyy4{kk} = lmd11{kk}.*uyy4{kk} + lmd22{kk}*u0 + lmd33{kk}*u1 + lmd44{kk}*u2;
    vyy4{kk} = lmd11{kk}.*vyy4{kk} + lmd22{kk}*v0 + lmd33{kk}*v1 + lmd44{kk}*v2;
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = lmd11{1}.*uy4{1} + lmd22{1}*u0 + lmd33{1}*u1 + lmd44{1}*u2;
    uy44{1} = uy4{1};
    vy4{1} = lmd11{1}.*vy4{1} + lmd22{1}*v0 + lmd33{1}*v1 + lmd44{1}*v2;
    vy44{1} = vy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = lmd11{2}.*uy4{2} + lmd22{2}*u0 + lmd33{2}*u1 + lmd44{2}*u2;
    vy4{2} = lmd11{2}.*vy4{2} + lmd22{2}*v0 + lmd33{2}*v1 + lmd44{2}*v2;
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
clear lmd11 lmd22 lmd33 lmd44; 
% ------------------------------------------------------
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);
    uy2{1} = lmd1{1}.*uy2{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2;
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2;
    for kk = 1:len
        uyy4{kk} = lmd1{kk}.*uyy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v0 + lmd3{kk}*v1 + lmd4{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd1{1}.*uy4{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd1{2}.*uy4{2} + lmd2{2}*u0 + lmd3{2}*u1 + lmd4{2}*u2;
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v0 + lmd3{2}*v1 + lmd4{2}*v2;
    end
    u0 = u1;  u1 = u2; v0 = v1; v1 = v2;
end
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*aw{1}.*uy2{1}));
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum(imag(exp(lambda{1}*(t(n+1)-btau(1))).*bw{1}.*vy2{1}));
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1); v2 = uu1(2);
    
    uy2{1} = lmd1{1}.*uy2{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2;
    vy2{1} = lmd1{1}.*vy2{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2;
    for kk = 1:len
        uyy4{kk} = lmd1{kk}.*uyy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        vyy4{kk} = lmd1{kk}.*vyy4{kk} + lmd2{kk}*v0 + lmd3{kk}*v1 + lmd4{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = lmd1{1}.*uy4{1} + lmd2{1}*u0 + lmd3{1}*u1 + lmd4{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = lmd1{1}.*vy4{1} + lmd2{1}*v0 + lmd3{1}*v1 + lmd4{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = lmd1{2}.*uy4{2} + lmd2{2}*u0 + lmd3{2}*u1 + lmd4{2}*u2;
        vy4{2} = lmd1{2}.*vy4{2} + lmd2{2}*v0 + lmd3{2}*v1 + lmd4{2}*v2;
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*aw{kk}.*uy2{kk}));
        RHS1 = RHS1 + taualf1*sum(imag(exp(lambda{kk}*(t(n+1)-btau(kk))).*bw{kk}.*vy2{kk}));
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1);  v2 = uu1(2);

%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,u1,u2);
%     [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
%         (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2......
       (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2_p2......
       (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,lmd4,v0,v1,v2);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
eeu;
%END

function [uh,vh,eeu,eev] ....
    = fode_fast_conv_p2_gauss(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv,BB)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% Correction terms are applied, Newton method is applied, quatratic
% inerpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t;
end 
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);  
taualf = tau^alf;        taualf1 = tau^alf1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%--------------------------------------------
no1 = length(sgm); no2 = length(sgmv);
no = max(no1,no2) + 1;
% use one correction term with smaller stepsize to get the starting values
if no > 1
    N = nt;
    [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
        (uv,vv,fv,gv,ut0,vt0,tau/N,mu,no*tau,alf0,Hv,sgm(1),sgmv(1));
    uh(1:no+1) = uh00(1:N:end);  vh(1:no+1) = vh00(1:N:end);
else
     N = nt; no = 2;
    [uh00,vh00] = ini_system_fode_nonlinear_newton_1.......
        (uv,vv,fv,gv,ut0,vt0,tau/N,mu,2*tau,alf0,Hv,[],[]);
    uh(1:3) = uh00(1:N:end);  vh(1:3) = vh00(1:N:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% [aw,lambda,lmd1,lmd2,lmd3,aphi1,aphi2] = ini_contour(tau,-alf,BB,nt);
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0] = ini_contour_p2_gauss(tau,-alf,BB,nt);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2(tau,-alf,BB,nt);
% aphi00 = tau^(-alf)*aphi00;  aphi11 = tau^(-alf)*aphi11; aphi22 = tau^(-alf)*aphi22;

% [bw,lambda,lmd1,lmd2,lmd3,bphi1,bphi2] = ini_contour(tau,-alf1,BB,nt);
[bw,blambda,blmd1,blmd2,blmd3,blmd4,bphi0,bphi1,bphi2] = ini_contour_p2_gauss(tau,-alf1,BB,nt);
bphi0 = tau^alf1*bphi0;  bphi1 = tau^alf1*bphi1;  bphi2 = tau^alf1*bphi2;
[~,~,blmd11,blmd22,blmd33,blmd44] = ini_contour_p2_gauss_2(tau,-alf1,BB,nt);
% bphi00 = tau^(-alf1)*bphi00;  bphi11 = tau^(-alf1)*bphi11;  bphi22 = tau^(-alf1)*bphi22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 

vy2 = yz0; vy3 = yz0;  vy4 = yz0; vyy4 = yz0;
vy33 = yz0; vy44 = yz0;  vy66 = yz0; vy77 = yz0;  vy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);  v0 = vh(1);  u1 = uh(2);  v1 = vh(2);
n = 2;    u2 = uh(n+1);    v2 = vh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
vy2{1} = blmd11{1}.*vy2{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2;
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    vyy4{kk} = blmd11{kk}.*vyy4{kk} + blmd22{kk}*v0 + blmd33{kk}*v1 + blmd44{kk}*v2;
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
    uy44{1} = uy4{1};
    vy4{1} = blmd11{1}.*vy4{1} + blmd22{1}*v0 + blmd33{1}*v1 + blmd44{1}*v2;
    vy44{1} = vy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2;
    vy4{2} = blmd11{2}.*vy4{2} + blmd22{2}*v0 + blmd33{2}*v1 + blmd44{2}*v2;
end
u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
clear almd11 almd22 almd33 almd44 blmd11 blmd22 blmd33 blmd44; 
% ------------------------------------------------------
for n = 3:no
    u2 = uh(n+1);    v2 = vh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2;
    end
    u0 = u1;  u1 = u2; v0 = v1; v1 = v2;
end
arrL = [aphi2 - a11*taualf, -a12*taualf; -a21*taualf1, bphi2 - a22*taualf1];
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum((exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1}));
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0.....
        + taualf1*sum((exp(blambda{1}*(t(n+1)-btau(1)-nT0(1))).*bw{1}.*vy2{1}));
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    
    uh(n+1) = uu1(1);    vh(n+1) = uu1(2);
    u2 = uu1(1); v2 = uu1(2);
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    vy2{1} = blmd1{1}.*vy2{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
        vyy4{kk} = blmd1{kk}.*vyy4{kk} + blmd2{kk}*v0 + blmd3{kk}*v1 + blmd4{kk}*v2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
        uy44{1} = uy4{1};
        vy4{1} = blmd1{1}.*vy4{1} + blmd2{1}*v0 + blmd3{1}*v1 + blmd4{1}*v2;
        vy44{1} = vy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
        vy4{2} = blmd1{2}.*vy4{2} + blmd2{2}*v0 + blmd3{2}*v1 + blmd4{2}*v2;
    end
    u0 = u1; u1 = u2; v0 = v1; v1 = v2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [vy2,vy3,vy4,b2,b3,d3,dd3] = update_1(vy2,vy3,vy4,vy6,vy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    RHS1 = bphi0*v0 + bphi1*v1 - n^(-alf1)/gamma(1-alf1)*vh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum((exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk}));
        RHS1 = RHS1 + taualf1*sum((exp(blambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*bw{kk}.*vy2{kk}));
    end
    uu1 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1(1),uu1(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(uu1(1),uu1(2),t(n+1));
        rhs = -(arrL*uu1 + [rhs0;rhs1]);
        J = arrL - E*Hv(uu1(1),uu1(2),t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1(1);  vh(n+1) = uu1(2);
    u2 = uu1(1);  v2 = uu1(2);

%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,u1,u2);
%     [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2......
%         (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,lmd1,lmd2,lmd3,v1,v2);
    
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2......
       (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2);
    [vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77] = contour_2_p2......
       (len,vy3,vy33,vy4,vy44,vyy4,vy6,vy66,vy77,b2,b3,d3,dd3,BB,n,blmd1,blmd2,blmd3,blmd4,v0,v1,v2);
    u0 = u1; u1 = u2;  v0 = v1; v1 = v2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
eeu;
%END

function [uh,vh,eeu,eev] .....
    = ini_system_fode_nonlinear_newton_1(uv,vv,fv,gv,ut0,vt0,tau,mu,T,alf0,Hv,sgm,sgmv)
% D_{0,t}^{alpha}u = mu11*u + mu12*v + f(u,v,t)
% D_{0,t}^{beta}v  = mu21*u + mu22*v + g(u,v,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
a11 = mu(1,1);  a12 = mu(1,2);  a21 = mu(2,1);  a22 = mu(2,2);
alf = alf0(1);  alf1 = alf0(2);
if isempty(sgm)
    [w0,B0,~,C] = my_weght_2(ne,alf,[]);
    [w10,B10,~,C1] = my_weght_2(ne,alf1,[]);
else
    [w0,B0,~,C] = my_weght_2(ne,alf,sgm(1));
    [w10,B10,~,C1] = my_weght_2(ne,alf1,sgmv(1));
end
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
taualf = tau^alf;  taualf1 = tau^alf1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = [w(1) - a11*taualf, -a12*taualf; -a21*taualf1, w1(1) - a22*taualf1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;
uh0 = ut0;             vh0 = vt0;
uh(1) = uh0;           vh(1) = vh0;  
E = diag([taualf,taualf1]);
%----------------  first step ---------------------
no = 1;  n = 1;     
RHS0 = -B(n+1)*uh0 + uh(1:n)*w(n+1:-1:2) - uh0*C(n+1);
RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) - vh0*C1(n+1);
u2 = [uh(n);vh(n)];
for kk = 1:400
    rhs0 = RHS0 - taualf*fv(u2(1),u2(2),t(n+1));
    rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),t(n+1));
    rhs = -((arrL+diag([C(n+1),C1(n+1)]))*u2 + [rhs0;rhs1]);
    J = arrL - E*Hv(u2(1),u2(2),t(n+1)) + diag([C(n+1),C1(n+1)]);
    uu2 = u2 + J\rhs;
    if norm(u2-uu2,inf)<1e-15
        u2 = uu2; break;
    end
    u2 = uu2;
end
uh(n+1) = u2(1);  vh(n+1) = u2(2);
uh00 = uh(2:no+1) - uh0;     vh00 = vh(2:no+1) - vh0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0 + uh(1:n)*w(n+1:-1:2) + uh00*C(:,n+1);
    RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) + vh00*C1(:,n+1);
    u2 = [uh(n);vh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2(1),u2(2),t(n+1));
        rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),t(n+1));
        rhs = -(arrL*u2 + [rhs0;rhs1]);
        J = arrL - E*Hv(u2(1),u2(2),t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    
end
%----------------------------------------------------- 
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));
eeu = abs(ue-uh);    eev = abs(ve-vh);  
%END

function [uh,vh,wh,eeu,eev,eew] .....
    = ini_system_fode_nonlinear_newton_system_3....
    (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T,alf0,Hv,sgm,sgmv,sgmw)
% D_{0,t}^{alpha1}u = mu11*u + mu12*v + mu13*w + f(u,v,w,t)
% D_{0,t}^{alpha2}v = mu21*u + mu22*v + mu23*w + g(u,v,w,t)
% D_{0,t}^{alpha3}w = mu31*u + mu32*v + mu33*w + h(u,v,w,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t; wv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
alf = alf0(1);  alf1 = alf0(2);  alf3 = alf0(3);
if isempty(sgm)
    [w0,B0,~,C] = my_weght_2(ne,alf,[]);
    [w10,B10,~,C1] = my_weght_2(ne,alf1,[]);
    [w30,B30,~,C3] = my_weght_2(ne,alf3,[]);
else
    [w0,B0,~,C] = my_weght_2(ne,alf,sgm(1));
    [w10,B10,~,C1] = my_weght_2(ne,alf1,sgmv(1));
    [w30,B30,~,C3] = my_weght_2(ne,alf3,sgmw(1));
end
w = w0;     B = B0;  
w1 = w10;   B1 = B10;  
w3 = w30;   B3 = B30; 
taualf = tau^alf;  taualf1 = tau^alf1;  taualf3 = tau^alf3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% arrL = [w(1) - a11*taualf, -a12*taualf, -a13*taualf;....
%     -a21*taualf1, w1(1) - a22*taualf1,- a23*taualf1;....
%     -a31*taualf3, -a32*taualf3, w3(1) -  a33*taualf3];
arrL = diag([w(1),w1(1),w3(1)]) - diag([taualf,taualf1,taualf3])*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    vh = uh;       wh = uh;
uh0 = ut0;             vh0 = vt0;     wh0 = wt0;
uh(1) = uh0;           vh(1) = vh0;   wh(1) = wh0; 
E = diag([taualf,taualf1,taualf3]);
%----------------  first step ---------------------
no = 1;  n = 1;     
RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  - uh0*C(n+1);
RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) - vh0*C1(n+1);
RHS3 = -B3(n+1)*wh0 + wh(1:n)*w3(n+1:-1:2) - wh0*C3(n+1);
u2 = [uh(n);vh(n);wh(n)];
for kk = 1:400
    rhs0 = RHS0  - taualf*fv(u2(1),u2(2),u2(3),t(n+1));
    rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),u2(3),t(n+1));
    rhs3 = RHS3 - taualf3*hv(u2(1),u2(2),u2(3),t(n+1));
    rhs = -((arrL+diag([C(n+1),C1(n+1),C3(n+1)]))*u2 + [rhs0;rhs1;rhs3]);
    J = arrL - E*Hv(u2(1),u2(2),u2(3),t(n+1)) + diag([C(n+1),C1(n+1),C3(n+1)]);
    uu2 = u2 + J\rhs;
    if norm(u2-uu2,inf)<1e-15
        u2 = uu2; break;
    end
    u2 = uu2;
end
uh(n+1) = u2(1);  vh(n+1) = u2(2);  wh(n+1) = u2(3);
uh00 = uh(2:no+1) - uh0;     vh00 = vh(2:no+1) - vh0;
wh00 = wh(2:no+1) - wh0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) + vh00*C1(:,n+1);
    RHS3 = -B3(n+1)*wh0 + wh(1:n)*w3(n+1:-1:2) + wh00*C3(:,n+1);
    u2 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2(1),u2(2),u2(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),u2(3),t(n+1));
        rhs3 = RHS3 - taualf3*hv(u2(1),u2(2),u2(3),t(n+1));
        rhs = -(arrL*u2 + [rhs0;rhs1;rhs3]);
        J = arrL - E*Hv(u2(1),u2(2),u2(3),t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    wh(n+1) = u2(3);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));  we = wv(t(1:nt+1)); 
eeu = abs(ue-uh);    eev = abs(ve-vh);    eew = abs(we-wh);  
%END

function [uh,vh,wh,eeu,eev,eew] = system_fode_nonlinear_newton_system_3...
    (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T,alf0,Hv,sgm,sgmv,sgmw)
% D_{0,t}^{alpha1}u = mu11*u + mu12*v + mu13*w + f(u,v,w,t)
% D_{0,t}^{alpha2}v = mu21*u + mu22*v + mu23*w + g(u,v,w,t)
% D_{0,t}^{alpha3}w = mu31*u + mu32*v + mu33*w + h(u,v,w,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ne = round(T/tau)+1;  
alf = alf0(1);  alf1 = alf0(2);  alf3 = alf0(3);
if isempty(sgm) && isempty(sgmv) && isempty(sgmw)
     [uh,vh,wh,eeu,eev,eew] = ini_system_fode_nonlinear_newton_system_3....
        (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau,mu,T,alf0,Hv,[],[],[]);
    return;
end
     
t = 0:tau:T+1e-10;     nt = length(t)-1;     
no1 = length(sgm);  no2 = length(sgmv);  no3 = length(sgmw);
no = max([no1,no2,no3]);
if isempty(uv)
[uh5,vh5,wh5] = ini_system_fode_nonlinear_newton_system_3....
        (uv,vv,wv,fv,gv,hv,ut0,vt0,wt0,tau/nt,mu,no*tau,alf0,Hv,sgm(1),sgmv(1),sgmw(1));
end
     
[w0,B0,~,C] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C1] = my_weght_2(ne,alf1,sgmv);
[w30,B30,~,C3] = my_weght_2(ne,alf3,sgmw);

w = w0;     B = B0;
w1 = w10;   B1 = B10;
w3 = w30;   B3 = B30; 
taualf = tau^alf;  taualf1 = tau^alf1;  taualf3 = tau^alf3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = diag([w(1),w1(1),w3(1)]) - diag([taualf,taualf1,taualf3])*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
uh = zeros(1,nt+1);    vh = uh;       wh = uh;
uh0 = ut0;             vh0 = vt0;     wh0 = wt0;
if isempty(uv)
    uh(1:no+1) = uh5(1:nt:end);
    vh(1:no+1) = vh5(1:nt:end);
    wh(1:no+1) = wh5(1:nt:end);
else
    uh(1:no+1) = uv(t(1:no+1));
    vh(1:no+1) = vv(t(1:no+1));
    wh(1:no+1) = wv(t(1:no+1));
end
E = diag([taualf,taualf1,taualf3]);
if isempty(uv)
    uv = @(t) 0.*t;  vv = @(t) 0.*t; wv = @(t) 0.*t;
end
%-------------------------------------
uh00 = uh(2:no1+1) - uh0;     vh00 = vh(2:no2+1) - vh0;
wh00 = wh(2:no3+1) - wh0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    RHS1 = -B1(n+1)*vh0 + vh(1:n)*w1(n+1:-1:2) + vh00*C1(:,n+1);
    RHS3 = -B3(n+1)*wh0 + wh(1:n)*w3(n+1:-1:2) + wh00*C3(:,n+1);
    u2 = [uh(n);vh(n);wh(n)];
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2(1),u2(2),u2(3),t(n+1));
        rhs1 = RHS1 - taualf1*gv(u2(1),u2(2),u2(3),t(n+1));
        rhs3 = RHS3 - taualf3*hv(u2(1),u2(2),u2(3),t(n+1));
        rhs = -(arrL*u2 + [rhs0;rhs1;rhs3]);
        J = arrL - E*Hv(u2(1),u2(2),u2(3),t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2(2);    wh(n+1) = u2(3);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));  ve = vv(t(1:nt+1));  we = wv(t(1:nt+1)); 
eeu = abs(ue-uh);    eev = abs(ve-vh);    eew = abs(we-wh);  
%END

function [uh,eeu] = system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = round(T/tau)+1;  alf = alf0(1); 
if isempty(sgm) 
     [uh,eeu] = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau,mu,T,alf0,Hv,[]);
    return;
end
     
t = 0:tau:T+1e-10;     nt = length(t)-1;     
no1 = length(sgm);     no = no1;
if ~isempty(uv)
    uh5 = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau/nt,mu,no*tau,alf0,Hv,sgm(1));
end 
[w0,B0,~,C] = my_weght_2(ne,alf,sgm);
w = w0;     B = B0;
taualf = tau^alf;  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
uh = zeros(1,nt+1);   
uh0 = ut0;            
if isempty(uv)
    uh(1:no+1) = uh5(1:nt:end);
else
    uh(1:no+1) = uv(t(1:no+1));
end
E = taualf;
if isempty(uv)
    uv = @(t) 0.*t;  
end
%-------------------------------------
uh00 = uh(2:no1+1) - uh0;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));    eeu = abs(ue-uh);      
%END

function [uh,eeu] = system_fode_nonlinear_newton_quadratic(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = round(T/tau)+1;  alf = alf0(1); 
if isempty(sgm) 
     [uh,eeu] = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau,mu,T,alf0,Hv,[]);
    return;
end
     
t = 0:tau:T+1e-10;     nt = length(t)-1;     
no1 = length(sgm);     no = no1;
if ~isempty(uv)
    uh5 = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau/nt,mu,no*tau,alf0,Hv,sgm(1));
end 
[w0,B0,~,C] = my_weght_2(ne,alf,sgm);
w = w0;     B = B0;
taualf = tau^alf;  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
uh = zeros(1,nt+1);   
uh0 = ut0;            
if isempty(uv)
    uh(1:no+1) = uh5(1:nt:end);
else
    uh(1:no+1) = uv(t(1:no+1));
end
E = taualf;
if isempty(uv)
    uv = @(t) 0.*t;  
end
%-------------------------------------
uh00 = uh(2:no1+1) - uh0;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));    eeu = abs(ue-uh);      
%END

function [uh,eeu] = fode_fast_conv_p2_gauss_correction_system_1......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
%--------------------------------------------
no1 = length(sgm);  no = max(no1,2);
% use one correction term with smaller stepsize to get the starting values
if isempty(uv)
    N = nt;
    if isempty(sgm)
        [uh00] = system_fode_nonlinear_newton_system_1.......
            (uv,fv,ut0,tau/N,mu,no*tau,alf0,Hv,[]);
        uh(1:no+1) = uh00(1:N:end); 
        uhm0 = 0; 
    else
        [uh00] = system_fode_nonlinear_newton_system_1.......
            (uv,v,fv,ut0,tau/N,mu,no*tau,alf0,Hv,sgm(1));
        uh(1:no+1) = uh00(1:N:end);  
        uhm0 = uh(2:no1+1) - uh(1);  
    end
else
    uh(1:no+1) = uv(t(1:no+1));  
    uhm0 = uh(2:no1+1) - uh(1);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isempty(uv) 
    uv = @(t) 0.*t;  
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
NN0 = 64;
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,tau,NN0);
aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,tau,NN0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);     acw2{kk} = acw1{kk};        
end

T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
        end
        acw2{kk} = wwww3;     acw1{kk} = zeros(nt,lenN);  
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:no
    u2 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = aphi2 - taualf*mu;
for n = no+1:2*BB-1
    b = tauL1(BB,n);   btau = tau*b;     
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0....
        + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;    
    u2 = uu1; 
    
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
    end
    u0 = u1; u1 = u2; 
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
%         if cond(J)>1e10
%             n,nt,uu1,J,return;
%         end
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;   u2 = uu1; 
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu] = fode_fast_conv_p2_gauss_correction_system_11......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,TT0,dTT0)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TT0 = 0.5;      
[~,id] = min(abs(t-TT0));       TT0 = t(id);
% dTT0 = 0.5;     
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);   dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0+dn0;   no = n0;  
% TT0,dTT0
% [aa,bb,cc,dd] = coefficients_p2(alf,n1);
% aa = aa(dn0-1:-1:1); bb = bb(dn0-1:-1:1); cc = cc(dn0-1:-1:1);
[aa2,bb2,cc2,dd2,ww2] = coefficients_p2_2(-alf,n1,sgm);
aa2 = aa2(dn0-1:-1:1); bb2 = bb2(dn0-1:-1:1); cc2 = cc2(dn0-1:-1:1);
dd2 = dd2(:,dn0);
%--------------------------------------------
no1 = length(sgm);   
% use one correction term with smaller stepsize to get the starting values
% if isempty(uv)
    if isempty(sgm)
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[]);
        [uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[]);
        uh(1:n1+1) = uh00(1:n1+1);
        uhm0 = 0;
    else
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
        [uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
        uh(1:n1+1) = uh00(1:n1+1);
        uhm0 = uh(2:no1+1) - uh(1);
    end

if isempty(sgm)
    uhm0 = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
NN0 = 64;
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
% aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,dTT0,NN0);



% % % [w,lambda,lmd11,lmd22,lmd33,lmd44,phi00,phi11,phi22]....
% % %     = initial_contour_gauss_p2_b(tau,alf,B,nN,varargin)
% % [aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
% %     = initial_contour_gauss_p2(tau,-alf,BB,nt,dTT0);
% % % aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
% % [~,~,almd11,almd22,almd33,almd44] = initial_contour_gauss_p2_b(tau,-alf,BB,nt,dTT0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);     acw2{kk} = acw1{kk};        
end
T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
        end
        acw2{kk} = wwww3;     acw1{kk} = zeros(nt,lenN);  
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:2*BB-1
    u2 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:n0   
    u2 = uh(n+1);    
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrL = aphi2 - taualf*mu;
% arrL = cc(end) + dd(3) - taualf*mu;
% cc = cc(1:end-1); dd = dd(1:2);  
% % % % % % % % % % % % % % % % % % % % 
arrL = cc2(end) - taualf*mu;   cc2 = cc2(1:end-1);
% % % % % % % % % % % % % % % % % % % 
uu = uh(n0+1:n1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = n0+1:nt-dn0+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
%      RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
%      RHS0 = sum(dd.*uu(end-1:end)) + sum(aa.*uu(1:end-1))....
%          + sum(bb.*uu(2:end)) + sum(cc.*uu(3:end))....
%          + uu(1)*(dn0)^(-alf)/gamma(1-alf) - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0;
     RHS0 = (dd2(1).*uu(1) + dd2(2).*uu(2) + dd2(3).*uu(3))...
         + sum(aa2.*uu(1:end-1)) + sum(bb2.*uu(2:end)) + sum(cc2.*uu(3:end))....
          - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1); 
    uu(1:end-1) = uu(2:end);   uu(end) = uu1;
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END


function [uh,eeu] = fode_fast_conv_p2_gauss_correction_s1......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,dTT0)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % %      
TT0 = 2*tau;   
[~,id] = min(abs(t-TT0));       TT0 = t(id);   
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);     dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0 + dn0;   
[aa2,bb2,cc2,dd2,ww2] = coefficients_p2_2(-alf,nt+1,sgm);
aa2 = aa2(dn0:-1:1); bb2 = bb2(dn0:-1:1); cc2 = cc2(dn0:-1:1);
% ww2(:,round(50/tau):end) = 0;
figure; semilogy(1:length(ww2(1,:)),abs(ww2));
%--------------------------------------------
% use one correction term with smaller stepsize to get the starting values
no1 = length(sgm);   
[uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
uh(1:n1+1) = uh00(1:n1+1);
uhm0 = 0;
if ~isempty(sgm)
    uhm0 = uh(2:no1+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
 


[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = initial_contour_gauss_p2(tau,-alf,BB,nt,dTT0);
% aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = initial_contour_gauss_p2_b(tau,-alf,BB,nt,dTT0); 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = cc2(end) - taualf*mu;   cc2 = cc2(1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
uu = uh(n0:n1); 
n = 2;     
u2 = uh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2; 
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2; 
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2;
end
u0 = u1; u1 = u2;  
%
for n = 3:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;
    RHS0 = sum(aa2.*uu(1:end-1)) + sum(bb2.*uu(2:end)) + sum(cc2.*uu(3:end))....
        - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dn0-1);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1);
    uu(1:end-1) = uu(2:end);   uu(end) = uu1;

    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
    end
    u0 = u1;  u1 = u2;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:(nt-dn0+1)      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 = sum(aa2.*uu(1:end-1)) + sum(bb2.*uu(2:end)) + sum(cc2.*uu(3:end))....
          - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dn0-1);  
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1); 
    uu(1:end-1) = uu(2:end);   uu(end) = uu1;
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu,cputime,nT] = fast_L1_method(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,dT)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % 
no = length(sgm);   
% dT = 10*tau;  
dnT = round(dT/tau);  dT = t(dnT+1);
no1 = max(no,dnT-1);  
uhm0 = 0;
uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no1*tau,alf0,Hv,sgm);
uh(1:no1+1) = uh00(1:no1+1); 
uu = uh(max(no-dnT+2,1):no1+1);
if no > 0
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% NN0 = 64;
% [aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
%     = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt,dT);
% [cc2,ww2] = coefficients_p1(-alf,nt+1,sgm);   
% [aa2,bb2,ww2] = coefficients_p1_1(-alf,nt+1,sgm);
% [aa2,bb2] = coefficients_p1_1(-alf,0,sgm);
[aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
aa3 = aa2'; bb3 = bb2';
bb20 = bb2(1);
aa2 = aa2(dnT:-1:1)';  bb2 = bb2(dnT:-1:2)';
% nt0 = round(30/tau); ww2(:,nt0:end) = 0;
s0 = 0;
for kk = 1:length(alambda)
    s0 = s0+length(alambda{kk});
end
fprintf('number of quadrature points = %-d\n',s0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
tic
d3 = tauL1(BB,2*BB);  u0 = uh(1);      
for n = 1:no-dnT+1
    u1 = uh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;
end
% % % % % % % % % % % % % % % % % % % % % % % % 
nT = T*(0.1:0.1:1)';
cputime = zeros(1,length(nT)); kkk = 1;
arrL = bb20 - taualf*mu; %
for n = max(0,no-dnT+1)+1:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;     
    RHS0 =  uu*aa2 + uu(2:end)*bb2.....
        - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-1);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
%     e = RHS0 - (uh(1:n+dnT-1)*aa3(n+dnT-1:-1:1) + uh(2:n+dnT-1)*bb3(n+dnT-1:-1:2)....
%         - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-1));
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dnT) = uu1; 
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    u1 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt-dnT+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 =  uu*aa2 + uu(2:end)*bb2.....
         - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-1);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    uh(n+dnT) = uu1;   u1 = uh(n+1); 
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
%        uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; 
  if abs(t(n+dnT)-nT(kkk))<1e-10
       cputime(kkk) = toc; kkk = kkk + 1;
   end
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
% cputime(2) = toc;
%END


function [uh,eeu,cputime,nT] = fast_L1_method_trap(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,dT)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % 
no = length(sgm);   
% dT = 10*tau;  
dnT = round(dT/tau);   
no1 = max(no,dnT-1);  
uhm0 = 0;
uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no1*tau,alf0,Hv,sgm);
uh(1:no1+1) = uh00(1:no1+1); 
uu = uh(max(no-dnT+2,1):no1+1);
if no > 0
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2] = ini_contour_trap(tau,-alf,nt);
aw = aw.*exp(dT*alambda);
% [cc2,ww2] = coefficients_p1(-alf,nt+1,sgm);   
% [aa2,bb2,ww2] = coefficients_p1_1(-alf,nt+1,sgm);
% [aa2,bb2] = coefficients_p1_1(-alf,0,sgm);
[aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
aa3 = aa2'; bb3 = bb2';  
bb20 = bb2(1);
aa2 = aa2(dnT:-1:1)';  bb2 = bb2(dnT:-1:2)';
% nt0 = round(30/tau); ww2(:,nt0:end) = 0;
s0 = length(aw);
fprintf('number of quadrature points = %-d\n',s0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uy2 = 0; 
% % % % % % % % % % % % % % % % % % % % % % % 
tic
u0 = uh(1);      
for n = 1:no-dnT+1
    u1 = uh(n+1);
    uy2 = almd1.*uy2 + almd2*u0 + almd3*u1;
    u0 = u1;
end
% % % % % % % % % % % % % % % % % % % % % % % % 
nT = T*(0.1:0.1:1)';
cputime = zeros(1,length(nT)); kkk = 1;
arrL = bb20 - taualf*mu; %
for n = max(0,no-dnT+1)+1:nt-dnT+1
    RHS0 =  uu*aa2 + uu(2:end)*bb2.....
        - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-1);
    RHS0 = RHS0 + taualf*sum(aw.*uy2);
%     e = RHS0 - (uh(1:n+dnT-1)*aa3(n+dnT-1:-1:1) + uh(2:n+dnT-1)*bb3(n+dnT-1:-1:2)....
%         - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-1));
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dnT) = uu1;  e = uu1 - uv(t(n+dnT));
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    u1 = uh(n+1);   
    uy2 = almd1.*uy2 + almd2*u0 + almd3*u1;
    u0 = u1;   
  if abs(t(n+dnT)-nT(kkk))<1e-10
       cputime(kkk) = toc; kkk = kkk + 1;
   end
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
% cputime(2) = toc;
%END


function [uh,eeu] = fast_L1_method_bak(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%--------------------------------------------
% use one correction term with smaller stepsize to get the starting values
no = length(sgm);   
uhm0 = 0;
if no > 0
    uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no*tau,alf0,Hv,sgm);
    uh(1:no+1) = uh00(1:no+1);
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% NN0 = 64;
% [aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
%     = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt);
% [cc2,ww2] = coefficients_p1(-alf,nt+1,sgm);   
% [aa2,bb2,ww2] = coefficients_p1_1(-alf,nt+1,sgm);
% [aa2,bb2] = coefficients_p1_1(-alf,0,sgm);
[aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);  u0 = uh(1);      
for n = 1:no
    u1 = uh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;
end
% % % % % % % % % % % % % % % % % % % % % % % % 
arrL = bb2(1) - taualf*mu; %
for n = no+1:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;     
    RHS0 =  aa2(1)*uh(n) - n^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;
    u1 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 =  aa2(1)*uh(n) - n^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;   u1 = uh(n+1); 
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
%        uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END


function [uh,eeu] = fast_L1_method_2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%--------------------------------------------
% use one correction term with smaller stepsize to get the starting values
no = length(sgm);   
uhm0 = 0;
if no > 0
    uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no*tau,alf0,Hv,sgm);
    uh(1:no+1) = uh00(1:no+1);
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% NN0 = 64;
% [aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
%     = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt);
% [cc2,ww2] = coefficients_p1(-alf,nt+1,sgm);   
[aa2,bb2,ww2] = coefficients_p1_1(-alf,nt+1,sgm);
[aa2,bb2] = coefficients_p1_1(-alf,0,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);  u0 = uh(1);      
for n = 1:no
    u1 = uh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;
end
% % % % % % % % % % % % % % % % % % % % % % % % 
arrL = aa2(1) - taualf*mu; %
for n = no+1:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;     
    RHS0 =  bb2*uh(n) - n^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;
    u1 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 =  bb2*uh(n) - n^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+1));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+1) = uu1;   u1 = uh(n+1); 
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
%        uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu] = fast_p2_method(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,varargin)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
% % % % % % % % % % % % % % % % % % % % % % % % % 
% varargin = (dT,vps,N)
novar = numel(varargin);   
dT = 1;  
if novar > 0
    dT = varargin{1};
end
if novar > 1
    vps = varargin{2};
end
if novar > 2
    N = varargin{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % 
no = length(sgm);   
% dT = 2*tau;  
dnT = round(dT/tau);  dT = t(dnT+1);
no1 = max(no,dnT);  
uhm0 = 0;
% uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no1*tau,alf0,Hv,sgm);

uh00 = fode_nonlinear_quadratic_p2_sisc(uv,fv,ut0,tau,mu,no1*tau,alf0,Hv,sgm);

uh(1:no1+1) = uh00(1:no1+1); 


uu = uh(1:dnT);
if no > 0
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end

% [aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt,dT);
% [aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
% bb20 = bb2(1);
% aa2 = aa2(dnT:-1:1)';  bb2 = bb2(dnT:-1:2)';
% % % % % % % % % % % % % % % % % % % % % % % 
% ini_sisc_LvXu(tau,alf,B,nN,varargin)
% varargin = (dT,vps,N)
if novar == 1
    [aw,alambda,almd1,almd2,almd3,almd4,aphi1,aphi1,aphi2,nT0]....
        = ini_sisc_LvXu(tau,-alf,BB,nt,dT);
end
if novar == 2
    [aw,alambda,almd1,almd2,almd3,almd4,aphi1,aphi1,aphi2,nT0]....
        = ini_sisc_LvXu(tau,-alf,BB,nt,dT,vps);
end
if novar == 3
    [aw,alambda,almd1,almd2,almd3,almd4,aphi1,aphi1,aphi2,nT0]....
        = ini_sisc_LvXu(tau,-alf,BB,nt,dT,vps,N);
end
[aa2,bb2,cc2,dd2,ww2] = LvXu_sisc_2016(-alf,nt+1,sgm);
cc20 = cc2(1);
% aa3 = aa2(dnT:end)';     bb3 = bb2(dnT:end)';     cc3 = cc2(dnT:end)';
aa2 = aa2(dnT-1:-1:1)';  bb2 = bb2(dnT-1:-1:1)';  cc2 = cc2(dnT-1:-1:2)';
% nt0 = round(30/tau); ww2(:,nt0:end) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0;  uy3 = yz0;   uy4 = yz0;  uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);  u0 = uh(1);     u1 = uh(2); 
% n = 2;
% u2 = uh(n+1);
% uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
% for kk = 1:len
%     uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
% end
% if n > d3(2) && n < d3(1)+1
%     uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
%     uy44{1} = uy4{1};
% end
% if n > d3(3) && n < d3(2)+1
%     uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
% end
% u0 = u1; u1 = u2;

% % % % % % % % % % % % % % % % % % % % % % % % 
arrL = cc20 + dd2(3) - taualf*mu; 
for n = 1:min(2*BB-1,nt-dnT+1)
    b = tauL1(BB,n);    btau = tau*b;     
    RHS0 = uu(1:end-1)*aa2 + uu(2:end)*bb2 + uu(3:end)*cc2....
        + dd2(1)*uu(end-1) + dd2(2)*uu(end).....
        - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-2);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
%     RHS0 = RHS0 + uh(n-1:-1:1)*aa3(1:n-1)....
%         + uh(n:-1:2)*bb3(1:n-1)+ uh(n+1:-1:3)*cc3(1:n-1); 
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-12
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dnT) = uu1; 
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    u2 = uh(n+2);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
    end
%     e = uv(t(n+dnT))-uu1;
    u0 = u1;    u1 = u2;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt-dnT+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 = uu(1:end-1)*aa2 + uu(2:end)*bb2 + uu(3:end)*cc2....
         + dd2(1)*uu(end-1) + dd2(2)*uu(end).....
         - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-2);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-12
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    uh(n+dnT) = uu1;   u2 = uh(n+2); 
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2);
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; u1 = u2;
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu] = fast_p2_method_bak(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,dT)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % 
no = length(sgm);   
% dT = 10*tau;  
dnT = round(dT/tau);  dT = t(dnT+1);
no1 = max(no,dnT);  
uhm0 = 0;
uh00 = L1_method_corection_m(uv,fv,ut0,tau,mu,no1*tau,alf0,Hv,sgm);
uh(1:no1+1) = uh00(1:no1+1); 
uu = uh(1:no1);
if no > 0
    uhm0 = uh(2:no+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end

% [aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0] = ini_contour_gauss(tau,-alf,BB,nt,dT);
% [aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
% bb20 = bb2(1);
% aa2 = aa2(dnT:-1:1)';  bb2 = bb2(dnT:-1:2)';

[aw,alambda,almd1,almd2,almd3,almd4,aphi1,aphi1,aphi2,nT0] = ini_sisc_LvXu(tau,-alf,BB,nt);
[aa2,bb2,cc2,dd2,ww2] = LvXu_sisc_2016(-alf,nt+1,sgm);
cc20 = cc2(1);
aa2 = aa2(dnT-1:-1:1)';  bb2 = bb2(dnT-1:-1:1)'; cc2 = cc2(dnT-1:-1:2)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);  u0 = uh(1);  u1 = uh(2); 
% for n = 2:no-dnT+1
%     u2 = uh(n+1); 
%     uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
%     for kk = 1:len
%         uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
%     end
%     if n > d3(2) && n < d3(1)+1
%         uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
%         uy44{1} = uy4{1};
%     end
%     if n > d3(3) && n < d3(2)+1
%         uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
%     end
%     u0 = u1; u1 = u2;
% end
% % % % % % % % % % % % % % % % % % % % % % % % 
% RHS0 = dd(1)*uh(n-1) + dd(2)*uh(n) + uh00*ww(:,n-1) - n^(-alf)/gamma(1-alf)*uh0;
% RHS0 = RHS0 + uh(1:n-1)*aa(n-1:-1:1) + uh(2:n)*bb(n-1:-1:1)...
%     + uh(3:n)*cc(n-1:-1:2);
arrL = cc20 + dd2(3) - taualf*mu; %
for n = 1:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;     
    RHS0 =  uu(1:end-1)*aa2 + uu(2:end)*bb2 + uu(3:end)*cc2...
        + dd2(1)*uu(end-1) + dd2(2)*uu(end).....
        - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-2);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dnT) = uu1; 
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    u2 = uh(n+2);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{1}*u2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{kk}*u2;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
    end
    u0 = u1;   u1 = u2; 
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt-dnT+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 =  uu(1:end-1)*aa2 + uu(2:end)*bb2 + uu(3:end)*cc2...
         + dd2(1)*uu(end-1) + dd2(2)*uu(end).....
         - (n+dnT-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dnT-2);
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uu(end);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dnT));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dnT));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uu1;
    uh(n+dnT) = uu1;   u2 = uh(n+2); 
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2);
%     [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
%         (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu] = fast_p2_method_2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,dTT0)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % %      
TT0 = 2*tau;   
[~,id] = min(abs(t-TT0));       TT0 = t(id);   
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);     dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0 + dn0;   
% [aa2,bb2,cc2,dd2,ww2] = coefficients_p2_2(-alf,nt+1,sgm);
% aa2 = aa2(dn0:-1:1); bb2 = bb2(dn0:-1:1); cc2 = cc2(dn0:-1:1);
% ww2(:,round(50/tau):end) = 0;
% figure; semilogy(1:length(ww2(1,:)),abs(ww2));
%--------------------------------------------
% use one correction term with smaller stepsize to get the starting values
no1 = length(sgm);   
[uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
uh(1:n1+1) = uh00(1:n1+1);
uhm0 = 0;
if ~isempty(sgm)
    uhm0 = uh(2:no1+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
 


[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0] = ini_sisc_LvXu(tau,-alf,BB,nt);
[aa2,bb2,cc2,dd2,ww2] = LvXu_sisc_2016(-alf,nt+1,sgm);
aa2 = aa2(dn0:-1:1); bb2 = bb2(dn0:-1:1); cc2 = cc2(dn0:-1:1);
figure; semilogy(1:length(ww2(1,:)),abs(ww2));

% [aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
%     = initial_contour_gauss_p2(tau,-alf,BB,nt,dTT0);
% % aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
% [~,~,almd11,almd22,almd33,almd44] = initial_contour_gauss_p2_b(tau,-alf,BB,nt,dTT0); 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
arrL = cc2(end) + dd2(3) - taualf*mu;   cc2 = cc2(1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
uu = uh(n0:n1); 
n = 2;     
u2 = uh(n+1);
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2; 
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2; 
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2;
end
u0 = u1; u1 = u2;  
%
for n = 3:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;
    RHS0 = sum(aa2.*uu(1:end-1)) + sum(bb2.*uu(2:end)) + sum(cc2.*uu(3:end))....
        - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dn0-1);
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1);
    uu(1:end-1) = uu(2:end);   uu(end) = uu1;

    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2;
    end
    u0 = u1;  u1 = u2;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:(nt-dn0+1)      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 = sum(aa2.*uu(1:end-1)) + sum(bb2.*uu(2:end)) + sum(cc2.*uu(3:end))....
          - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0 + uhm0*ww2(:,n+dn0-1);  
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1); 
    uu(1:end-1) = uu(2:end);   uu(end) = uu1;
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2_p2(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uu,eeu] = fast_conv_p2_gauss_correction_nonlinear_fpde......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,TT0,dTT0,uav,ubv,xa,xb,N)
% D_{0,t}^{alpha}u = mu*u_xx + f(u,x,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
% % % % % % % % % % % % % % % % % % % % % % % 
dx = (xb-xa)/N;   x = xa:dx:xb+1e-12;  
x = x';           x2 = x(1:end-1);
S = zeros(N,N+1);
for k = 1:N-1
    S(k+1,k:k+2) = [1 -2 1];
end
S(1,1:2) = [-2,1];  S(1,N) = 1; 
S(N,1) = 1;         S(:,end) = []; 
S = sparse(S);
EN1 = sparse(eye(N));
S = S/dx^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;     
t = 0:tau:T+1e-10;     nt = length(t)-1; 
% % % % % % % % % % % % % % % % % % % % % % % % % % %   
[~,id] = min(abs(t-TT0));       TT0 = t(id);   
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);   dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0+dn0;   no = n0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh = zeros(N,n1+1);    
uh0 = ut0(x2);    uh(:,1) = uh0;           
E = taualf;
% [aa,bb,cc,dd] = coefficients_p2(alf,n1);
% aa = aa(dn0-1:-1:1); bb = bb(dn0-1:-1:1); cc = cc(dn0-1:-1:1);
[aa2,bb2,cc2,dd2] = coefficients_p2_2(-alf,n1,sgm);
aa2 = aa2(dn0-1:-1:1); bb2 = bb2(dn0-1:-1:1); cc2 = cc2(dn0-1:-1:1);
dd2 = dd2(:,dn0);
%--------------------------------------------
no1 = length(sgm);   
% use one correction term with smaller stepsize to get the starting values
if isempty(sgm)
    uh00 = fpde_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[],uav,ubv,xa,xb,N);
    uh(:,1:n1+1) = uh00(:,1:n1+1);
    uhm0 = 0;
else
    uh00 = fpde_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm,uav,ubv,xa,xb,N);
    uh(:,1:n1+1) = uh00(:,1:n1+1);
    uhm0 = uh(:,2:no1+1) - uh(:,1)*ones(1,no1);
end

if isempty(sgm)
    uhm0 = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(x,t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,dTT0,NN0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(N,lenN,nt);     acw2{kk} = acw1{kk};        
end
T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n,:) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n,:) = uhm0*tmp;
        end
        acw2{kk} = wwww3;     acw1{kk} = zeros(N,lenN,nt);  
        acw1{kk}(:,:,1:N0-1) = wwww1(:,:,1:N0-1);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);     M0 = ones(N,1); 
for kk = 1:len
    yz0{kk} = 0;
    almd11{kk}=M0*almd11{kk};  almd1{kk}=M0*almd1{kk};
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(:,1);    u1 = uh(:,2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(:,n+1);    
uy2{1} = almd11{1}.*uy2{1} + u0*almd22{1} + u1*almd33{1} + u2*almd44{1};
uy2{1} = uy2{1} + acw2{1}(:,:,1);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + u0*almd22{kk} + u1*almd33{kk} + u2*almd44{kk};
    uyy4{kk} = uyy4{kk} + acw2{kk}(:,:,1);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + u0*almd22{1} + u1*almd33{1} + u2*almd44{1}...
        + acw2{1}(:,:,1);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + u0*almd22{2} + u1*almd33{2} + u2*almd44{2}....
        + acw2{2}(:,:,1);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:2*BB-1
    u2 = uh(:,n+1);   
    uy2{1} = almd1{1}.*uy2{1} + u0*almd2{1} + u1*almd3{1} + u2*almd4{1}...
        + acw1{1}(:,:,n-1);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + u0*almd2{kk} + u1*almd3{kk} + u2*almd4{kk}...
            + acw1{kk}(:,:,n-1);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + u0*almd2{1} + u1*almd3{1} + u2*almd4{1}...
            + acw1{1}(:,:,n-1);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + u0*almd2{2} + u1*almd3{2} + u2*almd4{2}....
            + acw1{2}(:,:,n-1);
    end
    u0 = u1;  u1 = u2;  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:n0   
    u2 = uh(:,n+1);    
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction_pde(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = cc2(end)*EN1 - taualf*mu*S;   cc2 = cc2(1:end-1);
aa2 = aa2'; bb2 = bb2'; cc2 = cc2';
% % % % % % % % % % % % % % % % % % % 
uu = uh(:,n0+1:n1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = n0+1:nt-dn0+1    
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 = (dd2(1).*uu(:,1) + dd2(2).*uu(:,2) + dd2(3).*uu(:,3))...
         + uu(:,1:end-1)*aa2 + uu(:,2:end)*bb2 + uu(:,3:end)*cc2....
          - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(M0*(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}).*uy2{kk},2);
    end
    uu1 = uu(:,1);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,x2,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - sparse(diag(E*Hv(uu1,x2,t(n+dn0))));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
%     uh(:,n+dn0) = uu1;    
    u2 = uu(:,2);
%     u2 = uh(:,n+1); 
    uu(:,1:end-1) = uu(:,2:end);   uu(:,end) = uu1;
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction_pde(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
% for k = 1:nt+1
%     ue = uv(x2,t(k)); 
%     eeu(k) = max(abs(ue-uh(:,k)));
% end
eeu = abs(uu1 - uv(x2,t(nt+1)));    
%END

function [uh,eeu] = fast_conv_p2_gauss_correction_nonlinear_fpde_2......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,TT0,dTT0,uav,ubv,xa,xb,N)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
% % % % % % % % % % % % % % % % % % % % % % % 
dx = (xb-xa)/N;   x = xa:dx:xb+1e-12;  
x = x';     x2 = x(2:end-1);
S = zeros(N-1,N+1);
for i=1:N-1
    S(i,i:i+2) = [1 -2 1];
end
S = sparse(S);
EN1 = sparse(eye(N-1));
S = S/dx^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(N-1,nt+1);    
uh0 = ut0(x2);    uh(:,1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TT0 = 0.5;      
[~,id] = min(abs(t-TT0));       TT0 = t(id);
% dTT0 = 0.5;     
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);   dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0+dn0;   no = n0;  
% TT0,dTT0
% [aa,bb,cc,dd] = coefficients_p2(alf,n1);
% aa = aa(dn0-1:-1:1); bb = bb(dn0-1:-1:1); cc = cc(dn0-1:-1:1);
[aa2,bb2,cc2,dd2,ww2] = coefficients_p2_2(-alf,n1,sgm);
aa2 = aa2(dn0-1:-1:1); bb2 = bb2(dn0-1:-1:1); cc2 = cc2(dn0-1:-1:1);
dd2 = dd2(:,dn0);
%--------------------------------------------
no1 = length(sgm);   
% use one correction term with smaller stepsize to get the starting values
% if isempty(uv)
    if isempty(sgm)
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[]);
        uh00 = fpde_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[],uav,ubv,xa,xb,N);
        uh(:,1:n1+1) = uh00(:,1:n1+1);
        uhm0 = 0;
    else
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
        uh00 = fpde_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm,uav,ubv,xa,xb,N);
        uh(:,1:n1+1) = uh00(:,1:n1+1);
        uhm0 = uh(:,2:no1+1) - uh(:,1)*ones(1,no1);
    end

if isempty(sgm)
    uhm0 = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(x,t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% NN0 = 64;
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
% aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,dTT0,NN0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(N-1,lenN,nt);     acw2{kk} = acw1{kk};        
end
T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n,:) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n,:) = uhm0*tmp;
        end
        acw2{kk} = wwww3;     acw1{kk} = zeros(N-1,lenN,nt);  
        acw1{kk}(:,:,1:N0-1) = wwww1(:,:,1:N0-1);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);     M0 = ones(N-1,1); 
for kk = 1:len
    yz0{kk} = 0;
    almd11{kk}=M0*almd11{kk};  almd1{kk}=M0*almd1{kk};
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(:,1);    u1 = uh(:,2);  

% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(:,n+1);    
uy2{1} = almd11{1}.*uy2{1} + u0*almd22{1} + u1*almd33{1} + u2*almd44{1};
uy2{1} = uy2{1} + acw2{1}(:,:,1);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + u0*almd22{kk} + u1*almd33{kk} + u2*almd44{kk};
    uyy4{kk} = uyy4{kk} + acw2{kk}(:,:,1);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + u0*almd22{1} + u1*almd33{1} + u2*almd44{1}...
        + acw2{1}(:,:,1);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + u0*almd22{2} + u1*almd33{2} + u2*almd44{2}....
        + acw2{2}(:,:,1);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:2*BB-1
    u2 = uh(:,n+1);   
    uy2{1} = almd1{1}.*uy2{1} + u0*almd2{1} + u1*almd3{1} + u2*almd4{1}...
        + acw1{1}(:,:,n-1);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + u0*almd2{kk} + u1*almd3{kk} + u2*almd4{kk}...
            + acw1{kk}(:,:,n-1);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + u0*almd2{1} + u1*almd3{1} + u2*almd4{1}...
            + acw1{1}(:,:,n-1);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + u0*almd2{2} + u1*almd3{2} + u2*almd4{2}....
            + acw1{2}(:,:,n-1);
    end
    u0 = u1;  u1 = u2;  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:n0   
    u2 = uh(:,n+1);    
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction_pde(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrL = aphi2 - taualf*mu;
% arrL = cc(end) + dd(3) - taualf*mu;
% cc = cc(1:end-1); dd = dd(1:2);  
% % % % % % % % % % % % % % % % % % % % 
arrL = cc2(end)*EN1 - taualf*mu*S(:,2:end-1);   cc2 = cc2(1:end-1);
bnd = taualf*mu*[S(:,1),S(:,end)];
aa2 = aa2'; bb2 = bb2'; cc2 = cc2';
% % % % % % % % % % % % % % % % % % % 
uu = uh(:,n0+1:n1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = n0+1:nt-dn0+1    
     ua = uav(t(n+dn0)); ub = ubv(t(n+dn0)); 
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
%      RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
%      RHS0 = sum(dd.*uu(end-1:end)) + sum(aa.*uu(1:end-1))....
%          + sum(bb.*uu(2:end)) + sum(cc.*uu(3:end))....
%          + uu(1)*(dn0)^(-alf)/gamma(1-alf) - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0;
     RHS0 = (dd2(1).*uu(:,1) + dd2(2).*uu(:,2) + dd2(3).*uu(:,3))...
         + uu(:,1:end-1)*aa2 + uu(:,2:end)*bb2 + uu(:,3:end)*cc2....
          - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0 - bnd*[ua;ub];
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(M0*(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}).*uy2{kk},2);
    end
    uu1 = uu(:,1);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,x2,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - sparse(diag(E*Hv(uu1,x2,t(n+dn0))));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(:,n+dn0) = uu1;    u2 = uu(:,2);
%     u2 = uh(:,n+1); 
    uu(:,1:end-1) = uu(:,2:end);   uu(:,end) = uu1;
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction_pde(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
for k = 1:nt+1
    ue = uv(x2,t(k)); 
    eeu(k) = max(abs(ue-uh(:,k)));
end
eeu;    
%END

function [uh,eeu] = fode_fast_conv_p2_gauss_correction_system_11_2......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB,NN0,TT0)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;           
E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TT0 = 0.5;      
[~,id] = min(abs(t-TT0));       TT0 = t(id);
dTT0 = tau;     
[~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);   dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);     
n1 = n0+dn0;   no = n0;  
% TT0,dTT0
% [aa,bb,cc,dd] = coefficients_p2(alf,n1);
% aa = aa(dn0-1:-1:1); bb = bb(dn0-1:-1:1); cc = cc(dn0-1:-1:1);
% [aa2,bb2,cc2,dd2,ww2] = coefficients_p2_2(-alf,n1,sgm);
% aa2 = aa2(dn0-1:-1:1); bb2 = bb2(dn0-1:-1:1); cc2 = cc2(dn0-1:-1:1);
% dd2 = dd2(:,dn0);
%--------------------------------------------
no1 = length(sgm);   
% use one correction term with smaller stepsize to get the starting values
% if isempty(uv)
    if isempty(sgm)
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[]);
        [uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,[]);
        uh(1:n1+1) = uh00(1:n1+1);
        uhm0 = 0;
    else
%         [uh00] = system_fode_nonlinear_newton_system_1.......
%             (uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
        [uh00] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
        uh(1:n1+1) = uh00(1:n1+1);
        uhm0 = uh(2:no1+1) - uh(1);
    end

if isempty(sgm)
    uhm0 = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
% NN0 = 64;
[aw,alambda,almd1,almd2,almd3,almd4,aa2,bb2,cc2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
aa2 = tau^alf*aa2;  bb2 = tau^alf*bb2; cc2 = tau^alf*cc2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,dTT0,NN0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);     acw2{kk} = acw1{kk};        
end
T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,10*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
        end
        acw2{kk} = wwww3;     acw1{kk} = zeros(nt,lenN);  
        acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:2*BB-1
    u2 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:n0   
    u2 = uh(n+1);    
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
arrL = cc2 - taualf*mu;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = n0+1:nt-dn0+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 = aa2*u0 + bb2*u1 - n^(-alf)/gamma(1-alf)*uh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(n+dn0));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(n+dn0));
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0) = uu1;   u2 = uh(n+1); 
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  
eeu = abs(ue-uh);   
eeu;
%END

function [uh,eeu] = fode_fast_conv_p2_gauss_correction_system_111......
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,BB)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf = alf0(1);  taualf = tau^alf;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;    uh(1) = uh0;     E = taualf;
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
TT0 = 1;      [~,id] = min(abs(t-TT0));       TT0 = t(id);
dTT0 = 1;     [~,id] = min(abs(t-TT0-dTT0));  
dTT0 = t(id);   dTT0 = dTT0 - TT0;
n0 = round(TT0/tau);  dn0 = round(dTT0/tau);  n1 = dn0 + n0;
[aa,bb,cc,dd] = coefficients_p2(alf,n1);
aa = aa(dn0-1:-1:1); bb = bb(dn0-1:-1:1); cc = cc(dn0-1:-1:1);
%--------------------------------------------
no1 = length(sgm);  no = max(no1,2);
% use the known method to derive first n1+1 values
if isempty(uv)
    uh00 = system_fode_nonlinear_newton_system_1.....
        (uv,v,fv,ut0,tau,mu,n1*tau,alf0,Hv,sgm);
    uh = uh00;   uhm0 = uh(2:no1+1) - uh(1);
else
    uh(1:n1+1) = uv(t(1:n1+1));    uhm0 = uh(2:no1+1) - uh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isempty(uv) 
    uv = @(t) 0.*t;  
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ------ initial parameters of fast convolution ------ % %
% BB = 5;
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
NN0 = 64;
[aw,alambda,almd1,almd2,almd3,almd4,aphi0,aphi1,aphi2,nT0]....
    = ini_contour_p2_gauss_b(tau,-alf,BB,nt,dTT0,NN0);
% aphi0 = tau^alf*aphi0;  aphi1 = tau^alf*aphi1; aphi2 = tau^alf*aphi2;
[~,~,almd11,almd22,almd33,almd44] = ini_contour_p2_gauss_2_b(tau,-alf,BB,nt,dTT0,NN0);
% [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT00] = ini_contour_p2_gauss_c(tau,-alf,BB,nt,NN0,dTT0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain starting weights
lenN = length(alambda{1});
for kk = 1:len
    acw1{kk} = zeros(nt,lenN);     acw2{kk} = acw1{kk};        
end

T0 = 1;   N0 = round(T0/tau)-1;
if N0 < no
    N0 = no;  
end
T0 = N0*tau;
for k = 1:no1
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    x0u{k} = x0; w0u{k} = w0;
end
[xx0,ww0] = RootsJacobiLobatto(0,0,32); 
lenN = length(alambda{1});
if ~isempty(sgm)
    for kk = 1:len
        for n = 1:lenN
            tmp = correction_weights_11(tau,T0,alambda{kk}(n),almd2{kk}(n),almd3{kk}(n),.....
                almd4{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww1(:,n) = uhm0*tmp;
            tmp = correction_weights_22(tau,11*tau,alambda{kk}(n),almd22{kk}(n),almd33{kk}(n),.....
                almd44{kk}(n),x0u,w0u,sgm,xx0,ww0);
            wwww3(:,n) = uhm0*tmp;
        end
        acw2{kk} = zeros(10,lenN); acw1{kk} = zeros(nt,lenN);
        acw2{kk} = wwww3(1:10,:);  acw1{kk}(1:N0-1,:) = wwww1(1:N0-1,:);  
    end
end
clear wwww1 wwww3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
d3 = tauL1(BB,2*BB);
u0 = uh(1);    u1 = uh(2);  
% --------------- the fist several steps  --------------
% u(k) and v(k) are known, k = 1,2,...,no+1
n = 2;    u2 = uh(n+1);    
uy2{1} = almd11{1}.*uy2{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2;
uy2{1} = uy2{1} + acw2{1}(1,:);   
for kk = 1:len
    uyy4{kk} = almd11{kk}.*uyy4{kk} + almd22{kk}*u0 + almd33{kk}*u1 + almd44{kk}*u2;
    uyy4{kk} = uyy4{kk} + acw2{kk}(1,:);   
end
if n > d3(2) && n < d3(1)+1
    uy4{1} = almd11{1}.*uy4{1} + almd22{1}*u0 + almd33{1}*u1 + almd44{1}*u2 + acw2{1}(1,:);
    uy44{1} = uy4{1};
end
if n > d3(3) && n < d3(2)+1
    uy4{2} = almd11{2}.*uy4{2} + almd22{2}*u0 + almd33{2}*u1 + almd44{2}*u2 + acw2{2}(1,:);
end
u0 = u1; u1 = u2;  
clear almd11 almd22 almd33 almd44 acw2;  
%
for n = 3:2*BB-1
    u2 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1 + almd4{kk}*u2 + acw1{kk}(n-1,:);
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1 + almd4{1}*u2 + acw1{1}(n-1,:);
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1 + almd4{2}*u2 + acw1{2}(n-1,:);
    end
    u0 = u1;  u1 = u2;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
for n = 2*BB:n0      
    u2 = uh(n+1);
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrL = cc(end) + dd(3) - taualf*mu;
cc = cc(1:end-1);   dd = dd(1:2);  uu = uh(n0+2:n1+1); 
for n = n0+1:nt-dn0      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
    [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
    
    RHS0 = sum(dd.*uu(end-1:end)) + sum(aa.*uu(1:end-1))....
        + sum(bb.*uu(2:end)) + sum(cc.*uu(3:end))....
        + uu(1)*(dn0)^(-alf)/gamma(1-alf) - (dn0+n-1)^(-alf)/gamma(1-alf)*uh0;
%     RHS0 = aphi0*u0 + aphi1*u1 - n^(-alf)/gamma(1-alf)*uh0;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    uu1 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(uu1,t(dn0+n+1));
        rhs = -(arrL*uu1 + rhs0);
        J = arrL - E*Hv(uu1,t(dn0+n+1));
%         if cond(J)>1e10
%             n,nt,uu1,J,return;
%         end
        uu2 = uu1 + J\rhs;
        if norm(uu1-uu2,inf) < 1e-15
            uu1 = uu2; break;
        end
        uu1 = uu2;
    end
    uh(n+dn0+1) = uu1;   uu(1:end-1) = uu(2:end);   uu(end) = uu1;
    u2 = uh(n+1);
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = update_p2_correction(len,uy3,uy33,uy4,....
       uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,almd4,u0,u1,u2,acw1);
    u0 = u1; u1 = u2; 
end
%--------------------------------------
ue = uv(t(1:nt+1));  eeu = abs(ue-uh);   
eeu;
%END
%END

function [uh,eeu] = ini_system_fpde_nonlinear_newton_system_1.....
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,uav,ubv,xa,xb,N)
% D_{0,t}^{alpha1}u = mu*u_xx + f(u,t)
% One correction term is applied, period boundary conditions
dx = (xb-xa)/N;   x = xa:dx:xb+1e-12;  
x = x'; x2 = x(1:end-1);
S = zeros(N,N+1);
for k = 1:N-1
    S(k+1,k:k+2) = [1 -2 1];
end
S(1,1:2) = [-2,1];  S(1,N) = 1; 
S(N,1) = 1;         S(:,end) = []; 
S = sparse(S);
EN1 = sparse(eye(N));
S = S/dx^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(x,t) 0.*x; 
end
ne = round(T/tau)+1;  
alf = alf0(1); 
if isempty(sgm)
    [w0,B0,~,C] = my_weght_2(ne+1,alf,[]);
else
    [w0,B0,~,C] = my_weght_2(ne+1,alf,sgm(1));
end
w = w0;     B = B0;  
taualf = tau^alf;  
arrL = w(1)*EN1 - mu*taualf*S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(N,nt+1);  uh0 = ut0(x2);             
uh(:,1) = uh0;         E = taualf;
%----------------  first step ---------------------
no = 1;  n = 1;    
RHS0 = -B(n+1)*uh0 + uh(:,1:n)*w(n+1:-1:2) - uh0*C(n+1);
u2 = uh(:,n);  
for kk = 1:400
    rhs0 = RHS0 - taualf*(fv(u2,x2,t(n+1)));
    rhs = -((arrL + C(n+1)*EN1)*u2 + rhs0);
    J = arrL - E*sparse(diag(Hv(u2,x2,t(n+1)))) + C(n+1)*EN1;
    uu2 = u2 + J\rhs;
    if norm(u2-uu2,inf)<1e-15
        u2 = uu2; break;
    end
    u2 = uu2;
end
uh(:,n+1) = u2; 
uh00 = uh(:,2:no+1) - uh0;    
e = uv(x2,t(2))-u2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt     
    RHS0 = -B(n+1)*uh0  + uh(:,1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    u2 = uh(:,n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,x2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*sparse(diag(Hv(u2,x2,t(n+1))));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(:,n+1) = u2;  
end
%-----------------------------------------------------
% for k = 1:nt+1
%     ue = uv(x2,t(k)); 
%     eeu(k) = max(abs(ue-uh(:,k)));
% end
eeu = abs(u2 - uv(x2,t(nt+1)));  
eeu;     
%END



function [uh,eeu] .....
    = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;  
alf = alf0(1); 
if isempty(sgm)
    [w0,B0,~,C] = my_weght_2(ne+1,alf,[]);
else
    [w0,B0,~,C] = my_weght_2(ne+1,alf,sgm(1));
end
w = w0;     B = B0;  
taualf = tau^alf;  
arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    
uh0 = ut0;             
uh(1) = uh0;           
E = taualf;
%----------------  first step ---------------------
no = 1;  n = 1;     
RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  - uh0*C(n+1);
u2 = uh(n);
for kk = 1:400
    rhs0 = RHS0  - taualf*fv(u2,t(n+1));
    rhs = -((arrL+C(n+1))*u2 + rhs0);
    J = arrL - E*Hv(u2,t(n+1)) + C(n+1);
    uu2 = u2 + J\rhs;
    if norm(u2-uu2,inf)<1e-15
        u2 = uu2; break;
    end
    u2 = uu2;
end
uh(n+1) = u2; 
uh00 = uh(2:no+1) - uh0;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = -B(n+1)*uh0  + uh(1:n)*w(n+1:-1:2)  + uh00*C(:,n+1);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);     
%END

function [uh,eeu] = L1_method_corection_1(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;   alf = alf0(1); 
if isempty(sgm)
    [w,C] = coefficients_p1(-alf,ne,sgm);
else
    [w,C] = coefficients_p1(-alf,ne,sgm(1));
end
w = w';
taualf = tau^alf;   arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    uh0 = ut0;             
uh(1) = uh0;           E = taualf;
%----------------  first step ---------------------
no = 1;  n = 1;     
RHS0 = - uh0*C(n) - w(1)*uh(n);
u2 = uh(n);
for kk = 1:400
    rhs0 = RHS0 - taualf*fv(u2,t(n+1));
    rhs = -((arrL+C(n))*u2 + rhs0);
    J = arrL - E*Hv(u2,t(n+1)) + C(n);
    uu2 = u2 + J\rhs;
    if norm(u2-uu2,inf)<1e-15
        u2 = uu2; break;
    end
    u2 = uu2;
end
uh(n+1) = u2; 
uh00 = uh(2:no+1) - uh0;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    RHS0 = (uh(2:n)-uh(1:n-1))*w(n:-1:2) + uh00*C(:,n) - w(1)*uh(n);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);     
%END

function [uh,eeu] = L1_method_quadratic(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;   alf = alf0(1); 
if isempty(sgm)
    [w,C] = coefficients_yan(-alf,ne,sgm);
else
    [w,C] = coefficients_yan(-alf,ne,sgm(1));
end
w = w';
taualf = tau^alf;   arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    uh0 = ut0;             
uh(1) = uh0;           E = taualf;
%----------------  first step ---------------------  
uh(1:3) = uv(t(1:3));    uh0 = uh(1)/gamma(1-alf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = 3:nt      
    RHS0 = uh(1:n)*w(n:-1:1) - n^(-alf)*uh0;
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);     
%END

function [uh,eeu,cputime] = L1_method_corection_m(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cputime = zeros(1,2);
tic
ne = round(T/tau)+1; alf = alf0(1);   M = round(1/tau);
no = length(sgm);
if no < 2
    [uh,eeu] = L1_method_corection_1(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm); 
    cputime(1) = toc;    return;
else
    vh0 = L1_method_corection_1(uv,fv,ut0,tau/M,mu,no*tau,alf0,Hv,sgm(1));
end
if isempty(uv)
    uv = @(t) 0.*t; 
end
[w,C] = coefficients_p1(-alf,ne,sgm);
w = w';
taualf = tau^alf;   arrL = w(1) - taualf*mu; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    uh0 = ut0;             
uh(1:no+1) = vh0(1:M:end);           
E = taualf;
%----------------  first step ---------------------
uh00 = uh(2:no+1) - uh0;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
cputime(1) = toc;
tic
for n = no+1:nt      
    RHS0 = (uh(2:n)-uh(1:n-1))*w(n:-1:2) + uh00*C(:,n) - w(1)*uh(n);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);  
cputime(2) = toc;
%END

function [uh,eeu] = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
alf = alf0(1);   nt0 = 4;
if isempty(sgm)
%     [w0,B0,~,C] = my_weght_2(ne,alf,[]);
    [aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,[]);
    uh5 = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau/nt0,mu,2*tau,alf0,Hv,[]);
else
%     [w0,B0,~,C] = my_weght_2(ne,alf,sgm(1));
    [aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,sgm(1));
    uh5 = ini_system_fode_nonlinear_newton_system_1(uv,fv,ut0,tau/nt0,mu,2*tau,alf0,Hv,sgm(1));
end
  
taualf = tau^alf;  
arrL = cc(1) - taualf*mu; 
%----------------  first step ---------------------  
uh = zeros(1,nt+1);           
E = taualf;  uh0 = ut0;      
uh(1:3) = uh5(1:nt0:end);
uh10 = uh(2) - uh0;       uh20 = uh(3) - uh0;  
vh = uh - uh(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa'; bb = bb'; cc = cc';
for n = 3:nt      
    RHS0 = -cc(1)*uh0 + uh10*dd(2,n) + uh20*dd(3,n) + ww(n)*uh10;
    RHS0 = RHS0 + vh(1:n-1)*aa(n-1:-1:1) + vh(2:n)*bb(n-1:-1:1)...
        + vh(3:n)*cc(n-1:-1:2);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2 - uh0;
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);     
%END

function [uh,eeu] = ini_fpde_nonlinear_quadratic_p2....
    (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,uav,ubv,xa,xb,N)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = (xb-xa)/N;   x = xa:dx:xb+1e-12;  
x = x'; x2 = x(1:end-1);
S = zeros(N,N+1);
for k = 1:N-1
    S(k+1,k:k+2) = [1 -2 1];
end
S(1,1:2) = [-2,1];  S(1,N) = 1; 
S(N,1) = 1;         S(:,end) = []; 
S = sparse(S);
EN1 = sparse(eye(N));
S = S/dx^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(x,t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
alf = alf0(1);   nt0 = 4;
if isempty(sgm)
%     [w0,B0,~,C] = my_weght_2(ne,alf,[]);
    [aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,[]);
    uh5 = ini_system_fpde_nonlinear_newton_system_1(uv,fv,ut0,tau/nt0,mu,2*tau,alf0,Hv,[],uav,ubv,xa,xb,N);
else
%     [w0,B0,~,C] = my_weght_2(ne,alf,sgm(1));
    [aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,sgm(1));
    uh5 = ini_system_fpde_nonlinear_newton_system_1(uv,fv,ut0,tau/nt0,mu,2*tau,alf0,Hv,sgm(1),uav,ubv,xa,xb,N);
end
  
taualf = tau^alf;  
arrL = cc(1)*EN1 - taualf*mu*S; 
%----------------  first step ---------------------  
uh = zeros(N,nt+1);           
E = taualf;  uh0 = ut0(x2);      
uh(:,1:3) = uh5(:,1:nt0:end);
uh10 = uh(:,2) - uh0;       uh20 = uh(:,3) - uh0; 
vh = uh - uh(:,1)*ones(1,nt+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa'; bb = bb'; cc = cc';
for n = 3:nt      
    RHS0 = -cc(1)*uh0 + uh10*dd(2,n) + uh20*dd(3,n) + ww(n)*uh10;
    RHS0 = RHS0 + vh(:,1:n-1)*aa(n-1:-1:1) + vh(:,2:n)*bb(n-1:-1:1)...
        + vh(:,3:n)*cc(n-1:-1:2);
    u2 = uh(:,n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,x2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*sparse(diag(Hv(u2,x2,t(n+1))));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(:,n+1) = u2;  vh(:,n+1) = u2 - uh0;
end
%-----------------------------------------------------
% for k = 1:nt+1
%     ue = uv(x2,t(k)); 
%     eeu(k) = max(abs(ue-uh(:,k)));
% end
eeu = abs(u2 - uv(x2,t(nt+1)));
eeu;   
%END

function [uh,eeu] = fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
alf = alf0(1);   nt0 = 1/tau;
if isempty(sgm)
    [uh,eeu] = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm);
    return;
end
no1 = length(sgm);
no = max(no1,2);

uh5 = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau/nt0,mu,no*tau,alf0,Hv,sgm(1));
[aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,sgm);
taualf = tau^alf;  
arrL = cc(1) - taualf*mu; 
%----------------  first step ---------------------  
uh = zeros(1,nt+1);           
E = taualf;  uh0 = ut0;      
uh(1:no+1) = uh5(1:nt0:end);
% % % % % % % % % % % % % % % % % % 
uh(1:no+1) = uv(t(1:no+1));
ww = ww*0;
% % % % % % % % % % % % % % % % % % 
% uh(1:no+1) = uv(t(1:no+1));  %%%%%%%%%%%%%%%
uh10 = uh(2) - uh0;       uh20 = uh(3) - uh0;  
uh00 = uh(2:no1+1) - uh0;
vh = uh - uh(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa'; bb = bb'; cc = cc';
for n = no+1:nt      
    RHS0 = -cc(1)*uh0 + uh10*dd(2,n) + uh20*dd(3,n) + uh00*ww(:,n);
    RHS0 = RHS0 + vh(1:n-1)*aa(n-1:-1:1) + vh(2:n)*bb(n-1:-1:1)...
        + vh(3:n)*cc(n-1:-1:2);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);  vh(n+1) = u2 - uh0;
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);        
%END


function [uh,eeu] = fode_nonlinear_quadratic_p2_sisc(uv,fv,ut0,tau,mu,T,alf,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
uh = zeros(1,nt+1);           
uh0 = ut0;      
no1 = length(sgm);  no = max(no1,1);
if isempty(sgm)
    nt0 = round(1/tau);
    uh5 = L1_method_corection_m(uv,fv,ut0,tau/nt0,mu,tau,alf,Hv,sgm);
    uh00 = 0;
    uh(1:no+1) = uh5(1:nt0:end);
    no1 = no1+1;
else
    nt0 = 1;
    uh5 = L1_method_corection_m(uv,fv,ut0,tau/nt0,mu,no*tau,alf,Hv,sgm);
    uh(1:no+1) = uh5(1:nt0:end);
    uh00 = uh(2:no+1) - uh0;
end
[aa,bb,cc,dd,ww] = LvXu_sisc_2016(-alf,ne,sgm);
taualf = tau^alf;    E = taualf;  
arrL = cc(1) + dd(3) - taualf*mu; 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% ww = ww*0;
uh(1:no+1) = uv(t(1:no+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa'; bb = bb'; cc = cc';
for n = no1+1:nt      
    RHS0 = dd(1)*uh(n-1) + dd(2)*uh(n) + uh00*ww(:,n-1) - n^(-alf)/gamma(1-alf)*uh0;
    RHS0 = RHS0 + uh(1:n-1)*aa(n-1:-1:1) + uh(2:n)*bb(n-1:-1:1)...
        + uh(3:n)*cc(n-1:-1:2);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-14
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1); 
end
%------------------------
ue = uv(t(1:nt+1));   eeu = abs(ue-uh);        
%END

function [uh,eeu] = fode_nonlinear_linear_interpolation_p1(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
alf = alf0(1);   nt0 = 1/tau;
if isempty(sgm)
    [uh,eeu] = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm);
    return;
end
no1 = length(sgm);
no = max(no1,2);

uh5 = ini_fode_nonlinear_quadratic_p2(uv,fv,ut0,tau/nt0,mu,no*tau,alf0,Hv,sgm(1));
% [aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,sgm);
[aa,ww] = coefficients_p1(-alf,ne,sgm);
taualf = tau^alf;  
arrL = aa(1) - taualf*mu; 
%----------------  first step ---------------------  
uh = zeros(1,nt+1);           
E = taualf;  uh0 = ut0;      
uh(1:no+1) = uh5(1:nt0:end);
% uh(1:no+1) = uv(t(1:no+1));  %%%%%%%%%%%%%%% 
uh00 = uh(2:no1+1) - uh0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa';  
for n = 2:nt      
    RHS0 = (uh(2:n-1)-uh(1:n-2))*aa(n-1:-1:2) + uh00*ww(:,n) - aa(1)*uh(n-1);
    u2 = uh(n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*Hv(u2,t(n+1));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(n+1) = u2(1);   
end
%-----------------------------------------------------
ue = uv(t(1:nt+1));   
eeu = abs(ue-uh);        
%END

function [uh,eeu] = fpde_nonlinear_quadratic_p2(uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,uav,ubv,xa,xb,N)
% D_{0,t}^{alpha1}u = mu*u + f(u,t)
% One correction term is applied, Newton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = (xb-xa)/N;   x = xa:dx:xb+1e-12;  
x = x'; x2 = x(1:end-1);
S = zeros(N,N+1);
for k = 1:N-1
    S(k+1,k:k+2) = [1 -2 1];
end
S(1,1:2) = [-2,1];  S(1,N) = 1; 
S(N,1) = 1;         S(:,end) = []; 
S = sparse(S);
EN1 = sparse(eye(N));
S = S/dx^2;
% % % % % % % % % % % % % % % % % % % % % 
if isempty(uv)
    uv = @(x,t) 0.*t; 
end
ne = round(T/tau)+1;
t = 0:tau:T+1e-10;     nt = length(t)-1;    
alf = alf0(1);   nt0 = 1/tau;
if isempty(sgm)
    [uh,eeu] = ini_fpde_nonlinear_quadratic_p2....
        (uv,fv,ut0,tau,mu,T,alf0,Hv,sgm,uav,ubv,xa,xb,N);
    return;
end
no1 = length(sgm);
no = max(no1,2);

uh5 = ini_fpde_nonlinear_quadratic_p2...
    (uv,fv,ut0,tau/nt0,mu,no*tau,alf0,Hv,sgm(1),uav,ubv,xa,xb,N);
[aa,bb,cc,dd,ww] = coefficients_p2_2(-alf,ne,sgm);
taualf = tau^alf;  
arrL = cc(1)*EN1 - taualf*mu*S; 
% arrL = cc(1)*EN1 - mu*taualf*S(:,2:end-1);
%----------------  first step ---------------------  
uh = zeros(N,nt+1);           
E = taualf;  uh0 = ut0(x2);      
uh(:,1:no+1) = uh5(:,1:nt0:end);
% uh(1:no+1) = uv(t(1:no+1));  %%%%%%%%%%%%%%%
uh10 = uh(:,2) - uh0;       uh20 = uh(:,3) - uh0;  
uh00 = uh(:,2:no1+1) - uh0*ones(1,no1);
vh = uh - uh(:,1)*ones(1,nt+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
aa = aa'; bb = bb'; cc = cc';
for n = no+1:nt     
    RHS0 = -cc(1)*uh0 + uh10*dd(2,n) + uh20*dd(3,n) + uh00*ww(:,n);
    RHS0 = RHS0 + vh(:,1:n-1)*aa(n-1:-1:1) + vh(:,2:n)*bb(n-1:-1:1)...
        + vh(:,3:n)*cc(n-1:-1:2);
    u2 = uh(:,n);
    for kk = 1:400
        rhs0 = RHS0 - taualf*fv(u2,x2,t(n+1));
        rhs = -(arrL*u2 + rhs0);
        J = arrL - E*sparse(diag(Hv(u2,t(n+1))));
        uu2 = u2 + J\rhs;
        if norm(u2-uu2,inf)<1e-15
            u2 = uu2; break;
        end
        u2 = uu2;
    end
    uh(:,n+1) = u2;  vh(:,n+1) = u2 - uh0;
end
%-----------------------------------------------------
% for k = 1:nt+1
%     ue = uv(x2,t(k)); 
%     eeu(k) = max(abs(ue-uh(:,k)));
% end
eeu = abs(u2 - uv(x2,t(nt+1)));
eeu; 
%END

function [uh,e,e2,ee] = fdm_multi_term_fode_2(uv,fv,ut0,tau,mu,T,alf0,nu,sgm,d0,D,d)
% D_{0,t}^{alf}u + nv*D_{0,t}^{beta}u = mu*u + f
% the method is exact for u = t^smg, alf=<smg<2 is not an integer. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
% [w,B,C,w1,w2] = weght_2(ne,alf,sgm);
alf = alf0(1); alf1 = alf0(2);
[w0,B0,~,C0] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgm);
w = nu*tau^(1-alf)*w0;  B = nu*tau^(1-alf)*B0;  C = C0*nu*tau^(1-alf);
w1 = tau^(1-alf1)*w10;  B1 = tau^(1-alf1)*B10;  C1 = C10*tau^(1-alf1);
w = w+w1; B = B+B1; C = C+C1;
% C(:,round(ne/10):end) = 0;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = w(1) - mu*tau;
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
uhin0 = fdm_multi_term_fode_adapt(uv,fv,ut0,tau,mu,d0,alf0,nu,sgm,D,d);
p = length(uhin0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh0 = ut0;  
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);  
no = length(sgm);
% for k = 1:no
%     uh(k+1) = uuh0(k) + uh0; 
% end
uh(1) = uh0;  rhs1 = 0;
uh(1:p) = uhin0;    uuh0 = 0;
if no>0
    uuh0 = uhin0(2:no+1) - uh0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = p:nt      
    rhs0 = B(n)*uh0 - uuh0*C(:,n+1);
    rhs = rhs0  + tau*fv(t(n+1));
    if n>1
        rhs1 =  uh(2:n)*w(n:-1:2); 
    end
    rhs = rhs - rhs1;
    u2 = rhs/arrL; 
    uh(n+1) = u2;    
end
%----------------------------------------------------- 
ue = uh;
for k = 1:nt+1
    ue(k) = uv(t(k));
end
ee = abs(ue-uh);
e2 = norm(ue-uh,inf);
e = abs(ue(nt)-uh(nt));
%END

function [uh,e,e2] = fdm_multi_term_fode_3(uv,fv,ut0,tau,mu,T,alf0,nu,sgm,T0,m0)
% D_{0,t}^{alf}u + nv*D_{0,t}^{beta}u = mu*u + f
% the method is exact for u = t^smg, alf=<smg<2 is not an integer. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(uv)
    uv = @(t) 0.*t;
end
ne = round(T/tau)+1;  
% [w,B,C,w1,w2] = weght_2(ne,alf,sgm);
alf = alf0(1); alf1 = alf0(2);
[w0,B0,~,C0] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgm);
w = nu*tau^(1-alf)*w0;  B = nu*tau^(1-alf)*B0;  C = C0*nu*tau^(1-alf);
w1 = tau^(1-alf1)*w10;  B1 = tau^(1-alf1)*B10;  C1 = C10*tau^(1-alf1);
w = w+w1; B = B+B1; C = C+C1;
% C(:,round(ne/10):end) = 0;
%%%%%%%%%%%%%%%% end of check f2phi %%%%%%%%%%%%%%%%%
arrL = w(1) - mu*tau;
% uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
% uhin0 = fdm_multi_term_fode_adapt(uv,fv,ut0,tau,mu,m0,alf0,nu,sgm,m);
tau0 = tau/(2^m0);
uhin0 = fdm_multi_term_fode(uv,fv,ut0,tau0,mu,T0,alf0,nu,sgm);
uhin0 = uhin0(1:2^m0:end);
p = length(uhin0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh0 = ut0;  
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);  
no = length(sgm);
% for k = 1:no
%     uh(k+1) = uuh0(k) + uh0; 
% end
uh(1) = uh0;  rhs1 = 0;
uh(1:p) = uhin0;    uuh0 = 0;
if no>0
    uuh0 = uhin0(2:no+1) - uh0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = p:nt      
    rhs0 = B(n)*uh0 - uuh0*C(:,n+1);
    rhs = rhs0  + tau*fv(t(n+1));
    if n>1
        rhs1 =  uh(2:n)*w(n:-1:2); 
    end
    rhs = rhs - rhs1;
    u2 = rhs/arrL; 
    uh(n+1) = u2;    
end
%----------------------------------------------------- 
ue = uh;
for k = 1:nt+1
    ue(k) = uv(t(k));
end
e2 = norm(ue-uh,inf);
e = abs(ue(nt)-uh(nt));
%END

function [uh,e,e2] = fdm_multi_term_fode_adapt(uv,fv,ut0,tau,mu,d0,alf0,nu,sgm,D,d)
% D_{0,t}^{alf}u + nv*D_{0,t}^{beta}u = mu*u + f
% the method is exact for u = t^smg, alf=<smg<2 is not an integer. 
% graded mesh is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% m0 = 5;    
if isempty(uv)
    uv = @(t) 0.*t;
end
h = 2^d0*tau;
mtau = zeros(1,d);  mt = cell(1,d);  
for k = 1:d
    mtau(k) = h/(2^(D-k+1));
    mt{k} = 0:mtau(k):k*h+1e-12;
end
nT = (1:d)*h;   k=1:d;   
SN = k.*(2.^(D-k+1));     
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
tau = mtau(1);  T = nT(1);
ne = round(T/tau)+1;  
% [w,B,C,w1,w2] = weght_2(ne,alf,sgm);
alf = alf0(1); alf1 = alf0(2);
[w0,B0,~,C0] = my_weght_2(ne,alf,sgm);
[w10,B10,~,C10] = my_weght_2(ne,alf1,sgm);
w = nu*tau^(1-alf)*w0;  B = nu*tau^(1-alf)*B0;  C = C0*nu*tau^(1-alf);
w1 = tau^(1-alf1)*w10;  B1 = tau^(1-alf1)*B10;  C1 = C10*tau^(1-alf1);
w = w+w1; B = B+B1; C = C+C1;
% C(:,round(ne/10):end) = 0;
arrL = w(1) - mu*tau;
uuh0 = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh0 = ut0;
t = 0:tau:T+1e-10;     nt = length(t)-1;      
% ua = uav(t(2:end));   ub = ubv(t(2:end));
uh = zeros(1,nt+1);  
no = length(sgm);
for k = 1:no
    uh(k+1) = uuh0(k) + uh0; 
end
uh(1) = uh0;  rhs1 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = no+1:nt      
    rhs0 = B(n)*uh0 - uuh0*C(:,n+1);
    rhs = rhs0  + tau*fv(t(n+1));
    if n>1
        rhs1 =  uh(2:n)*w(n:-1:2); 
    end
    rhs = rhs - rhs1;
    u2 = rhs/arrL; 
    uh(n+1) = u2;    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for kk = 2:d
    uh(2:2:end) = []; 
%     uuu = uh(1:2:kk*2^(D-kk+1)+1); uh(1:kk*2^(D-kk)+1) = uuu;
%     uuu = uh(1:2:end); uh = [uuu,zeros(1,length(uuu)-1)];
    if no>0
        uuh0 = uh(2:no+1) - uh0;
    end
    tau = mtau(kk);         T = nT(kk);
    t = 0:tau:T+1e-10;      nt = length(t)-1;
    w = nu*tau^(1-alf)*w0;  B = nu*tau^(1-alf)*B0;  C = C0*nu*tau^(1-alf);
    w1 = tau^(1-alf1)*w10;  B1 = tau^(1-alf1)*B10;  C1 = C10*tau^(1-alf1);
    w = w+w1; B = B+B1; C = C+C1;
    arrL = w(1) - mu*tau;
    for n = SN(kk-1)/2+1:nt  %%%%%%%%%%%% check
        rhs0 = B(n)*uh0 - uuh0*C(:,n+1);
        rhs = rhs0  + tau*fv(t(n+1));
        rhs1 =  uh(2:n)*w(n:-1:2);
        rhs = rhs - rhs1;
        u2 = rhs/arrL;
        uh(n+1) = u2;
    end
end
%----------------------------------------------------- 
ue = uh;
for k = 1:nt+1
    ue(k) = uv(t(k));
end
e2 = norm(ue-uh,inf);
e = abs(ue(nt+1)-uh(nt+1));
% uh = uh(1:2^(D-d0-d+1):end);
uh = uh(1:2^(D-d0-d+1):d*2^(D-d+1)+1);
%END


function uh = ini_sem_frac_sub_n(fv,ut0,tau,mu,T,alf,alf1,nu,sgm)
% the method is exact for u = t^smg, alf=<smg<2 is not an integer. 
if isempty(sgm)
    uh = 0; return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = round(T/tau)+1;  
[w,~,~,w2] = my_weght_2(ne+1,alf,sgm);
[ww,~,~,ww2] = my_weght_2(ne+1,alf1,sgm);
w12 = nu*tau^(1-alf)*w2 + tau^(1-alf1)*ww2;   
  g = nu*tau^(1-alf)*w  + tau^(1-alf1)*ww;  
m = length(sgm); AA = zeros(m);
for k=1:m
    AA(k,1:k) = g(k:-1:1);
    AA(k,:) = AA(k,:) + w12(:,k+1).';
end   
uh0 = ut0;
% % % % % % % % % % % % % 
arr = zeros(m);  rhs = zeros(m,1);
for i=1:m
    for j=1:m
        arr(i,j) = AA(i,j);
    end
    arr(i,i) = AA(i,i) - mu*tau;
    rhs(i,1) = mu*tau*uh0 + tau*(fv(i*tau));
end
uh = arr\rhs; uh = uh';
%END


function [uv,fv,tau,mu,T,ut0,alpha,sgm,Hv] = solu_system_nonlinear_system_1(sel)
% u_t  = mu u_xx + f
switch sel
    case 0
        mu = 0;
        T = 4;    alf = 0.8;  
        tau = 1./[32 64 128 256 512]';
        nv = [alf 2*alf 3*alf 4 5 6 7 8]+0.1;
        nv = [alf 2*alf 3*alf 4 5 6 7 8]*0;
%         nv = [1 2 3 4 5 6 7 8];
        c = [1 1 1 1 1 1 1 1];
        c = [1 1 1 1 0 0 0 0];
        c = [0 0 0 0 0 0 0 0];
        nv1 = nv(1);   nv2 = nv(2);  nv3 = nv(3);  nv4 = nv(4);
        nv5 = nv(5);   nv6 = nv(6);  nv7 = nv(7);  nv8 = nv(8);
        f = gamma(nv+1)./gamma(nv+1-alf).*c;        
        d = 1;
        uv = @(t) (c(1)*t.^nv1 + c(2)*t.^nv2 + c(3)*t.^nv3 + c(4)*t.^nv4....
            + c(5)*t.^nv5 + c(6)*t.^nv6 + c(7)*t.^nv7 + c(8)*t.^nv8 + d);
        fv = @(u,t) (f(1)*t.^(nv1-alf) + f(2)*t.^(nv2-alf) + f(4)*t.^(nv4-alf) ....
                + f(3)*t.^(nv3-alf) + f(5)*t.^(nv5-alf) + f(6)*t.^(nv6-alf).....
                + f(7)*t.^(nv7-alf) + f(8)*t.^(nv8-alf)).....
                - mu*uv(t);
        Hv = @(u,t) 0;
        ut0 = uv(0);  
        alpha = alf;
        sgm = nv; 
    case 1
        mu = -1;
        T = 40;    alf = 0.1;  
        tau = 1./[32 64 128 256 512]';
        nv = [alf 2*alf 3*alf 4 5 6 7 8]+0.1;
        nv = [1 2 3 4 5 6 7 8];
        c = [1 1 1 1 1 1 1 1];
        c = [1 1 1 1 0 0 0 0];
%         c = [1 1 0 0 0 0 0 0];
        nv1 = nv(1);   nv2 = nv(2);  nv3 = nv(3);  nv4 = nv(4);
        nv5 = nv(5);   nv6 = nv(6);  nv7 = nv(7);  nv8 = nv(8);
        f = gamma(nv+1)./gamma(nv+1-alf).*c;        
        d = 2;
        uv = @(t) (c(1)*t.^nv1 + c(2)*t.^nv2 + c(3)*t.^nv3 + c(4)*t.^nv4....
            + c(5)*t.^nv5 + c(6)*t.^nv6 + c(7)*t.^nv7 + c(8)*t.^nv8 + d);
        fv = @(u,t) (f(1)*t.^(nv1-alf) + f(2)*t.^(nv2-alf) + f(4)*t.^(nv4-alf) ....
                + f(3)*t.^(nv3-alf) + f(5)*t.^(nv5-alf) + f(6)*t.^(nv6-alf).....
                + f(7)*t.^(nv7-alf) + f(8)*t.^(nv8-alf)).....
                - mu*uv(t) - (u.^2 - uv(t).^2);
        Hv = @(u,t) -2*u;
        ut0 = uv(0);  
        alpha = alf;
        sgm = nv; 
    case 2 % non-smooth solutions, linear case   
        mu = -1;
        T = 40;    alf = 0.7;   
        tau = 1./[32 64 128 256 512]';  
%         tau = 1/64;
        uv = @(t) ml(-t.^alf,alf);  
        fv = @(u,t) 0.*t;
        Hv = @(u,t) 0;
        ut0 = uv(0); 
        alpha = alf;
        sgm = (1:10)*alf;  
    case 3 % non-smooth solutions   
        mu = -1;
        T = 5;    alf = 0.5;  
        tau = 1./[32 64 128 256 512]';  
%         tau = T/512;
        uv = @(t) ml(-t.^alf,alf);  
        fv = @(u,t) -(1 + mu)*uv(t) - (u.^2 - uv(t).^2);
        Hv = @(u,t) -2*u;
        ut0 = uv(0); 
        alpha = alf;
        sgm = (1:10)*alf;  
    case 4  % no analytical solutions 
        mu = -1;
        T = 10;    alf = 0.1;   
        tau = T./[32 64 128 256 512]';
        uv = [];  
        fv = @(u,t) -u.*(1 - u.^2);

        Hv = @(u,t) -1 + 3*u.^2;
        ut0 = 1;  
        alpha = alf;
        sgm = alf + (0:10)*0.1; 
    case 5  % no analytical solutions 
        mu = -1;
        T = 10;    alf = 0.1;   
        tau = T./[32 64 128 256 512]';
        uv = [];  
        fv = @(u,t) -u.*(1 - u.^2) + cos(2*pi*t);

        Hv = @(u,t) -1 + 3*u.^2;
        ut0 = 1;  
        alpha = alf;
        sgm = alf + (0:10)*0.1; 
    case 7  % no analytical solutions 
        mu = -4;
        T = 10;    alf = 0.1;   
        tau = T./[32 64 128 256 512]';
        uv = [];  
        fv = @(u,t) - u.^2 + 2*sin(pi*t);

        Hv = @(u,t) -2*u;
        ut0 = 1;  
        alpha = alf;
        sgm = alf + (0:10)*0.1; 
end
%END


function c0 = initial(u0,Nx,SN,Q,aa,bb)
% u0 = ut0(nxpt);  
c0 = zeros(SN(end),1);     
c = utouh(aa,bb,u0{1},2);    c = (Q{1}.')\c;    
c0(SN(1):SN(2)) = c;
for k=2:Nx
    c = utouh(aa,bb,u0{k},2);   c = (Q{k}')\c;    
    c0(SN(k):SN(k+1)) = c;
end
%END xL2L

function Q=xL2L(N)
% (x*L_{0},...,x*L_{N})' = Q * (L_{0},...,L_{N+1})'
Q=sparse(zeros(N+1,N+2));
Q(1,2)=1;
for k=1:N
   Q(k+1,k:k+2)=[k/(2*k+1),0,(k+1)/(2*k+1)]; 
end
%END xL2L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mm,sm,xsm,ms,QLL,Q]=basis00(N)
% base functions:  phi_i = L_i - L_{i+2}, i=0,...,N-2
% Let  phi_{-2} = (L_0 - L_1)/2, phi_{-1} = (L_0 + L_1)/2
% L = [L_0,...,L_N]'
% phi = [phi_0,...,phi_{N-2}]' = Q*L
% chi = [phi_{-2},phi_0,...,phi_{N-2},phi_{-1}] = Q1*L
%>>> [mm]  = [(phi_i,phi_j)] ================ Q*LL*Q'            
%>>> [sm]  = [(diff(phi_i), diff(phi_j)] ==== (Q*D1L)*LL*(Q*D1L)'
%>>> [xsm] = [(x*diff(phi_i), diff(phi_j))] = (Q*D1L)*P*L1L*(Q*D1L)' 
%>>> [ms]  = [(phi_i, diff(phi_j))] ========= Q*LL*(Q*D1L)'
%>>> (phileg) = (phi_i,L_j) ============= Q*LL;
Q=zeros(N+1,N+1);
for k=0:N-2
    Q(k+2,k+1:k+3)=([1,0,-1]);
end
Q(1,1:2)=[1/2,-1/2]; Q(N+1,1:2)=[1/2,1/2]; Q = sparse(Q);
LL=sparse(diag((2./(2*(0:N)+1))));
QLL=Q(2:end-1,:)*LL;
L1L=LL;        L1L(N+2,:)=0;
D1L=DLegendre(N,1);   P=xL2L(N);
mm=Q*LL*Q';       sm=(Q*D1L)*LL*(Q*D1L)';
xsm=(Q*D1L)*P*L1L*(Q*D1L)';    ms=Q*LL*(Q*D1L).';
%END basis00

function [mm,sm,xsm,ms,QLL,Q]=basis0(N)
Nx = length(N); mm=cell(1,Nx); 
sm=mm; xsm=mm; ms=mm; QLL=mm; Q=mm;
for k=1:Nx
    [mm0,sm0,xsm0,ms0,QLL0,Q0]=basis00(N(k));
    mm{k}=mm0;  sm{k}=sm0;    xsm{k}=xsm0;
    ms{k}=ms0;  QLL{k}=QLL0;  Q{k}=Q0;
end
%END

function [mmat,smat,smmat,xsmat]=basis11(N,nodes)
M = length(nodes)-1;
scl = (nodes(2:M+1)-nodes(1:M))/2;
sclb = (nodes(2:M+1)+nodes(1:M))/2;
[mm,sm,xsm,ms] = basis0(N);
SN = ones(1,M+1);
for k=1:M
    SN(k+1) = SN(k)+N(k);
end
mmat = zeros(SN(end));  smat = mmat; 
smmat = mmat; xsmat = mmat;     
for k=1:M
     mmat(SN(k):SN(k+1),SN(k):SN(k+1)) = mm{k}*scl(k);
     smat(SN(k):SN(k+1),SN(k):SN(k+1)) = sm{k}/scl(k);
    xsmat(SN(k):SN(k+1),SN(k):SN(k+1)) = xsm{k} + sclb(k)/scl(k)*sm{k};
    smmat(SN(k):SN(k+1),SN(k):SN(k+1)) = ms{k};
end
for k=1:M-1
    mmat(SN(k+1),SN(k+1)) = mm{k}(end,end)*scl(k)+mm{k+1}(1,1)*scl(k+1);
    smat(SN(k+1),SN(k+1)) = sm{k}(end,end)/scl(k)+sm{k+1}(1,1)/scl(k+1);
    xsmat(SN(k+1),SN(k+1)) = xsm{k}(end,end)+xsm{k+1}(1,1)...
        +sclb(k)/scl(k)*sm{k}(end,end) + sclb(k+1)/scl(k+1)*sm{k+1}(1,1);
    smmat(SN(k+1),SN(k+1)) = ms{k}(end,end)+ms{k+1}(1,1);
end
%END


function l2u=legm(n,x)  % l2u=legm(n, scl*(x-a)-1);
% scl=2/(xb-xa);  x <- scl*(x-a)-1;  a<x<b;     
% u=l2u*uh=(L_0(x),L_1(x),...,L_n(x))*uh 
l2u=[ones(size(x)),x];
for k=2:n
    l2u=[l2u,((2*k-1)*x.*l2u(:,k)-(k-1)*l2u(:,k-1))/k]; 
end
%END legm

function uh=ifftr(u)
% fast cosine transform ,u is real vector
% u(x_j) = \sum_{k=0}^{k=N} uh_{k}*T_{k}(x_j)    
% u(x_j) -----> uh_k
no=size(u,1); n=no-1; 
v=[u;flipud(u(2:n,:))];  %(column vector x= 1 to -1)
% v=[flipud(u);u(2:n)];    %(column vector x= -1 to 1)
%v(1:no)=u, v(no+1:2*n)=fliplr(u(2:n))  % (for row vector)
z=ifft(v); uh=2*real(z(1:no,:));   
uh(1,:)=uh(1,:)/2; uh(end,:)=uh(end,:)/2;
% END ifftr
function y=ifftri(u)
y=ifftr(real(u));
% END  ifftri =========
function D1T=DT2T(N)
% T = [T_0,...,T_N]'
% T_x = [(T_0)_x,...,(T_N)_x]'
% T_x = D1T*T;
if N==0
    D1T=0;
end
if N>0
    c=ones(N+1,1);c(1)=2;
    D1T=sparse(zeros(N+1));
    for n=1:N
        for k=0:n-1
            if mod(k+n,2)>0.1
                D1T(n+1,k+1)=2*n/c(k+1);
            end
        end
    end
end
% D1T=sparse(D1T);
function DL=DLegendre(N,od)
% [diff(L_{0}),...,diff(L_{N})]' = D1L*[L_{0},...,L_{N}]'
% [diff(L_{0},2),...,diff(L_{N},2)]' = D2L*[L_{0},...,L_{N}]'
DL=sparse(zeros(N+1,N+1));
for n=1:N
    for k=0:n-1
        if abs(mod(k+n,2)-1)<1e-10
            DL(n+1,k+1)=2*k+1;
        end
    end
end
DL=DL^od;
% END
%=========== Chebyshev-Legendre transform ======================
function [L2T,T2L]=cheb_lgd_trsfm(N)
% p(x) = \sum f_{k}*T_{k}(x) 
%      = \sum g_{k}*L_{k}(x),k=0,...,N 
% f=[f_{0},...,f_{N}]', g=[g_{0},...,g_{N}]';
% L2T = a_{k,r}=2/(c_{k}*pi)*( T_{k} , L_{r} )_w, c_{0}=2,c_{k}=1,k>0;
% T2L = b_{k,r}=(k+1/2)*( L_{k} , T_{r} )
% f = L2T*g,    g = T2L*f,    A*T2L = I
% (L0,L1,...,L_N) = (T0,T1,...,T_N) * L2T;
% (T0,T1,...,T_N) = (L0,L1,...,L_N) * T2L
L2T=sparse(zeros(N+1,N+1)); 
T2L=L2T;
c=ones(1,N+1);c(1)=2;
% f = L2T*g
L2T(1,1)=1;L2T(2,2)=1;
for k=1:N-1
    L2T(1,k+2)=(2*k+1)/(2*k+2)*L2T(2,k+1)-k/(k+1)*L2T(1,k);
    L2T(k+2,k+2)=(2*k+1)*L2T(k+1,k+1)/(2*k+2); 
    for r=1:k
        L2T(r+1,k+2)=(2*k+1)/(2*k+2)*(L2T(r+2,k+1)*c(r+2)+L2T(r,k+1)*c(r))/c(r+1)-k/(k+1)*L2T(r+1,k);
    end
end
% g = T2L*f
T2L(1,1)=1;T2L(2,2)=1;
for k=1:N-1
     T2L(1,k+2)=2*T2L(2,k+1)/3-T2L(1,k);
     T2L(k+2,k+2)=(2*k+2)*T2L(k+1,k+1)/(2*k+1);  
    for r=1:k
         T2L(r+1,k+2)=(2*r+2)/(2*r+3)*T2L(r+2,k+1)+2*r*T2L(r,k+1)/(2*r-1)-T2L(r+1,k);
    end
end
%================= end of cheb_lgd_trsfm =====================

% ====== LGL-nodes and weights ======
function [varargout]=legnw(n) 
%  The function x=legnw(n) computes n nodes of  the Legendre-Gauss-Lobatto quadrature
%  The function [x,w]= legnw(n) also returns the weights
%   Newton iteration  method is used for computing nodes.  n->n+1  nn->n
  n1=n+1;  k=[1:n];                  % indices                      
  thetak=(4*k-1)*pi/(4*n+2);
  sigmak=(1-(n-1)/(8*n^3)-(39-28./sin(thetak).^2)/(384*n^4)).*cos(thetak);
  ze=(sigmak(1:n-1)+sigmak(2:n))/2;    % Set the intitial approximation
  ep=eps*10;                           % Error threshold
  ze1=ze+ep+1;
 while max(abs(ze1-ze))>=ep,           % Newton'w iteration procedure
      ze1=ze;  [dy,y]=legpl(n,ze);
      ze=ze-(1-ze.*ze).*dy./(2*ze.*dy-n*n1*y); 
 end;                                   % Around 6 iterations are required for n=100.
   varargout{1}=flipud([1,ze,-1]');
 if nargout==1, return; end; 
   varargout{2}=flipud([2/(n*n1),2./(n*n1*y.^2),2/(n*n1)]');
%%END function [varargout]=legnw(n);
% ====== values of J_{a,b} ======
function [varargout]=legpl(n,x) 
    % The function y=legpl(n,alp,bet,x) computes the values of Jacobi polynomial
    %        of  degree n, and parameters (alp,bet) at x.
    % The function [dy,y]=legpl(n,x) also returns the values of 1st-order 
    %        derivatives (in the first output argument dy).
 if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=x; return; end;
     polylst=ones(size(x));   poly=x;
     for  k=2:n,  polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k;
        polylst=poly;   poly=polyn;   end; 
        varargout{1}=poly;
 end;
 if nargout==2,
     if n==0,       varargout{2}=ones(size(x)); 
      varargout{1}=zeros(size(x));  return;end;
     if n==1,  varargout{2}=x;
       varargout{1}=ones(size(x));  return; end;
     polylst=ones(size(x));   pderlst=zeros(size(x));
     poly=x;   pder=ones(size(x));
    for  k=2:n,
      polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k;
      pdern=pderlst+(2*k-1)*poly;
      polylst=poly;   poly=polyn;
      pderlst=pder;  pder=pdern;
    end;
      varargout{2}=poly;  varargout{1}=pder;
 end;
%%END function [varargout]=legpl(n,x)
function [r,w]=RootsJacobiLobatto(alf,bet,N)
[r0,w0] = RootsJacobi(alf+1,bet+1,N-2);
r = [-1;r0;1];
w(1,1)=2^(alf+bet+1)*gamma(bet+2)*gamma(alf+1)/gamma(alf+bet+2);
w(N,1)=2^(alf+bet+1)*gamma(bet+1)*gamma(alf+2)/gamma(alf+bet+2);
w(1,1)=w(1,1)/(N-1);  w(N,1)=w(N,1)/(N-1);
for k=1:N-1
    w(1,1)=w(1,1)*k/(k+1+alf+bet)*(alf+k)/(bet+k);
    w(N,1)=w(N,1)*k/(k+1+alf+bet)*(bet+k)/(alf+k);
end
w(2:N-1,1)=w0./(1+r0)./(1-r0);
%END
function [r,w]=RootsJacobiRadauL(alf,bet,N)
[r0,w0] = RootsJacobi(alf,bet+1,N-1);
r = [-1;r0];
w(1,1)=2^(alf+bet+1)*gamma(bet+1)*gamma(alf+1)/gamma(alf+bet+2);
for k=1:N-1
    w(1,1)=w(1,1)*k/(k+1+alf+bet)*(alf+k)/(bet+k+1);
end
w(2:N,1)=w0./(1+r0);
%END

function [r,w]=RootsJacobiRadauR(alf,bet,N)
[r0,w0] = RootsJacobi(alf+1,bet,N-1);
r = [r0;1];
w(N,1)=2^(alf+bet+1)*gamma(bet+1)*gamma(alf+1)/gamma(alf+bet+2);
for k=1:N-1
    w(N,1)=w(N,1)*k/(k+1+alf+bet)*(bet+k)/(alf+k+1);
end
w(1:N-1,1)=w0./(1-r0);
%END
function [y,w]=RootsJacobi(alf,bet,N)
dt=1/N^2;
a0=(N+alf+bet+1)/2; vps = 1e-15;
xx=-1:dt:1;  xx(end+1)=1;
yy=jacobi(alf,bet,N,xx);
yy1=yy>0; 
str=num2str(yy1);
str(find(isspace(str)))=[] ;
id1=findstr(str, '01');  id2=findstr(str, '10');
xc1=(xx(id1)+xx(id1+1))/2;
xc2=(xx(id2)+xx(id2+1))/2;
xc=[xc1,xc2]; xc=sort(xc).';
%-------------- Newton Method ---------------
for j=1:100
    xd = xc - jacobi(alf,bet,N,xc)./jacobi(alf+1,bet+1,N-1,xc)/a0;
    err=norm(xd-xc,inf);
    if err<vps
        xc=xd;
        break;
    else
        xc=xd;
    end
end
y=sort(xc);
%------------------------------------------
N = N-1;
GN = 2^(alf+bet)*(2*N+alf+bet+2)*gamma(alf+1)*gamma(bet+1)/gamma(alf+bet+2);
for k=1:N
    GN=GN*(alf+k)/(k+1+alf+bet)*(bet+k)/(k+1);
end
t1 = jacobi(alf,bet,N,y);  t2=djacobi(alf,bet,N+1,1,y);
w=GN./(t1.*t2);
%END
function y=jacobi(alf,bet,N,x)
if N==0
    y=ones(length(x),1); return;
end
if N==1
    y=((alf+bet+2)*x + alf-bet)/2;  return;
end
y0 = 1;
y1 = ((alf+bet+2)*x + alf-bet)/2;
n = 1:N;
a = (2*n+alf+bet+1)./(2*n+2).*((2*n+alf+bet+2)./(n+alf+bet+1));
b = (bet^2-alf^2)*(2*n+alf+bet+1)./(2*n+2)./(n+alf+bet+1)./(2*n+alf+bet);
c = (n+alf)./(n+1).*((n+bet)./(n+alf+bet+1)).*((2*n+alf+bet+2)./(2*n+alf+bet));
for n=1:N-1
    y = (a(n)*x-b(n)).*y1 - c(n)*y0;
    y0 = y1; y1 = y;
end
%END

function y = jacobi3(alf,bet,N,x)
if N==0
    y=ones(length(x),1); return;
end
if N==1
    y=((alf+bet+2)*x + alf-bet)/2; y = [ones(length(x),1),y]; return;
end

y0 = 1;  y1 = ((alf+bet+2)*x + alf-bet)/2;
z = ones(length(x),N+1);   z(:,1) = y0;  z(:,2) = y1; 
n = 1:N;
a = (2*n+alf+bet+1)./(2*n+2).*((2*n+alf+bet+2)./(n+alf+bet+1));
b = (bet^2-alf^2)*(2*n+alf+bet+1)./(2*n+2)./(n+alf+bet+1)./(2*n+alf+bet);
c = (n+alf)./(n+1).*((n+bet)./(n+alf+bet+1)).*((2*n+alf+bet+2)./(2*n+alf+bet));
for n=1:N-1
    y = (a(n)*x-b(n)).*y1 - c(n)*y0;
    y0 = y1; y1 = y;  z(:,n+2) = y;
end
y = z;
%END

function dy=djacobi(alf,bet,N,m,x)
dnk=1; dnk0 = 1;
for k=1:m
    dnk=dnk*(N+k+alf+bet)/2;
    dnk0=dnk0*2/(N+k+alf+bet);
end
dnk0 = 1/dnk0;
% dnk0=gamma(N+m+alf+bet+1)/gamma(N+alf+bet+1)/2^m;
dy=0;
if m<N
    dy=jacobi(alf+m,bet+m,N-m,x)*dnk;
end
dy;
%END
function y=normjacobi(alf,bet,n)
y=ones(n+1,1);
y(1) = 2^(alf+bet+1)*gamma(alf+1)*gamma(bet+1)/gamma(alf+bet+2);
y(2) = 2^(alf+bet+1)*gamma(alf+2)*gamma(bet+2)/gamma(alf+bet+2)/(3+alf+bet);
for k=1:n-1
    y(k+2)=(k+alf+1)/(k+1)*((k+bet+1)/(k+alf+bet+1)).....
        *((2*k+alf+bet+1)/(2*k+alf+bet+3)*y(k+1));
end
%END


function [l,x0,w0] = utouh(alf,bet,u,flg)
% u(x_j) ===> c_k : \sum_{k=0}^N c_k J_k(x)
N = length(u); l=zeros(N,1);
dlt = normjacobi(alf,bet,N-1);
switch flg
    case 1
        [x0,w0]=RootsJacobi(alf,bet,N); % Jacobi-Gauss
    case 2
        dlt(end)=(2+(alf+bet+1)/(N-1))*dlt(end);
        [x0,w0]=RootsJacobiLobatto(alf,bet,N);   % Jacobi-Gauss-Lobatto
    case 3
        [x0,w0]=RootsJacobiRadauL(alf,bet,N);   %Jacobi-Gauss-Radau  x0=-1
    case 4
        [x0,w0]=RootsJacobiRadauR(alf,bet,N);  % Jacobi-Gauss-Radau  xN=1
    otherwise
        l=[]; return;
end
for n=0:N-1
    J = jacobi(alf,bet,n,x0);
    l(n+1)=w0'*(u.*J);
end
l=l./dlt;
%END

function cc = lag2jacobi(a,b,N)
dlt = normjacobi(a,b,N);
dlt(end) = (2+(a+b+1)/N)*dlt(end);
[x0,w0] = RootsJacobiLobatto(a,b,N+1);   % Jacobi-Gauss-Lobatto
cc = zeros(N+1);
for j=0:N
    cc(j+1,:) = jacobi(a,b,j,x0).*w0./dlt(j+1);
end
%END
function y = fast_conv(u,v)
n1 = length(u); n2 = length(v);
n = n1+n2;    nm = max(n1,n2);
u1 = zeros(n,1); v1 = u1;
u1(1:n1,1) = u; v1(1:n2,1) = v;
y = ifft(fft(u1).*fft(v1));
y = y(1:nm);
%END

function y = my_conv(u,v)
n1 = length(u); n2 = length(v);  nm = max(n1,n2);
u1 = zeros(nm,1);  v1 = u1;  y = u1;
u1(1:n1,1) = u; v1(1:n2,1) = v;  u1 = u1.';
for n = 1:nm
    y(n) = u1(n:-1:1)*v1(1:n);
end
%END

function w2 = weights(alf,n0,p)
n = n0+1;       
switch p
    case 1
        w = (ones(n,1)); 
        for k=1:n-1
            w(k+1)=(k-1-alf)/k*w(k); 
        end
        w2 = w;
    case 2
        w = (ones(n,1)); 
        u=zeros(n,1);  u(1:3) = [3/2;  -2; 1/2];  
        w2=w; w2(1) = (3/2)^alf;
        ii=0:n-1;  u=sparse(u);
        for j=1:n-1
            w2(j+1) = (alf - ii(1:j)/j*(alf+1))*(w2(1:j).*u(j+1:-1:2))/u(1);
        end
    case 3
        u=zeros(n,1);  u(1:p+1) = [11/6; -3; 3/2;-1/3];  
        w2=ones(n,1);  w2(1) = (u(1))^alf;
        ii=0:n-1;  u=sparse(u);
        for j=1:n-1
            w2(j+1) = (alf - ii(1:j)/j*(alf+1))*(w2(1:j).*u(j+1:-1:2))/u(1);
        end
    case 4
        u=zeros(n,1);  u(1:p+1) = [25/12; -4; 3;-4/3;1/4];  
        w2=ones(n,1);  w2(1) = (u(1))^alf;
        ii=0:n-1;  u=sparse(u);
        for j=1:n-1
            w2(j+1) = (alf - ii(1:j)/j*(alf+1))*(w2(1:j).*u(j+1:-1:2))/u(1);
        end
    case 5
        u=zeros(n,1);  u(1:p+1) = [137/60; -5; 5; -10/3; 5/4; -1/5];  
        w2=ones(n,1);  w2(1) = (u(1))^alf;
        ii=0:n-1;  u=sparse(u);
        for j=1:n-1
            w2(j+1) = (alf - ii(1:j)/j*(alf+1))*(w2(1:j).*u(j+1:-1:2))/u(1);
        end
    case 6
        u=zeros(n,1);  u(1:p+1) = [147/60; -6; 15/2; -20/3; 15/4; -6/5; 1/6];  
        w2=ones(n,1);  w2(1) = (u(1))^alf;
        ii=0:n-1;  u=sparse(u);
        for j=1:n-1
            w2(j+1) = (alf - ii(1:j)/j*(alf+1))*(w2(1:j).*u(j+1:-1:2))/u(1);
        end
end

function output2(fmts,s,fmt,varargin)
fprintf('\n');
row = length(varargin{1}); 
% fmt0 = []; s = [];
col = length(varargin);

for j=1:col
    fprintf(fmts{j},s{j});
end
fprintf('\n');
for i=1:row
    for j=1:col
        fprintf(fmt{j},varargin{j}(i));
    end
    fprintf('\n');
end
fprintf('\n');
% END

function output3(fid,fmts,s,fmt,varargin)
fprintf(fid,'\n');
row = length(varargin{1}); 
col = length(varargin);

for j=1:col
    fprintf(fid,fmts{j},s{j});
end
fprintf(fid,'\n');
for i=1:row
    for j=1:col
        fprintf(fid,fmt{j},varargin{j}(i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% END

function output(s,A)
[row,col] = size(A);
fprintf('%-10s%-10s%-16s%-8s%-16s%-8s',s{1},s{2},s{3},s{4},s{5},s{6});
fprintf('\n');
for k=1:row
    fprintf('%-10s%-10s%-16.4e%-8.4f%-16.4e%-8.4f%',.....
        num2str(A(k,1)),num2str(A(k,2)),A(k,3),A(k,4),A(k,5),A(k,6));
    fprintf('\n');
end
fprintf('\n');
% END

function output1(s,A,fid)
[row,col] = size(A);
fprintf(fid,'%-10s%-10s%-16s%-8s%-16s%-8s',s{1},s{2},s{3},s{4},s{5},s{6});
fprintf(fid,'\n');
for k=1:row
    fprintf(fid,'%-10s%-10s%-16.4e%-8.4f%-16.4e%-8.4f',.....
        num2str(A(k,1)),num2str(A(k,2)),A(k,3),A(k,4),A(k,5),A(k,6));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% END
function odr=myorder(ntau,NN,err)
n=length(ntau); 
tau(1:n,1)=ntau; N(1:n,1)=NN(1:n); 
er(1:n,1)=err(1:n);
if n==1
    odr=0; return;
end
if n>1
    if abs(ntau(1)-ntau(2))>1e-10
        odr=(log(er(1:end-1)./er(2:end)))./log(tau(1:n-1,:)./tau(2:n,:));
    else 
        odr=(log(er(1:end-1)./er(2:end)))./log(N(1:n-1,:)./N(2:n,:));
    end
end
odr=[0;odr];
%END 


function [w,B,w1,w2] = my_weght_2(ne,alf,sgm)
% Weighted GL : exact for u = t^sigma
% w0 = weights(alf,ne+1,1);    w = w0;  
% w(1) = (2+alf)/2*w0(1); w(2:end) = (2+alf)/2*w0(2:end) - alf/2*w0(1:end-1);
w = weights(alf,ne+1,3);
B = w;   
for k = 1:ne
    B(k+1) = B(k) + w(k+1);
end
% no = numel(varargin); 
no = length(sgm);
n = 0:ne;  
arr2 = (zeros(no));  bb = (zeros(no,ne+1));
% x0 = RootsJacobi(0,0,no);   
if ~isempty(sgm)
    for k = 1:no
        s = (sgm(k));    arr2(k,:) = ((1:no).^s);
        t1 = gamma(1+s)/gamma(1+s-alf)*n.^(s-alf);
        if isinf(t1(1)) || isnan(t1(1))
            t1(1) = 0;
        end
        t2 = (my_conv(n.^s,w)');
        bb(k,:) = (t1(1:ne+1) - t2(1:ne+1));
    end
    w2 = arr2\bb;
%     e=max(max(abs(arr2*w2-bb)));
end
if isempty(sgm)
    w2 = zeros(1,ne+1);
end
w1 = (w2(:,1:end-1) + w2(:,2:end))/2;
%END

function myplot0(uh,tau,T,p0,varargin)
tt = 0:tau:T;   d = round(T/tau/p0);
tt = [tt(1:d:end),tt(end)]; 
u1 = [uh(1:d:end),uh(end)]; 

figure;
plot(tt,u1,'r*-');
xlabel('t'); ylabel('u');
% legend('u');
m = length(varargin);
if m >0
    saveas(gcf,varargin{1})
end
%END

function myplot1(uh,vh,tau,T,p0,varargin)
tt = 0:tau:T;   d = round(T/tau/p0);
tt = [tt(1:d:end),tt(end)]; 
u1 = [uh(1:d:end),uh(end)]; 
v1 = [vh(1:d:end),vh(end)];
figure;
plot(tt,u1,'r*-',tt,v1,'md-');
xlabel('t'); ylabel('numerical solution');
legend('X','Y');
m = length(varargin);
if m >0
    saveas(gcf,varargin{1})
end
%END



function myplot11(uh,vh,wh,tau,T,p0,varargin)
tt = 0:tau:T;   d = round(T/tau/p0);
tt = [tt(1:d:end),tt(end)]; 
u1 = [uh(1:d:end),uh(end)]; 
v1 = [vh(1:d:end),vh(end)];
w1 = [wh(1:d:end),wh(end)];
figure;
plot(tt,u1,'r*-',tt,v1,'md-',tt,w1,'k+-');
xlabel('t'); ylabel('Numerical solutions');
legend('X','Y','Z');
m = length(varargin);
if m >0
    saveas(gcf,varargin{1})
end
%END

function myplot2(erru,errv,tau,T,p,varargin)
clr = {'r+-','kd-','mo-','b>-','k*-','c*-','b-'};
p0 = p;
if min(T./tau) < p
    p0 = min(T./tau);
end  
len = length(erru);
figure;
for k = 1:len
    ee = erru{k};   d = round(T/tau(k)/p0); 
    tt = 0:tau(k):T;  
    tt = [tt(2:d:end),tt(end)];
    ee = [ee(2:d:end),ee(end)];
    semilogy(tt,ee,clr{k},'LineWidth',1.12);
    hold on;
end
% legend('\tau = 2^{-5}','\tau = 2^{-6}','\tau = 2^{-7}','\tau = 2^{-8}','\tau = 2^{-9}');
legend('k = 5','k = 6','k = 7','k = 8','k = 9');
xlabel('t'); ylabel('Error-X');
% title('\tau = T/2^k');
m = length(varargin);
if m >0
    saveas(gcf,varargin{1})
end
hold off;


len = length(errv);
figure;
for k = 1:len
    ee = errv{k};    d = round(T/tau(k)/p0); 
    tt = 0:tau(k):T;  
    tt = [tt(2:d:end),tt(end)];
    ee = [ee(2:d:end),ee(end)];
    semilogy(tt,ee,clr{k},'LineWidth',1.12);
    hold on;
end
% legend('\tau = 2^{-5}','\tau = 2^{-6}','\tau = 2^{-7}','\tau = 2^{-8}','\tau = 2^{-9}');
legend('k = 5','k = 6','k = 7','k = 8','k = 9');
xlabel('t'); ylabel('Error-Y');
title('\tau = T/2^k');
if m >1
    saveas(gcf,varargin{2})
end
hold off;
%END

function myplot22(erru,errv,errw,tau,T,p,varargin)
clr = {'r+-','kd-','mo-','b>-','k*-','c*-','b-'};
p0 = p;
if min(T./tau) < p
    p0 = min(T./tau);
end  
len = length(erru);
figure;
for k = 1:len
    ee = erru{k};   d = round(T/tau(k)/p0); 
    tt = 0:tau(k):T;  
    tt = [tt(2:d:end),tt(end)];
    ee = [ee(2:d:end),ee(end)];
    semilogy(tt,ee,clr{k},'LineWidth',1.12);
    hold on;
end
% legend('\tau = 2^{-5}','\tau = 2^{-6}','\tau = 2^{-7}','\tau = 2^{-8}','\tau = 2^{-9}');
legend('k = 5','k = 6','k = 7','k = 8','k = 9');
xlabel('t'); ylabel('Error-X');
title('\tau = T/2^k');
m = length(varargin);
if m >0
    saveas(gcf,varargin{1})
end
hold off;


len = length(errv);
figure;
for k = 1:len
    ee = errv{k};    d = round(T/tau(k)/p0); 
    tt = 0:tau(k):T;  
    tt = [tt(2:d:end),tt(end)];
    ee = [ee(2:d:end),ee(end)];
    semilogy(tt,ee,clr{k},'LineWidth',1.12);
    hold on;
end
% legend('\tau = 2^{-5}','\tau = 2^{-6}','\tau = 2^{-7}','\tau = 2^{-8}','\tau = 2^{-9}');
legend('k = 5','k = 6','k = 7','k = 8','k = 9');
xlabel('t'); ylabel('Error-Y');
title('\tau = T/2^k');
if m >1
    saveas(gcf,varargin{2})
end
hold off;

len = length(errw);
figure;
for k = 1:len
    ee = errw{k};    d = round(T/tau(k)/p0); 
    tt = 0:tau(k):T;  
    tt = [tt(2:d:end),tt(end)];
    ee = [ee(2:d:end),ee(end)];
    semilogy(tt,ee,clr{k},'LineWidth',1.12);
    hold on;
end
% legend('\tau = 2^{-5}','\tau = 2^{-6}','\tau = 2^{-7}','\tau = 2^{-8}','\tau = 2^{-9}');
legend('k = 5','k = 6','k = 7','k = 8','k = 9');
xlabel('t'); ylabel('Error-Z');
title('\tau = T/2^k');
if m >1
    saveas(gcf,varargin{3})
end
hold off;
%END

function [w,lambda,lmd1,lmd2,lmd3,phi1,phi2] = ini_contour(tau,alf,B,nN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 32;  % N can be changed
F = @(alf,s) s.^(-alf);
len = ceil(log(nN/2)/log(B))+1;
if abs(len - (ceil(log(nN/2)/log(B))+1)) < 1e-14
    len = len + 1;
end
k = 1:len;        nT0 = (2*B.^k-1)*tau;
k = -N:N-1;       
theta = (2*k+1)*pi/2/N;
s0 = 0.4814;  m0 = 0.6443;  n0 = 0.5653;
for kk = 1:len
    x0 = nT0(kk);
    sigma0 = -s0*N/x0;  mu0 = m0*N/x0;   nu0 = n0;
    lambda{kk} = sigma0 + mu0*(theta.*cot(theta) + i*nu0*theta);
    w{kk} = mu0*(cot(theta) - theta.*(csc(theta).^2) + nu0*i)/2/N;
    w{kk} = F(alf,lambda{kk}).*w{kk};
    lmd1{kk} = exp(tau*lambda{kk}); 
    lmd2{kk} = (exp(tau*lambda{kk})-1)./lambda{kk}; 
    lmd3{kk} = (exp(tau*lambda{kk})-1-tau*lambda{kk})./(tau*lambda{kk})./lambda{kk};  
    lmd2{kk} = lmd2{kk} - lmd3{kk};
end
phi1 = sum(imag(exp(lambda{1}*tau).*w{1}./lambda{1}));
phi2 = sum(imag(exp(lambda{1}*tau).*w{1}./(lambda{1}.^2)));
%END

function [w,lambda,lmd1,lmd2,lmd3,phi1,phi2,nT0,N,idn] = ini_contour_gauss(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff2 = exp(z) - 1 - (exp(z) - z - 1)/z;
ff3 = (exp(z) - z - 1)/z;

no = numel(varargin);   T0 = tau;
if no > 0
    T0 = varargin{1};
end
% % % % % % % % % % % % % % % % % % % % % 
vps0 = 1e-16;  vps = 1e-6;
if no > 1
    vps = varargin{2};
end
L = 1:len;  B = ones(1,len)*B;  
T_ell = (2*B-1-B.^(1-L))./(1+(T0/tau-1)*B.^(1-L));
N = round(log(vps)./log(T_ell./(T_ell+1))/2);
% % % % % % % % % % % % % % % % % % % % % 
for kk = 1:len
    x0 = nT0(kk) + T0 - tau;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N(kk),-alf,0,x0);
    lambda{kk} = lambda{kk}'; w{kk} = w{kk}';
% % % % % % % % % % % % % %     
    id = min(N(kk),ceil(2/pi*sqrt((N(kk)+1)*log(N(kk)^(-alf)/vps0))-1));
%     id = min(id,N(kk));
    idn(kk) = id;
    lambda{kk} = lambda{kk}(1:id);   w{kk} = w{kk}(1:id); 
% % % % % % % % % % % % % % 
    lambda{kk} = -lambda{kk};
    w{kk} = sin(alf*pi)/pi*w{kk};
    lmd1{kk} = exp(tau*lambda{kk}); 
%     lmd2{kk} = (exp(tau*lambda{kk})-1)./lambda{kk}; 
%     lmd3{kk} = (exp(tau*lambda{kk})-1-tau*lambda{kk})./(tau*lambda{kk})./lambda{kk};  
%     lmd2{kk} = lmd2{kk} - lmd3{kk};
    z0 = sym(lambda{kk})*sym(tau);
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk};
    lmd3{kk} = double(vpa(subs(ff3,z,z0),60))./lambda{kk};
end
phi1 = tau^alf/gamma(1+alf);
phi2 = tau^(alf+1)/gamma(2+alf);
%END

function [w,lambda,lmd1,lmd2,lmd3,phi1,phi2,N] = ini_contour_trap(tau,alf,nt,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
no = numel(varargin);    
vps0 = 1e-20;  vps = 1e-20;
if no > 1
    vps = varargin{1};
end
T = nt*tau;  vps0 = vps0*T^(alf-1);
xa = log(vps0)/(1-alf);   xb = log(-log(vps0)/tau); 
N = round(1.2*(xb-xa)*log(2^(1-alf)/vps+1)/2/pi^2);
h = (xb-xa)/(N-1);   x = xa:h:xb; x = x';
w = h*sin(alf*pi)/pi*exp((1-alf)*x);
lambda = exp(x);

syms z;
ff2 = exp(z) - 1 - (exp(z) - z - 1)/z;
ff3 = (exp(z) - z - 1)/z;
% % % % % % % % % % % % % % % % % % % % %
lambda = -lambda;
lmd1 = exp(tau*lambda);
% w = w.*lmd1;
z0 = sym(lambda)*sym(tau);
lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda;
lmd3 = double(vpa(subs(ff3,z,z0),60))./lambda;

phi1 = tau^alf/gamma(1+alf);
phi2 = tau^(alf+1)/gamma(2+alf);
%END

function [y2,y3,y4,b2,b3,d3,dd3] = update_1(y2,y3,y4,y6,y66,b2,b3,d3,dd3,yz0,B,n)
% function [y2,y3,y4,b2,b3,d3,dd3] = update_1(y2,y3,y4,y6,y66,b2,b3,d3,dd3,yz0,B,n)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);
if L0 == L
    for kk = 1:L-1
        if b(kk) == b0(kk) && b(kk+1) == b0(kk+1)
        elseif b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            if kk == L-1 && L > 2
                y2{kk} = y6;
            else
                y2{kk} = y66{kk};
            end
        else
            y2{kk} = y3{kk};
            if n == b(kk+1) + B^kk
                y3{kk} = 0;  b2(kk,:) = b(kk:kk+1) + B^kk;
                b3(kk) = b2(kk,1) + B^(kk-1);
            end
        end
    end
else
    d3 = tauL1(B,n*B);
    for kk = 1:L-1
        y2{kk} = y4{kk};  y3{kk} = 0;
        b2(kk,:) = b(kk:kk+1) + B^kk;
        b3(kk) = b2(kk,1) + B^(kk-1);
        dd3(kk) = d3(kk) + B^(kk-1);
    end
    dd3(L) = d3(L) + B^(L-1);    y4 = yz0;
end
%END

function [y3,y33,y4,y44,yy4,y6,y66,y77] = contour_2......
    (len,y3,y33,y4,y44,yy4,y6,y66,y77,b2,b3,d3,dd3,B,n,lmd1,lmd2,lmd3,u1,u2)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);
for kk = 1:L-1
    if n > d3(kk+1) && n < d3(kk)+1
        y4{kk} = lmd1{kk}.*y4{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
        y44{kk} = y4{kk};
    end
    if n > d3(kk) && n < dd3(kk)
        y44{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
    end
end
% % -------------------------------------------
if L0 == L
    for kk = 1:L-1
        y77{kk} = lmd1{kk}.*y77{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
        if n > b2(kk,2) && n < b2(kk,1) + 1
            y3{kk} = lmd1{kk}.*y3{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
            y33{kk} = y3{kk};
        end
        if n > b2(kk,1) && n < b3(kk)
            y33{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
        end
        if b(kk) ~= b0(kk) && b(kk+1) ~= b0(kk+1)
            y66{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
            y77{kk} = y66{kk};
        end
        if b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            y66{kk} = y77{kk};
        end
    end
end
% % ---------------------------------------------
for kk = L-1:len
    yy4{kk} = lmd1{kk}.*yy4{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
end
if n == B^(L-1)
    y4{L} = yy4{L};
end
if L ~= L0
    for kk = 1:L-1
        y66{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u1 + lmd3{kk}*u2;
        y77{kk} = y66{kk};
    end
    y6 = yy4{L-1};
elseif b(L-1) ~= b0(L-1)
    y6 = yy4{L-1};
end
%END

function [y3,y33,y4,y44,yy4,y6,y66,y77] = update_2_p1......
    (len,y3,y33,y4,y44,yy4,y6,y66,y77,b2,b3,d3,dd3,B,n,lmd1,lmd2,lmd3,u1,u2,cw)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);
for kk = 1:L-1
    if n > d3(kk+1) && n < d3(kk)+1
        y4{kk} = lmd1{kk}.*y4{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
        y44{kk} = y4{kk};
    end
    if n > d3(kk) && n < dd3(kk)
        y44{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
    end
end
% % -------------------------------------------
if L0 == L
    for kk = 1:L-1
        y77{kk} = lmd1{kk}.*y77{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
        if n > b2(kk,2) && n < b2(kk,1) + 1
            y3{kk} = lmd1{kk}.*y3{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
            y33{kk} = y3{kk};
        end
        if n > b2(kk,1) && n < b3(kk)
            y33{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
        end
        if b(kk) ~= b0(kk) && b(kk+1) ~= b0(kk+1)
            y66{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
            y77{kk} = y66{kk};
        end
        if b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            y66{kk} = y77{kk};
        end
    end
end
% % ---------------------------------------------
for kk = L-1:len
    yy4{kk} = lmd1{kk}.*yy4{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
end
if n == B^(L-1)
    y4{L} = yy4{L};
end
if L ~= L0
    for kk = 1:L-1
        y66{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u1 + lmd3{kk}*u2 + cw{kk}(n,:);
        y77{kk} = y66{kk};
    end
    y6 = yy4{L-1};
elseif b(L-1) ~= b0(L-1)
    y6 = yy4{L-1};
end
%END

function b = tauL1(B,n)
L = ceil(log(n/2)/log(B));
if n == 2*B^L
    L = L + 1;
end
b = zeros(1,L+1);
b(1) = n-1;  vps = 1e-14;
if L>1
    BK = B.^(1:L-1);
    q1 = (n+1)./BK-2 - vps;  q2 = n./BK - 1 + vps;
    q = round((ceil(q1) + floor(q2))/2);
%     q = floor(q2);
    b(2:L) = q.*BK;
end
if n ==1
    b = [0,0];
end
%END

function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2] = ini_contour_p2(tau,alf,B,nN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 32;  % N can be changed
F = @(alf,s) s.^(-alf);
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
len = ceil(log(nN/2)/log(B))+2;    
k = 1:len;        nT0 = (2*B.^k-1)*tau;
k = -N:N-1;       
theta = (2*k+1)*pi/2/N;
s0 = 0.4814;  m0 = 0.6443;  n0 = 0.5653; 
for kk = 1:len
    x0 = nT0(kk);
    sigma0 = -s0*N/x0;  mu0 = m0*N/x0;   nu0 = n0;
    lambda{kk} = sigma0 + mu0*(theta.*cot(theta) + i*nu0*theta);
    w{kk} = mu0*(cot(theta) - theta.*(csc(theta).^2) + nu0*i)/2/N;
    w{kk} = F(alf,lambda{kk}).*w{kk};
    z0 = sym(lambda{kk})*sym(tau);
    lmd1{kk} = exp(tau*lambda{kk});    
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
    lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
    lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END

function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT0] = ini_contour_p2_gauss(tau,alf,B,nN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 80;  % N can be changed
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;

for kk = 1:len
    x0 = nT0(kk);
    [lambda{kk},w{kk}] = gen_laguerre_rule(N,-alf,0,x0);
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    z0 = sym(lambda{kk})*sym(tau);
    lmd1{kk} = exp(tau*lambda{kk});    
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
    lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
    lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END


function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT0] = ini_sisc_LvXu(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% varargin = (dT,vps,N)
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
% syms z;
% ff2 = ((2-3*z+2*z^2)*exp(z) + z - 2)/z^2;
% ff3 = ((2-2*z)*exp(z) + z^2 - 2)/z^2;
% ff4 = ((2-z)*exp(z) - z - 2)/z^2;

g2 = @(z) ((2-3*z+2*z.^2).*exp(z) + z - 2)./z.^3;
g3 = @(z) ((2-2*z).*exp(z) + z.^2 - 2)./z.^3;
g4 = @(z) ((2-z).*exp(z) - z - 2)./z.^3;

f2 = @(z) (2-3*z+2*z.^2).*exp_n(z) + z + 1/2;
f3 = @(z) (2-2*z).*exp_n(z) - 1;
f4 = @(z) (2-z).*exp_n(z) - 1/2;

% z0 = -0.39:0.05:0; 
% e2 = (f2(z0) - g2(z0))./g2(z0); 
% e3 = (f3(z0) - g3(z0))./g3(z0); 
% e4 = (f4(z0) - g4(z0))./g4(z0);
vps0 = 1e-16;  vps = 1e-12; T0 = tau;
no = numel(varargin);   
if no > 0
    T0 = varargin{1};
end
% % % % % % % % % % % % % % % % % % % % % 
if no > 1
    vps = varargin{2};
end
L = 1:len;  B = ones(1,len)*B;  
T_ell = (2*B-1-B.^(1-L))./(1+(T0/tau-1)*B.^(1-L));
N = round(log(vps)./log(T_ell./(T_ell+1))/2);
% % % % % % % % % % % % % % % % % % % % % 
if no > 2
    N0 = varargin{3};
    if length(N0)==1
        N = ones(1,length(N))*N0;
    else
        N = N0;
    end
end

for kk = 1:len
    x0 = nT0(kk) + T0 - tau;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N(kk),-alf,0,x0);
    lambda{kk} = -lambda{kk}';
% % % % % % % % % % % %     
    id = min(N(kk),ceil(2/pi*sqrt((N(kk)+1)*log(N(kk)^(-alf)/vps0))-1));
    lambda{kk} = lambda{kk}(1:id);   w{kk} = w{kk}(1:id); 
% % % % % % % % % % % % % %     
    w{kk} = sin(alf*pi)/pi*w{kk}';
    lmd1{kk} = exp(tau*lambda{kk});   
    z0 = (lambda{kk})*(tau);
    id = abs(z0)<0.4;    z1 = z0(id);
    lmd2{kk} = tau*g2(z0)/2;  
    lmd3{kk} = -tau*g3(z0);   
    lmd4{kk} = tau*g4(z0)/2;  
    lmd2{kk}(id) = tau*f2(z1)/2;
    lmd3{kk}(id) = -tau*f3(z1);
    lmd4{kk}(id) = tau*f4(z1)/2;
    
%     z0 = sym(lambda{kk})*sym(tau);
%     lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
%     lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
%     lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END

function [w,lambda,lmd11,lmd22,lmd33,lmd44,phi00,phi11,phi22]....
    = ini_contour_p2_gauss_2(tau,alf,B,nN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 80;  % N can be changed
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff22 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
ff33 = ((2-z)*exp(z) - 2 - z)/z^2;
ff44 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;

for kk = 1:len
    x0 = nT0(kk);
    [lambda{kk},w{kk}] = gen_laguerre_rule(N,-alf,0,x0);
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    zz0 = sym(2*lambda{kk})*sym(tau);
    lmd11{kk} = exp(2*tau*lambda{kk});
    lmd22{kk} = double(vpa(subs(ff22,z,zz0),60))./lambda{kk};
    lmd33{kk} = -double(vpa(subs(ff33,z,zz0),60))./lambda{kk}*4;
    lmd44{kk} = double(vpa(subs(ff44,z,zz0),60))./lambda{kk};
end
phi00 = alf^2*2^alf*tau^(alf)/gamma(3+alf);
phi11 = 4*alf*2^alf*tau^(alf)/gamma(3+alf);
phi22=(2-alf)*2^alf*tau^(alf)/gamma(3+alf);
%END

function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT0]....
    = ini_contour_p2_gauss_b(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 40;  % N can be changed
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
no = numel(varargin);
T0 = 0;
if no>0
    T0 = varargin{1} - tau;
end
if no>1
    N = varargin{2};
end
for kk = 1:len
    x0 = nT0(kk) + T0;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N,-alf,0,x0);
    id = x0^(1-alf)*w{kk} > 1e-16;  
    lambda{kk} = lambda{kk}(id);   w{kk} = w{kk}(id); 
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    z0 = sym(lambda{kk})*sym(tau);
    lmd1{kk} = exp(tau*lambda{kk});    
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
    lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
    lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END

function [w,lambda,lmd11,lmd22,lmd33,lmd44,phi00,phi11,phi22]....
    = ini_contour_p2_gauss_2_b(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 40;  % N can be changed
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff22 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
ff33 = ((2-z)*exp(z) - 2 - z)/z^2;
ff44 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;

no = numel(varargin);
T0 = 0;
if no>0
    T0 = varargin{1} - tau;
end
if no>1
    N = varargin{2};
end
for kk = 1:len
    x0 = nT0(kk) + T0;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N,-alf,0,x0);
    id = x0^(1-alf)*w{kk} > 1e-16;   
    lambda{kk} = lambda{kk}(id);   w{kk} = w{kk}(id); 
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    zz0 = sym(2*lambda{kk})*sym(tau);
    lmd11{kk} = exp(2*tau*lambda{kk});
    lmd22{kk} = double(vpa(subs(ff22,z,zz0),60))./lambda{kk};
    lmd33{kk} = -double(vpa(subs(ff33,z,zz0),60))./lambda{kk}*4;
    lmd44{kk} = double(vpa(subs(ff44,z,zz0),60))./lambda{kk};
end
phi00 = alf^2*2^alf*tau^(alf)/gamma(3+alf);
phi11 = 4*alf*2^alf*tau^(alf)/gamma(3+alf);
phi22=(2-alf)*2^alf*tau^(alf)/gamma(3+alf);
%END

function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT0]....
    = initial_contour_gauss_p2(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
vps0 = 1e-16;  vps = 1e-12;
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
no = numel(varargin);
T0 = tau;
if no>0
    T0 = varargin{1};
end
if no>1
    vps = varargin{2};
end

L = 1:len;  B = ones(1,len)*B;  
T_ell = (2*B-1-B.^(1-L))./(1+(T0/tau-1)*B.^(1-L));
N = round(log(vps)./log(T_ell./(T_ell+1))/2);
% N = ones(1,length(N))*128;

for kk = 1:len
    x0 = nT0(kk) + T0 - tau;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N(kk),-alf,0,x0);
    id = min(N(kk),ceil(2/pi*sqrt((N(kk)+1)*log(N(kk)^(-alf)/vps0))-1));
    lambda{kk} = lambda{kk}(1:id);   w{kk} = w{kk}(1:id); 
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    z0 = sym(lambda{kk})*sym(tau);
    lmd1{kk} = exp(tau*lambda{kk});    
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
    lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
    lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END

function [w,lambda,lmd11,lmd22,lmd33,lmd44,phi00,phi11,phi22,nT0]....
    = initial_contour_gauss_p2_b(tau,alf,B,nN,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
vps0 = 1e-16;  vps = 1e-12;
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff22 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
ff33 = ((2-z)*exp(z) - 2 - z)/z^2;
ff44 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;

no = numel(varargin);
T0 = tau;
if no>0
    T0 = varargin{1};
end
if no>1
    vps = varargin{2};
end

L = 1:len;  B = ones(1,len)*B;  
T_ell = (2*B-1-B.^(1-L))./(1+(T0/tau-1)*B.^(1-L));
N = round(log(vps)./log(T_ell./(T_ell+1))/2);
% N = ones(1,length(N))*128;

for kk = 1:len
    x0 = nT0(kk) + T0 - tau;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N(kk),-alf,0,x0);
    id = min(N(kk),ceil(2/pi*sqrt((N(kk)+1)*log(N(kk)^(-alf)/vps0))-1));
%     id = x0^(1-alf)*w{kk} > 1e-16;   
    lambda{kk} = lambda{kk}(1:id);   w{kk} = w{kk}(1:id); 
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    zz0 = sym(2*lambda{kk})*sym(tau);
    lmd11{kk} = exp(2*tau*lambda{kk});
    lmd22{kk} = double(vpa(subs(ff22,z,zz0),60))./lambda{kk};
    lmd33{kk} = -double(vpa(subs(ff33,z,zz0),60))./lambda{kk}*4;
    lmd44{kk} = double(vpa(subs(ff44,z,zz0),60))./lambda{kk};
end
phi00 = alf^2*2^alf*tau^(alf)/gamma(3+alf);
phi11 = 4*alf*2^alf*tau^(alf)/gamma(3+alf);
phi22=(2-alf)*2^alf*tau^(alf)/gamma(3+alf);
%END


function [w,lambda,lmd11,lmd22,lmd33,lmd44,phi00,phi11,phi22]....
    = ini_contour_p2_2(tau,alf,B,nN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
N = 32;  % N can be changed
F = @(alf,s) s.^(-alf);
syms z;
ff22 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
ff33 = ((2-z)*exp(z) - 2 - z)/z^2;
ff44 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;

len = ceil(log(nN/2)/log(B))+2;    
k = 1:len;        nT0 = (2*B.^k-1)*tau;
k = -N:N-1;       
theta = (2*k+1)*pi/2/N;
s0 = 0.4814;  m0 = 0.6443;  n0 = 0.5653;
for kk = 1:len
    x0 = nT0(kk);
    sigma0 = -s0*N/x0;  mu0 = m0*N/x0;   nu0 = n0;
    lambda{kk} = sigma0 + mu0*(theta.*cot(theta) + i*nu0*theta);
    w{kk} = mu0*(cot(theta) - theta.*(csc(theta).^2) + nu0*i)/2/N;
    w{kk} = F(alf,lambda{kk}).*w{kk};
    zz0 = sym(2*lambda{kk})*sym(tau);
    lmd11{kk} = exp(2*tau*lambda{kk});
    lmd22{kk} = double(vpa(subs(ff22,z,zz0),60))./lambda{kk};
    lmd33{kk} = -double(vpa(subs(ff33,z,zz0),60))./lambda{kk}*4;
    lmd44{kk} = double(vpa(subs(ff44,z,zz0),60))./lambda{kk};
end
phi00 = alf^2*2^alf*tau^(alf)/gamma(3+alf);
phi11 = 4*alf*2^alf*tau^(alf)/gamma(3+alf);
phi22=(2-alf)*2^alf*tau^(alf)/gamma(3+alf);
%END

function [y3,y33,y4,y44,yy4,y6,y66,y77] = contour_2_p2......
    (len,y3,y33,y4,y44,yy4,y6,y66,y77,b2,b3,d3,dd3,B,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);
for kk = 1:L-1
    if n > d3(kk+1) && n < d3(kk)+1
        y4{kk} = lmd1{kk}.*y4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        y44{kk} = y4{kk};
    end
    if n > d3(kk) && n < dd3(kk)
        y44{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
    end
end
% % -------------------------------------------
if L0 == L
    for kk = 1:L-1
        y77{kk} = lmd1{kk}.*y77{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        if n > b2(kk,2) && n < b2(kk,1) + 1
            y3{kk} = lmd1{kk}.*y3{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
            y33{kk} = y3{kk};
        end
        if n > b2(kk,1) && n < b3(kk)
            y33{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        end
        if b(kk) ~= b0(kk) && b(kk+1) ~= b0(kk+1)
            y66{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
            y77{kk} = y66{kk};
        end
        if b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            y66{kk} = y77{kk};
        end
    end
end
% % ---------------------------------------------
for kk = L-1:len
    yy4{kk} = lmd1{kk}.*yy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
end
if n == B^(L-1)
    y4{L} = yy4{L};
end
if L ~= L0
    for kk = 1:L-1
        y66{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2;
        y77{kk} = y66{kk};
    end
    y6 = yy4{L-1};
elseif b(L-1) ~= b0(L-1)
    y6 = yy4{L-1};
end
%END

function [y3,y33,y4,y44,yy4,y6,y66,y77] = update_p2_correction......
    (len,y3,y33,y4,y44,yy4,y6,y66,y77,b2,b3,d3,dd3,B,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2,cw)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);   
for kk = 1:L-1
    if n > d3(kk+1) && n < d3(kk)+1
        y4{kk} = lmd1{kk}.*y4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
        y44{kk} = y4{kk};
    end
    if n > d3(kk) && n < dd3(kk)
        y44{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
    end
end
% % -------------------------------------------
if L0 == L
    for kk = 1:L-1
        y77{kk} = lmd1{kk}.*y77{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
        if n > b2(kk,2) && n < b2(kk,1) + 1
            y3{kk} = lmd1{kk}.*y3{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
            y33{kk} = y3{kk};
        end
        if n > b2(kk,1) && n < b3(kk)
            y33{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
        end
        if b(kk) ~= b0(kk) && b(kk+1) ~= b0(kk+1)
            y66{kk} = lmd1{kk}.*y33{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
            y77{kk} = y66{kk};
        end
        if b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            y66{kk} = y77{kk};
        end
    end
end
% % ---------------------------------------------
for kk = L-1:len
    yy4{kk} = lmd1{kk}.*yy4{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
end
if n == B^(L-1)
    y4{L} = yy4{L};
end
if L ~= L0
    for kk = 1:L-1
        y66{kk} = lmd1{kk}.*y44{kk} + lmd2{kk}*u0 + lmd3{kk}*u1 + lmd4{kk}*u2 + cw{kk}(n-1,:);
        y77{kk} = y66{kk};
    end
    y6 = yy4{L-1};
elseif b(L-1) ~= b0(L-1)
    y6 = yy4{L-1};
end
%END


function [y3,y33,y4,y44,yy4,y6,y66,y77] = update_p2_correction_pde......
    (len,y3,y33,y4,y44,yy4,y6,y66,y77,b2,b3,d3,dd3,B,n,lmd1,lmd2,lmd3,lmd4,u0,u1,u2,cw)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b0 = tauL1(B,n-1);   L0 = length(b0);
b = tauL1(B,n);      L = length(b);   
for kk = 1:L-1
    if n > d3(kk+1) && n < d3(kk)+1
        y4{kk} = (lmd1{kk}.*y4{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}...
            + cw{kk}(:,:,n-1);
        y44{kk} = y4{kk};
    end
    if n > d3(kk) && n < dd3(kk)
        y44{kk} = (lmd1{kk}.*y44{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}....
            + cw{kk}(:,:,n-1);
    end
end
% % -------------------------------------------
if L0 == L
    for kk = 1:L-1
        y77{kk} = (lmd1{kk}.*y77{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}....
            + cw{kk}(:,:,n-1);
        if n > b2(kk,2) && n < b2(kk,1) + 1
            y3{kk} = (lmd1{kk}.*y3{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}....
                + cw{kk}(:,:,n-1);
            y33{kk} = y3{kk};
        end
        if n > b2(kk,1) && n < b3(kk)
            y33{kk} = (lmd1{kk}.*y33{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}...
                + cw{kk}(:,:,n-1);
        end
        if b(kk) ~= b0(kk) && b(kk+1) ~= b0(kk+1)
            y66{kk} = (lmd1{kk}.*y33{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}....
                + cw{kk}(:,:,n-1);
            y77{kk} = y66{kk};
        end
        if b(kk) ~= b0(kk) && b(kk+1) == b0(kk+1)
            y66{kk} = y77{kk};
        end
    end
end
% % ---------------------------------------------
for kk = L-1:len
    yy4{kk} = (lmd1{kk}.*yy4{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}....
        + cw{kk}(:,:,n-1);
end
if n == B^(L-1)
    y4{L} = yy4{L};
end
if L ~= L0
    for kk = 1:L-1
        y66{kk} = (lmd1{kk}.*y44{kk}) + u0*lmd2{kk} + u1*lmd3{kk} + u2*lmd4{kk}...
            + cw{kk}(:,:,n-1);
        y77{kk} = y66{kk};
    end
    y6 = yy4{L-1};
elseif b(L-1) ~= b0(L-1)
    y6 = yy4{L-1};
end
%END

function [w,lmd2,lmd3,lmd4] = correction_weights_1(tau,T,lambda,sgm) 
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
t = 0:tau:T+1e-14;
z0 = sym(lambda)*sym(tau); 
lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda/2;
lmd3 = -double(vpa(subs(ff3,z,z0),60))./lambda;
lmd4 = double(vpa(subs(ff4,z,z0),60))./lambda/2;

N = round(T/tau) - 1;
m = length(sgm); w = zeros(1,N);
if m == 0
    return;
end
[xx0,ww0]=RootsJacobiLobatto(0,0,16); ww0 = ww0.';
T0 = 10*tau;   N0 = round(T0/tau);  

b = zeros(m,N);  arr = zeros(m); 
for k = 1:m
    s = sgm(k);   tt = t.^s;    tt1 = t.^(s+1);
    arr(k,:) = tt(2:m+1);
    r = 2^(-1-s);
    [x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
    for n = 1:N0-1
        xx1 = t(n+1)*(1+x0)/2;   f1 = exp(-lambda*(xx1 - t(n+2)));
        xx2 = t(n+2)*(1+x0)/2;   f2 = exp(-lambda*(xx2 - t(n+2)));
        b(k,n) = r*(w0*f2*tt1(n+2) - w0*f1*tt1(n+1))......
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
    for n = N0:N
        xx = (tau*xx0 + t(n+1) + t(n+2))/2;
        ff = exp(-lambda*(xx - t(n+2))).*xx.^s;
        b(k,n) = tau/2*sum(ww0*ff).....
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
end
w = arr\b;
%END

function w = correction_weights_11(tau,T,lambda,lmd2,lmd3,lmd4,xx0,ww0,sgm,r0,wt0) 
% syms z;
% ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
% ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
% ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
% z0 = sym(lambda)*sym(tau); 
% lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda/2;
% lmd3 = -double(vpa(subs(ff3,z,z0),60))./lambda;
% lmd4 = double(vpa(subs(ff4,z,z0),60))./lambda/2;

N = round(T/tau) - 1;
m = length(sgm);  w = zeros(1,N);
if m == 0
    return;
end
b = zeros(m,N);  arr = zeros(m); 

T0 = 10*tau;   N0 = round(T0/tau);
if N0 == 0
    N0 = 1;
end
t = 0:tau:max((N+2),N0+1)*tau+1e-14;

for k = 1:m
    s = sgm(k);   tt = t.^s;    tt1 = t.^(s+1);
    arr(k,:) = tt(2:m+1);
    r = exp(lambda*tau)*2^(-1-s);
%     [x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
    x0 = xx0{k}; w0 = ww0{k}.';
    for n = 1:N0
        xx1 = t(n+1)*(1+x0)/2;   f1 = exp(-lambda*(xx1 - t(n+1)));
        xx2 = t(n+2)*(1+x0)/2;   f2 = exp(-lambda*(xx2 - t(n+1)));
        b(k,n) = r*(w0*f2*tt1(n+2) - w0*f1*tt1(n+1))......
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
    for n = N0+1:N-1
        xx = (tau*r0 + t(n+1) + t(n+2))/2;
        ff = exp(-lambda*(xx - t(n+2))).*xx.^s;
        b(k,n) = tau/2*sum(wt0'*ff).....
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
end
w = arr\b;
%END

function w = correction_weights_p1(tau,T,lambda,lmd3,lmd4,xx0,ww0,sgm,r0,wt0) 
% syms z;
% ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
% ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
% ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
% z0 = sym(lambda)*sym(tau); 
% lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda/2;
% lmd3 = -double(vpa(subs(ff3,z,z0),60))./lambda;
% lmd4 = double(vpa(subs(ff4,z,z0),60))./lambda/2;

N = round(T/tau)-1;
m = length(sgm); w = zeros(1,N);
if m == 0
    return;
end
b = zeros(m,N+1);  arr = zeros(m); 
t = 0:tau:(N+1)*tau+1e-14;
T0 = 10*tau;   N0 = round(T0/tau);
if N0 == 0
    N0 = 1;
end
for k = 1:m
    s = sgm(k);   tt = t.^s;    tt1 = t.^(s+1);
    arr(k,:) = tt(2:m+1);
    r = exp(lambda*tau)*2^(-1-s);
%     [x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
    x0 = xx0{k}; w0 = ww0{k}.';
    for n = 0:N0-1
        xx1 = t(n+1)*(1+x0)/2;   f1 = exp(-lambda*(xx1 - t(n+1)));
        xx2 = t(n+2)*(1+x0)/2;   f2 = exp(-lambda*(xx2 - t(n+1)));
        b(k,n+1) = r*(w0*f2*tt1(n+2) - w0*f1*tt1(n+1)) - (lmd3*tt(n+1) + lmd4*tt(n+2));
    end
    for n = N0:N
        xx = (tau*r0 + t(n+1) + t(n+2))/2;
        ff = exp(-lambda*(xx - t(n+2))).*xx.^s;
        b(k,n+1) = tau/2*sum(wt0'*ff) - (lmd3*tt(n+1) + lmd4*tt(n+2));
    end
end
w = arr\b;
%END

function [w,lmd2,lmd3,lmd4] = correction_weights_2(tau,T,lambda,sgm) 
syms z;
ff2 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
ff3 = ((2-z)*exp(z) - 2 - z)/z^2;
ff4 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;
t = 0:tau:T+1e-14;
z0 = sym(lambda)*sym(2*tau); 
lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda;
lmd3 = -double(vpa(subs(ff3,z,z0),60))./lambda*4;
lmd4 = double(vpa(subs(ff4,z,z0),60))./lambda;

N = round(T/tau) - 1;
m = length(sgm); w = zeros(1,N);
if m == 0
    return;
end
b = zeros(m,N);  arr = zeros(m); 
for k = 1:m
    s = sgm(k);   tt = t.^s;    tt1 = t.^(s+1);
    arr(k,:) = tt(2:m+1);
    r = exp(2*lambda*tau)*2^(-1-s);
    [x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
    for n = 1:N-1
        xx1 = t(n)*(1+x0)/2;     f1 = exp(-lambda*(xx1 - t(n)));
        xx2 = t(n+2)*(1+x0)/2;   f2 = exp(-lambda*(xx2 - t(n)));
        b(k,n) = r*(w0*f2*tt1(n+2) - w0*f1*tt1(n))......
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
end
w = arr\b;
%END
function w = correction_weights_22(tau,T,lambda,lmd2,lmd3,lmd4,xx0,ww0,sgm,r0,wt0) 
% syms z;
% ff2 = ((4-3*z+z^2)*exp(z) - 4 - z)/z^2;
% ff3 = ((2-z)*exp(z) - 2 - z)/z^2;
% ff4 = ((4-z)*exp(z) - 4 - 3*z - z^2)/z^2;
% 
% z0 = sym(lambda)*sym(2*tau); 
% lmd2 = double(vpa(subs(ff2,z,z0),60))./lambda;
% lmd3 = -double(vpa(subs(ff3,z,z0),60))./lambda*4;
% lmd4 = double(vpa(subs(ff4,z,z0),60))./lambda;

N = round(T/tau) - 1;
t = 0:tau:(N+1)*tau+1e-14;
T0 = T/10;
if T0 > 0.1
    T0 = 0.1;
end
N0 = round(T0/tau);  

if N0 == 0
    N0 = 1;
end
m = length(sgm); w = zeros(1,N);
if m == 0
    return;
end
b = zeros(m,N);  arr = zeros(m); 
for k = 1:m
    s = sgm(k);   tt = t.^s;    tt1 = t.^(s+1);
    arr(k,:) = tt(2:m+1);
    r = exp(2*lambda*tau)*2^(-1-s);
%     [x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
    x0 = xx0{k}; w0 = ww0{k}.';
    for n = 1:N0
        xx1 = t(n)*(1+x0)/2;     f1 = exp(-lambda*(xx1 - t(n)));
        xx2 = t(n+2)*(1+x0)/2;   f2 = exp(-lambda*(xx2 - t(n)));
        b(k,n) = r*(w0*f2*tt1(n+2) - w0*f1*tt1(n))......
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
    for n = N0+1:N-1
        xx = (2*tau*r0 + t(n) + t(n+2))/2;
        ff = exp(-lambda*(xx - t(n+2))).*xx.^s;
        b(k,n) = tau*sum(wt0'*ff).....
            - (lmd2*tt(n) + lmd3*tt(n+1) + lmd4*tt(n+2));
    end
end
w = arr\b;
%END

function test_correction_weights_p1
tau = 1e-4;
T = 1;
lambda = 1.5;
m = 4;
sigma = [0.1 0.2 0.3 0.4 1 2 3 4];
a = [1 1 1 1 0 0 0 0];
len = length(sigma);
sgm = sigma(1:m);
for k = 1:m
    s = sgm(k); 
    [x0,w0] = RootsJacobiLobatto(0,s,64);  
    xx0{k} = x0; ww0{k} = w0;
end
[r0,wt0] = RootsJacobiLobatto(0,0,32); 
z = lambda*tau;
lmd2 = ((z-1)*exp(z)+1)/lambda/z;
lmd3 = exp(z) - 1 - z;  lmd3 = lmd3/lambda/z;
w = correction_weights_p1(tau,T,lambda,lmd2,lmd3,xx0,ww0,sgm,r0,wt0);

t = 0:tau:T+1e-14;
d = 1.4;
uv = @(t) d + a(1)*t.^sigma(1) + a(2)*t.^sigma(2) + a(3)*t.^sigma(3) + a(4)*t.^sigma(4)....
     + a(5)*t.^sigma(5) + a(6)*t.^sigma(6) + a(7)*t.^sigma(7) + a(8)*t.^sigma(8);
for k = 1:len
    s = sigma(k);
    bb2(k,:) = lambda^(-1-s)*exp(lambda*t(2:end))....
        .*(gammainc(lambda*t(2:end),s+1) - gammainc(lambda*t(1:end-1),s+1));
    bb2(k,:) = a(k)*bb2(k,:)*gamma(s+1);
end
f1 = sum(bb2) + d*(exp(lambda*tau)-1)/lambda;
u = uv(t);  uk0 = u(2:m+1) - u(1);
if m == 0
    uk0 = 0;
end
% f2 = lmd2*u(1:end-1) + lmd3*u(2:end) + uk0*w....
%     + u(1)*(exp(lambda*tau)-1)/lambda - (lmd2 + lmd3)*u(1); 
f2 = lmd2*u(1:end-1) + lmd3*u(2:end) + uk0*w;
e = abs(f1 - f2);
figure(2); semilogy(t(2:end),e);
ee = norm(e,inf)
%END

function test_correction_weights_1
tau = 1e-2;
T = 40*tau;
lambda = -20;
m = 4;
sigma = [0.1 0.2 0.3 0.4 1 2 3 4];
a = [1 1 1 1 0 0 0 0];
len = length(sigma);
sgm = sigma(1:m);
[w,lmd2,lmd3,lmd4] = correction_weights_1(tau,T,lambda,sgm); 
t = 0:tau:T+1e-14;
d = 2.1;
uv = @(t) d + a(1)*t.^sigma(1) + a(2)*t.^sigma(2) + a(3)*t.^sigma(3) + a(4)*t.^sigma(4)....
     + a(5)*t.^sigma(5) + a(6)*t.^sigma(6) + a(7)*t.^sigma(7) + a(8)*t.^sigma(8);
for k = 1:len
    s = sigma(k);
    bb2(k,:) = lambda^(-1-s)*exp(lambda*t(3:end))....
        .*(gammainc(lambda*t(3:end),s+1) - gammainc(lambda*t(2:end-1),s+1));
    bb2(k,:) = a(k)*bb2(k,:)*gamma(s+1);
end
f1 = sum(bb2) + d*(exp(lambda*tau)-1)/lambda;
u = uv(t);  uk0 = u(2:m+1) - u(1);
if m == 0
    uk0 = 0;
end
% f2 = lmd2*u(1:end-2) + lmd3*u(2:end-1) + lmd4*u(3:end) + uk0*w....
%     + u(1)*(exp(lambda*tau)-1)/lambda - (lmd2 + lmd3 + lmd4)*u(1); 
f2 = lmd2*u(1:end-2) + lmd3*u(2:end-1) + lmd4*u(3:end) + uk0*w;
% (lmd2 + lmd3 + lmd4)-(exp(lambda*tau)-1)/lambda,
e = abs(f1 - f2);
figure(2); semilogy(t(1:end-2),e);
ee = norm(e,inf)/norm(f2,inf)
%END

function test_correction_weights_2
tau = 1e-4;
T = 1;
lambda = 2;
m = 5;
sigma = [0.1 0.2 0.3 0.4 1 2 3 4];
len = length(sigma);
sgm = sigma(1:m);
[w,lmd2,lmd3,lmd4] = correction_weights_2(tau,T,lambda,sgm); 
t = 0:tau:T+1e-14;
d = 2.1;
uv = @(t) d + t.^sigma(1) + t.^sigma(2) + t.^sigma(3) + t.^sigma(4)....
     + t.^sigma(5) + t.^sigma(6) + t.^sigma(7) + t.^sigma(8);
for k = 1:len
    s = sigma(k);
    bb2(k,:) = lambda^(-1-s)*exp(lambda*t(3:end))....
        .*(gammainc(lambda*t(3:end),s+1) - gammainc(lambda*t(1:end-2),s+1));
    bb2(k,:) = bb2(k,:)*gamma(s+1);
end
f1 = sum(bb2) + d*(exp(2*lambda*tau)-1)/lambda;
u = uv(t);  uk0 = u(2:m+1) - u(1);
if m == 0
    uk0 = 0;
end
% f2 = lmd2*u(1:end-2) + lmd3*u(2:end-1) + lmd4*u(3:end) + uk0*w....
%     + (exp(2*lambda*tau)-1)/lambda*u(1) - (lmd2 + lmd3 + lmd4)*u(1);
f2 = lmd2*u(1:end-2) + lmd3*u(2:end-1) + lmd4*u(3:end) + uk0*w;
e = f1 - f2;
ee = norm(e,inf)
%END
function test_3
tau = 1e-2;
T = 1;
sgm = 0.1;
t = 0:tau:T+1e-14;
N = round(T/tau);
s = sgm;
[x0,w0]=RootsJacobiLobatto(0,s,64); w0 = w0.';
for n = 1:N
    xx1 = (tau*x0 + t(n) + t(n+1))/2;   f1 = xx1.^s;
    bb1(1,n) = (w0*f1)*tau/2;
end
bb2 = (t(2:end).^(s+1) - t(1:end-1).^(s+1))/(1+s);
e = abs(bb1 - bb2);
semilogy(t(3:end),e(2:end));
ee = norm(e,inf);
ee;
%END
function test_coefficients_p2_2
tau = 1e-3; alf = -0.4;
T = 1;
N = round(T/tau);
t = 0:tau:T+1e-12;
s1 = 0.5; s2 = 0.6; s3 = 3;
sgm = [s1 s2];
no = length(sgm);
[a,b,c,d,w] = coefficients_p2_2(alf,N+1,sgm);
lambda = 1;  lambda0 = 0;
uv = @(t) lambda0 + t.^s1 + t.^s2 + lambda*t.^s3;
fv = @(t) lambda0*t.^(alf)/gamma(1+alf) + t.^(s1+alf)*gamma(s1+1)/gamma(s1+1+alf)....
    + t.^(s2+alf)*gamma(s2+1)/gamma(s2+1+alf)...
    + lambda*gamma(s3+1)/gamma(s3+1+alf)*t.^(s3+alf);
u = uv(t); du = zeros(1,N);
k = 1; du(k) = d(1,k)*u(1) + d(2,k)*u(2) + d(3,k)*u(3) + u(2:no+1)*w(:,k);
for k = 2:N
    du(k) = d(1,k)*u(1) + d(2,k)*u(2) + d(3,k)*u(3);
    du(k) = du(k) + sum(a(1:k-1).*u(k-1:-1:1)) + sum(b(1:k-1).*u(k:-1:2))....
        + sum(c(1:k-1).*u(k+1:-1:3)) + u(2:no+1)*w(:,k); 
end
du = tau^alf*du;
du2 = fv(t(2:end));
e = du - du2;
ee = norm(e,inf)
ee;
%END

function test_LvXu_sisc_2016
tau = 1e-3; alf = -0.4;
T = 1;
N = round(T/tau);
t = 0:tau:T+1e-12;
s1 = 0.5; s2 = 0.6; s3 = 3;
sgm = [s1 s2];
sgm = [s1];
no = length(sgm);
% [a,b,c,d,w] = coefficients_p2_2(alf,N+1,sgm);
[a,b,c,d,w] = LvXu_sisc_2016(alf,N+1,sgm);

lambda = 1;  lambda0 = 1;
uv = @(t) lambda0 + t.^s1 + t.^s2 + lambda*t.^s3;
fv = @(t) lambda0*t.^(alf)/gamma(1+alf) + t.^(s1+alf)*gamma(s1+1)/gamma(s1+1+alf)....
    + t.^(s2+alf)*gamma(s2+1)/gamma(s2+1+alf)...
    + lambda*gamma(s3+1)/gamma(s3+1+alf)*t.^(s3+alf);
u = uv(t); du = zeros(1,N-1);
uhm0 = u(2:no+1) - u(1);
for k = 1:N-1
    du(k) = d(1)*u(k) + d(2)*u(k+1) + d(3)*u(k+2)....
        + sum(a(1:k).*u(k:-1:1)) + sum(b(1:k).*u(k+1:-1:2))....
        + sum(c(1:k).*u(k+2:-1:3)); 
    du(k) = du(k) + uhm0*w(:,k);
end
du = tau^alf*du;
du2 = fv(t(3:end));
e = du - du2;
ee = norm(e,inf),e(end)
ee;
%END

function [a,b,c,d,w] = coefficients_p2(alf,nN,varargin)
% fractional derivative of the quadratic interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
% Lv and Xu, SISC, 2016
% ini_sisc_LvXu

j = 1:nN;  j1 = j+1;
a = -3/2*(2-alf)*j1.^(1-alf) + 1/2*(2-alf)*j.^(1-alf) + j1.^(2-alf) - j.^(2-alf);
b = 2*(2-alf)*j1.^(1-alf) - 2*j1.^(2-alf) + 2*j.^(2-alf);
c = -1/2*(2-alf)*(j1.^(1-alf) + j.^(1-alf)) + j1.^(2-alf) - j.^(2-alf);
a = a/gamma(3-alf);  b = b/gamma(3-alf); c = c/gamma(3-alf);
d = [alf/2,-2,(4-alf)/2]/gamma(3-alf);
% % % % % % % % % % % % % % % % % % % % % 
N = nN;
no1 = numel(varargin);   w = zeros(1,N);
if no1<1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s-alf)*n(2:N+1).^(s-alf);
        t2(1,1) = d(1,1)*ns(1) + d(2,1)*ns(2) + d(3,1)*ns(3);
        
        t2(1,2:N) = my_conv(ns(1:N-1),a(1:N-1))' + my_conv(ns(2:N),b(1:N-1))'....
            + my_conv(ns(3:N+1),c(1:N-1))' + d(1,2:N)*ns(1)....
            + d(2,2:N)*ns(2) + d(3,2:N)*ns(3);
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
%END

function [a,b,c,d,w] = LvXu_sisc_2016_2(alf,nt,varargin)
% fractional derivative of the quadratic interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
% Lv and Xu, SISC, 2016
j = 0:nt;   alf1 = alf+1; alf2 = alf+2;
fv = @(alf,j) (j+1).^alf - j.^alf;
f0 = fv(alf,j);    f1 = fv(alf1,j);  f2 = fv(alf2,j); 
f0(1) = 1;         f1(1) = 1;        f2(1) = 1;
a = alf*alf1*f2 - alf*alf2*(2*j-1).*f1 + alf1*alf2*j.*(j-1).*f0;
b = alf*alf1*f2 - alf*alf2*(2*j).*f1   + alf1*alf2*(j+1).*(j-1).*f0;
c = alf*alf1*f2 - alf*alf2*(2*j+1).*f1 + alf1*alf2*j.*(j+1).*f0;
a = a/gamma(3+alf)/2;  b = -b/gamma(3+alf); c = c/gamma(3+alf)/2;
a = a(2:end); b = b(2:end); c = c(2:end);
d = [-alf/2,alf*(3+alf),(4+alf)/2]/gamma(3+alf);
% % % % % % % % % % % % % % % % % % % % % 
N = nt;
no1 = numel(varargin);   w = zeros(1,N);
if no1<1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+2;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(3:N+2).^(s+alf);
        tmp = conv(ns(1:N),a(1:N)) + conv(ns(2:N+1),b(1:N))....
            + conv(ns(3:N+2),c(1:N));
        t2 = d(1)*ns(1:N) + d(2)*ns(2:N+1) + d(3)*ns(3:N+2);
        t2 = t2 + tmp(1:N);
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
%END


function [a,b,c,d,w] = LvXu_sisc_2016(alf0,nt,varargin)
% fractional derivative of the quadratic interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
% Lv and Xu, SISC, 2016
tic
alf = alf0;
j = 1:nt;   alf1 = alf+1; alf2 = alf+2;
fv = @(alf,j) (j+1).^alf - j.^alf;
f0 = fv(alf,j);    f1 = fv(alf1,j);  f2 = fv(alf2,j); 
a = alf*alf1*f2 - alf*alf2*(2*j-1).*f1 + alf1*alf2*j.*(j-1).*f0;
b = alf*alf1*f2 - alf*alf2*(2*j).*f1   + alf1*alf2*(j+1).*(j-1).*f0;
c = alf*alf1*f2 - alf*alf2*(2*j+1).*f1 + alf1*alf2*j.*(j+1).*f0;


% y = power_expansion(n,sigma,n0)
j0 = 200; j = j(j0:end);
f0 = power_expansion(j,alf,3);  
f1 = power_expansion(j,alf1,3);   
f2 = power_expansion(j,alf2,3);  
a(j0:end) = alf*alf1*f2 - alf*alf2*(2*j-1).*f1...
    + alf1*alf2*j.*(j-1).*f0 + alf*alf1*alf2/2*j.^(alf-1);
b(j0:end) = alf*alf1*f2 - alf*alf2*(2*j).*f1...
    + alf1*alf2*(j+1).*(j-1).*f0...
    - alf*alf1*alf2*(j.^(alf-1) + (alf-1)/2*j.^(alf-2));
c(j0:end) = alf*alf1*f2 - alf*alf2*(2*j+1).*f1...
    + alf1*alf2*j.*(j+1).*f0 - alf*alf1*alf2/2*j.^(alf-1);

a = a/gamma(3+alf)/2;  b = -b/gamma(3+alf); c = c/gamma(3+alf)/2;
d = [-alf/2,alf*(3+alf),(4+alf)/2]/gamma(3+alf);
% % % % % % % % % % % % % % % % % % % % % 
N = nt;
no1 = numel(varargin);   w = zeros(1,N);
if no1<1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+2;
arr2 = (zeros(no));  bb = (zeros(no,N)); 


tic
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(3:N+2).^(s+alf);
        tmp = conv(ns(1:N),a(1:N)) + conv(ns(2:N+1),b(1:N))....
            + conv(ns(3:N+2),c(1:N));
        
%         tmp = fast_conv(ns(1:N),a(1:N)) + fast_conv(ns(2:N+1),b(1:N))....
%             + fast_conv(ns(3:N+2),c(1:N));
        
        t2 = d(1)*ns(1:N) + d(2)*ns(2:N+1) + d(3)*ns(3:N+2);
        t2 = t2 + tmp(1:N);
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
cputime = toc;
cputime;
%END

function [a,b,c,d,w] = coefficients_p2_2(alf,N,varargin)
% fractional derivative of the quadratic interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
alf1 = alf+1; alf2 = alf+2;
n = 1:N+1; n1 = 0:N;
fv = @(alf,m,n) n.^(m+alf) - (n-1).^(m+alf);
d(1,:) = fv(alf,2,n) - alf2/2*(3*n.^alf1 - n1.^alf1) + alf1*alf2*n.^alf;
d(2,:) = -(2*fv(alf,2,n) - 2*alf2*n.^alf1 + alf1*alf2*n1.^alf);
d(3,:) = fv(alf,2,n) - alf2/2*(n.^alf1 + n1.^alf1);
d(1,1) = 1 - alf2/2*3 + alf1*alf2;
d(2,1) =  -(2 - 2*alf2);  d(3,1) = 1 - alf2/2;
% d = d*0;
% % % % % % % % % % % % % % % % % % % % % % % % % % 
a = fv(alf,2,n) - alf2/2*(n.^alf1 + n1.^alf1);
b = -(2*fv(alf,2,n) - 2*alf2*n1.^alf1 - alf1*alf2*n.^alf);
c = fv(alf,2,n) + alf2/2*n.^alf1 - 3*alf2/2*n1.^alf1 - alf1*alf2*n1.^alf;
a(1) = -alf/2;  b(1) = alf*(3+alf);  c(1) = (4+alf)/2;
a = a/gamma(3+alf);  b = b/gamma(3+alf); 
c = c/gamma(3+alf);  d = d/gamma(3+alf); 
% % % % % % % % % % % % % % % % % % % % % % % 
no1 = numel(varargin);   w = zeros(1,N);
if no1<1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(2:N+1).^(s+alf);
        t2(1,1) = d(1,1)*ns(1) + d(2,1)*ns(2) + d(3,1)*ns(3);
        
        tp = conv(ns(1:N-1),a(1:N-1)) + conv(ns(2:N),b(1:N-1))....
            + conv(ns(3:N+1),c(1:N-1));
        t2(1,2:N) = tp(1:N-1) + d(1,2:N)*ns(1) + d(2,2:N)*ns(2) + d(3,2:N)*ns(3);
        
%         t2(1,2:N) = my_conv(ns(1:N-1),a(1:N-1))' + my_conv(ns(2:N),b(1:N-1))'....
%             + my_conv(ns(3:N+1),c(1:N-1))'...
%             + d(1,2:N)*ns(1) + d(2,2:N)*ns(2) + d(3,2:N)*ns(3);
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
%END


function [a,b,w] = coefficients_p1_2(alf,N,varargin)
% linear interpolation
% D_{t_{j},t}^alf(I_hu)(t) 
n = 0:N;
a = (-(n+1).^(alf+1) + n.^(alf+1) + (1+alf)*(n+1).^alf)/gamma(2+alf); 
b = ((n+1).^(alf+1) - n.^(alf+1) - (1+alf)*n.^alf)/gamma(2+alf);
b(1) = 1/gamma(2+alf); 
a(1) = alf/gamma(2+alf); 
% % % % % % % % % % % % % % % % % % % % % % % 
no1 = numel(varargin);   w = zeros(1,N);
if no1 < 1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(2:N+1).^(s+alf);
        t2 = conv(ns(1:N),a(1:N)) + conv(ns(2:N+1),b(1:N));
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
%END


function [a,w] = coefficients_p1(alf,N,varargin)
% linear interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
n = 0:N;
a = ((n+1).^(alf+1) - n.^(alf+1))/gamma(2+alf);  
n = 200:N;
a(201:end) = pow_n(n,alf+1)/gamma(2+alf); 

% % % % % % % % % % % % % % % % % % % % % % % 
no1 = numel(varargin);   w = zeros(1,N);
if no1 < 1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(2:N+1).^(s+alf);
        t2 = conv(ns(2:N+1)-ns(1:N),a(1:N));
        bb(k,1:N) = t1(1:N) - t2(1:N);
    end
    w = arr2\bb;
end
%END


function [a,w] = coefficients_yan(alf,N,varargin)
% linear interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
n = 3:N;
n2 = (n+2).^(2+alf);  n1 = (n+1).^(2+alf); 
nm1 = (n-1).^(2+alf); n0 = n.^(2+alf);
a = n2 - (1+alf/2)*n2 - 3*n1 + (3+3*alf/2)*n1....
    + 3*n0 - (3+3/2*alf)*n0 - nm1 + (1+alf/2)*nm1;  
a = [(1-alf/2)*2^(1+alf),(2-alf/2)*3^(1+alf)-(3-3*alf/2)*2^(1+alf),.....
    (3-alf/2)*4^(1+alf)-(6-3*alf/2)*3^(1+alf)+(3-3*alf/2)*2^(1+alf),a];
a = a/gamma(3+alf);
% % % % % % % % % % % % % % % % % % % % % % % 
no1 = numel(varargin);   w = zeros(1,N+1);
if no1 < 1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N+1)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n(2:N+2).^(s+alf);
        t2 = conv(ns(1:N+1),a(1:N+1));
        bb(k,1:N+1) = t1(1:N+1) - t2(1:N+1);
    end
    w = arr2\bb;
end
%END

function [a,b,w] = coefficients_p1_1(alf,N,varargin)
% linear interpolation
% D_{t_{j},t}^alf(I_hu)(t) at t = t_{j+1},  [t_{j-1},t_j,t_{j+1}] 
n = 0:N;
a = [1,((n+2).^(alf+1) - 2*(n+1).^(alf+1) + n.^(alf+1))]; 
b = n.^(alf+1) - (n-alf).*(n+1).^alf;
a = a/gamma(2+alf);   b = b/gamma(2+alf);
% % % % % % % % % % % % % % % % % % % % % % % 
no1 = numel(varargin);   w = zeros(1,N);
if no1 < 1
    return;
else
    sgm = varargin{1};
end
no = length(sgm);    n = 0:N+1;
arr2 = (zeros(no));  bb = (zeros(no,N)); 
if ~isempty(sgm)
    for k = 1:no
        s = sgm(k);    arr2(k,:) = (1:no).^s;
        ns = n.^s;
        t1 = gamma(1+s)/gamma(1+s+alf)*n.^(s+alf);
        t2 = conv(ns(2:N+1),a(1:N));
        bb(k,1:N) = t1(2:N+1) - t2(1:N);
    end
    w = arr2\bb;
end
%END


function [w,lambda,lmd1,lmd2,lmd3,lmd4,phi0,phi1,phi2,nT0].....
    = ini_contour_p2_gauss_c(tau,alf,B,nN,N,varargin)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% N = 128*8;  % N can be changed
len = ceil(log(nN/2)/log(B))+2;    
k = 0:len;        nT0 = B.^k*tau;
syms z;
ff2 = ((2-z)*exp(z) - 2 - z)/z^2;
ff3 = ((2-z^2)*exp(z) - 2 - 2*z)/z^2;
ff4 = ((2+z)*exp(z) - 2 - 3*z - 2*z^2)/z^2;
no = numel(varargin);
T0 = 0;
if no>0
    T0 = varargin{1} - tau;
end
if no>1
    N = varargin{2};
end
for kk = 1:len
    x0 = nT0(kk) + T0;
    [lambda{kk},w{kk}] = gen_laguerre_rule(N,-alf,0,x0);
%     id = w{kk} > 1e-20;  
    id = x0^(1-alf)*w{kk} > 1e-16;   
    lambda{kk} = lambda{kk}(id);   w{kk} = w{kk}(id); 
    lambda{kk} = -lambda{kk}';
    w{kk} = sin(alf*pi)/pi*w{kk}';
    z0 = sym(lambda{kk})*sym(tau);
    lmd1{kk} = exp(tau*lambda{kk});    
    lmd2{kk} = double(vpa(subs(ff2,z,z0),60))./lambda{kk}/2;
    lmd3{kk} = -double(vpa(subs(ff3,z,z0),60))./lambda{kk};
    lmd4{kk} = double(vpa(subs(ff4,z,z0),60))./lambda{kk}/2;
end
phi0 = -alf*tau^(alf)/gamma(3+alf)/2;
phi1 = alf*(alf+3)*tau^(alf)/gamma(3+alf);
phi2 = (4+alf)*tau^(alf)/gamma(3+alf)/2;
%END


function  test_compact_fdm
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
xa = -1.5; xb = 1;
p = 2.5;
uv = @(x) sin(p*x);
fv = @(x) -p^2*sin(p*x);
NN = [16 32 64 128]';

for kk = 1:length(NN)
    N = NN(kk);
    h = (xb-xa)/N;
    x = (xa:h:xb+1e-12)';
    
    S = zeros(N-1,N+1);
    for i=1:N-1
        S(i,i:i+2) = [1 -2 1];
    end
    E = eye(N+1);  E(1,:) = []; E(end,:) = [];
    H = S/h^2;  A = E + S/12;
    bnd = [H(:,1),H(:,end)];
    H(:,1) = []; H(:,end) = [];
    f = fv(x);
    ua = uv(xa); ub = uv(xb);
    rhs = A*f - bnd*[ua;ub];
    uu = H\rhs;
    u = uv(x(2:end-1));
    e = abs(uu - u);
    ee(kk,1) = norm(e,inf);
    plot(x(2:end-1),e);
end
fprintf('%-12.4e\n',ee);
rate=-log(ee(2:end)./ee(1:end-1))./(log(NN(2:end)./NN(1:end-1)));
rate


function  y = exp_n(z0)
% (exp(z)-1-z)/z^2
M = 32; N = length(z0);
z = zeros(N,1); z(1:N) = z0;
c = zeros(N,M);
c(:,1) = 1/6;
for n = 0:M-2
    c(:,n+2) = z.*c(:,n+1)/(n+4);
end
c = c(:,end:-1:1);
c1 = c(:,1:2:end); c2 = c(:,2:2:end);
y = sum(c1,2) + sum(c2,2);
% yy = sum(c,2); e = (y - yy)./yy;
if (size(z0,2)>1)
    y = y';
end
%END

function  y = pow_n(n,sigma)
% (n+1)^sigma - n^sigma
m = 20;
w = weights(sigma,m,1);
w(2:2:end) = -w(2:2:end);
w(1) = [];
mm = length(w);
y = 0;
for k = mm:-1:1
    y = y + n.^(-k)*w(k);
end
y = y.*n.^(sigma);
%END

function  y = power_expansion(n,sigma,n0)
% (n+1)^sigma - n^sigma
m = 20;
w = weights(sigma,m,1);
w(2:2:end) = -w(2:2:end);
w(1:n0) = [];
mm = length(w);
y = 0;
for k = mm:-1:1
    y = y + n.^(-k+1)*w(k);
end
y = y.*n.^(sigma-n0);
%END

function  correction_weights
% (exp(z)-1-z)/z^2
alf = -0.1;
m = 5;
nt = 1024*256;
sgm = (1:m)*0.1;
[aa2,bb2,cc2,dd2,ww2] = LvXu_sisc_2016(-alf,nt+1,sgm);
figure;
semilogy(1:length(ww2(1,:)),abs(ww2));

[aa2,bb2,ww2] = coefficients_p1_2(-alf,nt+1,sgm);
figure;
semilogy(1:length(ww2(1,:)),abs(ww2));

function approximation_of_kernel_2
close all;
vps = 1e-10;
B = [2 3 4 5 8 10 15 20];
B = [30 40 50 60 70 80 90 100];
B = [2 3 4 5 8 10 15 20 30 40 50 60 70 80 90 100];

% B = [5 5 5 5];
% nvps = [1e-12 1e-10 1e-8 1e-6];
% B = [5 5];  nvps = [1e-5 1e-4];

for k = 1:length(B)
%     vps = nvps(k);
    [uh,eeu,cputime] = approximation_of_kernel(B(k),vps);
end
%END

function [uh,eeu,cputime] = approximation_of_kernel(BB,vps)
% D_{0,t}^{alpha}u = mu*u + f(u,t)
% Correction terms are applied, Newton method is applied, 
% quatratic inerpolation, gauss quadrature
alf = 0.5;
% dT = 1; 
% BB = 5;
% vps = 1e-10;
T = 1000;
tau = 0.1;
dT = 1;
c0 = 1;
uv = @(t) 1 + c0*t;
duv = @(t) t.^(-alf)/gamma(1-alf) + c0*t.^(1-alf)/gamma(2-alf);
sgm = [];
tic
taualf = tau^alf;        
t = 0:tau:T+1e-10;     nt = length(t)-1;      
uh = zeros(1,nt+1);    duh = uh;
uh0 = uv(0);    uh(1) = uh0;    
% % % % % % % % % % % % % % % 
no = length(sgm);    
dnT = round(dT/tau);  dT = t(dnT+1);
no1 = max(no,dnT-1);  
uh(1:no1+2) = uv(t(1:no1+2)); 
uu = uh(max(no-dnT+2,1):no1+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ initial parameters of fast convolution ------ % %
len = ceil(log(nt/2)/log(BB))+1;  
if abs(len - (ceil(log(nt/2)/log(BB))+1)) < 1e-14
    len = len + 1;
end
[aw,alambda,almd1,almd2,almd3,aphi1,aphi2,nT0,N,idn]....
    = ini_contour_gauss(tau,-alf,BB,nt,dT,vps);
[aa2,bb2] = coefficients_p1_2(-alf,dnT+1,sgm);
aa2 = aa2(dnT:-1:1)';  bb2 = bb2(dnT:-1:1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yz0 = cell(1,len);      
for kk = 1:len
    yz0{kk} = 0;
end
uy2 = yz0; uy3 = yz0;  uy4 = yz0; uyy4 = yz0;
uy33 = yz0; uy44 = yz0;  uy66 = yz0; uy77 = yz0;  uy6 = []; 
% % % % % % % % % % % % % % % % % % % % % % % 
cputime(1) = toc;
tic
d3 = tauL1(BB,2*BB);  u0 = uh(1);      
for n = 1:no-dnT+1
    u1 = uh(n+1);
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;
end
% % % % % % % % % % % % % % % % % % % % % % % % 
for n = max(0,no-dnT+1)+1:2*BB-1
    b = tauL1(BB,n);    btau = tau*b;    
    uh(n+dnT) = uv(t(n+dnT));
    RHS0 =  uu(1:end-1)*aa2 + uu(2:end)*bb2;
    RHS0 = RHS0 + taualf*sum(exp(alambda{1}*(t(n+1)-btau(1)-nT0(1))).*aw{1}.*uy2{1});
    duh(n+dnT) = RHS0;
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uh(n+dnT);
    u1 = uh(n+1);   
    uy2{1} = almd1{1}.*uy2{1} + almd2{1}*u0 + almd3{1}*u1;
    for kk = 1:len
        uyy4{kk} = almd1{kk}.*uyy4{kk} + almd2{kk}*u0 + almd3{kk}*u1;
    end
    if n > d3(2) && n < d3(1)+1
        uy4{1} = almd1{1}.*uy4{1} + almd2{1}*u0 + almd3{1}*u1;
        uy44{1} = uy4{1};
    end
    if n > d3(3) && n < d3(2)+1
        uy4{2} = almd1{2}.*uy4{2} + almd2{2}*u0 + almd3{2}*u1;
    end
    u0 = u1;    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
b2 = zeros(len,2);     b3 = zeros(1,len);   dd3 = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2*BB:nt-dnT+1      
     b = tauL1(BB,n);    btau = tau*b;    L = length(b); 
     uh(n+dnT) = uv(t(n+dnT));   
     [uy2,uy3,uy4,b2,b3,d3,dd3] = update_1(uy2,uy3,uy4,uy6,uy66,b2,b3,d3,dd3,yz0,BB,n);
     RHS0 =  uu(1:end-1)*aa2 + uu(2:end)*bb2;
    for kk = 1:L-1
        RHS0 = RHS0 + taualf*sum(exp(alambda{kk}*(t(n+1)-btau(kk)-nT0(kk))).*aw{kk}.*uy2{kk});
    end
    duh(n+dnT) = RHS0;
    uu(:,1:end-1) = uu(:,2:end);     uu(:,end) = uh(n+dnT);
    u1 = uh(n+1); 
    [uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77] = contour_2......
        (len,uy3,uy33,uy4,uy44,uyy4,uy6,uy66,uy77,b2,b3,d3,dd3,BB,n,almd1,almd2,almd3,u0,u1);
    u0 = u1; 
end
%--------------------------------------
duh = duh/tau^alf;
u = duv(t);
% eeu = abs(duh-u);  eeu(1:dnT+1) = 0;   
eeu = abs((duh-u)./u); eeu(1:dnT+1) = 0;
cputime(2) = toc;
fprintf('alpha = %-s    tau = %-s   dT = %-s   T = %-s   B = %-d\n',...
    num2str(alf),num2str(tau),num2str(dT),num2str(T),BB);
fprintf('max-error = %-10.4e\n',norm(eeu,inf));
fprintf('       N = '); fprintf('%-5d',N); fprintf('\n');
fprintf('kappa(N) = '); fprintf('%-5d',idn); fprintf('\n');
fprintf('sum(N) = '); fprintf('%-5d',sum(N)); fprintf('\n');
fprintf('sum(kappa(N)) = %-5d',sum(idn)); fprintf('\n\n');
figure; loglog(t,eeu);
xlabel('t'); ylabel('log(Error)');
%END

