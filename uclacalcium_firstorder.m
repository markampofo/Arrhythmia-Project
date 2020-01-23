set(0,'defaultlinelinewidth',3)
set(0,'defaultaxesfontsize',20)
%set(0,'defaultfontsize',20) %will this affect legends? Previous line
%doesn't. 
% UCLA calcium map
% from PRE 75, 011927 (2007)
% variables:
% a: APD
% d: DI
% cp: peak Ca
% r: Ca released; difference between peak Ca and prev ca (end of prev)
% c: ca just before next stimulus
% l: SR calcium load just before next stimulus
% T: period of stimulation

%clear variables


fixedptfolder = 'fixedptfolder/'; % folder where fixed point values will be saved
Jacfolder = 'Jacfolder/'; % folder where jacobians will be saved
Eigfolder = 'Eigfolder/'; % folder where eigenvalues will be saved


% number of iterations
niters = 1000;

% initial values

% % set parameter values
% % Seth's values; vary nu 0.1, 0.28, 0.4, 0.7
gamma = 0.001;% can be 0 or negative
A0 = 220; D0 = 40; tau0 = 30; A1 = 0; D1 = 60; tau1 = 500; % steeper
%%A0 = 220; D0 = 40; tau0 = 60; A1 = 0; D1 = 60; tau1 = 500; % shallower
sigma = 0.5; tauq = 80; 
alpha = 0.036; lc = 93.5; beta = 5;
%rho = 0.15; tauu = 200;
rho=0.5; tauu=200;
nu = 0.1;%orig
%nu = 0.25;
c0h = 50; theta = 20;
kappa = 0.1; % or 0.5, or lower changes alternans in r from period-doubling bif.
eta = 0.1;%original % between -1 and 1
c0c = 28; eps = 2; tauc = 300;

parameterflag = 4; 

if parameterflag == 0
    parametersetname = ('Shallower (Restitution Curve)');
    tau0 = 60;    
elseif parameterflag == 1 %steep parameters
    parametersetname = 'Voltage-Driven Alternans';%'Steep Restitution Curve';% voltage driven alternans
    
elseif parameterflag == 2 %quasiperiodicity 
    parametersetname = 'Quasiperiodicity (Restitution Curve)';
    eta = -0.5;
elseif parameterflag == 3 %And, of course, the biphasic restitution curve wreaks much havoc.
    parametersetname = 'Biphasic (Restitution Curve)';
    tau0 = 60;
    A1 = 60;
elseif parameterflag == 4
    parametersetname = 'Calcium-Driven Alternans';%'('Shallow Restitution, Ca-induced');
    tau0 = 60;
    nu = 0.4;
elseif parameterflag == 5 
    parametersetname = ('Shallow, Ca-induced, Neg Coup');
    tau0 = 60;    
    nu = 0.4;
    gamma = -0.001;
elseif parameterflag == 6 %steep restitution + Ca-induced
    parametersetname = 'Steep + Ca';
    nu = 0.4;
elseif parameterflag == 7
    parametersetname = ('Shallow Restitution, Ca-induced 2');
    tau0 = 60;
    gamma = 0.005;
elseif parameterflag == 8
    parametersetname = ('Shallow Restitution, Ca-induced 3');
    tau0 = 60;
    nu = 0.7;
elseif parameterflag == 9
    parametersetname = ('\nu = 0.35');
    nu=0.35;
elseif parameterflag == 10
    parametersetname = ('\tau_0 = 20');
    tau0=20;

else
    disp('flag not recognized')
end



% period of stimulation
%T = 220;
%Tarray=[500:-50:300 290:-10:200];
Tarray=[600:-10:80];
%Tarray = 200;

nT = length(Tarray);
cumcpsave=zeros(nT,niters+1);
cumasave=zeros(nT,niters+1);
cumbsave=zeros(nT,niters+1);
cumrsave=zeros(nT,niters+1);
cumlsave=zeros(nT,niters+1);
cumdsave=zeros(nT,niters+1);%new

% define functions
f = @(d) A0 * (1 - 1./(1+exp((d-D0)/tau0))) + A1 * exp(-((d-D1).^2)./tau1);
%q = @(d) 1-sigma*exp(-d/tauq);
%g = @(l) l*(1 - (1-alpha)/(1+exp((l-lc)/beta)));
%h = @(cp) nu * cp * (1-1/(1+exp((cp-c0h)/theta)));
%u_func = @(T) (1-rho*exp(-T/tauu));
%c_func = @(T) c0c*(1+eps*exp(-T/tauc));

for k=1:nT
    T=Tarray(k);
    % initial values: trying some 
%     a_old = 250; % APD
%     b_old = 120; % total Ca in cell
%     r_old = 5;  % SR Ca release
%     l_old = 80;  % SR Ca load in
    % initial values: from a no-alt case before alt
    a_old = 197.57; % APD
    b_old = 117.09; % total Ca in cell
    r_old = 2.3762;  % SR Ca release
    l_old = 67.791;  % SR Ca load in

    c_old = b_old - l_old;
    cp_old = c_old + r_old;
    asave=zeros(niters+1,1);
    cpsave=asave;
    bsave=cpsave; rsave=cpsave; lsave=cpsave; dsave = cpsave;
    cpsave(1,1)=cp_old;
    asave(1,1)=a_old;
    bsave(1,1)=b_old;
    rsave(1,1)=r_old;
    lsave(1,1)=l_old;
    %dsave(1,1)=d_old; %new
    
    for j=1:niters
         d_old = T - a_old; %new
         dsave(1,1)=d_old; %new
    
%        c_old = b_old - l_old;
%        cp_new = (b_old - l_old) + r_old;
%        u_new = (1-rho*exp(-T/tauu)) * (nu * ((b_old - l_old) + r_old) * (1-1/(1+exp((((b_old - l_old) + r_old)-c0h)/theta))));
        r_new = (1-sigma*exp(-(T - a_old)/tauq)) * (l_old*(1 - (1-alpha)/(1+exp((l_old-lc)/beta))));
        %    a_new = f(d_old) + p(cp_new) * a_new % will need to rearrange
        a_new = (A0 * (1 - 1./(1+exp(((T - a_old)-D0)/tau0))) + A1 * exp(-(((T - a_old)-D1).^2)./tau1)) / (1 - gamma*((b_old - l_old) + r_new));
        b_new = b_old - kappa*((b_old - l_old) - (c0c*(1+eps*exp(-T/tauc)))) + eta*(((A0 * (1 - 1./(1+exp(((T - a_old)-D0)/tau0))) + A1 * exp(-(((T - a_old)-D1).^2)./tau1)) / (1 - gamma*((b_old - l_old) + r_new))) - a_old);
        l_new = l_old - ((1-sigma*exp(-(T - a_old)/tauq)) * (l_old*(1 - (1-alpha)/(1+exp((l_old-lc)/beta))))) + ((1-rho*exp(-T/tauu)) * (nu * ((b_old - l_old) + r_new) * (1-1/(1+exp((((b_old - l_old) + r_new)-c0h)/theta)))));
        
        % save updated values
        cpsave(j+1,1) = (b_old - l_old) + r_old;
        asave(j+1,1) = a_new;
        bsave(j+1,1) = b_new;
        rsave(j+1,1) = r_new;
        lsave(j+1,1) = l_new;
        dsave(j+1,1) = d_old;%new
        
        % overwrite "old"
        %    cp_old = cp_new;
        a_old = a_new;
        b_old = b_new;
        r_old = r_new;
        l_old = l_new;
    end
    cumcpsave(k,:)=cpsave;
    cumasave(k,:)=asave;
    cumbsave(k,:)=bsave;
    cumlsave(k,:)=lsave;
    cumrsave(k,:)=rsave;
    cumdsave(k,:)=dsave;%new
end
    

iters=0:niters;
set(0,'defaultaxesfontsize',20)


p1 = figure(1);
ib=niters-8;ie=niters+1;
subplot(3,1,1)
plot(Tarray,cumasave(:,ib:ie),'ko')
ylabel('APD')
ylim([floor(min(min(cumasave(:,ib:ie))/10))*10 ceil(max(max(cumasave(:,ib:ie)))/10)*10])
subplot(3,1,2)
plot(Tarray,cumcpsave(:,ib:ie),'ro')
ylabel('Peak Ca')
ylim([floor(min(min(cumcpsave(:,ib:ie))/10))*10 ceil(max(max(cumcpsave(:,ib:ie)))/10)*10])
subplot(3,1,3)
plot(Tarray,cumlsave(:,ib:ie),'bo')
xlabel('Period'),ylabel('SR Ca')
ylim([floor(min(min(cumlsave(:,ib:ie))/10))*10 ceil(max(max(cumlsave(:,ib:ie)))/10)*10])
%saveas(gcf,'fig1cumsave.png')


% Plot of 2:1 regime
mm = 0:600;
nn = 0:600;
figure()
hold on
plot(Tarray,cumasave(:,ib:ie),'ko')
plot(mm,nn)
%ylim([floor(min(min(cumasave(:,ib:ie))/10))*10 ceil(max(max(cumasave(:,ib:ie)))/10)*10])
%xlim([floor(min(min(cumasave(:,ib:ie))/10))*10 ceil(max(max(cumasave(:,ib:ie)))/10)*10])


iiter = 31;
p2 = figure();
subplot(2,1,1)
plot(iters,cumasave(iiter,:),'ko','markerfacecolor','k')
%legend('APD','Total Ca'),legend boxoff
xlabel('Iteration i'),ylabel('a_i')
titlestring = ['T = ' num2str(Tarray(iiter))]
title(titlestring)
xlim([0 50])

subplot(2,1,2)
plot(iters,cumcpsave(iiter,:),'ko','markerfacecolor','k')
%legend('Peak Ca','SR Ca','Ca rel'),legend boxoff
xlabel('Iteration i'),ylabel('c^p_{i}')
title(titlestring)
xlim([0 50])
saveas(gcf,'fig2cumsave.png')
% 
% p3 = figure();
% subplot(2,1,1)
% plot(iters,cumasave(iiter,:),'k.')
% ylabel('APD')
% xlabel('Iteration')%,ylabel('APD')
% titlestring = ['T = ' num2str(Tarray(iiter))]
% title(titlestring)
% 
% subplot(2,1,2)
% plot(iters,cumcpsave(iiter,:),'k.')
% ylabel('Peak Ca')
% xlabel('Iteration')%,ylabel('APD')
% 

%%
%p4 = figure();
% subplot(2,1,1)
% plot(iters(1:end-1),abs(cumasave(iiter,2:end)-cumasave(iiter,1:end-1)),'k.')
% ylabel('\DeltaAPD')
% xlabel('Iteration')%,ylabel('APD')
% titlestring = ['T = ' num2str(T)]
% title(titlestring)
% 
% subplot(2,1,2)
% plot(iters(1:end-1),abs(cumcpsave(iiter,2:end)-cumcpsave(iiter,1:end-1)),'k.')
% ylabel('\Deltac^p')
% xlabel('Iteration')%,ylabel('APD')
%%



set(0,'defaultaxesfontsize',20)
figure()
DIarray = [0:400];
plot(DIarray,f(DIarray),'k-',[0 140],[75 215],'r--')
ylabel('APD (ms)')
xlabel('DI (ms)')
daspect([1 1 1])
ylim([0 250])
title(parametersetname)
%saveas(gcf,'figdisave.png')
%%

% y = [r; a; b; l];


%T = 250;
Tarray= [600:-10:80];

Tlength     = length(Tarray);
fixedptsave = zeros(4,Tlength);

RK1save     = zeros(1,Tlength);
RK2save     = zeros(1,Tlength);
RK3save     = zeros(1,Tlength);
RK4save     = zeros(1,Tlength);
RK5save     = zeros(1,Tlength);
RK6save     = zeros(1,Tlength);
RK7save     = zeros(1,Tlength);
RK8save     = zeros(1,Tlength);
RK9save     = zeros(1,Tlength);
RK10save    = zeros(1,Tlength);
RK11save    = zeros(1,Tlength);
RK12save    = zeros(1,Tlength);
RK13save    = zeros(1,Tlength);
RK14save    = zeros(1,Tlength);
RK15save    = zeros(1,Tlength);

MSV1save    = zeros(1,Tlength);
MSV2save    = zeros(1,Tlength);
MSV3save    = zeros(1,Tlength);
MSV4save    = zeros(1,Tlength);
MSV5save    = zeros(1,Tlength);
MSV6save    = zeros(1,Tlength);
MSV7save    = zeros(1,Tlength);
MSV8save    = zeros(1,Tlength);
MSV9save    = zeros(1,Tlength);
MSV10save   = zeros(1,Tlength);
MSV11save   = zeros(1,Tlength);
MSV12save   = zeros(1,Tlength);
MSV13save   = zeros(1,Tlength);
MSV14save   = zeros(1,Tlength);
MSV15save   = zeros(1,Tlength);

fperrornorm = zeros(1,Tlength); % store fixed-point error norms

fperrormax  = 10^-4; % maxiumum acceptable size for fixed-point error norm
epsln       = 10^-5; % Jacobian step size

% jacfwdsave  = cell(4,Tlength);
% jacbacksave = cell(4,Tlength);
% jacsave     = cell(4,Tlength);  
jacfwdsave  = cell(1,Tlength);
jacbacksave = cell(1,Tlength);
jacsave     = cell(1,Tlength); 

Esave       = zeros(4,Tlength); 
ReEsave     = zeros(4,Tlength);
ImEsave     = zeros(4,Tlength);
EModsave    = zeros(4,Tlength);

Vsave      = cell(Tlength);
alleigsave = cell(Tlength);
alleigsabs = cell(Tlength);
Wsave      = cell(Tlength);

mjjsave    = zeros(4,Tlength);
%nT = length(Tarray);
%Csave = cell(nT);

%R1save(counter1,counter2) = R1;
for counter1 = 1:length(Tarray)
    T = Tarray(counter1)
% Vary additional parameters
% eta = 0:0.1:1;
% nu  = 0:0.1:1;
% 
% etalength = length(eta);
% nulength = length(nu);
% Csave = cell(etalength,nulength);
% R1save = zeros(etalength,nulength);
% R2save = zeros(etalength,nulength);
% R3save = zeros(etalength,nulength);
% R4save = zeros(etalength,nulength);
% R5save = zeros(etalength,nulength);
% fperrornorm = zeros(etalength,nulength); % store fixed-point error norms



% standard iterating function
% ff = @(y) [(1-sigma*exp(-(T - y(2))/tauq)) * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))));
%           (A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / (1 - gamma*((y(3) - y(4)) + y(1)));
%           y(3) - kappa*((y(3) - y(4)) - (c0c*(1+eps*exp(-T/tauc)))) + eta*(((A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / (1 - gamma*((y(3) - y(4)) + y(1)))) - y(2));
%           y(4) - ((1-sigma*exp(-(T - y(2))/tauq)) * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) + ((1-rho*exp(-T/tauu)) * (nu * ((y(3) - y(4)) + y(1)) * (1-1/(1+exp((((y(3) - y(4)) + y(1))-c0h)/theta)))))];

% updated version
ff = @(y) [
(1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))));
(A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) ...
    + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / ...
    (1 - gamma*((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))));
y(3) - kappa*((y(3) - y(4)) ...
    - (c0c*(1+eps*exp(-T/tauc)))) ...
    + eta*(((A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) ...
    + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / ...
    (1 - gamma*((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))))) - y(2));
y(4) - ((1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) ...
    + ((1-rho*exp(-T/tauu)) * (nu * ((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) ...
    * (1-1/(1+exp((((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta)))))-c0h)/theta)))))];
            
% f=y -> f - y for finding fixed points
% gg = @(y) [((1-sigma*exp(-(T - y(2))/tauq)) * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta)))))-y(1);
%           ((A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / (1 - gamma*((y(3) - y(4)) + y(1)))) - y(2);
%           (y(3) - kappa*((y(3) - y(4)) - (c0c*(1+eps*exp(-T/tauc)))) + eta*(((A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / (1 - gamma*((y(3) - y(4)) + y(1)))) - y(2))) - y(3);
%           (y(4) - ((1-sigma*exp(-(T - y(2))/tauq)) * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) + ((1-rho*exp(-T/tauu)) * (nu * ((y(3) - y(4)) + y(1)) * (1-1/(1+exp((((y(3) - y(4)) + y(1))-c0h)/theta)))))) - y(4); ]

% updated version
gg = @(y) [
(1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))));
(A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) ...
    + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / ...
    (1 - gamma*((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))));
y(3) - kappa*((y(3) - y(4)) ...
    - (c0c*(1+eps*exp(-T/tauc)))) ...
    + eta*(((A0 * (1 - 1./(1+exp(((T - y(2))-D0)/tau0))) ...
    + A1 * exp(-(((T - y(2))-D1).^2)./tau1)) / ...
    (1 - gamma*((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))))) - y(2));
y(4) - ((1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) ...
    + ((1-rho*exp(-T/tauu)) * (nu * ((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta))))) ...
    * (1-1/(1+exp((((y(3) - y(4)) + (1-sigma*exp(-(T - y(2))/tauq)) ...
    * (y(4)*(1 - (1-alpha)/(1+exp((y(4)-lc)/beta)))))-c0h)/theta)))))];
      
fsolve(gg,[50;50;50;50])

hh=@(y) ff(y)-y
%hh=@(y) ff(ff(y))-y   % this one will solve alternans states
if(counter1==1)
%    initialguess = [50;50;50;50];
    initialguess = [4;100;150;100];
else
    initialguess = fixedptsave(:,counter1-1);
end
options = optimoptions('fsolve','Display','iter');
%[rr,qq1]=fsolve(hh,initialguess,options)
qq1=fsolve(hh,initialguess);
qq2=qq1;

%[x,fval] = fsolve(@myfun,x0,)

%qq2=ff(qq1) % this would find the other alternans value by iterating once
% ff(qq2)
%qq1-ff(qq1) % check the difference, this should be ~0 (for period-1)


% B Matrix
B1  = [1;1;1;1];
B2  = [1;0;0;0];
B3  = [0;1;0;0];
B4  = [0;0;1;0];
B5  = [0;0;0;1];
B6  = [1;1;0;0];
B7  = [0;0;1;1];
B8  = [0;1;1;0];
B9  = [1;0;1;0];
B10 = [0;1;0;1];
B11 = [1;0;0;1];
B12 = [0;1;1;1];
B13 = [1;0;1;1];
B14 = [1;1;0;1];
B15 = [1;1;1;0];

% fperrornorm(counter1) = norm(ff(y)); %error
%         if fperrornorm(counter1) > fperrormax % stop if error is too large
%             disp(['Fixed point error of ' num2str(fperrornorm(counter1)) ' exceeds threshold of ' num2str(fperrormax) '. Program will exit now.']);
%             return;
%         end
        

fixedpt = [qq2(1); qq2(2); qq2(3); qq2(4)];
fixedptsave(:,counter1) = fixedpt;



 myfilename=[fixedptfolder 'fixedpt' num2str(T) '_pflag' num2str(parameterflag)];
        save([myfilename '.mat'],'fixedpt','T','counter1','parameterflag')%%%



        diffeq_Conly = @(qq2) ff(qq2);
        
        % Forward-difference numerical Jacobian: 
        %jacfwd = diffjac_mod(fixedpt,diffeq_Conly,feval(diffeq_Conly,fixedpt),epsln); % Solve for empirical bcl-to-bcl Jacobian.
        jacfwd = diffjac_mod(fixedptsave(:,counter1),diffeq_Conly,feval(diffeq_Conly,fixedptsave(:,counter1)),epsln); % Solve for empirical bcl-to-bcl Jacobian.
        jacfwdsave{counter1} = jacfwd;
        %jacfwdsave = jacfwd;

        
        % Backward-difference numerical Jacobian: 
        %jacback = diffjac_mod(fixedpt,diffeq_Conly,feval(diffeq_Conly,fixedpt),-epsln); % Solve for empirical bcl-to-bcl Jacobian.
        jacback = diffjac_mod(fixedptsave(:,counter1),diffeq_Conly,feval(diffeq_Conly,fixedptsave(:,counter1)),-epsln); % Solve for empirical bcl-to-bcl Jacobian.
        jacbacksave{counter1} = jacback;
        %jacbacksave = jacback;
        
        % Central-difference Jacobian
        jac               = (jacfwd+jacback)/2;
        %jacsave{counter1} = jac;
        jacsave = jac; 

        myfilename=[Jacfolder 'jac' num2str(T) '_pflag' num2str(parameterflag)];
        save([myfilename '.mat'],'jac','T','counter1','parameterflag')%%%


        
        
        %Eigenvalues
        E                    = sort(eig(jacsave));
        Esave(:,counter1)    = E ;
        Ereal                = real(E);
        ReEsave(:,counter1)  = Ereal ;
        ImE                  = imag(E);
        ImEsave(:,counter1)  = ImE ;
        EMod                 = abs(E);
        EModsave(:,counter1) = EMod;
        
        %Eigenvectors
%        [V, D]              = eig(jac);
        %Calculate the right eigenvectors, V, the eigenvalues, D, and the left eigenvectors, W.
        [V, D, W]            = eig(jac,'nobalance');
        alleigsave{counter1} = diag(D);
        alleigsabs{counter1} = abs(alleigsave{counter1});
        Wsave{counter1}      = W;
        Vsave{counter1}      = V;

        myfilename = [Eigfolder 'alleigsave' num2str(T) '_pflag' num2str(parameterflag)];
        save([myfilename '.mat'],'alleigsave','Wsave','Vsave','alleigsabs','T','counter1','parameterflag')%%%


         % Controllability Matrices
        K1  = [B1 jac*B1 jac^2*B1 jac^3*B1];
        K2  = [B2 jac*B2 jac^2*B2 jac^3*B2];%r
        K3  = [B3 jac*B3 jac^2*B3 jac^3*B3];%a
        K4  = [B4 jac*B4 jac^2*B4 jac^3*B4];%b
        K5  = [B5 jac*B5 jac^2*B5 jac^3*B5];%l
        K6  = [B6 jac*B6 jac^2*B6 jac^3*B6];
        K7  = [B7 jac*B7 jac^2*B7 jac^3*B7];
        K8  = [B8 jac*B8 jac^2*B8 jac^3*B8];
        K9  = [B9 jac*B9 jac^2*B9 jac^3*B9];
        K10 = [B10 jac*B10 jac^2*B10 jac^3*B10];
        K11 = [B11 jac*B11 jac^2*B11 jac^3*B11];
        K12 = [B12 jac*B12 jac^2*B12 jac^3*B12];
        K13 = [B13 jac*B13 jac^2*B13 jac^3*B13];
        K14 = [B14 jac*B14 jac^2*B14 jac^3*B14];
        K15 = [B15 jac*B15 jac^2*B15 jac^3*B15];

        % Minimum of the Singular Values of K
        MSV1                = min(svd(K1'*K1));
        MSV1save(counter1)  = MSV1;
        MSV2                = min(svd(K2'*K2));
        MSV2save(counter1)  = MSV2;
        MSV3                = min(svd(K3'*K3));
        MSV3save(counter1)  = MSV3;
        MSV4                = min(svd(K4'*K4));
        MSV4save(counter1)  = MSV4;
        MSV5                = min(svd(K5'*K5));  
        MSV5save(counter1)  = MSV5;
        MSV6                = min(svd(K6'*K6));
        MSV6save(counter1)  = MSV6;
        MSV7                = min(svd(K7'*K7));
        MSV7save(counter1)  = MSV7;
        MSV8                = min(svd(K8'*K8));
        MSV8save(counter1)  = MSV8;
        MSV9                = min(svd(K9'*K9));
        MSV9save(counter1)  = MSV9;
        MSV10               = min(svd(K10'*K10));
        MSV10save(counter1) = MSV10;
        MSV11               = min(svd(K11'*K11));
        MSV11save(counter1) = MSV11;
        MSV12               = min(svd(K12'*K12));
        MSV12save(counter1) = MSV12;
        MSV13               = min(svd(K13'*K13));
        MSV13save(counter1) = MSV13;
        MSV14               = min(svd(K14'*K14));
        MSV14save(counter1) = MSV14;
        MSV15               = min(svd(K14'*K15));
        MSV15save(counter1)  = MSV15;
        
        % Ranks
        RK1               = rank(K1);
        RK1save(counter1) = RK1;
        RK2               = rank(K2);%r
        RK2save(counter1) = RK2;
        RK3               = rank(K3);%a
        RK3save(counter1) = RK3;
        RK4               = rank(K4);%b
        RK4save(counter1) = RK4;
        RK5               = rank(K5);%l
        RK5save(counter1) = RK5;
        RK6               = rank(K6);
        RK6save(counter1) = RK6;
        RK7               = rank(K7);
        RK7save(counter1) = RK7;
        RK8               = rank(K8);
        RK8save(counter1) = RK8;
        RK9               = rank(K9);
        RK9save(counter1) = RK9;
        RK10              = rank(K10);
        RK10save(counter1)= RK10;
        RK11              = rank(K11);
        RK11save(counter1)= RK11;
        RK12              = rank(K12);
        RK12save(counter1)= RK12;
        RK13              = rank(K13);
        RK13save(counter1)= RK13;
        RK14              = rank(K14);
        RK14save(counter1)= RK14;
        RK15              = rank(K15);
        RK15save(counter1)= RK15;
        
        save('jacsave.mat','jacsave');
        %load('jacsave')
        example = matfile('jacsave.mat');
         A = [     0   -0.0003       0  0.3966;
              0.2387   -0.0001  0.2387 -0.2387;
              0.0239   -0.1000  0.9239  0.0761;
              0.2479    0.0003  0.2479  0.3555]; 
        %A = example.jacsave;
        %
      %  numstep = 51;
        numstep = 53;
        x       = [1;1;1;1];
        for k = 1:numstep-1
          %  x(:,k+1) = jacsave*x(:,k); %correct one
            x(:,k+1) = A*x(:,k);
            %xsave = x
        end
        
end

%%
set(0,'defaultaxesfontsize',20);
set(0,'defaultlinelinewidth',2);

% figure()
% hold on
% for counter1 = 1:length(Tarray)
%     T = Tarray(counter1);
% plot(T,fixedptsave{counter1},'.','markerfacecolor','k')
% xlabel('Tarray'),ylabel('Fixed Points')
% legend('r','a','b','l','Location', 'Best')
% end

% figure()
% hold on
% plot(Tarray,fixedptsave(1,:),'b-o' )
% plot(Tarray,fixedptsave(2,:),'r-x')
% plot(Tarray,fixedptsave(3,:),'g-^')
% plot(Tarray,fixedptsave(4,:),'k-.')
% xlabel('Period','FontSize',14),ylabel('Fixed Points','FontSize',14)
% legend('r','a','b','l','Location', 'Best')
% grid on
% 
% figure()
% hold on
% plot(Tarray,Esave(1,:),'b-o' )
% plot(Tarray,Esave(2,:),'r-x')
% plot(Tarray,Esave(3,:),'g-^')
% plot(Tarray,Esave(4,:),'k-.')
% ylim([-1 1])
% xlabel('Period','FontSize',14),ylabel('Eigenvalues','FontSize',14)
% legend('Eig1','Eig2','Eig3','Eig4','Location', 'Best')
% grid on


% figure()
% plot(Tarray,Esave{counter1}(1,:), 'k.','markerfacecolor','k')
% xlabel('Tarray','FontSize',14),ylabel('Eigenvalues','FontSize',14)
% %legend({'MSV4','MSV5'},'Location', 'Best','FontSize',4)

% figure()
% %plot(Tarray,R1save,'k','markerfacecolor','k')
% plot(Tarray,MSV1save,Tarray,MSV2save,Tarray,MSV3save,Tarray,MSV4save,Tarray,...
%     MSV5save,Tarray,MSV6save,Tarray,MSV7save,Tarray,MSV8save,Tarray,MSV9save,...
%     Tarray,MSV10save,Tarray,MSV11save,Tarray,MSV12save,Tarray,MSV13save,...
%     Tarray,MSV14save,Tarray,MSV15save, 'k.','markerfacecolor','k')
% xlabel('Tarray','FontSize',14),ylabel('Minimum Singular Values','FontSize',14)
% legend({'MSV1','MSV2','MSV3','MSV4','MSV5','MSV6','MSV7','MSV8','MSV9','MSV10',...
%     'MSV11','MSV12','MSV13','MSV14','MSV15'},'Location', 'Best','FontSize',4)

% figure()
% hold on
% %plot(Tarray,R1save,'k','markerfacecolor','k')
% % plot(Tarray,MSV2save,Tarray,MSV3save,Tarray,MSV4save,Tarray,...
% %     MSV5save, 'k.','markerfacecolor','k')
% plot(Tarray,MSV2save,'b-o' )
% plot(Tarray,MSV3save,'r-x')
% plot(Tarray,MSV4save,'g-^')
% plot(Tarray,MSV5save,'k-.')
% set(gca, 'yscale', 'log')
% xlabel('Period','FontSize',14),ylabel('Minimum Singular Values','FontSize',14)
% legend({'MSV2','MSV3','MSV4','MSV5'},'Location', 'Best','FontSize',4)
% grid on
% 
% figure()
% subplot(3,1,1)
% hold on
% plot(Tarray,fixedptsave(1,:),'b-o' )
% plot(Tarray,fixedptsave(2,:),'r-x')
% plot(Tarray,fixedptsave(3,:),'g-^')
% plot(Tarray,fixedptsave(4,:),'k-.')
% %xlabel('Tarray','FontSize',14)
% ylabel('Fixed Points','FontSize',10)
% legend('r','a','b','l','Location', 'Best')
% grid on
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(Tarray,Esave(1,:),'b-o' )
% plot(Tarray,Esave(2,:),'r-x')
% plot(Tarray,Esave(3,:),'g-^')
% plot(Tarray,Esave(4,:),'k-.')
% ylim([-1 1])
% %xlabel('Tarray','FontSize',14),
% ylabel('Eigenvalues','FontSize',10)
% legend('Eig1','Eig2','Eig3','Eig4','Location', 'Best')
% grid on
% hold off

% subplot(3,1,3)
% hold on
% plot(Tarray,MSV2save,'b-o' )
% plot(Tarray,MSV3save,'r-x')
% plot(Tarray,MSV4save,'g-^')
% plot(Tarray,MSV5save,'k-.')
% set(gca, 'yscale', 'log')
% xlabel('Period','FontSize',14)
% ylabel('Minimum Singular Values','FontSize',10)
% legend({'MSV2','MSV3','MSV4','MSV5'},'Location', 'Best')
% grid on
%hold off

%figure()
figure1 = figure('Position', [100, 10, 1024, 1200]);
ax1 = subplot(3,1,1);
hold on
plot(Tarray,fixedptsave(2,:),'ko','markerfacecolor','k')
plot(Tarray,fixedptsave(4,:),'rd','markerfacecolor','r')
plot(Tarray,fixedptsave(1,:),'bx' )
plot(Tarray,fixedptsave(3,:),'g^','markerfacecolor','g')
%xlabel('Tarray','FontSize',14)
ylabel('Fixed Points')
%%%legend('r','a','b','l','Location', 'Best')
%set(gca,'xtickMode', 'auto')
%set(gca,'XTick',[]);
%grid on
set(gca,'XTicklabel',[]);
%title(parametersetname)
hold off

% ax2=subplot(4,1,2);
% hold on
% plot(Tarray,ReEsave(1,:),'bo' )
% plot(Tarray,ReEsave(2,:),'bo')
% plot(Tarray,ReEsave(3,:),'bo')
% plot(Tarray,ReEsave(4,:),'bo')
% %ylim([-1 1])
% %xlabel('Tarray','FontSize',14),
% ylabel('Re(\lambda_{i})')
% %legend('Re1','Re2','Re3','Re4','Location', 'Best')
% grid on
% set(gca,'XTicklabel',[]);
% hold off

% subplot(3,1,2)
% hold on
% plot(Tarray,ImEsave(1,:),'b-o' )
% plot(Tarray,ImEsave(2,:),'r-x')
% plot(Tarray,ImEsave(3,:),'g-^')
% plot(Tarray,ImEsave(4,:),'k-.')
% set(gca, 'yscale', 'log')
% xlabel('Period','FontSize',14)
% ylabel('Im(Eigenvalue)','FontSize',10)
% legend({'Im1','Im2','Im3','Im4'},'Location', 'Best')
% grid on
% hold off

% ax3=subplot(4,1,3);
% hold on
% plot(Tarray,EModsave(1,:),'bo' )
% plot(Tarray,EModsave(2,:),'bo')
% plot(Tarray,EModsave(3,:),'bo')
% plot(Tarray,EModsave(4,:),'bo')
% %set(gca, 'yscale', 'log')
% %xlabel('Period','FontSize',14)
% ylabel('|\lambda_{i}|')
% %legend({'Mod1','Mod2','Mod3','Mod4'},'Location', 'Best')
% grid on
% set(gca,'XTicklabel',[]);
% %set(ax3,'xgrid','on');
% % ax3.XGrid = 'on';
% % ax3.YGrid = 'on';
% hold off


ax3 = subplot(3,1,2);
hold on
posreal1=find(ReEsave(1,:)>=0);
posreal2=find(ReEsave(2,:)>=0);
posreal3=find(ReEsave(3,:)>=0);
posreal4=find(ReEsave(4,:)>=0);
negreal1=find(ReEsave(1,:)<0);
negreal2=find(ReEsave(2,:)<0);
negreal3=find(ReEsave(3,:)<0);
negreal4=find(ReEsave(4,:)<0);
if (size(posreal1,2)>1),plot(Tarray(posreal1),EModsave(1,posreal1),'ko','markerfacecolor','k'),end
if (size(posreal2,2)>1),plot(Tarray(posreal2),EModsave(2,posreal2),'ko','markerfacecolor','k'),end
if (size(posreal3,2)>1),plot(Tarray(posreal3),EModsave(3,posreal3),'ko','markerfacecolor','k'),end
if (size(posreal4,2)>1),plot(Tarray(posreal4),EModsave(4,posreal4),'ko','markerfacecolor','k'),end
if (size(negreal1,2)>1),plot(Tarray(negreal1),EModsave(1,negreal1),'mo','markerfacecolor','r'),end
if (size(negreal2,2)>1),plot(Tarray(negreal2),EModsave(2,negreal2),'mo','markerfacecolor','r'),end
if (size(negreal3,2)>1),plot(Tarray(negreal3),EModsave(3,negreal3),'mo','markerfacecolor','r'),end
if (size(negreal4,2)>1),plot(Tarray(negreal4),EModsave(4,negreal4),'mo','markerfacecolor','r'),end
plot([ 0 600],[1 1],'k--')
%set(gca, 'yscale', 'log')
xlabel('Period (ms)')
ylabel('|\lambda_{i}|')
%legend({'Mod1','Mod2','Mod3','Mod4'},'Location', 'Best')
%grid on
%set(gca,'XTicklabel',[]);
%set(ax3,'xgrid','on');
% ax3.XGrid = 'on';
% ax3.YGrid = 'on';
hold off


ax4 = subplot(3,1,3);
hold on
semilogy(Tarray,MSV2save,'b-o' )
semilogy(Tarray,MSV3save,'r-x')
semilogy(Tarray,MSV4save,'g-^')
semilogy(Tarray,MSV5save,'k-.')
% plot(Tarray,MSV2save,'b-o' )
% plot(Tarray,MSV3save,'r-x')
% plot(Tarray,MSV4save,'g-^')
% plot(Tarray,MSV5save,'k-.')
% set(gca, 'yscale', 'log')
xlabel('Period (ms)','FontSize',14)
ylabel('\sigma_{min}')
%legend({'\sigma_{min,r}','\sigma_{min,a}','\sigma_{min,b}','\sigma_{min,l}'},'Location', 'Best')
legend({'r','a','b','l'},'Location', 'Best')
legend boxoff
hold off
linkaxes([ax1,ax3,ax4],'x')
set(gca,'XTickLabel');
%saveas(gcf,'figfxdpt,mod,eig.png')
%set(gcf,'position',[1.0000    1.0000  500  1000.0000])%

set(0,'defaultaxesfontsize',20)
figure(10)
mygreen=[0 0.7 0];
myxlims = [80 540];
set(gcf,'Position', [100, 10, 1024, 1200]);
ib=niters-8;ie=niters+1;
subplot(3,1,1)
plot(Tarray,cumasave(:,ib:ie),'ko','markerfacecolor','k')
ylabel('a_i (ms)')
xlim(myxlims)
ylim([floor(min(min(cumasave(:,ib:ie))/10))*10 ceil(max(max(cumasave(:,ib:ie)))/10)*10])
subplot(3,1,2)
plot(Tarray,cumcpsave(:,ib:ie),'ko','markerfacecolor','k')
ylabel('Ca_i^p (\muM)')
xlim(myxlims)
ylim([floor(min(min(cumcpsave(:,ib:ie))/10))*10 ceil(max(max(cumcpsave(:,ib:ie)))/10)*10])
subplot(3,1,3)
semilogy(Tarray,MSV3save,'ko','markerfacecolor','k')
hold on
semilogy(Tarray,MSV5save,'rd','markerfacecolor','r')
semilogy(Tarray,MSV2save,'bx','markerfacecolor','b')
semilogy(Tarray,MSV4save,'^','color',mygreen,'markerfacecolor',mygreen)
hold off
xlabel('Period (ms)','FontSize',14)
ylabel('\sigma_{min}')
xlim(myxlims)
ylim([1e-10 1e0])
%ylim([1e-9 1e0])
%legend({'\sigma_{min,r}','\sigma_{min,a}','\sigma_{min,b}','\sigma_{min,l}'},'Location', 'Best')
%legend({'a','l','r','b'},'Location', 'southoutside')
legend boxoff
%saveas(gcf,'figminsin.png')

set(0,'defaultaxesfontsize',20)
figure(6)
%plot(Tarray,R1save,'k','markerfacecolor','k')
plot(Tarray,RK1save,Tarray,RK2save,Tarray,RK3save,Tarray,RK4save,Tarray,...
    RK5save,Tarray,RK6save,Tarray,RK7save,Tarray,RK8save,Tarray,RK9save,...
    Tarray,RK10save,Tarray,RK11save,Tarray,RK12save,Tarray,RK13save,...
    Tarray,RK14save,Tarray,RK15save, 'k.','markerfacecolor','k')
xlabel('Period'),ylabel('Ranks')
legend({'RK1','RK2','RK3','RK4','RK5','RK6','RK7','RK8','RK9','RK10',...
    'RK11','RK12','RK13','RK14','RK15'},'Location', 'eastoutside')
grid on
%saveas(gcf,'figranks.png')


set(0,'defaultaxesfontsize',32)
set(0,'defaultlinelinewidth',3)
figure(7)
%MarkerSize = 80;
% ax = gca;
% ax.FontSize = 20;
plot(Tarray,RK2save,'b-o','markerfacecolor','b','MarkerSize',8)
hold on
plot(Tarray,RK3save,'r-o','markerfacecolor','r','MarkerSize',8)
plot(Tarray,RK4save,'g-x','markerfacecolor','g','MarkerSize',8)
plot(Tarray,RK5save,'k-.','markerfacecolor','k','MarkerSize',8)
xlabel('Period'),ylabel('Ranks')
legend({'r','a','b','l'},'Location', 'eastoutside')
%title(['Ranks for ' parametersetname, ], 'FontSize', 15);
title(['Ranks for ' parametersetname, ]);
grid on
hold off

set(0,'defaultaxesfontsize',20)
figure(8)
plot(Tarray,RK3save, 'r.','markerfacecolor','r')
xlabel('Period'),ylabel('Ranks')
legend({'RK3'},'Location', 'eastoutside')
grid on

set(0,'defaultaxesfontsize',20)
figure(9)
plot(Tarray,RK3save, 'k.','markerfacecolor','k')
xlabel('Period'),ylabel('Ranks')
legend({'RK3'},'Location', 'eastoutside')
grid on

%figure()
%hold on
%plot(Tarray, Vsave(:,1:end))
%cellfun(@plot,Vsave);
%plot(Tarray, Vsave(2,:))
%plot(Tarray, Vsave(3,:))
%plot(Tarray, Vsave(4,:))

%plot(Tarray,Vsave(1,:),'b-o' )
%plot(Tarray,Vsave(2,:),'r-x')
%plot(Tarray,Vsave(3,:),'g-^')
%plot(Tarray,Vsave(4,:),'k-.')
%ylim([-1 1])
%%legend('Eig1','Eig2','Eig3','Eig4','Location', 'Best')
grid on

% figure()
% plot(Tarray,R3save,Tarray,R4save,Tarray,R5save,'k','markerfacecolor','k')
% xlabel('Tarray'),ylabel('Ranks (R2)')


% figure()
% surf(Tarray,R1save(:,Tlength),R2save(:,Tlength));%shading interp
% xlabel('Tarray'),ylabel('Ranks (R1)'),zlabel('Ranks (R2)')
% 
% figure()
% surf(Tarray,R2save(:,Tlength),R1save(:,Tlength));%shading interp
% xlabel('Tarray'),ylabel('Ranks (R2)'),zlabel('Ranks (R1)')


% figure()
% surf(b1values,b2values,R2save);%shading interp
% xlabel('b1 Values'),ylabel('b2 Values'),zlabel('Ranks (R2)')
% 
% figure()
% surf(b1values,b2values,R3save);%shading interp
% xlabel('b1 Values'),ylabel('b2 Values'),zlabel('Ranks (R3)')
% 
% b2plotindex = 1; % Fix b2 value and plot ranks vs. b1
% figure()
% %plot(b1values,R1save,'m-s',b1values,R2save,'r-^',b1values,R3save,'b-o')
% plot(b1values,R1save(:,b2plotindex),'m-s',b1values,R2save(:,b2plotindex),'r-^',b1values,R3save(:,b2plotindex),'b-o')
% xlabel('b1 values '),ylabel('Rank Value')%,legend('R1save','R2save','R3save','location','best'),legend boxoff
% title(['Controllability ranks for b2 = ' num2str(b2values(b2plotindex))])
% legend('R1','R2','R3')
% T = 600;
% figure()
% hold on
% plot((1:numstep)*T, x(1,:),'bo')
% plot((1:numstep)*T, x(2,:),'ro')
% plot((1:numstep)*T, x(3,:),'go')
% plot((1:numstep)*T, x(4,:),'ko')
% xlabel('Time (ms)'),ylabel('x(k+1)') %miliseconds
% legend('r','a','b','l','Location', 'Best')

% figure()
% hold on
% plot((1:numstep)*600, x(1,:),'b')
% plot((1:numstep)*600, x(2,:),'r')
% plot((1:numstep)*600, x(3,:),'g')
% plot((1:numstep)*600, x(4,:),'k')
% xlabel('Time (ms)'),ylabel('x(k+1)') %miliseconds
% legend('r','a','b','l','Location', 'Best')
set(0,'defaultaxesfontsize',20)%%

%figure()
set(0,'defaultaxesfontsize',20)
figure1=figure('Position', [100, 10, 1024, 1200]);
ax1=subplot(4,1,1)
hold on
plot(Tarray,fixedptsave(1,:),'b-o' )
plot(Tarray,fixedptsave(2,:),'r-x')
plot(Tarray,fixedptsave(3,:),'g-^')
plot(Tarray,fixedptsave(4,:),'k-.')
%xlabel('Tarray','FontSize',14)
ylabel('Fixed Points','FontSize',10)
legend('r','a','b','l','Location', 'Best')
%set(gca,'xtickMode', 'auto')
%set(gca,'XTick',[]);
grid on
set(gca,'XTicklabel',[]);
hold off


ax2=subplot(4,1,2);
hold on
plot(Tarray,ReEsave(1,:),'bo' )
plot(Tarray,ReEsave(2,:),'bo')
plot(Tarray,ReEsave(3,:),'bo')
plot(Tarray,ReEsave(4,:),'bo')
%ylim([-1 1])
%xlabel('Tarray','FontSize',14),
ylabel('Re(\lambda_{i})','FontSize',10)
%legend('Re1','Re2','Re3','Re4','Location', 'Best')
grid on
set(gca,'XTicklabel',[]);
hold off

% subplot(3,1,2)
% hold on
% plot(Tarray,ImEsave(1,:),'b-o' )
% plot(Tarray,ImEsave(2,:),'r-x')
% plot(Tarray,ImEsave(3,:),'g-^')
% plot(Tarray,ImEsave(4,:),'k-.')
% set(gca, 'yscale', 'log')
% xlabel('Period','FontSize',14)
% ylabel('Im(Eigenvalue)','FontSize',10)
% legend({'Im1','Im2','Im3','Im4'},'Location', 'Best')
% grid on
% hold off

ax3=subplot(3,1,2);
hold on
posreal1=find(ReEsave(1,:)>=0);
posreal2=find(ReEsave(2,:)>=0);
posreal3=find(ReEsave(3,:)>=0);
posreal4=find(ReEsave(4,:)>=0);
negreal1=find(ReEsave(1,:)<0);
negreal2=find(ReEsave(2,:)<0);
negreal3=find(ReEsave(3,:)<0);
negreal4=find(ReEsave(4,:)<0);
if (size(posreal1,2)>1),plot(Tarray(posreal1),EModsave(1,posreal1),'ko','markerfacecolor','k'),end
if (size(posreal2,2)>1),plot(Tarray(posreal2),EModsave(2,posreal2),'ko','markerfacecolor','k'),end
if (size(posreal3,2)>1),plot(Tarray(posreal3),EModsave(3,posreal3),'ko','markerfacecolor','k'),end
if (size(posreal4,2)>1),plot(Tarray(posreal4),EModsave(4,posreal4),'ko','markerfacecolor','k'),end
if (size(negreal1,2)>1),plot(Tarray(negreal1),EModsave(1,negreal1),'mo','markerfacecolor','r'),end
if (size(negreal2,2)>1),plot(Tarray(negreal2),EModsave(2,negreal2),'mo','markerfacecolor','r'),end
if (size(negreal3,2)>1),plot(Tarray(negreal3),EModsave(3,negreal3),'mo','markerfacecolor','r'),end
if (size(negreal4,2)>1),plot(Tarray(negreal4),EModsave(4,negreal4),'mo','markerfacecolor','r'),end
plot([ 0 600],[1 1],'k--')
%set(gca, 'yscale', 'log')
xlabel('Period (ms)')
ylabel('|\lambda_{i}|')
%legend({'Mod1','Mod2','Mod3','Mod4'},'Location', 'Best')
%grid on
%set(gca,'XTicklabel',[]);
%set(ax3,'xgrid','on');
% ax3.XGrid = 'on';
% ax3.YGrid = 'on';
hold off

ax4=subplot(4,1,4);
hold on
plot(Tarray,MSV2save,'b-o' )
plot(Tarray,MSV3save,'r-x')
plot(Tarray,MSV4save,'g-^')
plot(Tarray,MSV5save,'k-.')
set(gca, 'yscale', 'log')
xlabel('Period','FontSize',14)
ylabel('\sigma_{min}','FontSize',10)
%legend({'\sigma_{min,r}','\sigma_{min,a}','\sigma_{min,b}','\sigma_{min,l}'},'Location', 'Best')
legend({'r','a','b','l'},'Location', 'Best')
grid on
hold off
linkaxes([ax1,ax2,ax3,ax4],'x')
set(gca,'XTickLabel');

set(0,'defaultaxesfontsize',32)
set(0,'defaultlinelinewidth',3)
figure()
%figure1=figure('Position', [100, 10, 1024, 1200]);
%ax1=subplot(4,1,1)
hold on
plot(Tarray,fixedptsave(1,:),'b-o' )
plot(Tarray,fixedptsave(2,:),'r-x')
plot(Tarray,fixedptsave(3,:),'g-^')
plot(Tarray,fixedptsave(4,:),'k-.')
%xlabel('Tarray','FontSize',14)
%ylabel('Fixed Points','FontSize',18)
ylabel('Fixed Points')
legend('r','a','b','l','Location', 'Best')
%set(gca,'xtickMode', 'auto')
%set(gca,'XTick',[]);
%plot([ 0 600],[1 1],'k--')
xlabel('Period (ms)')
grid on
%set(gca,'XTicklabel',[]);
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontsize',24)
set(0,'defaultlinelinewidth',3)
figure()
mygreen=[0 0.7 0];
myxlims = [80 540];
set(gcf,'Position', [100, 10, 1024, 1200]);
ib=niters-8;ie=niters+1;
ax1=subplot(3,1,1);
plot(Tarray,cumasave(:,ib:ie),'ko','markerfacecolor','k')
ylabel('a_i (ms)')
xlim(myxlims)
ylim([floor(min(min(cumasave(:,ib:ie))/10))*10 ceil(max(max(cumasave(:,ib:ie)))/10)*10])
ax2=subplot(3,1,2);
plot(Tarray,cumcpsave(:,ib:ie),'ro','markerfacecolor','r')
ylabel('Ca_i^p (\muM)')
xlim(myxlims)
ylim([floor(min(min(cumcpsave(:,ib:ie))/10))*10 ceil(max(max(cumcpsave(:,ib:ie)))/10)*10])
hold off
ax3=subplot(3,1,3);
hold on
plot(Tarray,MSV2save,'b-o' )
plot(Tarray,MSV3save,'r-x')
plot(Tarray,MSV4save,'g-^')
plot(Tarray,MSV5save,'k-.')
set(gca, 'yscale', 'log')
xlabel('Period')
ylabel('\sigma_{min}')
h = legend({'r','a','b','l'},'Location', 'Best')
legend boxoff
set(h,'FontSize',18);
grid on
hold off
ylim([1e-10 1e0])
linkaxes([ax1,ax2,ax3],'x')
set(gca,'XTickLabel');
