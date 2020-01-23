%This script takes in jacobians (from jacfolder), 
%control strategies (Bmatrix), penalty matrix for state vector (Q) and 
%penalty matrix (or scaler) for an input(R) and eigenvalues (eigfolder). 
%The code then produces optimum gains, the solution to the Ricatti equation 
%and the closed loop eigenvalues using matlabs in-built discrete linear 
%quadratic regulator. Plots of the ff are also produced 
%a)Absolute closed loop eigenvalues (|eig(A - Bk_{lqr})| vs. BCL)  
%b)Deviations from fixed point (x(n+1) = (A - B*K)*x(n) vrs Iterations)
clear variables
set(0,'defaultlinelinewidth',5)
set(0,'defaultaxesfontsize',23)

selected_bcls = [600:-10:80]; %periods
parameterflag = 1; % 
B             = 4; %

%B matrix (control strategies)
if (B == 1)
    bflag   = 'B_{r}';%'B1'; 
    Bmatrix = [1;0;0;0];
elseif (B == 2)
    bflag   = 'B_{a}';%'B2'; 
    Bmatrix = [0;1;0;0];
elseif (B == 3)
    bflag   = 'B_{b}';%'B3'; 
    Bmatrix = [0;0;1;0];
elseif (B == 4)
    bflag   = 'B_{l}';%'B4'; 
    Bmatrix = [0;0;0;1];
end

 if (parameterflag == 0)
 %   param= 'mixed mechanism 1';
 elseif (parameterflag == 1)
    param = 'voltage-driven alternans';
 elseif (parameterflag == 2)
  %   param = 'calcium-driven alternans';
 elseif (parameterflag == 3)
  %   param = 'calcium-driven alternans';
 elseif (parameterflag == 4)
     param = 'calcium-driven alternans';

 end

jacfolder  = 'jacfolder/'; % folder where jacobians are stored
eigfolder  = 'Eigenvalues/'; %folder where eigenvalues will be saved. 

% alleigs    = cell(1,length(selected_bcls)); % Store eigenvalues here.
alleigsabs = cell(1,length(selected_bcls)); % Store eigenvalue magnitudes here.
sorteigs   = cell(1,length(selected_bcls)); % Store eigenvalue magnitudes here.
% allv       = cell(1,length(selected_bcls)); % Store right eigenvectors here.
% allw       = cell(1,length(selected_bcls)); % Store left eigenvectors here.
% 
% eigw       = cell(length(selected_bcls)); 
% eigv       = cell(length(selected_bcls));
% eigr       = cell(length(selected_bcls));

% Tsave    = cell(length(selected_bcls));
% xn       = cell(1,length(selected_bcls));
% normxn   = cell(1,length(selected_bcls));
%xneig    = cell(1,length(selected_bcls));
%absxneig = cell(1,length(selected_bcls));
lqreig     = cell(1,length(selected_bcls));
Klqr       = cell(1,length(selected_bcls));
normlqreig = cell(1,length(selected_bcls));
normKlqr   = cell(1,length(selected_bcls));
maxabslqreig = cell(1,length(selected_bcls));
normmaxabslqreig = cell(1,length(selected_bcls));
for i = 1:length(selected_bcls)
    
    eval(['load ' jacfolder 'jac' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
    eval(['load ' eigfolder 'alleigs' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from eigenvalue folder
    bcl = selected_bcls(i);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])

%  [sorteigs,sorteigsind] = sort(alleigs{:,i},'descend','ComparisonMethod','abs');
%    eigenvalues = [sorteigs(1); sorteigs(2); sorteigs(3); sorteigs(4)]


    n = size(jac); %this size works because the matrix is an m*m
    I = eye(n); %identity matrix
   %Q = diag(1./varamp.^2); % penalty matrix for state vector
    Q = I; % penalty matrix for state vector
    Qn= norm(Q);
   % penalty matrix (or scalar) for input. R should be more
   %like 1/||umax,i||^2, where i is the index of the input channel
%r = 2^(0); %1;
%r = 2^(-1); %0.5;
%r = 2^(-2); %0.25;
%r = 2^(-3); %0.1250;
%r = 2^(-4); %0.0625;
%r = 2^(-5); %0.0313;
%r = 2^(-6); %0.0157;
%r = 2^(-7); %0.0078;
%r = 2^(-8); %0.0039;
r = 2^(-9);%0.0019;
%r is the ratios of RQ
   R = r*Qn; 
            
   %B = eye(n);

   %[sorteigs,sorteigsind] = sort(alleigs{:,i},'descend','ComparisonMethod','abs');

    % Compute LQR gain (= Klq) system 
    %dlqr for discrete and lqr for continuous system.
    %This function produces the optimum gains, the solution to the Ricatti
    %equation and the closed loop eigenvalues.
    [Klqr{i},lqrsoln,lqreig{i}] = dlqr(jac,Bmatrix,Q,R);
    Klqr{i} = Klqr{i}';
    normKlqr{i} = norm(Klqr{i});
    maxabslqreig{i} = max(abs(lqreig{i}));
    normmaxabslqreig{i} = norm(maxabslqreig{i});
    normlqreig{i} = norm(lqreig{i});
    
end

%Table
%for i=1:length(selected_bcls)
%for m=37; %32:37 %periods not in the 2:1 regime for voltage driven alternans 
%%% for m=34:37 %periods not in the 2:1 regime for voltage driven alternans 
%%%for m=39:44  %periods not in the 2:1 regime for calcium driven alternans   
for m=37:40  %periods not in the 2:1 regime for calcium driven alternans (use this)   
   bcl = selected_bcls(m);
   T1 = table(bcl(:),parameterflag,B, normKlqr{m},maxabslqreig{m},normlqreig{m},normmaxabslqreig{m},...
           'VariableNames',{'BCL','parameterflag','B','norm_gain','maximum_closed_loop_eig','norm_closed_loop_eig','norm_max_abs_lqreig'});
   T2 = table(Bmatrix,Klqr{m},lqreig{m},alleigs{:,m},...
           'VariableNames',{'Bmatrix','lqr_gain','closed_loop_eig','original_eigenvalues'});
end
nsteps = 100;
x0 = [0.1 0.1 0.1 0.1]';
if parameter_flag == 1 % Laura added
    bcl = selected_bcls(37) %choose bcl=240 for voltage-driven alt
elseif parameter_flag == 4
    bcl = selected_bcls(40) %choose bcl=210 for voltage-driven alt
end

[xlqr, lqr_bcl] = lqrsf(x0,jac,Bmatrix,Klqr{m},nsteps,bcl);

%Plots of the absolute closed loop eigenvalues with respect to period
figure()
title(['Eigenvalues of (A - Bk_{lqr}) for ' param ' parameters, ' bflag '']);
ylabel('|eig(A - Bk_{lqr})|');
xlabel('Period (ms)');

grid on;
hold on
for h=1:length(selected_bcls)
%     for m = 1:4
%        scatter(selected_bcls(h), abs(lqreig{h}(m)), 'b*');
%     end

%this eigenvalue plot has absolute
    scatter(selected_bcls(h), abs(lqreig{h}(1)), 'b*','LineWidth',8);%,'LineWidth',2
    scatter(selected_bcls(h), abs(lqreig{h}(2)), 'b*','LineWidth',8);
    scatter(selected_bcls(h), abs(lqreig{h}(3)), 'b*','LineWidth',8);
    scatter(selected_bcls(h), abs(lqreig{h}(4)), 'b*','LineWidth',8);

%plot([x x],[0 y_max])        % Vertical Line
%plot([o x_max],[y y])        % Horizontal line
% Illustrating 2:1 regime
    if (parameterflag == 1)
    plot([selected_bcls(32) selected_bcls(32)],[0 1],'k--') % Vertical Line Voltage driven alternans
    plot([selected_bcls(37) selected_bcls(37)],[0 1],'k--') % Vertical Line Voltage driven alternans
  %  plot([selected_bcls(33) selected_bcls(33)],[0 1],'g--') % Vertical Line Voltage driven alternans
    elseif (parameterflag == 4) 
    plot([selected_bcls(37) selected_bcls(37)],[0 1],'k--') % Vertical Line Calcium driven alternans
%    plot([selected_bcls(40) selected_bcls(40)],[0 1],'r--') % Vertical Line Calcium driven alternans
    plot([selected_bcls(44) selected_bcls(44)],[0 1],'k--') % Vertical Line Calcium driven alternans
    end
    %this eigenvalue plot has no absolute
%     scatter(selected_bcls(h), (lqreig{h}(1)), 'r*');
%     scatter(selected_bcls(h), (lqreig{h}(2)), 'm*');
%     scatter(selected_bcls(h), (lqreig{h}(3)), 'b*');
%     scatter(selected_bcls(h), (lqreig{h}(4)), 'g*');

end
%daspect([1 2 1])
%pbaspect([2 2 1])
%legend('r','a','b','l')
%pbaspect([0.5 0.5 0.5])
set(gcf, 'Units', 'pixels', 'Position', [10, 100, 1000, 400]);


% Plots for lqr state feedback
figure()
title(['Deviations from fixed point for ' param ' parameters, ' bflag '']);
ylabel('xlqr_{n+1} = (A - Bklqr)xlqr_n (lqr design)');
xlabel('Iterations');
grid on;
hold on
for k=1:nsteps
    scatter(k, xlqr(1,k), 'r*','LineWidth',8);
    scatter(k, xlqr(2,k), 'm*','LineWidth',8);
    scatter(k, xlqr(3,k), 'b*','LineWidth',8);
    scatter(k, xlqr(4,k), 'g*','LineWidth',8);
end
legend('r','a','b','l')
set(gcf, 'Units', 'pixels', 'Position', [10, 100, 1000, 400]);
%pbaspect([1 1 1])
%pbaspect([0.5 0.5 0.5])
%daspect([1 1 1])

function [xlqr,bcl] = lqrsf(x0,jac,Bmatrix,Klqr,nsteps,bcl)
xlqr       = zeros(4,length(nsteps));
%xxlqreig    = zeros(4,length(nsteps));
%absxlqreig = zeros(4,length(nsteps));
%normxlqr   = zeros(1,length(nsteps));
xlqr(:,1) = x0;
%xlqreig = eig(jac - Bmatrix*Klqr);  
%maxabsxlqreig = max(abs(xlqreig));
%absxeig = abs(xeig);
%normx = norm(x); %norm at each time step of xn
%maxabsxeig = absxeig; We computed absxeig for purposes of plotting and
%observing the eigenvalues with respect to the iterations.

for j = 1:nsteps
    xlqr(:,j+1)    = jac*xlqr(:,j) - Bmatrix*Klqr'*xlqr(:,j); %x(n+1) = (A - B*K)*x(n)
   % xxlqreig(:,j+1)   = eig(jac - Bmatrix*Klqr);
   % absxlqreig(:,j+1) = abs(xxlqreig(:,j+1));
   % normxlqr(:,j+1)   = norm(xlqr(:,j+1)); %norm at each time step of xn
end
end