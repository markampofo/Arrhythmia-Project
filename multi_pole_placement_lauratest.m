% "lauratest" version: I'm just testing out some scenarios to check
% numbers, e.g. what happens if we change the position of the nonzero Ksc entry, 
% to compare values with those in the thesis presentation. 
% I also used this to produce some Kpp gains for cases where the largest
% eigenvalue was moved to 0.95. 
% I found a bug in the code: alleigs{} indexing doesn't work when you just
% load one Jacobian. eigenvalues should be recomputed. 
% Laura: an improvement to the code would be to automatically switch
% nonzero Ksc entry to the correct location based on B=B1, B=B2, etc.,
% instead of manually editing the vector. 
% 11/04/2018
%This code compares the multi-variable and single variable state feedback. 
%It uses the matlab place function to compute for optimal gains for 
%a multi-variable state feedback. It obtains eigenvalues and jacobians from
%the jacfolder and eigfolder. Eigenvalues greater than one are selected 
%using a sort matlab inbuilt function and replaced with a desired 
%eigenvalue. This is used by the place function in matlab which
%then obtains an optimal gain to be used for a multi-variable state 
%feedback. The gains for the single variable state feedback are actually 
%chosen by myself using Ksc. We also plot the deviations from the fixed 
%points for both the multi-variable and state variable state feedback. It 
%also produces the ff plots  
%a)Absolute closed loop eigenvalues for the single variable state feedback
%(|eig(A - Bksv)| vs. BCL)  
%b)Deviations from fixed point for the single variable state feedback
%((xsv_{n+1} = (A - Bksv)xsv_n) vrs Iterations)
%c)Norm of deviations from fixed point for the single variable state 
%feedback (Norm(xsv_{n+1}) vs Iterations);
%d)Absolute closed loop eigenvalues for the multi variable state feedback
%(|eig(A - Bk)| vs. BCL)  
%e)Deviations from fixed point for the multi variable state feedback
%((x_{n+1} = (A - Bk)x_n) vrs Iterations)
%f)Norm of deviations from fixed point for the multi variable state 
%feedback (Norm(X_{n+1}) vs Iterations);


clear variables

%selected_bcls = [600:-10:80];
%selected_bcls = [600:-10:270];
selected_bcls = 240;%240;
parameterflag = 1 % 1=voltage-driven alternans,4=calcium-driven alt
B             = 2

%B matrix (control strategies)
if (B == 1)
    bflag   = 'B1'; 
    Bmatrix = [1;0;0;0];
elseif (B == 2)
    bflag   = 'B2'; 
    Bmatrix = [0;1;0;0];
elseif (B == 3)
    bflag   = 'B3'; 
    Bmatrix = [0;0;1;0];
elseif (B == 4)
    bflag   = 'B4'; 
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

Tsave    = cell(length(selected_bcls));
xn       = cell(1,length(selected_bcls));
normxn   = cell(1,length(selected_bcls));
xneig    = cell(1,length(selected_bcls));
absxneig = cell(1,length(selected_bcls));
for i = 1:length(selected_bcls)
    
    eval(['load ' jacfolder 'jac' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
%    eval(['load ' eigfolder 'alleigs' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from eigenvalue folder
    bcl = selected_bcls(i);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    [v,d,w] = eig(jac);     
    % Test: w'A = lambda w'
 %   w = allw{i}; %left eigenvectors
%    G = (w'*jac) - (eigv{i}*w'); %w'A - Lamda*w' = 0
        
    % sorted eigenvalues
%    [sorteigs,sorteigsind] = sort(alleigs{:,i},'descend','ComparisonMethod','abs');
    [sorteigs,sorteigsind] = sort(eig(jac),'descend','ComparisonMethod','abs');
    %sorted left eigenvectors
%    sortedlefteigs  = allw{:,i}(:, sorteigsind);
    sortedlefteigs = w(:,sorteigsind); 
    sortedlefteigs1 = sortedlefteigs(:,1);
    sortedlefteigs2 = sortedlefteigs(:,2);
    sortedlefteigs3 = sortedlefteigs(:,3);
    sortedlefteigs4 = sortedlefteigs(:,4);
    
    C1 = abs(sortedlefteigs1'*Bmatrix)/(norm(sortedlefteigs1)*norm(Bmatrix));
    C2 = abs(sortedlefteigs2'*Bmatrix)/(norm(sortedlefteigs2)*norm(Bmatrix));
    C3 = abs(sortedlefteigs3'*Bmatrix)/(norm(sortedlefteigs3)*norm(Bmatrix));
    C4 = abs(sortedlefteigs4'*Bmatrix)/(norm(sortedlefteigs4)*norm(Bmatrix));
    C  = [C1; C2; C3; C4];
    Cnorm = norm(C);
    
    eigenvalues = [sorteigs(1); sorteigs(2); sorteigs(3); sorteigs(4)];
    
    % replacing eigenvalues not in the unit circle
    % if (sorteigs{:,i}(1) > .90)
    %sorteigs(1) = (sorteigs(1))/2;

    %desired eigenvalues in the unit circle
%    desiredeigenvalues = [sorteigs(1); sorteigs(2); sorteigs(3); sorteigs(4)];
    desiredeigenvalues = [0.95; sorteigs(2); sorteigs(3); sorteigs(4)]
   
    if (desiredeigenvalues(:,:) == eigenvalues(:,:))
        disp(['Eigenvalues in unit circle for BCL = ' num2str(bcl) ' ms']) 
    else 
        [Kpl,placeprecision,placemessage] = place(jac,Bmatrix,desiredeigenvalues);
        Kplnorm = norm(Kpl)
               
        %Table
%         T1 = table(bcl(:),parameterflag,B,placeprecision, Kplnorm, Cnorm,...
%             'VariableNames',{'BCL','parameterflag','B','placeprecision','Norm_of_Kpl','Norm_of_C'})
%         T2 = table(Bmatrix,eigenvalues(:,:),desiredeigenvalues(:,:),Kpl',C,...
%             'VariableNames',{'Bmatrix','eigenvalues','desired_eigenvalues','Kpl','C'})
             
       


    end
end

%% functions for single and multiple variable state feedback 
nsteps = 100;
x0 = [0.1 0.1 0.1 0.1]';
%Ksc = 2.10;
%Ksc = -2.35;
Ksc = 1.10;
[xsv,xsveig,maxabsxsveig,normxsv,absxsveig] = singlevarsf(x0,jac,Bmatrix,Ksc,nsteps); %single variable state feedback
maxabsxsveig
normallxsv = norm(xsv);
normallxsveig = norm(xsveig);
[x,xeig,maxabsxeig,normx,absxeig] = multvarsf(x0,jac,Bmatrix,Kpl,nsteps); %multi-variable state feedback
maxabsxeig
%maxabsxeig = absxeig; We computed absxeig for purposes of plotting and
%observing the eigenvalues with respect to the iterations.
normallx = norm(x);
normallxeig = norm(xeig);

%% Plots for single variable state feedback
% figure()
% title(['Deviations from fixed point for ' param ' parameters, ' bflag '']);
% ylabel('xsv_{n+1} = (A - Bksv)xsv_n (single varaiable)');
% xlabel('Iterations');
% grid on;
% hold on
% for k=1:nsteps
%     scatter(k, xsv(1,k), 'r*');
%     scatter(k, xsv(2,k), 'm*');
%     scatter(k, xsv(3,k), 'b*');
%     scatter(k, xsv(4,k), 'g*');
% end
% legend('r','a','b','l')
% 
% figure()
% title(['Norm of deviations from fixed point for ' param ' parameters, ' bflag '']);
% ylabel('Norm(xsv_{n+1})');
% xlabel('Iterations');
% grid on;
% hold on
% for k=1:nsteps
%    scatter(k,normxsv(k),'b*')
% end
%  
% % figure()
% % title(['Eigenvalues of (A - Bksv) for ' param ' parameters, ' bflag '']);
% % ylabel('eig(A - Bksv)');
% % xlabel('Iterations');
% % grid on;
% % hold on
% % for h=1:nsteps
% % %     for m = 1:4
% % %     scatter(h, xeig(m,h), 'b*');
% % %     end  
% %     scatter(h, xsveig(1,h), 'r*');
% %     scatter(h, xsveig(2,h), 'm*');
% %     scatter(h, xsveig(3,h), 'b*');
% %     scatter(h, xsveig(4,h), 'g*');
% % end
% % legend('r','a','b','l')
% 
% figure()
% title(['Eigenvalues of (A - Bksv) for ' param ' parameters, ' bflag '']);
% ylabel('|eig(A - Bksv)|');
% xlabel('Iterations');
% grid on;
% hold on
% for f=1:nsteps
% %     for n = 1:4
% %     scatter(f, absxeig(n,f), 'b*');
% %     end     
%     scatter(f, absxsveig(1,f), 'r*');
%     scatter(f, absxsveig(2,f), 'm*');
%     scatter(f, absxsveig(3,f), 'b*');
%     scatter(f, absxsveig(4,f), 'g*');
% end
% legend('r','a','b','l')
% 
% 
% %% Plots for multiple variable state feedback
% figure()
% title(['Deviations from fixed point for ' param ' parameters, ' bflag '']);
% ylabel('x_{n+1} = (A - Bk)x_n (multiple variable)');
% xlabel('Iterations');
% grid on;
% hold on
% for k=1:nsteps
%     scatter(k, x(1,k), 'r*');
%     scatter(k, x(2,k), 'm*');
%     scatter(k, x(3,k), 'b*');
%     scatter(k, x(4,k), 'g*');
% end
% legend('r','a','b','l')
% 
% figure()
% title(['Norm of deviations from fixed point for ' param ' parameters, ' bflag '']);
% ylabel('Norm(x_{n+1})');
% xlabel('Iterations');
% grid on;
% hold on
% for k=1:nsteps
%    scatter(k,normx(k),'b*')
% end
%  
% % figure()
% % title(['Eigenvalues of (A - Bk) for ' param ' parameters, ' bflag '']);
% % ylabel('eig(A - Bk)');
% % xlabel('Iterations');
% % grid on;
% % hold on
% % for h=1:nsteps
% % %     for m = 1:4
% % %     scatter(h, xeig(m,h), 'b*');
% % %     end  
% %     scatter(h, xeig(1,h), 'r*');
% %     scatter(h, xeig(2,h), 'm*');
% %     scatter(h, xeig(3,h), 'b*');
% %     scatter(h, xeig(4,h), 'g*');
% % end
% % legend('r','a','b','l')
% 
% figure()
% title(['Eigenvalues of (A - Bk) for ' param ' parameters, ' bflag '']);
% ylabel('|eig(A - Bk)|');
% xlabel('Iterations');
% grid on;
% hold on
% for f=1:nsteps
% %     for n = 1:4
% %     scatter(f, absxeig(n,f), 'b*');
% %     end     
%     scatter(f, absxeig(1,f), 'r*');
%     scatter(f, absxeig(2,f), 'm*');
%     scatter(f, absxeig(3,f), 'b*');
%     scatter(f, absxeig(4,f), 'g*');
% end
% legend('r','a','b','l')


%% Single variable state feedback function
function [xsv,xsveig,maxabsxsveig,normxsv,absxsveig] = singlevarsf(x0,jac,Bmatrix,Ksc,nsteps)
    
xsv       = zeros(4,length(nsteps));
xxsveig    = zeros(4,length(nsteps));
absxsveig = zeros(4,length(nsteps));
normxsv   = zeros(1,length(nsteps));
xsv(:,1)  = x0;
%Ksv = [0 0 0 -Ksc]; %The position of -Ksc depends on where control is being applied hence can change
Ksv = [0 0 0 0]; %The position of -Ksc depends on where control is being applied hence can change
Ksv(Bmatrix==1) = -Ksc; 
xsveig    = eig(jac - Bmatrix*Ksv);
maxabsxsveig = max(abs(xsveig));
for j = 1:nsteps
    
    %u   = Ksv*xsv;
    xsv(:,j+1)       = jac*xsv(:,j) - Bmatrix*Ksv*xsv(:,j); %x(n+1) = (A - B*K)*x(n)
    xxsveig(:,j+1)    = eig(jac - Bmatrix*Ksv);
    absxsveig(:,j+1) = abs(xxsveig(:,j+1));
    normxsv(:,j+1)   = norm(xsv(:,j+1)); %norm at each time step of xn
end
end


%% multi variable state feedback function
function [x,xeig,maxabsxeig,normx,absxeig] = multvarsf(x0,jac,Bmatrix,Kpl,nsteps)
    
x       = zeros(4,length(nsteps));
xxeig    = zeros(4,length(nsteps));
absxeig = zeros(4,length(nsteps));
normx   = zeros(1,length(nsteps));
x(:,1) = x0;
xeig = eig(jac - Bmatrix*Kpl);  
maxabsxeig = max(abs(xeig));
%absxeig = abs(xeig);
%normx = norm(x); %norm at each time step of xn
%maxabsxeig = absxeig; We computed absxeig for purposes of plotting and
%observing the eigenvalues with respect to the iterations.

for j = 1:nsteps
    x(:,j+1)    = jac*x(:,j) - Bmatrix*Kpl*x(:,j); %x(n+1) = (A - B*K)*x(n)
    xxeig(:,j+1)   = eig(jac - Bmatrix*Kpl);
    absxeig(:,j+1) = abs(xxeig(:,j+1));
    normx(:,j+1)   = norm(x(:,j+1)); %norm at each time step of xn
end
end
