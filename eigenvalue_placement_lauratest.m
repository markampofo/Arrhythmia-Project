%This code uses the matlab place function to compute for optimal gains for 
%a multi-variable state feedback. It obtains eigenvalues and jacobians from
%the jacfolder and eigfolder. Eigenvalues greater than one are selected 
%using a sort matlab inbuilt function and replaced with a desired 
%eigenvalue. This is used by the place function in matlab which
%then obtains an optimal gain to be used for a multi-variable state 
%feedback. 

%Given the single- or multi-input system x'=Ax+Bu and a vector p of desired 
%self-conjugate closed-loop pole locations, place computes a gain matrix K 
%such that the state feedback u = Kx places the closed-loop poles at the 
%locations p. In other words, the eigenvalues of A  BK match the entries 
%of p (up to the ordering). K = place(A,B,p) places the desired closed-loop 
%poles p by computing a state-feedback gain matrix K. All the inputs of the 
%plant are assumed to be control inputs. The length of p must match the row 
%size of A. place works for multi-input systems. The algorithm place relies 
%on uses the extra degrees of freedom to find a solution that minimizes the 
%sensitivity of the closed-loop poles to perturbations in A or B.
%[K,prec,message] = place(A,B,p) returns prec, an estimate of how closely 
%the eigenvalues of A  BK match the specified locations p (prec measures 
%the number of accurate decimal digits in the actual closed-loop poles). 
%If some nonzero closed-loop pole is more than 10% off from the desired 
%location, message contains a warning message.


clear variables

%selected_bcls = [600:-10:80];
%selected_bcls = [600:-10:200];
selected_bcls = 200;
parameterflag = 4;
B             = 4;

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

jacfolder  = 'jacfolder/'; % folder where jacobians are stored
eigfolder  = 'Eigenvalues/'; %folder where eigenvalues will be saved. 

% alleigs    = cell(1,length(selected_bcls)); % Store eigenvalues here.
alleigsabs = cell(1,length(selected_bcls)); % Store eigenvalue magnitudes here.
sorteigs = cell(1,length(selected_bcls)); % Store eigenvalue magnitudes here.
% allv       = cell(1,length(selected_bcls)); % Store right eigenvectors here.
% allw       = cell(1,length(selected_bcls)); % Store left eigenvectors here.
% 
% eigw       = cell(length(selected_bcls)); 
% eigv       = cell(length(selected_bcls));
% eigr       = cell(length(selected_bcls));

Tsave = cell(length(selected_bcls));

for i = 1:length(selected_bcls)
    
    eval(['load ' jacfolder 'jac' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
    eval(['load ' eigfolder 'alleigs' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from eigenvalue folder
    bcl = selected_bcls(i);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
        
    % Test: w'A = lambda w'
    w = allw{i}; %left eigenvectors
    G = (w'*jac) - (eigv{i}*w'); %w'A - ?w' = 0
        
    % sorted eigenvalues
    [sorteigs,sorteigsind] = sort(alleigs{:,i},'descend','ComparisonMethod','abs');
    %sorted left eigenvectors
    sortedlefteigs  = allw{:,i}(:, sorteigsind);
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
    %eigenvalues = [alleigsabs{:,i}(1); alleigsabs{:,i}(2); alleigsabs{:,i}(3); alleigsabs{:,i}(4)];
    %replacing eigenvalues not in the unit circle
%     if (sorteigs{:,i}(1) > 1)
%         sorteigs{:,i}(1) = 0.8;
%     elseif (sorteigs{:,i}(1) < -1)
%         sorteigs{:,i}(1) = -0.8;
%     elseif (sorteigs{:,i}(2) > 1)
%         sorteigs{:,i}(2) = 0.8;
%     elseif (sorteigs{:,i}(2) < -1)
%         sorteigs{:,i}(2) = -0.8;
%     elseif (sorteigs{:,i}(3) > 1)
%         sorteigs{:,i}(3) = 0.8;
%     elseif (sorteigs{:,i}(3) < -1)
%         sorteigs{:,i}(3) = -0.8;
%     elseif (sorteigs{:,i}(4) > 1)
%         sorteigs{:,i}(4) = 0.8;
%     elseif (sorteigs{:,i}(4) < -1)
%         sorteigs{:,i}(4) = -0.8;
%     end
%  if (sorteigs{:,i}(1) > .90)
        sorteigs(1) = (sorteigs(1))/2;
%     elseif (sorteigs{:,i}(1) < -.95)
%         sorteigs{:,i}(1) = (sorteigs{:,i}(1))/2;
%     elseif (sorteigs{:,i}(2) > .95)
%         sorteigs{:,i}(2) = (sorteigs{:,i}(2))/2;
%     elseif (sorteigs{:,i}(2) < -.95)
%         sorteigs{:,i}(2) = (sorteigs{:,i}(2))/2;
%     elseif (sorteigs{:,i}(3) > .95)
%         sorteigs{:,i}(3) = (sorteigs{:,i}(3))/2;
%     elseif (sorteigs{:,i}(3) < -.95)
%         sorteigs{:,i}(3) = (sorteigs{:,i}(3))/2;
%     elseif (sorteigs{:,i}(4) > .95)
%         sorteigs{:,i}(4) = (sorteigs{:,i}(4))/2;
%     elseif (sorteigs{:,i}(4) < -.95)
%         sorteigs{:,i}(4) = (sorteigs{:,i}(4))/2;
%  end
    %desired eigenvalues in the unit circle
%     desiredeigenvalues = [sorteigs{:,i}(1); sorteigs{:,i}(2); sorteigs{:,i}(3); sorteigs{:,i}(4)];
  desiredeigenvalues = [sorteigs(1); sorteigs(2); sorteigs(3); sorteigs(4)];
   
    if (desiredeigenvalues(:,:) == eigenvalues(:,:))
        disp(['Eigenvalues in unit circle for BCL = ' num2str(bcl) ' ms']) 
    else 
        [Kpl,placeprecision,placemessage] = place(jac,Bmatrix,desiredeigenvalues);
        Kplnorm = norm(Kpl);
        
        %Table
        T1 = table(bcl(:),parameterflag,B,placeprecision, Kplnorm, Cnorm,...
            'VariableNames',{'BCL','parameterflag','B','placeprecision','Norm_of_Kpl','Norm_of_C'})
        T2 = table(Bmatrix,eigenvalues(:,:),desiredeigenvalues(:,:),Kpl',C,...
            'VariableNames',{'Bmatrix','eigenvalues','desired_eigenvalues','Kpl','C'})
        
%         writetable(T1,'myData.txt','Delimiter',' ')  
%type 'myData.txt'
 %Tsave{i} = [T1 T2]; 
    end


    
%     myfilename=[eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename '.mat'],'alleigs','alleigsabs','allw','allv','bcl','parameterflag')%%%

   % eval(['save ' eigfolder 'eigfile' num2str(logepsln) ' *'])
end

