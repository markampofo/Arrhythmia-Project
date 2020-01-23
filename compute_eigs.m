% This code computes the eigenvalues by loading the jacobians obtained from
% Jacfolder (The jacobians in this folder was obtained from the 
%uclacalcium_firstorder.m). This code also saves the eigenvalue magnitudes
%(absolute egenvaues), and both the left and right eigenvalues in the 
%folder "jacfolder". 

clear variables;

selected_bcls = [600:-10:80];
parameterflag = 4;

jacfolder  = 'jacfolder/'; % folder where jacobians are stored
eigfolder  = 'Eigenvalues/'; %folder where eigenvalues will be saved. 

alleigs    = cell(1,length(selected_bcls)); % Store eigenvalues here.
alleigsabs = cell(1,length(selected_bcls)); % Store eigenvalue magnitudes here.
allv       = cell(1,length(selected_bcls)); % Store right eigenvectors here.
allw       = cell(1,length(selected_bcls)); % Store left eigenvectors here.

eigw       = cell(length(selected_bcls)); 
eigv       = cell(length(selected_bcls));
eigr       = cell(length(selected_bcls));

for i = 1:length(selected_bcls)
    eval(['load ' jacfolder 'jac' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
    bcl = selected_bcls(i);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
        
    [allv{i}, eigv{i}, allw{i}] = eig(jac,'nobalance');
    alleigs{i}    = diag(eigv{i});
    alleigsabs{i} = abs(alleigs{i});
    
    
    myfilename=[eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag)];
    save([myfilename '.mat'],'alleigs','alleigsabs','allw','allv','eigv','bcl','parameterflag')%%%

   % eval(['save ' eigfolder 'eigfile' num2str(logepsln) ' *'])
end

