%This code helps with the eigenvector test (cos? = (???)/????b?). The code
%computes ??? and normalizes it (???)/????B?). The B represents the 
%control strategies (Bmatrix).The left eigenvectors(w) is obtained from 
%eigfolder (The eigenvalues in this folder were obtained from the 
%compute_eigs.m). The eigfolder also contains the eigenvalue 
%magnitudes(absolute egenvaues), and both the left and right eigenvalues.  
%The code also computes the inverse cosine to obtain angles (this is the 
%controllability angles). cos? = (???)/????b? and ? (controllability 
%angles) are then stored in wdotbfolder.

clear variables;

%selected_bcls = [600:-10:80];
selected_bcls = [600:-10:240];
parameterflag = 1;
B             = 2;

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

% % B Matrix
% if (B == 1)
%     bflag   = 'B1'; 
%     Bmatrix = [1;1;1;1];
% elseif (B == 2)
%     bflag   = 'B2'; 
%     Bmatrix = [1;0;0;0];
% elseif (B == 3)
%     bflag   = 'B3'; 
%     Bmatrix = [0;1;0;0];
% elseif (B == 4)
%     bflag   = 'B4'; 
%     Bmatrix = [0;0;1;0];
% elseif (B == 5)
%     bflag   = 'B5'; 
%     Bmatrix = [0;0;0;1];
% elseif (B == 6)
%     bflag   = 'B6'; 
%     Bmatrix = [1;1;0;0];
% elseif (B == 7)
%     bflag   = 'B7'; 
%     Bmatrix = [0;0;1;1];
% elseif (B == 8)
%     bflag   = 'B8'; 
%     Bmatrix = [0;1;1;0];
% elseif (B == 9)
%     bflag   = 'B9'; 
%     Bmatrix = [1;0;1;0];
% elseif (B == 10)
%     bflag   = 'B10'; 
%     Bmatrix = [0;1;0;1];
% elseif (B == 11)
%     bflag   = 'B11'; 
%     Bmatrix = [1;0;0;1];
% elseif (B == 12)
%     bflag   = 'B12'; 
%     Bmatrix = [0;1;1;1];
% elseif (B == 13)
%     bflag   = 'B13'; 
%     Bmatrix = [1;0;1;1];
% elseif (B == 14)
%     bflag   = 'B14'; 
%     Bmatrix = [1;1;0;1];   
% elseif (B == 15)
%     bflag   = 'B15'; 
%     Bmatrix = [1;1;1;0];
% end

%folders
eigfolder    = 'Eigenvalues/'; %folder where eigenvalues are stored. 
wdotbfolder  = 'wdotbfolder/'; %folder where wdotb will be saved. 
wdotbfolder1 = 'wdotbfolder/wdotb1/';
wdotbfolder2 = 'wdotbfolder/wdotb2/';
wdotbfolder3 = 'wdotbfolder/wdotb3/';
wdotbfolder4 = 'wdotbfolder/wdotb4/';
wdotbfolder5 = 'wdotbfolder/wdotb/';
wdotb1       = cell(1,length(selected_bcls)); 
wdotb2       = cell(1,length(selected_bcls));
wdotb3       = cell(1,length(selected_bcls));
wdotb4       = cell(1,length(selected_bcls));
wdotb        = cell(1,length(selected_bcls));
theta        = cell(1,length(selected_bcls));
%wdotbscaled  = cell(1,length(selected_bcls));
for i = 1:length(selected_bcls)
    eval(['load ' eigfolder 'alleigs' num2str(selected_bcls(i)) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
    bcl = selected_bcls(i);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
        
    Lefteig  = allw{:,i};
    Lefteig1 = allw{:,i}(:,1);
    Lefteig2 = allw{:,i}(:,2);
    Lefteig3 = allw{:,i}(:,3);
    Lefteig4 = allw{:,i}(:,4);

%   C1 = dot(Lefteig1,B3);
%   C2 = dot(Lefteig2,B3);
%   C3 = dot(Lefteig3,B3);
%   C4 = dot(Lefteig4,B3);
% computing for w.b and normalizing it (w.b)/(||w|| ||b||)
    C1 = (Lefteig1'*Bmatrix)/(norm(Lefteig1)*norm(Bmatrix));
    C2 = (Lefteig2'*Bmatrix)/(norm(Lefteig2)*norm(Bmatrix));
    C3 = (Lefteig3'*Bmatrix)/(norm(Lefteig3)*norm(Bmatrix));
    C4 = (Lefteig4'*Bmatrix)/(norm(Lefteig4)*norm(Bmatrix));
    
    P  = (Lefteig'*Bmatrix)/(norm(Lefteig)*norm(Bmatrix));
% inverse cosine to obtain angles
theta1 = acosd(C1);
theta2 = acosd(C2);
theta3 = acosd(C3);
theta4 = acosd(C4);

% Putting the normalized C's in an array     
    C  = [C1; C2; C3; C4];
    Cnorm = norm(C);
    Pnorm = norm(P);
% Putting the angles in an array
theta{i} = [theta1; theta2; theta3; theta4];
%   C  = sort(C);
%   C  = C;

% Taking the absolute of the C's inorder not to have negative values
    wdotb{i} = abs(C);
%   wdotb1{i} = C(1);
%   wdotb2{i} = C(2);
%   wdotb3{i} = C(3);
%   wdotb4{i} = C(4);
    
%scaling
%     wdotbmax    = max(wdotb{i});
%     wdotbscaled{i} = (1/wdotbmax)*wdotb{i};

        T1 = table(bcl(:),parameterflag,B, Cnorm,...
            'VariableNames',{'BCL','parameterflag','B','Norm_cos_theta'})
        T2 = table(Bmatrix,wdotb{i},...
            'VariableNames',{'Bmatrix','Cos_theta'})


    
    myfilename5=[wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)];
    save([myfilename5 '.mat'],'wdotb','bcl','parameterflag', 'Bmatrix','B','theta')%%%

%    myfilename=[wdotbfolder 'wdotb' num2str(bcl) '_pflag' num2str(parameterflag)];
%    save([myfilename '.mat'],'wdotb1','wdotb2','wdotb3','wdotb4','bcl','parameterflag')%%%
    
%     myfilename=[wdotbfolder 'wdotb' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename '.mat'],'wdotb1','wdotb2','wdotb3','wdotb4','bcl','parameterflag', 'Bmatrix','B')%%%
%     
%     myfilename1=[wdotbfolder1 'wdotb1_' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename1 '.mat'],'wdotb1','bcl','parameterflag', 'Bmatrix','B')%%%
% 
%     myfilename2=[wdotbfolder2 'wdotb2_' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename2 '.mat'],'wdotb2','bcl','parameterflag', 'Bmatrix','B')%%%
% 
%     myfilename3=[wdotbfolder3 'wdotb3_' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename3 '.mat'],'wdotb3','bcl','parameterflag', 'Bmatrix','B')%%%
% 
%     myfilename4=[wdotbfolder4 'wdotb4_' num2str(bcl) '_pflag' num2str(parameterflag)];
%     save([myfilename4 '.mat'],'wdotb4','bcl','parameterflag', 'Bmatrix','B')%%%

end

