%This code takes eigenvalues (from eigfolder) and the normalized ???
%(ie. cos? = (???)/????b? from wdotbfolder). It then plots
%controllability magnitudes and angles of controllability magnitude 
set(0,'defaultlinelinewidth',3)
set(0,'defaultaxesfontsize',32)


selected_bcls = [600:-10:80];
parameterflag = 4;
B             = 4;

%folders
eigfolder    = 'Eigenvalues/';
wdotbfolder  = 'wdotbfolder/'; %folder where wdotb will be saved. 
wdotbfolder1 = 'wdotbfolder/wdotb1/';
wdotbfolder2 = 'wdotbfolder/wdotb2/';
wdotbfolder3 = 'wdotbfolder/wdotb3/';
wdotbfolder4 = 'wdotbfolder/wdotb4/';
wdotbfolder5 = 'wdotbfolder/wdotb/';

% if (parameterflag == 0)
%     param = 'mixed mechanism 1';
% elseif (parameterflag == 1)
%     param = 'mixed mechanism 2';
% elseif (parameterflag == 2)
%     param = 'calcium-driven alternans';
% elseif (parameterflag == 3)
%     param = 'voltage-driven alternans';
% end

if parameterflag == 0
    param = ('Shallower (Restitution Curve)');    
elseif parameterflag == 1 %steep parameters
   % param = 'Steep Restitution Curve';% voltage driven alternans
   param = ('Voltage-Driven Alternans ');
elseif parameterflag == 2 %quasiperiodicity 
    param = 'Quasiperiodicity (Restitution Curve)';
elseif parameterflag == 3 %And, of course, the biphasic restitution curve wreaks much havoc.
    param = 'Biphasic (Restitution Curve)'; 
elseif parameterflag == 4
    %param = ('Shallow Restitution, Ca-induced');
    param = ('Calcium-Driven Alternans ');
% elseif parameterflag == 5 
%     param = ('Shallow, Ca-induced, Neg Coup');
% elseif parameterflag == 6 %steep restitution + Ca-induced
%     param = 'Steep + Ca';
% elseif parameterflag == 7
%     param = ('Shallow Restitution, Ca-induced 2');
% elseif parameterflag == 8
%     param = ('Shallow Restitution, Ca-induced 3');
elseif parameterflag == 9
    param = ('\nu = 0.35');
elseif parameterflag == 10
    param = ('\tau_0 = 20');
else
    disp('flag not recognized')
end



if (B == 1)
    bflag   = 'B_{r}'; 
    Bmatrix = [1;0;0;0];
elseif (B == 2)
    bflag   = 'B_{a}'; 
    Bmatrix = [0;1;0;0];
elseif (B == 3)
    bflag   = 'B_{b}'; 
    Bmatrix = [0;0;1;0];
elseif (B == 4)
    bflag   = 'B_{l}'; 
    Bmatrix = [0;0;0;1];
else
    disp('flag not recognized')
end


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


% figure()
% %title(['Eigenvalue moduli for ' param ' parameters, epsilon = 10^{' num2str(logepsln) '}']);
% ylabel('w.B');
% xlabel('Period (ms)');
% grid on;
% hold on
% for j=1:length(selected_bcls)
%     bcl = selected_bcls(j);
%     eval(['load ' wdotbfolder1 'wdotb1_' num2str(bcl) '_pflag' num2str(parameterflag)])
%     eval(['load ' wdotbfolder2 'wdotb2_' num2str(bcl) '_pflag' num2str(parameterflag)])
%     eval(['load ' wdotbfolder3 'wdotb3_' num2str(bcl) '_pflag' num2str(parameterflag)])
%     eval(['load ' wdotbfolder4 'wdotb4_' num2str(bcl) '_pflag' num2str(parameterflag)])
% 
%     %eval(['load ' eigfolder 'alleigs' num2str(selected_bcls_for_fps(i)) '_pflag' num2str(parameterflag) '_epsln' num2str(log10(epsln))]) %Load data from jacobians
%     scatter(selected_bcls(j), wdotb1{j}, 'b*');
%     scatter(selected_bcls(j), wdotb2{j}, 'r*');
%     scatter(selected_bcls(j), wdotb3{j}, 'y*');
%     scatter(selected_bcls(j), wdotb4{j}, 'g*');
% end
% %end
% hold off;
% 
% figure
% %title(['Eigenvalue moduli for ' param ' parameters, epsilon = 10^{' num2str(logepsln) '}']);
% ylabel('w.B');
% xlabel('Period (ms)');
% grid on;
% hold on
% 
% for k=1:length(selected_bcls)
%     bcl = selected_bcls(k);
%     eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
%     for r = 1:4
%     scatter(selected_bcls(k), wdotb{k}(r), 'b*');
%     end  
% end
% 
% 
% 
% figure
 markersize = 80;
% title(['Controllability Magnitude for ' param, ], 'FontSize', 8);
% ax = gca;
% ax.FontSize = 8;
% ylabel('|\lambda|');
% xlabel('Period (ms)');
% grid on;
% hold on
% for ii = 1:length(selected_bcls)
%     bcl = selected_bcls(ii);
% %scatter(bcl(ii), magnitudes of eigenvalues for iith BCL, markersize, magnitudes of w*B products for iith BCL);
%     eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
%     eval(['load ' eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
% 
%     for s = 1:4
%     scatter(selected_bcls(ii), alleigs{ii}(s), markersize, wdotb{ii}(s));
%     colormap(parula)
% %    colorbar
%     c = colorbar;
%     c.Label.String = '|cos \theta_{ki}|';
% 
%     end 
% end
% hold off;

set(0,'defaultlinelinewidth',3)
set(0,'defaultaxesfontsize',32)
figure
%markersize = 80;
%MarkerFaceColor = [.49 1 .63];
title(['Controllability Magnitude For ' param,' Using ' bflag ''], 'FontSize', 21);
%title(['Controllability Magnitude For ' param, ' Using ' bflag '']);
%ax = gca;
%ax.FontSize = 20;
ylabel('|\lambda|');
xlabel('Period (ms)');
grid on;
hold on
for jj = 1:length(selected_bcls)
    bcl = selected_bcls(jj);
    eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
    eval(['load ' eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag) ]) %Load data from jacobians

    for g = 1:4
  s = scatter(selected_bcls(jj), alleigsabs{jj}(g), markersize, wdotb{jj}(g),'filled');
    %...
   %'MarkerFaceColor','g');
    colormap(copper)
    c = colorbar;
    c.Label.String = '|cos \theta_{ij}|';
    s.MarkerEdgeColor = 'k';
    caxis([0 1])
    end 
end
hold off;
%saveas(gcf,'shallow(ca_induced)_ctrlmag_period_B1.png')


figure
markersize = 80;
title(['Controllability Magnitude For ' param, ' Using ' bflag ''], 'FontSize', 16);
ax = gca;
ax.FontSize = 20;
ylabel('|\lambda|');
xlabel('Period (ms)');
grid on;
hold on
for rr = 1:length(selected_bcls)
    bcl = selected_bcls(rr);
    eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
    eval(['load ' eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag) ]) %Load data from jacobians

    for g = 1:4
    d = scatter(selected_bcls(rr), alleigsabs{rr}(g), markersize, theta{rr}(g), 'filled');
    colormap(winter)
    d.MarkerEdgeColor = 'k';
    c = colorbar;
    c.Label.String = '\theta^{\circ}';
    %set(c,'FontSize',16);
    end 
end
hold off;
%saveas(gcf,'shallow(ca_induced)_ctrlmag_period_angle_B1.png')


figure
markersize = 80;
title(['Angles Of Controllability Magnitude For ' param, ' Using ' bflag ''], 'FontSize', 16);
ax = gca;
ax.FontSize = 20;
ylabel('\theta^{\circ}');
xlabel('Period (ms)');
grid on;
hold on
for uu = 1:length(selected_bcls)
    bcl = selected_bcls(uu);
    eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
    eval(['load ' eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag) ]) %Load data from jacobians

for g = 1:4
    p = scatter(selected_bcls(uu),theta{uu}(g),'b');
    p.LineWidth = 2;
    p.MarkerEdgeColor = 'k';
    p.MarkerFaceColor = [0 0 1];
end 
end
hold off;
%saveas(gcf,'shallow(ca_induced)_period_angle_B1.png')



% figure
% markersize = 80;
% title(['Controllability Magnitude for ' param, ' using ' bflag ', scaled'], 'FontSize', 8);
% ax = gca;
% ax.FontSize = 12;
% ylabel('|\lambda|');
% xlabel('Period (ms)');
% grid on;
% hold on
% for k = 1:length(selected_bcls)
%     bcl = selected_bcls(k);
%     eval(['load ' wdotbfolder5 'wdotb_' num2str(bcl) '_pflag' num2str(parameterflag) '_B' num2str(B)])
%     eval(['load ' eigfolder 'alleigs' num2str(bcl) '_pflag' num2str(parameterflag) ]) %Load data from jacobians
% 
%     for g = 1:4
%     scatter(selected_bcls(k), alleigsabs{k}(g), markersize, wdotbscaled{k}(g));
%     colormap(parula)
%     c = colorbar;
%     c.Label.String = '|cos \theta_{ki}|';
%     end 
% end
% hold off;


%end