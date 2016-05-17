clear

flux = load('flux/test_flux_n3_A10_Mt2m2.dat');

figure(3)
clf
[ax,p1,p2] = plotyy(flux(:,1), flux(:,4)/10, flux(:,1), flux(:,5));
axpos = get(ax(1), 'Position');
set(ax(1),'Position', axpos + [0 .02 -.03 -.04]);
set(ax(2),'Position', get(ax(1), 'Position'));
set(gca,'XMinorTick','on','YMinorTick','on');
xlabel('u');
ylabel(ax(1),'D_{11} subintegrand [a.u.]');
ylabel(ax(2),'$\Delta \bar{\eta}$', 'Interpreter', 'LaTex') 
set(ax(2), 'YTick', [0,1]);
title('A=10, m=0, n=3');