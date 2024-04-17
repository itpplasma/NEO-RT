clear

% colormap
c = [     0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840 ];

flux1 = load('srange/test_srange_n3_Mt1m5_dqdr.dat');
flux2 = load('srange/test_srange_n3_Mtm1m5_dqdr.dat');
flux3 = load('srange/test_srange_n3_Mt0.dat');
shaing1 = load('srange/D11DpShaing_n3_Mt1m5');
shaing2 = load('srange/D11DpShaing_n3_Mtm1m5');
shaing3 = load('srange/D11DpShaing_n3_Mt0');

figure(1)
clf
line(flux2(flux2(:,2)<.18,2), flux2(flux2(:,2)<.18,4), 'Color', c(1,:))
line(flux1(:,2), flux1(:,4), 'Color', c(2,:))
line(shaing2.epsrShaing, shaing2.D11DpShaing, 'Color', c(1,:), ...
     'LineStyle', '--')
line(shaing1.epsrShaing, shaing1.D11DpShaing, 'Color', c(2,:), ...
     'LineStyle', '--')
xlabel('1/A');
ylabel('D_{11}/D_p');
ylim([0,2e7]);
set(gca,'XMinorTick','on','YMinorTick','on')
l = legend('\Omega_{tE} = -\Omega_{\rm ref}', '\Omega_{tE} = +\Omega_{\rm ref}', 'Location', 'northwest')
lpos = get(l,'Position');
set(l, 'Position', lpos + [0.05, 0, 0, 0])
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'XAxisLocation','top',...
    'Color','none', 'Xtick', []);
%linkaxes([ax1, ax2], 'x')
line(flux1(:,2), flux1(:,3), 'Parent', ax2, 'Color', 'k', 'LineStyle', ...
     '-.');

%plot(flux3(:,2), flux3(:,4), 'g')
%plot(shaing3.epsrShaing, shaing3.D11DpShaing.*(1d5./(pi*shaing3.qShaing.')), 'g--')
%legend('M_t = -10^{-5}', 'M_t = 10^{-5}', 'q profile')
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('q');
ylim([1,2]);
title('m=0, n=3')
