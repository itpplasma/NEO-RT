% colormap
c = [     0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840 ];

flux1 = load('machrange/test_machrange2_n18.dat');
flux2 = load('machrange/test_machrange2_n18_A4.dat');
flux3 = load('machrange/test_machrange2_n3.dat');
flux4 = load('machrange/test_machrange2_n3_A4.dat');
ntvout = load('machrange/ntv_out_all_nustar3m5_lag6.dat');
ntvout2 = load('machrange/ntv_out_all_nustar3m5_lag7.dat');
ntvout3 = load('machrange/ntv_out_all_nustar3m5_lag6_n18.dat');
ntvout4 = load('machrange/ntv_out_all_nustar3m5_lag7_n18.dat');
shaing = load('shaing/shaing_joint_n3.mat');

n0 = 3;
R0 = ntvout(1,19)
bmod0 = ntvout(1,23)*1e4   % this is innermost flux surface
%bmod0 = 17799             % this is manual from current flux surface
Bthcov = ntvout(1,17)
Bphcov = ntvout(1,16)
q = 1./ntvout(1,10)
sqrtgBth = ntvout(1,14)
epsm2 = ntvout(1,22)

%ntvDp = pi*vth**3*q/(16d0*R0*(qi*bmod0/(mi*c))**2)
Drp3 = 4/sqrt(pi)*R0/q^2*bmod0*n0*(Bthcov+q*Bphcov)/sqrtgBth^2
Drp18 = 4/sqrt(pi)*R0/q^2*bmod0*18*(Bthcov+q*Bphcov)/sqrtgBth^2
ntvMt = ntvout(1:end-1, 2);  % D11/Dp from Neo2
ntvD11 = ntvout(1:end-1, 8)*1e6;  % D11/Dp from Neo2
ntvMt2 = ntvout2(1:end-1, 2);  % D11/Dp from Neo2
ntvD112 = ntvout2(1:end-1, 8)*1e6;  % D11/Dp from Neo2
Mtcond = ntvMt<0.045;


figure(22)
clf
plot(flux1(:,1), flux1(:,3))
hold on
%plot(flux2(:,1), flux2(:,3))
plot([1e-5,0.1],[Drp18,Drp18], '-.')
%plot(ntvout3(1:end-1,2), ntvout3(1:end-1,8)*1e6)
%plot(ntvout4(1:end-1,2), ntvout4(1:end-1,8)*1e6)
set(gca,'XMinorTick','on','YMinorTick','on')
legend('Hamiltonian', 'Ripple plateau', 'Location', 'east')
xlabel('M_t');
ylabel('D_{11}/D_p');
title('A=10, m=0, n=18')

figure(21)
clf
plot(flux3(:,1), flux3(:,3))
hold on
%plot(flux4(:,1), flux4(:,3))
plot([1e-5,0.1],[Drp3,Drp3], '-.')
%plot(ntvMt(Mtcond), ntvD11(Mtcond))
plot(ntvMt(Mtcond), ntvD112(Mtcond), 's-')
plot(ntvMt, shaing.flux, '.-')
flux3i = interp1(flux3(:,1), flux3(:,3), ntvMt);
plot(ntvMt, shaing.flux + flux3i, 'x-');
set(gca,'XMinorTick','on','YMinorTick','on')
l = legend({'Hamiltonian', 'Ripple plateau', 'NEO-2', 'Shaing',...
        'Hamiltonian+Shaing'}, 'Location', 'northwest');
lpos = get(l,'Position');
set(l, 'Position', lpos + [0.05, -0.02, 0, 0])

xlabel('M_t');
ylabel('D_{11}/D_p');
ylim([0,740]);
title('A=10, m=0, n=3')