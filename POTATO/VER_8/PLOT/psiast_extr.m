load psiast_vsR_p.dat
load psiast_vsR_m.dat
load xpoints.dat
load opoints.dat
load seps.dat
load sepcr.dat
load rhop_dom.dat
plot(psiast_vsR_m(:,1),psiast_vsR_m(:,2),psiast_vsR_p(:,1),psiast_vsR_p(:,2),xpoints(:,1),xpoints(:,2),'kx',opoints(:,1),opoints(:,2),'ko',sepcr(:,1),sepcr(:,2),'k*',rhop_dom(:,1),rhop_dom(:,2),'k',seps(:,1),seps(:,2),'k--')
xlim([140 202])
xlabel('R_c [cm]')
legend('\psi^*,  v_{||}<0','\psi^*,  v_{||}>0','X-points','O-points','sep. cr.','\rho_{pol} domain','location','northwest')
title('\psi^*=\psi^*(H_0,J_\perp,\sigma,R_c)')
text(151.2,-4e6,'1')
text(149,-4e6,'2')
text(167,-5e5,'5')
text(176.5,-1.2e6,'6')
text(187,-4.55e6,'7')
text(172,5.9e5,'8')
text(186.5,-1.2e6,'9')
text(194.7,-4.55e6,'10')
print -dpng psiastofRc.png
%
load boundaries_beg.dat
load boundaries_end.dat
%
plot(xpoints(:,1),xpoints(:,2),'kx',boundaries_end(:,1),boundaries_end(:,2),'ks',psiast_vsR_m(:,1),psiast_vsR_m(:,2),psiast_vsR_p(:,1),psiast_vsR_p(:,2),seps(:,1),seps(:,2),'k--')
xlim([158 158.4])
ylim([-1.83e6 -1.74e6])
legend('X-point','v_{||}=0','location','northwest')
xlabel('R_c [cm]')
text(158.245,-1.797e6,'1')
text(158.12,-1.775e6,'2')
text(158.28,-1.7635e6,'3')
print -dpng psiastofRc_X1.png
%
%
plot(xpoints(:,1),xpoints(:,2),'kx',boundaries_beg(:,1),boundaries_beg(:,2),'ks',psiast_vsR_m(:,1),psiast_vsR_m(:,2),psiast_vsR_p(:,1),psiast_vsR_p(:,2),seps(:,1),seps(:,2),'k--')
xlim([163.2 163.6])
ylim([-7.5e5 -6.9e5])
legend('X-point','v_{||}=0','location','northeast')
xlabel('R_c [cm]')
text(163.307,-7.282e5,'4')
text(163.5,-7.282e5,'5')
text(163.307,-6.99e5,'8')
print -dpng psiastofRc_X2.png
