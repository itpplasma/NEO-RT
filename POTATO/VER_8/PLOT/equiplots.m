load compeqpar_noE.dat
load testeqpar_noE.dat
load compeqpar_E.dat
load testeqpar_E.dat
load compeqpar_mE.dat
load testeqpar_mE.dat
%
rho_c=sqrt(compeqpar_noE(:,1));
rho_t=sqrt(testeqpar_noE(:,1));
plot(rho_t,testeqpar_noE(:,2),'k:','linewidth',2,rho_c,compeqpar_noE(:,2),'b',rho_c,compeqpar_E(:,2),'r',rho_c,compeqpar_mE(:,2),'m')
xlim([0 0.7])
ylim([0 1.1])
legend('reference','no E_r','E_r > 0','E_r < 0')
xlabel('\rho_{pol}')
ylabel('n_\alpha')
print -dpng dens.png
%
plot(rho_c,compeqpar_noE(:,3),'b',rho_c,compeqpar_E(:,3),'r',rho_c,compeqpar_mE(:,3),'m',rho_c,compeqpar_E(:,3)-compeqpar_noE(:,3),'g',rho_t,testeqpar_noE(:,3),'b--',rho_t,testeqpar_E(:,3),'r--',rho_t,testeqpar_mE(:,3),'m--',rho_t,testeqpar_E(:,3)-testeqpar_noE(:,3),'g--')
xlim([0 0.7])
legend('(1) no E_r','(2) E_r > 0','(3) E_r < 0','(3) - (2)')
xlabel('\rho_{pol}')
ylabel('<n_\alpha V^\phi_{g\alpha}>')
print -dpng torvel.png

