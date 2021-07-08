load_freqs
setcolors
%
% columns: x,psi^*,tau_b,delta_phi,...
%
%
plot(fi1(:,2),fi1(:,4),'color',clrs(1,:))
hold on
plot(fi2(:,2),fi2(:,4),'color',clrs(2,:))
plot(fi3(:,2),fi3(:,4),'color',clrs(3,:))
plot(fi4(:,2),fi4(:,4),'color',clrs(4,:))
plot(fi5(:,2),fi5(:,4),'color',clrs(5,:))
plot(fi6(:,2),fi6(:,4),'color',clrs(6,:))
plot(fi7(:,2),fi7(:,4),'color',clrs(7,:))
plot(fi8(:,2),fi8(:,4),'color',clrs(8,:))
plot(fi9(:,2),fi9(:,4),'color',clrs(9,:))
plot(fi10(:,2),fi10(:,4),'color',clrs(10,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('all classes')
print -dpng delphi_all.png
%
%
plot(fi1(:,2),fi1(:,4),'color',clrs(1,:))
hold on
plot(fg1(:,2),fg1(:,4),'o','color',clrs(1,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 1')
print -dpng delphi1.png
plot(fi1(:,1),fi1(:,4),'color',clrs(1,:))
hold on
plot(fg1(:,1),fg1(:,4),'o','color',clrs(1,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 1')
print -dpng delphi_x1.png
%
plot(fi2(:,2),fi2(:,4),'color',clrs(2,:))
hold on
plot(fg2(:,2),fg2(:,4),'o','color',clrs(2,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 2')
print -dpng delphi2.png
plot(fi2(:,1),fi2(:,4),'color',clrs(2,:))
hold on
plot(fg2(:,1),fg2(:,4),'o','color',clrs(2,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 2')
print -dpng delphi_x2.png
%
plot(fi3(:,2),fi3(:,4),'color',clrs(3,:))
hold on
plot(fg3(:,2),fg3(:,4),'o','color',clrs(3,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 3')
print -dpng delphi3.png
plot(fi3(:,1),fi3(:,4),'color',clrs(3,:))
hold on
plot(fg3(:,1),fg3(:,4),'o','color',clrs(3,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 3')
print -dpng delphi_x3.png
%    
plot(fi4(:,2),fi4(:,4),'color',clrs(4,:))
hold on
plot(fg4(:,2),fg4(:,4),'o','color',clrs(4,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 4')
print -dpng delphi4.png
plot(fi4(:,1),fi4(:,4),'color',clrs(4,:))
hold on
plot(fg4(:,1),fg4(:,4),'o','color',clrs(4,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 4')
print -dpng delphi_x4.png
%    
plot(fi5(:,2),fi5(:,4),'color',clrs(5,:))
hold on
plot(fg5(:,2),fg5(:,4),'o','color',clrs(5,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 5')
print -dpng delphi5.png
plot(fi5(:,1),fi5(:,4),'color',clrs(5,:))
hold on
plot(fg5(:,1),fg5(:,4),'o','color',clrs(5,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 5')
print -dpng delphi_x5.png
%    
plot(fi6(:,2),fi6(:,4),'color',clrs(6,:))
hold on
plot(fg6(:,2),fg6(:,4),'o','color',clrs(6,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 6')
print -dpng delphi6.png
plot(fi6(:,1),fi6(:,4),'color',clrs(6,:))
hold on
plot(fg6(:,1),fg6(:,4),'o','color',clrs(6,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 6')
print -dpng delphi_x6.png
%    
plot(fi7(:,2),fi7(:,4),'color',clrs(7,:))
hold on
plot(fg7(:,2),fg7(:,4),'o','color',clrs(7,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 7')
print -dpng delphi7.png
plot(fi7(:,1),fi7(:,4),'color',clrs(7,:))
hold on
plot(fg7(:,1),fg7(:,4),'o','color',clrs(7,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 7')
print -dpng delphi_x7.png
%    
plot(fi8(:,2),fi8(:,4),'color',clrs(8,:))
hold on
plot(fg8(:,2),fg8(:,4),'o','color',clrs(8,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 8')
print -dpng delphi8.png
plot(fi8(:,1),fi8(:,4),'color',clrs(8,:))
hold on
plot(fg8(:,1),fg8(:,4),'o','color',clrs(8,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 8')
print -dpng delphi_x8.png
%      
plot(fi9(:,2),fi9(:,4),'color',clrs(9,:))
hold on
plot(fg9(:,2),fg9(:,4),'o','color',clrs(9,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 9')
print -dpng delphi9.png
plot(fi9(:,1),fi9(:,4),'color',clrs(9,:))
hold on
plot(fg9(:,1),fg9(:,4),'o','color',clrs(9,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 9')
print -dpng delphi_x9.png
%
plot(fi10(:,2),fi10(:,4),'color',clrs(10,:))
hold on
plot(fg10(:,2),fg10(:,4),'o','color',clrs(10,:))
hold off
xlim([-8e6 1e6])
ylim([-12 12])
xlabel('\psi^*')
ylabel('\Delta\phi_b')
title('class 10')
print -dpng delphi10.png
plot(fi10(:,1),fi10(:,4),'color',clrs(10,:))
hold on
plot(fg10(:,1),fg10(:,4),'o','color',clrs(10,:))
hold off
xlabel('x')
ylabel('\Delta\phi_b')
title('class 10')
print -dpng delphi_x10.png
%

