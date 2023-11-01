load_freqs
setcolors
%
% columns: x,psi^*,tau_b,delta_phi,...
%
plot(fi1(:,2),fi1(:,3),'color',clrs(1,:))
hold on
plot(fi2(:,2),fi2(:,3),'color',clrs(2,:))
plot(fi3(:,2),fi3(:,3),'color',clrs(3,:))
plot(fi4(:,2),fi4(:,3),'color',clrs(4,:))
plot(fi5(:,2),fi5(:,3),'color',clrs(5,:))
plot(fi6(:,2),fi6(:,3),'color',clrs(6,:))
plot(fi7(:,2),fi7(:,3),'color',clrs(7,:))
plot(fi8(:,2),fi8(:,3),'color',clrs(8,:))
plot(fi9(:,2),fi9(:,3),'color',clrs(9,:))
plot(fi10(:,2),fi10(:,3),'color',clrs(10,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('all classes')
print -dpng taub_all.png
%
%
plot(fi1(:,2),fi1(:,3),'color',clrs(1,:))
hold on
plot(fg1(:,2),fg1(:,3),'o','color',clrs(1,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 1')
print -dpng taub1.png
plot(fi1(:,1),fi1(:,3),'color',clrs(1,:))
hold on
plot(fg1(:,1),fg1(:,3),'o','color',clrs(1,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 1')
print -dpng taub_x1.png
%
plot(fi2(:,2),fi2(:,3),'color',clrs(2,:))
hold on
plot(fg2(:,2),fg2(:,3),'o','color',clrs(2,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 2')
print -dpng taub2.png
plot(fi2(:,1),fi2(:,3),'color',clrs(2,:))
hold on
plot(fg2(:,1),fg2(:,3),'o','color',clrs(2,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 2')
print -dpng taub_x2.png
%
plot(fi3(:,2),fi3(:,3),'color',clrs(3,:))
hold on
plot(fg3(:,2),fg3(:,3),'o','color',clrs(3,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 3')
print -dpng taub3.png
plot(fi3(:,1),fi3(:,3),'color',clrs(3,:))
hold on
plot(fg3(:,1),fg3(:,3),'o','color',clrs(3,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 3')
print -dpng taub_x3.png
%    
plot(fi4(:,2),fi4(:,3),'color',clrs(4,:))
hold on
plot(fg4(:,2),fg4(:,3),'o','color',clrs(4,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 4')
print -dpng taub4.png
plot(fi4(:,1),fi4(:,3),'color',clrs(4,:))
hold on
plot(fg4(:,1),fg4(:,3),'o','color',clrs(4,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 4')
print -dpng taub_x4.png
%    
plot(fi5(:,2),fi5(:,3),'color',clrs(5,:))
hold on
plot(fg5(:,2),fg5(:,3),'o','color',clrs(5,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 5')
print -dpng taub5.png
plot(fi5(:,1),fi5(:,3),'color',clrs(5,:))
hold on
plot(fg5(:,1),fg5(:,3),'o','color',clrs(5,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 5')
print -dpng taub_x5.png
%    
plot(fi6(:,2),fi6(:,3),'color',clrs(6,:))
hold on
plot(fg6(:,2),fg6(:,3),'o','color',clrs(6,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 6')
print -dpng taub6.png
plot(fi6(:,1),fi6(:,3),'color',clrs(6,:))
hold on
plot(fg6(:,1),fg6(:,3),'o','color',clrs(6,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 6')
print -dpng taub_x6.png
%    
plot(fi7(:,2),fi7(:,3),'color',clrs(7,:))
hold on
plot(fg7(:,2),fg7(:,3),'o','color',clrs(7,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 7')
print -dpng taub7.png
plot(fi7(:,1),fi7(:,3),'color',clrs(7,:))
hold on
plot(fg7(:,1),fg7(:,3),'o','color',clrs(7,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 7')
print -dpng taub_x7.png
%    
plot(fi8(:,2),fi8(:,3),'color',clrs(8,:))
hold on
plot(fg8(:,2),fg8(:,3),'o','color',clrs(8,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 8')
print -dpng taub8.png
plot(fi8(:,1),fi8(:,3),'color',clrs(8,:))
hold on
plot(fg8(:,1),fg8(:,3),'o','color',clrs(8,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 8')
print -dpng taub_x8.png
%      
plot(fi9(:,2),fi9(:,3),'color',clrs(9,:))
hold on
plot(fg9(:,2),fg9(:,3),'o','color',clrs(9,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 9')
print -dpng taub9.png
plot(fi9(:,1),fi9(:,3),'color',clrs(9,:))
hold on
plot(fg9(:,1),fg9(:,3),'o','color',clrs(9,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 9')
print -dpng taub_x9.png
%
plot(fi10(:,2),fi10(:,3),'color',clrs(10,:))
hold on
plot(fg10(:,2),fg10(:,3),'o','color',clrs(10,:))
hold off
xlim([-8e6 1e6])
ylim([0 8e4])
xlabel('\psi^*')
ylabel('\tau_b')
title('class 10')
print -dpng taub10.png
plot(fi10(:,1),fi10(:,3),'color',clrs(10,:))
hold on
plot(fg10(:,1),fg10(:,3),'o','color',clrs(10,:))
hold off
xlabel('x')
ylabel('\tau_b')
title('class 10')
print -dpng taub_x10.png
%

