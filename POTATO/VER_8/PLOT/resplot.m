noE10=load('subintegrand94_vsJperp.noE10');
noE20=load('subintegrand94_vsJperp.noE20');
E10=load('subintegrand94_vsJperp.E10');
mE10=load('subintegrand94_vsJperp.mE10');
%
%
plot(noE10(:,1),noE10(:,2),'k:','linewidth',2,noE10(:,1),noE10(:,9-2),'b',noE10(:,1),noE10(:,9-1),'g',noE10(:,1),noE10(:,9),'r')
xlim([4.6e-5 5.8e-5])
ylim([-5e4 5e5])
legend('total','m_2 = -2','m_2 = -1','m_2 = 0','location','northeast')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('E_r = 0,  m_3 = 2')
print -dpng res_noE_modes.png
%
xlim([5.25e-5 5.32e-5])
xlabel('normalized J_\perp')
ylabel('integral over x')
title('E_r = 0,  m_3 = 2')
print -dpng res_noE_modes_zoom.png
%
%
plot(E10(:,1),E10(:,2),'k:','linewidth',2,E10(:,1),E10(:,9-3),'g',E10(:,1),E10(:,9-2),'b',E10(:,1),E10(:,9-1),'m',E10(:,1),E10(:,9+1),'r')
xlim([4.e-5 5.5e-5])
ylim([-2.5e5 2.7e6])
legend('total','m_2 = -3','m_2 = -2','m_2 = -1','m_2 = 1','location','northeast')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('E_r > 0,  m_3 = 2')
print -dpng res_E_modes.png
%
xlim([4.03e-5 4.17e-5])
ylim([-5e4 3e5])
xlabel('normalized J_\perp')
ylabel('integral over x')
legend('total','m_2 = -3','m_2 = -2','m_2 = -1','m_2 = 1','location','northwest')
title('E_r > 0,  m_3 = 2')
print -dpng res_E_modes_zoom.png
%
plot(mE10(:,1),mE10(:,2),'k:','linewidth',2,mE10(:,1),mE10(:,9-2),'b',mE10(:,1),mE10(:,9-1),'m',mE10(:,1),mE10(:,9),'g',mE10(:,1),mE10(:,9+1),'r',mE10(:,1),mE10(:,9+2),'k',mE10(:,1),mE10(:,9+3),'c',mE10(:,1),mE10(:,9+4),'y')
xlim([3.e-5 6.e-5])
ylim([-3.5e4 1.5e4])
legend('total','m_2 = -2','m_2 = -1','m_2 = 0','m_2 = 1','m_2 = 2','m_2 = 3','m_2 = 4','location','southwest')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('E_r < 0,  m_3 = 2')
print -dpng res_mE_modes.png
%
xlim([5.5e-5 5.9e-5])
legend('total','m_2 = -2','m_2 = -1','m_2 = 0','m_2 = 1','m_2 = 2','m_2 = 3','m_2 = 4','location','southeast')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('E_r < 0,  m_3 = 2')
print -dpng res_mE_modes_zoom.png
%
%
plot(noE10(:,1),noE10(:,2),'b','linewidth',2,E10(:,1),E10(:,2),'r','linewidth',2,mE10(:,1),mE10(:,2),'g','linewidth',2)
xlim([0 6.5e-5])
ylim([-2.5e5 2.7e6])
legend('E_r = 0','E_r > 0','E_r < 0','location','northwest')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('total,  m_3 = 2')
print -dpng res_total.png
%
xlim([4e-5 5.8e-5])
ylim([-2.5e5 2.7e6])
legend('E_r = 0','E_r > 0','E_r < 0','location','northwest')
xlabel('normalized J_\perp')
ylabel('integral over x')
title('total,  m_3 = 2')
print -dpng res_total_zoom.png
