load_orbs
setcolors
%
plot(o1(:,1),o1(:,3),'color',clrs(1,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 1')
print -dpng o1.png
%
plot(o2(:,1),o2(:,3),'color',clrs(2,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 2')
print -dpng o2.png
%
plot(o3(:,1),o3(:,3),'color',clrs(3,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 3')
print -dpng o3.png
%
plot(o4(:,1),o4(:,3),'color',clrs(4,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 4')
print -dpng o4.png
%
plot(o5(:,1),o5(:,3),'color',clrs(5,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 5')
print -dpng o5.png
%
plot(o6(:,1),o6(:,3),'color',clrs(6,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 6')
print -dpng o6.png
%
plot(o7(:,1),o7(:,3),'color',clrs(7,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 7')
print -dpng o7.png
%
plot(o8(:,1),o8(:,3),'color',clrs(8,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 8')
print -dpng o8.png
%
plot(o9(:,1),o9(:,3),'color',clrs(9,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 9')
print -dpng o9.png
%
plot(o10(:,1),o10(:,3),'color',clrs(10,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('class 10')
print -dpng o10.png
%
plot(o1(:,1),o1(:,3),'color',clrs(1,:))
hold on
plot(o2(:,1),o2(:,3),'color',clrs(2,:))
plot(o3(:,1),o3(:,3),'color',clrs(3,:))
plot(o4(:,1),o4(:,3),'color',clrs(4,:))
plot(o5(:,1),o5(:,3),'color',clrs(5,:))
plot(o6(:,1),o6(:,3),'color',clrs(6,:))
plot(o7(:,1),o7(:,3),'color',clrs(7,:))
plot(o8(:,1),o8(:,3),'color',clrs(8,:))
plot(o9(:,1),o9(:,3),'color',clrs(9,:))
plot(o10(:,1),o10(:,3),'color',clrs(10,:))
xlim([140 202])
ylim([-42 49])
xlabel('R [cm]')
ylabel('Z [cm]')
title('all classes')
hold off
print -dpng oall.png
