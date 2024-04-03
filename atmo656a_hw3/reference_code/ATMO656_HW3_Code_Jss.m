cd('D:\work-school\school\FALL_2019\ATMO 656\HW3')
clear variables
close all

tau = (1:0.1:10)';
muo = [0.1 0.5 0.9];
wo = [0.9 1];
g = 0.75;

g1 = 0.25*(7-wo.*(4+3*g));
g2 = -0.25*(1-wo.*(4-3*g));
g3 = 0.25*(2-3*g*muo);
g4 = 1 - g3;
k = (g1.^2-g2.^2).^0.5;

a1 = g1.*g4'+g2.*g3';
a2 = g1.*g3'+g2.*g4';

R1 = zeros(length(tau),length(wo)*length(muo));
T1 = zeros(length(tau),length(wo)*length(muo));
str1 = string(zeros(6,1));
z = 1;
for i1 = 1:length(wo)
    
    for i2 = 1:length(muo)
        C1 = (wo(i1)*(((1-(k(i1)^2)*muo(i2)^2)*(k(i1)+g1(i1))*exp(k(i1)*tau)+(k(i1)-g1(i1))*exp(-k(i1)*tau))).^(-1));
        if i1 == 2
            R1(:,z) = (1./(1+g1(i1)*tau)).*(g1(i1)*tau+(g3(i2)-g1(i1)*muo(i2))*(1-exp(-tau/muo(i2))));
            T1(:,z) = 1-R1(:,z);
        else
            T1(:,z) = exp(-tau/muo(i2)).*(1-(C1).*(((1+k(i1)*muo(i2)).*(a1(i2,i1)+k(i1)*g4(i2)))*(exp(k(i1)*tau))-(1-k(i1)*muo(i2))*(a1(i2,i1)-k(i1)*g4(i2))*exp(-k(i1)*tau)-2*k(i1)*(g4(i2)+a1(i2,i1)*muo(i2)).*exp(tau/muo(i2))));
            R1(:,z) = C1.*(((1-k(i1)*muo(i2)).*(a2(i2,i1)+k(i1)*g3(i2)))*(exp(k(i1)*tau))-(1+k(i1)*muo(i2))*(a2(i2,i1)-k(i1)*g3(i2))*exp(-k(i1)*tau)-2*k(i1)*(g3(i2)-a2(i2,i1)*muo(i2)).*exp(-tau./muo(i2)));
        end
        
        str1(z,1) = "Eddington; ("+sprintf('%.1f',wo(i1))+","+sprintf('%.1f',muo(i2))+","+sprintf('%.2f',g)+")";
        z = z + 1;
    end
    
end

g1 = ((3^0.5)/2)*(2-wo.*(1+g));
g2 = ((3^0.5)/2)*(1-g).*wo;
g3 = 0.5*(1-(3^0.5)*g*muo);
g4 = 1 - g3;
k = (g1.^2-g2.^2).^0.5;

a1 = g1.*g4'+g2.*g3';
a2 = g1.*g3'+g2.*g4';

R2 = zeros(length(tau),length(wo)*length(muo));
T2 = zeros(length(tau),length(wo)*length(muo));
str2 = string(zeros(6,1));
z = 1;
for i1 = 1:length(wo)
    
    for i2 = 1:length(muo)
        C1 = (wo(i1)*(((1-(k(i1)^2)*muo(i2)^2)*(k(i1)+g1(i1))*exp(k(i1)*tau)+(k(i1)-g1(i1))*exp(-k(i1)*tau))).^(-1));
        if i1 == 2
            R2(:,z) = (1./(1+g1(i1)*tau)).*(g1(i1)*tau+(g3(i2)-g1(i1)*muo(i2))*(1-exp(-tau/muo(i2))));
            T2(:,z) = 1-R2(:,z);
        else
            T2(:,z) = exp(-tau/muo(i2)).*(1-(C1).*(((1+k(i1)*muo(i2)).*(a1(i2,i1)+k(i1)*g4(i2)))*(exp(k(i1)*tau))-(1-k(i1)*muo(i2))*(a1(i2,i1)-k(i1)*g4(i2))*exp(-k(i1)*tau)-2*k(i1)*(g4(i2)+a1(i2,i1)*muo(i2)).*exp(tau/muo(i2))));
            R2(:,z) = C1.*(((1-k(i1)*muo(i2)).*(a2(i2,i1)+k(i1)*g3(i2)))*(exp(k(i1)*tau))-(1+k(i1)*muo(i2))*(a2(i2,i1)-k(i1)*g3(i2))*exp(-k(i1)*tau)-2*k(i1)*(g3(i2)-a2(i2,i1)*muo(i2)).*exp(-tau./muo(i2)));
        end
        str2(z,1) = "Quadrature; ("+sprintf('%.1f',wo(i1))+","+sprintf('%.1f',muo(i2))+","+sprintf('%.2f',g)+")";
        z = z + 1;
    end
    
end

Fo = 1;
wo = [0.9 0.99999];
mu1 = [(1/3)^(1/2) -(1/3)^(1/2)];
w1 = 0;
k1 = (1/mu1(1)).*((1-wo).*(1-w1.*mu1(1).^2)).^(0.5);
z = 1;
R3 = zeros(length(tau),length(wo)*length(muo));
T3 = zeros(length(tau),length(wo)*length(muo));
str3 = string(zeros(6,1));
for i1 = 1:length(wo)
    
    for i2 = 1:length(muo)
        Z = zeros(1,2);
        W1 = zeros(1,2);
        for i3 = 1:size(mu1,2)
%              k1 = (1/mu1(i3)).*((1-wo)*(1-w1*mu1(i3).^2)).^(0.5);
            W1(1,i3) = ((1+mu1(i3)*k1(i1)).^(-1))*(wo(i1)-w1*(1-wo(i1)*mu1(i3)/k1(i1)));
            Z(1,i3) = muo(i2)*(mu1(i3)-muo(i2))*(wo(i1)-w1*(1-wo(i1))*mu1(i3)*muo(i2))*(4*(mu1(i3)^2)*(1-(k1(i1)^2)*muo(i2).^2)).^(-1);
        end
        
%         k1 = (1/mu1(1)).*((1-wo)*(1-w1*mu1(1).^2)).^(0.5);
        C(:,1) = -(Z(1,2)+Z(1,1)*exp(-tau/muo(i2)));
        C(:,2) = -(Z(1,2)-Z(1,1)*exp(-tau/muo(i2)));
        
        A(:,1) = W1(1,2)+W1(1,1)*exp(-k1(i1)*tau);
        A(:,2) = W1(1,2)-W1(1,1)*exp(-k1(i1)*tau);
        
        L(:,1) = 0.5*(C(:,1)./A(:,1)+C(:,2)./A(:,2));
        L(:,2) = 0.5*(C(:,1)./A(:,1)-C(:,2)./A(:,2)).*exp(-k1(i1)*tau);
        
        R3(:,z) = 2*(mu1(1)/muo(i2))*(L(:,1).*W1(1,1)+L(:,2).*W1(1,2)+Z(1,1));
        T3(:,z) = -(-2*(mu1(1)/muo(i2))*(L(:,1).*W1(1,2).*exp(-k1(i1)*tau)+L(:,2).*W1(1,1).*exp(k1(i1)*tau)+Z(1,2)*exp(-tau/muo(i2)))-exp(-tau/muo(i2)));
        
        str3(z,1) = "Liou; ("+sprintf('%.1f',wo(i1))+","+sprintf('%.1f',muo(i2))+","+sprintf('%.2f',g)+")";
        z =z + 1;
    end
    
end
ldgtxt = cat(1,str1,str2,str3);

figure;
set(gcf, 'Position',  [100, 0, 755, 700])
graph1 = subplot(2,1,1);
hold on
clrs = [0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0 0.4470 0.7410;...
    0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330];
for i = 1:size(R1,2)
    plot(tau,R1(:,i),'-','LineWidth',1,'Color',clrs(i,:))
end
for i = 1:size(R1,2)
    plot(tau,R2(:,i),'--','LineWidth',1,'Color',clrs(i,:))
end
for i = 1:size(R3,2)
    plot(tau,R3(:,i),'-.','LineWidth',1,'Color',clrs(i,:))
end
ylabel('R');
xlim([1 10])
set(gca,'XMinorTick','on','YMinorTick','on','Xticklabel',{''});
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1;
lgd = legend(ldgtxt);
lgd.Orientation = 'horizontal';
lgd.NumColumns = 3;
lgd.Position = [0.5 0.88 0.03 0.03];
lgd.FontSize = 12;
hold off

graph2 = subplot(2,1,2);
hold on
for i = 1:size(R1,2)
    plot(tau,T1(:,i),'-','LineWidth',1,'Color',clrs(i,:))
end
for i = 1:size(R1,2)
    plot(tau,T2(:,i),'--','LineWidth',1,'Color',clrs(i,:))
end
for i = 1:size(R3,2)
    plot(tau,T3(:,i),'-.','LineWidth',1,'Color',clrs(i,:))
end
xlim([1 10])
xlabel('\tau (m)');
ylabel('T');
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1;
hold off

graph1.Position(4) = graph1.Position(4)-0.02;
graph2.Position(4) = graph1.Position(4);
graph1.Position(3) = graph1.Position(3);
graph2.Position(3) = graph1.Position(3);
graph1.Position(2) = graph1.Position(2)-.12;
graph2.Position(2) = graph2.Position(2)-.02;

savefig("ATMO656_HW3_plots")
print("ATMO656_HW3_plots",'-djpeg')
close all


% 
% figure;
% set(gcf, 'Position',  [100, 0, 750, 700])
% graph1 = subplot(2,1,1);
% hold on
% clrs = [0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0 0.4470 0.7410;...
%     0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330];
% for i = 1:size(R3,2)
%     plot(tau,R3(:,i),'-','LineWidth',1,'Color',clrs(i,:))
% end
% ylabel('R');
% xlim([1 10])
% set(gca,'XMinorTick','on','YMinorTick','on','Xticklabel',{''});
% ax = gca;
% ax.FontSize = 14;
% ax.LineWidth = 1;
% lgd = legend(ldgtxt);
% lgd.Orientation = 'horizontal';
% lgd.NumColumns = 3;
% lgd.Position = [0.5 0.9 0.03 0.03];
% lgd.FontSize = 12;
% hold off
% 
% graph2 = subplot(2,1,2);
% hold on
% 
% for i = 1:size(R3,2)
%     plot(tau,T3(:,i),'-','LineWidth',1,'Color',clrs(i,:))
% end
% xlim([1 10])
% xlabel('\tau (m)');
% ylabel('T');
% set(gca,'XMinorTick','on','YMinorTick','on');
% ax = gca;
% ax.FontSize = 14;
% ax.LineWidth = 1;
% hold off
% 
% graph1.Position(4) = graph1.Position(4)-0.02;
% graph2.Position(4) = graph1.Position(4);
% graph1.Position(3) = graph1.Position(3);
% graph2.Position(3) = graph1.Position(3);
% graph1.Position(2) = graph1.Position(2)-.1;
% graph2.Position(2) = graph2.Position(2)-.02;
% 
% savefig("ATMO656_HW3_2_plots")
% print("ATMO656_HW3_2_plots",'-djpeg')
% close all