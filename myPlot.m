function myPlot(DATA, net, figName)

w = net.weight;
edge = net.edge;
[N,D] = size(w);

% connection = graph(edge ~= 0); % If matlab version is 2015b above, this line can be used.
% label = conncomp(connection);

label = net.LabelCluster;


if D==2
    plot(DATA(:,1),DATA(:,2),'cy.','MarkerSize',1);
elseif D==3
    scatter3(DATA(:,1),DATA(:,2),DATA(:,3),'cy.','MarkerFaceColor','cy');
end

hold on;

for i=1:N-1
    for j=i:N
        if edge(i,j)==1
            if D==2
                plot([w(i,1) w(j,1)],[w(i,2) w(j,2)],'w','LineWidth',1.5);
            elseif D==3
                plot3([w(i,1) w(j,1)],[w(i,2) w(j,2)],[w(i,3) w(j,3)],'w','LineWidth',1.5);
            end
        end
    end
end



% Change Node color based on LebelCluster
color = [
    [1 0 0]; 
    [0 1 0]; 
    [0 0 1]; 
%     [0 1 1]; 
    [1 0 1];
%     [1 1 0];
%     [0 0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.9290 0.6940 0.1250];
    [0.4940 0.1840 0.5560];
    [0.4660 0.6740 0.1880];
    [0.3010 0.7450 0.9330];
    [0.6350 0.0780 0.1840];
%     [1 1 1];
];
m = length(color);

for k = 1:N
    if D==2
        plot(w(k,1),w(k,2),'.','Color',color(mod(label(1,k)-1,m)+1,:),'MarkerSize',35);
    elseif D==3
        plot3(w(k,1),w(k,2),w(k,3),'.','Color',color(mod(label(1,k)-1,m)+1,:),'MarkerSize',35);
    end
end




for i=1:N
    str = num2str(i); % Node Number
%     str = num2str(LebelCluster(1,i)); % Node Label
%     str = num2str(net.CountCluster(i)/max(TKBAnet.CountCluster));  % Norm Count Cluster
%     str = num2str(net.CountCluster(i));  % Count Cluster
%     str = num2str(net.adaptiveSig(1,i));
%     text(w(i,1)+0.01,w(i,2)+0.01, str,'Color','y','FontSize',14);
end




axis equal
grid on
hold off
axis([0 1 0 1]);
ytickformat('%.1f')
xtickformat('%.1f')
% axis([0 10 0 10]);
title(figName,'fontsize',12);
pause(0.01);

end