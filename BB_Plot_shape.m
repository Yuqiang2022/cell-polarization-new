
f11 = figure;
hold on

% r_C=r_C*[cos(-theta_xy),-sin(-theta_xy);sin(-theta_xy),cos(-theta_xy)];
% r_N=r_N*[cos(-theta_xy),-sin(-theta_xy);sin(-theta_xy),cos(-theta_xy)];

% scatter(r_C(:,1),r_C(:,2),'filled');
% scatter(r_N(:,1),r_N(:,2),'filled');
%plot interaction points
% for i = 1:1:Nv
%     if r_contact_doubleCell(i)==1
%         plot(r_C(i,1),r_C(i,2),'.k','MarkerSize',13);
%     end
% end
%plot fix points
for i = 1:1:Nv
    inter_plot = plot([r_C(i,1),r_N(i,1)],[r_C(i,2),r_N(i,2)],'r');
    set(inter_plot, 'LineWidth',1); 
end 

% for i = 1:1:Nv 
%     if r_fix_x(i)==0
%         plot(r_C(i,1),r_C(i,2),'.k','MarkerSize',20);
%     end
% end
a1 = Actin(:,1);
b1 = Actin(:,2);
V_actin=zeros(size(Actin,1),2); %vector of AFs
u_V_actin=zeros(size(Actin,1),2); %vector of AFs
for i = 1:1:size(Actin,1)
%     if rem(i,size(Actin,1)/N_cell) <= ActinNumber && rem(i,size(Actin,1)/N_cell)>0
    if i <= ActinNumber  %single cell    
        Actin_plot = plot([r_C(a1(i),1),r_C(b1(i),1)],[r_C(a1(i),2),r_C(b1(i),2)],'r');
        set(Actin_plot, 'LineWidth',10);
        V_actin(i,:)=r_C(a1(i),:)-r_C(b1(i),:);
        u_V_actin(i,:)=V_actin(i,:)./sqrt(V_actin(i,1).*V_actin(i,1)+V_actin(i,2).*V_actin(i,2));
    else
        Actin_plot = plot([r_N(a1(i),1),r_C(b1(i),1)],[r_N(a1(i),2),r_C(b1(i),2)],'r');
        set(Actin_plot, 'LineWidth',10);
        V_actin(i,:)=r_N(a1(i),:)-r_C(b1(i),:);
        u_V_actin(i,:)=V_actin(i,:)./sqrt(V_actin(i,1).*V_actin(i,1)+V_actin(i,2).*V_actin(i,2));        
    end
end  
% 
% a_MT=find(r_MT_n(:,1)>0);
% for i = 1:1:size(a_MT,1)
%     MT_plot = plot([r_N(ind_MT_N,1),r_C(a_MT(i),1)],[r_N(ind_MT_N,2),r_C(a_MT(i),2)],'g');
%     set(MT_plot,'LineWidth',2)
% end

%link strain color
    
cell_r_C0 = mat2cell(r_1,(Nv/N_cell.*ones(N_cell,1)));
cellfun(@(x) plot([x(:,1);x(1,1)],[x(:,2);x(1,2)],'r','LineWidth',1,'LineStyle',':'),cell_r_C0,'UniformOutput',false);
cell_r_N0 = mat2cell(rN_1,(Nv/N_cell.*ones(N_cell,1)));
cellfun(@(x) plot([x(:,1);x(1,1)],[x(:,2);x(1,2)],'r','LineWidth',1,'LineStyle',':'),cell_r_N0,'UniformOutput',false);
cell_r_C = mat2cell(r_C,(Nv/N_cell.*ones(N_cell,1)));
cellfun(@(x) plot([x(:,1);x(1,1)],[x(:,2);x(1,2)],'b','LineWidth',5),cell_r_C,'UniformOutput',false);
cell_r_N = mat2cell(r_N,(Nv/N_cell.*ones(N_cell,1)));
cellfun(@(x) plot([x(:,1);x(1,1)],[x(:,2);x(1,2)],'b','LineWidth',5),cell_r_N,'UniformOutput',false);

link_strain_0 =X;%X;%sum(Ft,2);%polarity;%./max(X);
link_strain_1 = link_strain_0;
link_strain_1(1:Nv-1,:) = link_strain_0(2:Nv,:);
link_strain_1(Nv,:) = link_strain_0(1,:);
link_strain = (link_strain_0+link_strain_1)./2;

MAP_LS = transpose(linspace(min(link_strain),1.0*max(link_strain),256)); 
colormap(jet) 
cmp = colormap;
    for i = 1:1:Nv 
        for j = 1:1:255 
            if link_strain(i,1) >= MAP_LS(j,1) && link_strain(i,1) <= MAP_LS(j + 1,1) 
                COL(i,:) = cmp(j,:); 
            end 
        end 
    end 

for i=1:Nv
    line([r_C(Links(i,1),1) r_C(Links(i,2),1)],[r_C(Links(i,1),2)
          r_C(Links(i,2),2)],'Color',COL(i,:),'LineWidth',5) ;
end
% for i=1:Nv
%     line([r_0(Links(i,1),1).*polarity(Links(i,1),1) r_0(Links(i,2),1).*polarity(Links(i,2),1)],[r_0(Links(i,1),2).*polarity(Links(i,1),1)
%           r_0(Links(i,2),2).*polarity(Links(i,2),1)],'Color',COL(i,:),'LineWidth',5) ;
% end
axis equal tight
if max(link_strain) > min(link_strain)
    caxis([min(link_strain) max(link_strain)]) ;
end
c = colorbar('position',[0.85 0.7 0.03 0.2]);
set(c,'fontsize',15);
set(gca,'fontname','Times New Roman','fontsize',20); box on;

x = r_C(:,1);
y = r_C(:,2);

Vector =F_trac./ sqrt(F_trac(:,1).*F_trac(:,1)+F_trac(:,2).*F_trac(:,2));
% u2 = 1.*S_polar.*Vector(:,1);%.*r_index;%.*r_fix_z;
% v2 = 1.*S_polar.*Vector(:,2);%.*r_index;%.*r_fix_z;
F_CC = F_substrate.*[vel_head vel_head];
u2 = 1.*F_CC(:,1);%.*r_index;%.*r_fix_z;
v2 = 1.*F_CC(:,2);%.*r_index;%.*r_fix_z;
quiver(x,y,u2,v2,'k')

hold off
axis equal on tight
xlabel('X')
ylabel('Y')

%%

% % 
DATE_hRBC_Parameters_simulation = datestr(now,'mm_dd_yyyy__HH_MM_SS'); %date and time of completion of parameter derivation 
fname = sprintf('Cell_%s',DATE_hRBC_Parameters_simulation);  
%print('-dtiff', fname, '-r600')% Export to TIFF file at 300 dpi.

% figure
% plot(U_TOT+U_TOT_N)
% figure
% plot(V_TOT+V_TOT_N)
