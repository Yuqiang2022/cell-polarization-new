clear all
close all
clc
%% cell and nucleus membrane
radius_C = 30e-6; %cell redius [m]
radius_N = radius_C/3; %nucleus redius [m]
r_C0 = radius_C;
r_N0 = radius_N;
Nv0 = 34;
Vertex_scale = 1;
Nv=Nv0*Vertex_scale; %particle number
Nv_single = Nv;
theta = (1/Nv:1/Nv:1)*2*pi;
[x_C,y_C] = pol2cart(theta,radius_C);
x_N = x_C./3;
y_N = y_C./3;
r_C = [x_C;y_C].';
r_N = [x_N;y_N].';
r_0 = r_C;
rN_0 = r_N;
%% spreading search
radius_spreading = 5e-9; %cell redius [m]
Nv_spreading = 20; %particle number
theta_spreading = (1/Nv_spreading:1/Nv_spreading:1)*2*pi;
[x_spreading,y_spreading] = pol2cart(theta_spreading,radius_spreading);
r_spreading = [x_spreading;y_spreading].';
%% Links
L_index = [1:1:Nv;2:1:Nv+1]';
L_index(Nv,2) = 1;
Links = L_index;
Ns = size(Links,1); %link number
r_edge = r_C(Links(:,1),:)-r_C(Links(:,2),:);
r_edge_N = r_N(Links(:,1),:)-r_N(Links(:,2),:);
LEN_0 = sqrt(sum(r_edge.*r_edge,2));
d0_C = mean(LEN_0); %mean length of cell edge [m]
avg_len = d0_C; % average length of all edges
LEN_N = sqrt(sum(r_edge_N.*r_edge_N,2));
d0_N = mean(LEN_N); %mean length of nucleus edge [m]
%% Element
for i = 1:1:Ns
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]
    a_S(i,:) = r_C(a1,:)-mean(r_C);
    b_S(i,:) = r_C(b1,:)-mean(r_C);
    S(i,:)=0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0]);
    A0_L(i) = S(i,3); %local area of cell
end
A0_G = sum(A0_L); %global area of cell

for i = 1:1:Ns
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]
    a_S(i,:) = r_N(a1,:)-mean(r_N);
    b_S(i,:) = r_N(b1,:)-mean(r_N);
    S(i,:)=0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0]);
    A0_L_N(i) = S(i,3); %local area of cell
end
A0_G_N = sum(A0_L_N); %global area of cell

%  DT = delaunayTriangulation(r_C(:,1),r_C(:,2));
%  Tri_element = DT.Points;
%  Tri_connect = DT.ConnectivityList;
%  N_tri = size(Tri_connect,1);
%  for i = 1:1:N_tri
%     p1 = Tri_connect(i,1);
%     p2 = Tri_connect(i,2);
%     p3 = Tri_connect(i,3);    
%     r_p1 = Tri_element(p1,:);
%     r_p2 = Tri_element(p2,:);
%     r_p3 = Tri_element(p3,:);   
%     a_S(i,:) = r_p1 - r_p2;
%     b_S(i,:) = r_p3 - r_p2;
%     S(i,:)=abs(0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0])); 
%     A0_L(i) = S(i,3); %local area of cell
% end
%  A0_G = sum(A0_L); %global area of cell

%% Bending
bending_pts = [[Nv,1:1:Nv-1]',[1:1:Nv]',[2:1:Nv,1]']; %bending points index
%% Actin
%ActinNumber = ceil(0.1*size(r_C,1)); %define generated number of actin filaments
ActinNumber = 0;
ActinNumber_N = ceil(0.0*size(r_C,1)); %define generated number of actin filaments connected with nucleus
Actin = zeros(ActinNumber+ActinNumber_N,2); %actin filaments ends index
ActinMaxLength_N = 1.1*d0_C;%sqrt(norm(r_C(1,:))^2 - norm(r_N(1,:))^2);
ActinMinLength_N = 0.8*d0_C;
% ActinMaxLength_N = 2.1*r_C0;
% ActinMinLength_N = 1.999*r_C0;
d0_actin = zeros(size(Actin,1),1);
N_actin = zeros(size(r_C,1),1); %record the nodes of actin using in the leading edge
N_actin_N = zeros(size(r_C,1),1); %record the nodes of actin using in the leading edge
i = 1;
while i <= ActinNumber+ActinNumber_N
    if i <= ActinNumber
        RandomIndex_1 = randi([1,size(r_C,1)],1,1); %Create a 1-by-1 array of random integer values drawn from a discrete uniform distribution on the set of numbers
        RandomIndex_2 = randi([1,size(r_C,1)],1,1);
        pts = [r_C(RandomIndex_1,:);r_C(RandomIndex_2,:)];
        CalDistence = sqrt((pts(1,1)-pts(2,1)).^2+(pts(1,2)-pts(2,2)).^2); % Calculate the length of the acitin filament
        if CalDistence<ActinMaxLength_N && CalDistence>ActinMinLength_N
            Actin(i,1) = RandomIndex_1;%Store actin end vertice
            Actin(i,2) = RandomIndex_2;%Store actin end vertice
            N_actin(RandomIndex_1)=N_actin(RandomIndex_1)+1;
            N_actin(RandomIndex_2)=N_actin(RandomIndex_2)+1;
            d0_actin(i) = CalDistence;
            i = i + 1;
        end
    else
        RandomIndex_1 = randi([1,size(r_N,1)],1,1); %Create a 1-by-1 array of random integer values drawn from a discrete uniform distribution on the set of numbers
        RandomIndex_2 = randi([1,size(r_C,1)],1,1);
        pts = [r_N(RandomIndex_1,:);r_C(RandomIndex_2,:)];
        CalDistence = sqrt((pts(1,1)-pts(2,1)).^2+(pts(1,2)-pts(2,2)).^2); % Calculate the length of the acitin filament
        if CalDistence<=ActinMaxLength_N % && pts(1,1)>0 && pts(2,1)>0
            Actin(i,1) = RandomIndex_1;%Store actin end vertice
            Actin(i,2) = RandomIndex_2;%Store actin end vertice
            N_actin_N(RandomIndex_1)=N_actin_N(RandomIndex_1)+1;
            N_actin(RandomIndex_2)=N_actin(RandomIndex_2)+1;
            d0_actin(i) = CalDistence;
            i = i + 1;
        end
    end
end

%manually set the actin index
 % Actin=[34 1; 2 3; 17 18;19 20]; %#3
 % Actin=[34 1;5 6;17 18;22 23]; %#2
 % Actin=[34 1; 8 9; 17 18; 25 26]; %#1

 Actin=[2 3]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1); %#0
 %Actin=[34 1];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
 %Actin=[34 1;17 18];ActinNumber = 2;d0_actin=LEN_0(1).*ones(ActinNumber,1);%double actins
 %Actin=[34 1;11 12;23 24];ActinNumber = 3;d0_actin=LEN_0(1).*ones(ActinNumber,1);%trible actins
% Actin=[34 1; 8 9; 17 18; 25 26]; ActinNumber = 4;d0_actin=LEN_0(1).*ones(ActinNumber,1);%four actins
 %Actin=[34 1;4 5;8 9;12 13;17 18;21 22;25 26;29 30]; ActinNumber = 8;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 

% the direction of actin alignment
% Actin=[33 1;4 5;8 9;12 13;16 18;21 22;25 26;29 30]; ActinNumber = 8;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 
% Actin=[33 1]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 
% Actin=[4 5]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 
% Actin=[8 9]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 
% Actin=[12 13]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 
% Actin=[16 18]; ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%eight actins 

% Actin=[34 1];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[1 2];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[2 3];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[3 4];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[4 5];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[5 6];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[6 7];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[7 8];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[8 9];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[9 10];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[10 11];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[11 12];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[12 13];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[13 14];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[14 15];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[15 16];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[16 17];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[17 18];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[18 19];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[19 20];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[20 21];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[21 22];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[22 23];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[23 24];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[24 25];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[25 26];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[26 27];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[27 28];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[28 29];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[29 30];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[30 31];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[31 32];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[32 33];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin
% Actin=[33 34];ActinNumber = 1;d0_actin=LEN_0(1).*ones(ActinNumber,1);%single actin

% %%load('Fig2A-','Actin')
% %isotropic alignment of cortical actin
% Actin_1 = (1:2:Nv)';
% Actin_2 = Actin_1;
% Actin_2(1:Nv/2-1)=Actin_1(2:Nv/2,1);
% Actin_2(Nv/2)=Actin_1(1);
% Actin_11=[Actin_1,Actin_2];
% Actin_22=[Actin_1+1,Actin_2+1];
% Actin=[Actin_11;Actin_22];
% d0_actin(1:ActinNumber)=sqrt((r_C(1,1)-r_C(3,1))^2+(r_C(1,2)-r_C(3,2))^2);
% N_actin(1:ActinNumber)=2;


%% Intermediate filaments
d0_inter = sqrt(sum((r_C-r_N).*(r_C-r_N),2)); %individual intermediate filaments [m]
%% Microtubule
MT_N = ceil(0*size(r_C,1)); %define generated number of microtubules
ind_MT_N = 1;
r_MT = randi([1,size(r_C,1)],1,MT_N)';
v0_MT = ones(Nv,1)*r_N(ind_MT_N,:)-r_C;
d0_MT = sqrt(sum(v0_MT.*v0_MT,2)); % initial length of MT
%% doubling cells
N_cell = 1;
    bending_pts = repmat(bending_pts,N_cell,1) + kron([0:N_cell-1]',bending_pts./bending_pts).*Nv_single;
    Links = repmat(Links,N_cell,1) + kron([0:N_cell-1]',Links./Links).*Nv_single;
    LEN =repmat(LEN_0,N_cell,1);
    LEN_N = repmat(LEN_N,N_cell,1);
    d0_inter = repmat(d0_inter,N_cell,1);
    Actin = repmat(Actin,N_cell,1) + kron((0:N_cell-1)',Actin./Actin).*Nv_single;
    d0_actin = repmat(d0_actin,N_cell,1);
    N_actin = repmat(N_actin,N_cell,1,1);
    ActinNumber = ActinNumber;
    avg_len = avg_len;
    r_C0 = r_C0;
    r_N0 = r_N0;
    r_0 = repmat(r_0,N_cell,1);
    rN_0 = repmat(rN_0,N_cell,1);
    A0_L = repmat(A0_L',N_cell,1);
    A0_L_N = repmat(A0_L_N',N_cell,1);
    r_spreading = r_spreading;
    bias_r = zeros(N_cell,2);
    ii=2;
    lowx=1;upx=40; lowy=1;upy=40; %[um]
    while ii <= N_cell
        tempx=randi([lowx,upx],1)*1e-6;
        tempy=randi([lowy,upy],1)*1e-6;
        if sum(sqrt((tempx-bias_r(:,1)).^2+(tempy-bias_r(:,2)).^2)<2.2*r_C0)==0 %&&sum(sqrt((tempx-bias_r(:,1)).^2+(tempy-bias_r(:,2)).^2)>3.2*r_C0)==0
            bias_r(ii,1)=tempx;
            bias_r(ii,2)=tempy;
            ii=ii+1;
        end
    end 
    N_rr = mat2cell(bias_r,(ones(N_cell,1)));
    NN_rr = cell2mat(cellfun(@(x) repmat(x,Nv_single,1),N_rr,'UniformOutput',false));       
    rr0 = repmat(r_C,N_cell,1)+NN_rr;
    rr0_N = repmat(r_N,N_cell,1)+NN_rr;
    r_C = rr0;
    r_N = rr0_N;
    Ns = N_cell*Ns;
    Nv = N_cell*Nv_single;
%% Save framework:
DATE_hRBC_Triangulation = datestr(now,'mm_dd_yyyy__HH_MM_SS'); % date and time of simulation completion
save('CELL_FRAMEWORK.mat')
%% Plot
scatter(r_C(:,1),r_C(:,2),'filled');
hold on
scatter(r_N(:,1),r_N(:,2),'filled');

for i = 1:1:size(r_C,1)
    inter_plot = plot([r_C(i,1),r_N(i,1)],[r_C(i,2),r_N(i,2)],'b');
    set(inter_plot, 'LineWidth',0.05);
end

a1 = Actin(:,1);
b1 = Actin(:,2);

for i = 1:1:size(Actin,1)
%     if rem(i,size(Actin,1)/N_cell) <= ActinNumber && rem(i,size(Actin,1)/N_cell)>0
    if i <= ActinNumber  %single cell    
        Actin_plot = plot([r_C(a1(i),1),r_C(b1(i),1)],[r_C(a1(i),2),r_C(b1(i),2)],'r');
        set(Actin_plot, 'LineWidth',1);
    else
        Actin_plot = plot([r_N(a1(i),1),r_C(b1(i),1)],[r_N(a1(i),2),r_C(b1(i),2)],'r');
        set(Actin_plot, 'LineWidth',1);
    end
end

for i = 1:1:size(r_C,1)
    inter_plot = plot([r_C(i,1),r_N(i,1)],[r_C(i,2),r_N(i,2)],'b');
    set(inter_plot, 'LineWidth',0.05);
end

for i = 1:1:size(r_MT,1)
    MT_plot = plot([r_N(ind_MT_N,1),r_C(r_MT(i),1)],[r_N(ind_MT_N,2),r_C(r_MT(i),2)],'g');
    set(MT_plot,'LineWidth',2)
end

hold off
axis equal tight
