function [F_C,f_links_tensile,f_actin,f_area_g,U_link,U_bending,U_area_loc,U_area_glo,U_inter,U_potential,U_actin,LEN,f_intermediate_C,f_intermediate_N,f_potential_N,f_actin_N]...
    = f_cons_forces_mex(Links,bending_pts,r_C,r_N,Actin,ActinNumber,r_p0,rN_p0,Vertex_scale)

% Declare Global Variables and Initialize Terms:
global A0 A0t A0t_N kB ka kb kp kp_inter Lm Lp_actin m Ns Nv Nv_single N_cell T theta0 Lm_inter Lp_inter Lm_actinfilamet Lp_actinfilamet

f_links = zeros(Nv,2); % initialize in-plane force array [N]
f_links_tensile = zeros(Nv,2); % initialize in-plane force array [N]
f_area_loc = zeros(Nv,3); % initialize local area force array [N]
f_area_g = zeros(Nv,3); % initialize global area force array [N]
f_bending = zeros(Nv,2); % initialize bending force array [N]
LEN = zeros(Ns,1); % initialize link length vector [m]
S = zeros(Ns,3); % initialize element-wise area vector [m^2]
A_L = zeros(Ns,1);
A_L_N = zeros(Ns,1);
U_link = 0; % initialize links energy term [J]
U_bending = 0; % initialize bending energy term [J]
U_area_loc = 0; % initialize local area energy term [J]
U_inter = 0; % initialize the interspace potential energy term [J]
U_potential = 0;
U_actin = 0;
a_S = zeros(Nv,2);
b_S = zeros(Nv,2);
c_S = zeros(Nv,2);
LEN_inter = zeros(Nv,2); %initialize length of intermediate filaments
f_intermediate_C = zeros(Nv,2); %initialize force of intermediate filaments on cytomembrane
f_intermediate_N = zeros(Nv,2); %initialize force of intermediate filaments on nucleus
f_potential_C = zeros(Nv,2);
f_potential_N = zeros(Nv,2);
LEN_actin = zeros(length(Actin),2); %initialize length of actin filaments
f_actin_1 = zeros(Nv,2); %initialize force of actin filaments on one end
f_actin_2 = zeros(Nv,2); %initialize force of actin filaments on one end
f_actin_3 = zeros(Nv,2); %initialize force of actin filaments on one end
f_actin_4 = zeros(Nv,2); %initialize force of actin filaments on one end
% Calculate Link, Bending and local Area Forces and Energies:
Center_v = mat2cell(r_C,(Nv_single.*ones(N_cell,1)));
Center_mean = cellfun(@(x) mean(x),Center_v,'UniformOutput',false);
for i = 1:1:Ns
    % In-Plane (2-point interaction):
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]
    x_a1 = r_C(a1,1); % x-coordinate of node 1 [m]
    y_a1 = r_C(a1,2); % y-coordinate of node 1 [m]
    x_b1 = r_C(b1,1); % x-coordinate of node 2 [m]
    y_b1 = r_C(b1,2); % y-coordinate of node 2 [m]
    r_21_1 = [(x_b1 - x_a1) (y_b1 - y_a1)]; % relative position vector (from node 1 to node 2) [m]
    len = sqrt(r_21_1(1,1)*r_21_1(1,1) + r_21_1(1,2)*r_21_1(1,2)); % distance between node 1 and node 2 [m]
    e_21_1 = r_21_1./len; % unit relative position vector (from node 1 to node 2) [m]
    LEN(i,1) = len; % save length of current link [m]
    xn = len/Lm; % normalized length of spring extension [unitless]
    wlc_len =(-kB*T/Lp_actin)*(1/(4*(1 - xn)^2) - (1/4) + xn);   
    f_wlc_1 = wlc_len.*(-e_21_1); % WLC contribution to force on node 1 [N]
    f_wlc_2 = -f_wlc_1; % WLC contribution to force on node 2 [N]
    kp_len = kp/(len^m);
    f_pow_1 =  kp_len.*(-e_21_1); % POW contribution to force on node 1 [N]
    f_pow_2 = -f_pow_1; % POW contribution to force on node 2 [N]
    f_links(a1,:) = f_pow_1 + f_wlc_1+ f_links(a1,:);
    f_links(b1,:) = f_pow_2 + f_wlc_2+ f_links(b1,:);
    f_links_tensile(a1,:) = f_wlc_1 + f_links_tensile(a1,:);
    f_links_tensile(b1,:) = f_wlc_1 + f_links_tensile(b1,:);

    U_wlc = (kB*T*Lm/(4*Lp_actin))*((3*xn^2 - 2*xn^3)/(1 - xn)); % WLC contribution to in-plane potential energy [J]
    U_pow = kp/((m - 1)*len^(m - 1)); % POW contribution to in-plane potential energy [J]
    U_link = U_link + U_wlc + U_pow;% + U_lj*h2*CON(i,1); % total in-plane potential energy [J]
    % Bending (2-point interaction):    
    p1 = bending_pts(i,1);
    p2 = bending_pts(i,2);
    p3 = bending_pts(i,3);
    r_p1 = r_C(p1,:);
    r_p2 = r_C(p2,:);    
    r_p3 = r_C(p3,:);    
    a_v = r_p1 - r_p2;
    b_v = r_p3 - r_p2;
    a_m = norm(a_v);
    b_m = norm(b_v);
    cosine = sum((a_v/a_m).*(b_v/b_m));
    vector_1 = a_v*[cos(-pi/2) sin(-pi/2);-sin(-pi/2) cos(-pi/2)]; %rotation matrix;
    vector_2 = b_v*[cos(pi/2) sin(pi/2);-sin(pi/2) cos(pi/2)]; %rotation matrix;       
    tc_1 = (r_p1 + r_p2)./2;
    tc_2 = (r_p2 + r_p3)./2;
    check = sum((vector_1 - vector_2).*(tc_1 - tc_2));    
    if cosine > 0.999999500000042 %0.999999500000042
        cosine = 0.999999500000042;
    elseif cosine < -0.999999500000042        
        cosine = -0.999999500000042;
    end   
    sine = sqrt(1 - cosine^2);    
    if check > 0
        sine=-sine;
    end
    
%     if abs(cosine) > 999999500000042 && cosine < 0
%         sine = -sqrt(1 - cosine^2);
%     end       
    beta_bending = kb*(sine*cos(theta0) - cosine*sin(theta0))/sine; % bending force coefficient     *(1-cosine*cos(theta0)-sine*sin(theta0))
    %    beta_bending = kb*(cosine-1); % polygon bending model
    f_bending_1 = -beta_bending*(1/a_m/b_m.*b_v-cosine/a_m^2.*a_v); % bending force 1 [N]
    f_bending_2 = -beta_bending*((cosine/b_m^2 - 1/a_m/b_m).*b_v + (cosine/a_m^2 - 1/a_m/b_m).*a_v); % bending force 2 [N]
    f_bending_3 = -beta_bending*(1/a_m/b_m.*a_v-cosine/b_m^2.*b_v); %bending force 3 [N]    
    f_bending(p1,:) = -f_bending_1 + f_bending(p1,:); % bending force on node 1 [N]
    f_bending(p2,:) = -f_bending_2 + f_bending(p2,:); % bending force on node 2 [N]
    f_bending(p3,:) = -f_bending_3 + f_bending(p3,:); % bending force on node 3 [N]    
    v_1 = cosine; % sum(ksi_1.*zeta_1)/(ksi_1_ns*zeta_1_ns);
    if v_1 <= - 1
        v_2 = -1;
    elseif v_1 >= 1
        v_2 = 1;
    else
        v_2 = v_1;
    end
    theta = acos(v_2)/2; %angle between elements [rad]
    U_bending = U_bending + kb*(1 - cos(theta - theta0)); % total bending potential energy [J]        
%     %polygon form of bending force (Feng's paper)  
%     i_2 = i-2;
%     if i_2 == -1
%         i_2 = Ns - 1;
%     end
%     if i_2 == 0
%         i_2 = Ns;
%     end    
%     i_1 = i-1;
%     if i_1 == 0
%         i_1 = Ns;
%     end
%     i_11 = i+1;
%     if i_11 == Ns+1
%         i_11 = 1;
%     end
%     i_22 = i+2;
%     if i_22 == Ns+1
%         i_22 = 1;
%     end    
%     if i_22 == Ns+2
%         i_22 = 2;
%     end   
%     cos_i = sum((r_C(i_11,:)-r_C(i,:)).*(r_C(i,:)-r_C(i_1,:)))./norm(r_C(i_11,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:));
%     cos_i_1 = sum((r_C(i,:)-r_C(i_1,:)).*(r_C(i_1,:)-r_C(i_2,:)))./norm(r_C(i,:)-r_C(i_1,:))./norm(r_C(i_2,:)-r_C(i_1,:));
%     cos_i_11 = sum((r_C(i_22,:)-r_C(i_11,:)).*(r_C(i_11,:)-r_C(i,:)))./norm(r_C(i_22,:)-r_C(i_11,:))./norm(r_C(i,:)-r_C(i_11,:));   
%     f_bending(i,:)= kb.*(cos_i-1).*((r_C(i_11,:)+r_C(i_1,:)-2.*r_C(i,:))./norm(r_C(i_11,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:))...
%                     + cos_i.*((r_C(i_11,:)-r_C(i,:))./norm(r_C(i_11,:)-r_C(i,:))./norm(r_C(i_11,:)-r_C(i,:))...
%                     +(r_C(i_1,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:))))...
%                     + kb.*(cos_i_1-1).*((r_C(i_1,:)-r_C(i_2,:))./norm(r_C(i,:)-r_C(i_1,:))./norm(r_C(i_2,:)-r_C(i_1,:))...
%                     + cos_i_1.*(r_C(i_1,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:))./norm(r_C(i_1,:)-r_C(i,:)))...
%                     + kb.*(cos_i_11-1).*((r_C(i_11,:)-r_C(i_22,:))./norm(r_C(i_22,:)-r_C(i_11,:))./norm(r_C(i_11,:)-r_C(i,:))...
%                     + cos_i_11.*(r_C(i_11,:)-r_C(i,:))./norm(r_C(i_11,:)-r_C(i,:))./norm(r_C(i_11,:)-r_C(i,:)));
    
    %Area (2-point interaction);
    a_S(i,:) = r_C(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    b_S(i,:) = r_C(b1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    S(i,:)=0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0]);
    A_L(i) = S(i,3); %local area of cell
    %Area (2-point interaction);
    a_S(i,:) = r_N(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    b_S(i,:) = r_N(b1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    S(i,:)=0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0]);
    A_L_N(i) = S(i,3); %local area of cell
end
% global Area:
% AL_cell = mat2cell(A_L,(Ns/N_cell.*ones(N_cell,1)));
% AL_sum = cellfun(@(x) sum(x(:,1)),AL_cell,'UniformOutput',false);
% A_G = [AL_sum{:,1}]';
A_G = sum(A_L);
A_G_N = sum(A_L_N);
beta_area = -ka*(A_G-A_G_N - A0t + A0t_N)/(A0t-A0t_N); % coefficient in global area force calculation
for i = 1:1:Ns
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]
    a_S(i,:) = r_C(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    b_S(i,:) = -r_C(b1,:)+Center_mean{ceil(i/(Ns/N_cell))};
    c_S(i,:) = r_C(a1,:)-r_C(b1,:);
    ksi_2 = cross([a_S(i,:) 0],[b_S(i,:) 0]);        
    alpha_area_g = beta_area(ceil(i/(Ns/N_cell)))/(4*A_L(i,1)); % global area force coefficient   
    f_area_g_1 = alpha_area_g.*cross(ksi_2,[b_S(i,:) 0]); % global area force 1 [N]
    f_area_g_2 = alpha_area_g.*cross(ksi_2,[a_S(i,:) 0]); % global area force 2 [N]    
    f_area_g_3 = -alpha_area_g.*cross(ksi_2,[c_S(i,:) 0]); % global area force 2 [N]   
    f_area_g(a1,:) = -f_area_g_1 + f_area_g(a1,:);% + f_area_g_3./2; % global area force on node 1 [N]
    f_area_g(b1,:) = -f_area_g_2 + f_area_g(b1,:);% + f_area_g_3./2; % global area force on node 2 [N]
end
U_area_glo = sum(ka.*((A_G - A0t).^2)./(2.*A0t)); %area induced energy change

% % local Area:
% for i = 1:1:Ns
%     a1 = Links(i,1); % index of node 1 [unitless]
%     b1 = Links(i,2); % index of node 2 [unitless]
%     a_S(i,:) = r_C(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
%     b_S(i,:) = -r_C(b1,:)+Center_mean{ceil(i/(Ns/N_cell))};
%     c_S(i,:) = r_C(a1,:)-r_C(b1,:);
%     ksi_2 = cross([a_S(i,:) 0],[b_S(i,:) 0]);        
%     beta_area = -ka*(A_L - A0(i))/A0(i); % coefficient in global area force calculation    
%     alpha_area_l = beta_area(ceil(i/(Ns/N_cell)))/(4*A_L(i,1)); % global area force coefficient   
%     f_area_l_1 = alpha_area_l.*cross(ksi_2,[b_S(i,:) 0]); % global area force 1 [N]
%     f_area_l_2 = alpha_area_l.*cross(ksi_2,[a_S(i,:) 0]); % global area force 2 [N]    
%     f_area_l_3 = -alpha_area_l.*cross(ksi_2,[c_S(i,:) 0]); % global area force 2 [N]   
%     f_area_loc(a1,:) = f_area_g_3./2; %-f_area_l_1 + f_area_loc(a1,:);% + f_area_g_3./2; % global area force on node 1 [N]
%     f_area_loc(b1,:) = f_area_g_3./2; %-f_area_l_2 + f_area_loc(b1,:);% + f_area_g_3./2; % global area force on node 2 [N]   
% end
% U_area_loc = sum(ka.*((A_L - A0).^2)./(2.*A0)); %area induced energy change

% cytommembrane and nucleus interaction (intermediate filaments) (2-point interaction):
for i = 1:1:Nv
    a1 = r_C(i,:); % index of node 1 on cytomembrane [unitless]
    b1 = r_N(i,:); % index of node 2 on nucleus [unitless]
    x_a1 = a1(1); % x-coordinate of node 1 [m]
    y_a1 = a1(2); % y-coordinate of node 1 [m]
    x_b1 = b1(1); % x-coordinate of node 2 [m]
    y_b1 = b1(2); % y-coordinate of node 2 [m]
    r_21_1 = [(x_b1 - x_a1) (y_b1 - y_a1)]; % relative position vector (from node 1 to node 2) [m]
    len = sqrt(r_21_1(1,1)*r_21_1(1,1) + r_21_1(1,2)*r_21_1(1,2)); % distance between node 1 and node 2 [m]
    e_21_1 = r_21_1./len; % unit relative position vector (from node 1 to node 2) [m]
    LEN_inter(i,1) = len; % save length of current link [m]
    xn = len/Lm_inter(i,1); % normalized length of spring extension [unitless]
    f_wlc_1 = (-kB*T/Lp_inter(i,1))*(1/(4*(1 - xn)^2) - (1/4) + xn).*(-e_21_1); % WLC contribution to force on node 1 [N]
    f_wlc_2 = -f_wlc_1; % WLC contribution to force on node 2 [N]
    f_intermediate_C(i,:) = f_wlc_1 + f_intermediate_C(i,:); % total in-plane force on node 1 [N]
    f_intermediate_N(i,:) = f_wlc_2 + f_intermediate_N(i,:); % total in-plane force on node 2 [N]
    % calculate potential energy between cytomembrane and nucleus
    % cytomembrane    
%     cell_r_N = mat2cell(r_N,(Nv/N_cell.*ones(N_cell,1)));    
%     cell_r_N_0 = mat2cell(rN_p0,(Nv/N_cell.*ones(N_cell,1)));       
%    len_potential_1 = (r_C(i,:) - cell_r_N{ceil(i/(Nv/N_cell))});
    len_potential_1 = (r_C(i,:) - r_N);
    len_p_1 = sqrt(sum(len_potential_1.*len_potential_1,2));
%    len_potential_1_0 = (r_p0(i,:) - cell_r_N_0{ceil(i/(Nv/N_cell))});
    len_potential_1_0 = (r_p0(i,:) - rN_p0);     
    len_p_1_0 = sqrt(sum(len_potential_1_0.*len_potential_1_0,2));
    e_p_1 = len_potential_1./(len_p_1*ones(1,2));
%    f_pow_1_p = kp_inter.*(12.*len_p_1_0.^12./len_p_1.^13-6.*len_p_1_0.^6./len_p_1.^7)*ones(1,2).*e_p_1; % 【L-J potential】
    len_p_1_v = (len_p_1 - len_p_1_0)./len_p_1_0;
    p_matrix_1 = (len_p_1_v(:,1)<0);
    f_pow_1_p = ((pi/2*kp_inter.*(1 - cos(pi/2.*len_p_1_v).^-2))./len_p_1_0.*p_matrix_1).*(-e_p_1); %【tan】
    f_pow_1(i,:) = sum(f_pow_1_p);
    f_potential_N((ceil(i/(Nv/N_cell))*Nv_single-Nv_single+1):(ceil(i/(Nv/N_cell))*Nv_single),:) = -f_pow_1_p + f_potential_N((ceil(i/(Nv/N_cell))*Nv_single-Nv_single+1):(ceil(i/(Nv/N_cell))*Nv_single),:); 
    f_potential_C(i,:) = f_pow_1(i,:) + f_potential_C(i,:);% the nucleus contribute no force to the contact points
%    U_pow_inter = sum(kp_inter.*(len_p_1_0.^12./len_p_1.^12-len_p_1_0.^6./len_p_1.^6));
    U_pow_inter = kp_inter*sum((pi/2.*len_p_1_v - tan(pi/2.*len_p_1_v)).*p_matrix_1); %potential energy (there is no potential energy between the nucleus and the contact points)    
    U_wlc_inter = (kB*T*Lm_inter(i,1)/(4*Lp_inter(i,1)))*((3*xn^2 - 2*xn^3)/(1 - xn)); % WLC contribution to in-plane potential energy [J]
    U_inter = U_inter + U_wlc_inter;% intermediate filaments energy [J]
    U_potential = U_potential + U_pow_inter;% potential energy [J]
end

% Actin filaments force (2-point interaction):
for i = 1:1:size(Actin,1)
%    if rem(i,size(Actin,1)/N_cell) <= ActinNumber && rem(i,size(Actin,1)/N_cell)>0
    if i <= ActinNumber  %single cell        
        a1 = Actin(i,1); % index of node 1 [unitless]
        b1 = Actin(i,2); % index of node 2 [unitless]
        x_a1 = r_C(a1,1); % x-coordinate of node 1 [m]
        y_a1 = r_C(a1,2); % y-coordinate of node 1 [m]
        x_b1 = r_C(b1,1); % x-coordinate of node 2 [m]
        y_b1 = r_C(b1,2); % y-coordinate of node 2 [m]
        r_21_1 = [(x_b1 - x_a1) (y_b1 - y_a1)]; % relative position vector (from node 1 to node 2) [m]
        len = sqrt(r_21_1(1,1)*r_21_1(1,1) + r_21_1(1,2)*r_21_1(1,2)); % distance between node 1 and node 2 [m]
        e_21_1 = r_21_1./len; % unit relative position vector (from node 1 to node 2) [m]
        LEN_actin(i,1) = len; % save length of current link [m]
        xn = len/Lm_actinfilamet(i); % normalized length of spring extension [unitless]
        f_wlc_1 = ((-kB*T/Lp_actinfilamet(i,1))*(1/(4*(1 - xn)^2) - (1/4) + xn)).*(-e_21_1); %+ kp_actinfilament(i)/len^m).*(-e_21_1); % WLC contribution to force on node 1 [N]
        f_wlc_2 = -f_wlc_1; % WLC contribution to force on node 2 [N]
        f_actin_1(Actin(i,1),:) = f_wlc_1 + f_actin_1(Actin(i,1),:); % total actin filaments force on node 1 [N]
        f_actin_2(Actin(i,2),:) = f_wlc_2 + f_actin_2(Actin(i,2),:); % total actin filaments force on node 2 [N]
        U_wlc_actin = (kB*T*Lm_actinfilamet(i)/(4*Lp_actinfilamet(i,1)))*(3*xn^2 - 2*xn^3)/(1 - xn);% + kp_actinfilament(i)/((m - 1)*len^(m - 1));% WLC contribution to in-plane potential energy [J]
        U_actin = U_actin + U_wlc_actin; % total actin filaments potential energy [J]
    else
        a1 = Actin(i,1); % index of node 1 [unitless]
        b1 = Actin(i,2); % index of node 2 [unitless]
        x_a1 = r_N(a1,1); % x-coordinate of node 1 [m]
        y_a1 = r_N(a1,2); % y-coordinate of node 1 [m]
        x_b1 = r_C(b1,1); % x-coordinate of node 2 [m]
        y_b1 = r_C(b1,2); % y-coordinate of node 2 [m]
        r_21_1 = [(x_b1 - x_a1) (y_b1 - y_a1)]; % relative position vector (from node 1 to node 2) [m]
        len = sqrt(r_21_1(1,1)*r_21_1(1,1) + r_21_1(1,2)*r_21_1(1,2)); % distance between node 1 and node 2 [m]
        e_21_1 = r_21_1./len; % unit relative position vector (from node 1 to node 2) [m]
        LEN_actin(i,1) = len; % save length of current link [m]
        xn = len/Lm_actinfilamet(i); % normalized length of spring extension [unitless]
        f_wlc_3 = (((-kB*T/Lp_actinfilamet(i,1))*(1/(4*(1 - xn)^2) - (1/4) + xn))).*(-e_21_1);%+ kp_actinfilament(i)/len^m).*(-e_21_1); % WLC contribution to force on node 1 [N]
        f_wlc_4 = -f_wlc_3; % WLC contribution to force on node 2 [N]
        f_actin_3(Actin(i,1),:) = f_wlc_3 + f_actin_3(Actin(i,1),:); % total actin filaments force on node 1 [N]
        f_actin_4(Actin(i,2),:) = f_wlc_4 + f_actin_4(Actin(i,2),:); % total actin filaments force on node 2 [N]
        U_wlc_actin = (kB*T*Lm_actinfilamet(i)/(4*Lp_actinfilamet(i,1)))*(3*xn^2 - 2*xn^3)/(1 - xn);% + kp_actinfilament(i)/((m - 1)*len^(m - 1));% WLC contribution to in-plane potential energy [J]
        U_actin = U_actin + U_wlc_actin; % total actin filaments potential energy [J]
    end
end
f_actin = f_actin_1 + f_actin_2 + f_actin_4;
f_actin_N = f_actin_3;
f_intermediate_C=f_intermediate_C./Vertex_scale; %coarse-grained force of intermediate fimalent 
f_intermediate_N=f_intermediate_N./Vertex_scale; %coarse-grained force of intermediate fimalent 
F_C = f_links  +  f_area_g(:,1:2) +  f_area_loc(:,1:2) + f_bending + f_intermediate_C + f_potential_C + f_actin; %total conservative force on all nodes [N]
end