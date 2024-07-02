function [F_C_N,U_link_N,U_bending_N,U_area_loc_N,U_area_glo_N,LEN_N,A_G_N] = f_cons_forces_mex_N(Links,bending_pts,r_N,f_intermediate_N,f_potential_N,f_actin_N)
% Declare Global Variables and Initialize Terms:
global A0_N A0t_N kB kb kp_N Lm Lm_N Lp_nucleus m Ns Nv Nv_single N_cell T theta0 ka_N
f_links = zeros(Nv,2); % initialize in-plane force array [N]
f_area_loc = zeros(Nv,3); % initialize local area force array [N]
f_area_g = zeros(Nv,3); % initialize global area force array [N]
f_bending = zeros(Nv,2); % initialize bending force array [N]
LEN_N = zeros(Ns,1); % initialize link length vector [m]
S = zeros(Ns,3); % initialize element-wise area vector [m^2]
A_L = zeros(Ns,1);
U_link_N = 0; % initialize links energy term [J]
U_bending_N = 0; % initialize bending energy term [J]
U_area_loc_N = 0; % initialize local area energy term [J]
kb_nucleus = kb; %the bending stiffness of nucleus is 2 times larger than that of cytomembrane
a_S = zeros(Nv,2);
b_S = zeros(Nv,2);
c_S = zeros(Nv,2);
ka_nucleus = ka_N; %the global area constrain of nucleus
%kd_nucleus = 3*kd_N; %the local area force coefficient of nucleus
% Calculate In-Plane and Bending Forces and Energies:
Center_v = mat2cell(r_N,(Nv_single.*ones(N_cell,1)));
Center_mean = cellfun(@(x) mean(x),Center_v,'UniformOutput',false);
for i = 1:1:Ns
    % In-Plane (2-point interaction):
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]
    x_a1 = r_N(a1,1); % x-coordinate of node 1 [m]
    y_a1 = r_N(a1,2); % y-coordinate of node 1 [m]
    x_b1 = r_N(b1,1); % x-coordinate of node 2 [m]
    y_b1 = r_N(b1,2); % y-coordinate of node 2 [m]
    r_21_1 = [(x_b1 - x_a1) (y_b1 - y_a1)]; % relative position vector (from node 1 to node 2) [m]
    len = sqrt(r_21_1(1,1)*r_21_1(1,1) + r_21_1(1,2)*r_21_1(1,2)); % distance between node 1 and node 2 [m]
    e_21_1 = r_21_1./len; % unit relative position vector (from node 1 to node 2) [m]
    LEN_N(i,1) = len; % save length of current link [m]
    xn = len/Lm_N; % normalized length of spring extension [unitless]
%     if xn > 0.9
%         xn=0.9;
%     end        
    wlc_len =(-kB*T/Lp_nucleus)*(1/(4*(1 - xn)^2) - (1/4) + xn);
    f_wlc_1 = wlc_len.*(-e_21_1); % WLC contribution to force on node 1 [N]
    f_wlc_2 = -f_wlc_1; % WLC contribution to force on node 2 [N]
    kp_len_N = kp_N/(len^m);
    f_pow_1 =  kp_len_N.*(-e_21_1); % POW contribution to force on node 1 [N]
    f_pow_2 = -f_pow_1; % POW contribution to force on node 2 [N]
    f_links(a1,:) = f_pow_1 + f_wlc_1+ f_links(a1,:);
    f_links(b1,:) = f_pow_2 + f_wlc_2+ f_links(b1,:);
    U_wlc = (kB*T*Lm/(4*Lp_nucleus))*((3*xn^2 - 2*xn^3)/(1 - xn)); % WLC contribution to in-plane potential energy [J]
    U_pow = kp_N/((m - 1)*len^(m - 1)); % POW contribution to in-plane potential energy [J]
    U_link_N = U_link_N + U_wlc + U_pow;% + U_lj*h2*CON(i,1); % total in-plane potential energy [J]
    % Bending (2-point interaction):    
    p1 = bending_pts(i,1);
    p2 = bending_pts(i,2);
    p3 = bending_pts(i,3);
    r_p1 = r_N(p1,:);
    r_p2 = r_N(p2,:);    
    r_p3 = r_N(p3,:);     
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

    beta_bending = kb_nucleus*(sine*cos(theta0) - cosine*sin(theta0))/sine; % bending force coefficient    
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
    U_bending_N = U_bending_N + kb*(1 - cos(theta - theta0)); % total bending potential energy [J]      
    %Area (2-point interaction);
    a_S(i,:) = r_N(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    b_S(i,:) = r_N(b1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    S(i,:)=0.5.*cross([a_S(i,:) 0],[b_S(i,:) 0]);
    A_L(i) = S(i,3); %local area of cell
end
% global Area:
%AL_cell = mat2cell(A_L,(Ns/N_cell.*ones(N_cell,1)));
%AL_sum = cellfun(@(x) sum(x(:,1)),AL_cell,'UniformOutput',false);
%A_G_N = [AL_sum{:,1}]';
A_G_N = sum(A_L);
beta_area = -ka_nucleus*(A_G_N - A0t_N)/A0t_N; % coefficient in global area force calculation
for i = 1:1:Ns
    a1 = Links(i,1); % index of node 1 [unitless]
    b1 = Links(i,2); % index of node 2 [unitless]    
    a_S(i,:) = r_N(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
    b_S(i,:) = -r_N(b1,:)+Center_mean{ceil(i/(Ns/N_cell))};
    c_S(i,:) = r_N(a1,:)-r_N(b1,:);
    ksi_2 = cross([a_S(i,:) 0],[b_S(i,:) 0]);
    alpha_area_g = beta_area(ceil(i/(Ns/N_cell)))/(4*A_L(i,1)); % global area force coefficient
    f_area_g_1 = alpha_area_g.*cross(ksi_2,[b_S(i,:) 0]); % global area force 1 [N]
    f_area_g_2 = alpha_area_g.*cross(ksi_2,[a_S(i,:) 0]); % global area force 2 [N]
    f_area_g_3 = -alpha_area_g.*cross(ksi_2,[c_S(i,:) 0]); % global area force 3 [N]           
    f_area_g(a1,:) = -f_area_g_1 + f_area_g(a1,:);% - f_area_g_3./2; % global area force on node 1 [N]
    f_area_g(b1,:) = -f_area_g_2 + f_area_g(b1,:);% - f_area_g_3./2; % global area force on node 2 [N]    
end
U_area_glo_N = sum(ka_nucleus.*((A_G_N - A0t_N).^2)./(2.*A0t_N)); %area induced energy change
% % local Area:
% for i = 1:1:Ns
%     a1 = Links(i,1); % index of node 1 [unitless]
%     b1 = Links(i,2); % index of node 2 [unitless]
%     a_S(i,:) = r_N(a1,:)-Center_mean{ceil(i/(Ns/N_cell))};
%     b_S(i,:) = -r_N(b1,:)+Center_mean{ceil(i/(Ns/N_cell))};
%     c_S(i,:) = r_N(a1,:)-r_N(b1,:);
%     ksi_2 = cross([a_S(i,:) 0],[b_S(i,:) 0]);        
%     beta_area = -ka_N*(A_L - A0_N(i))/A0_N(i); % coefficient in global area force calculation    
%     alpha_area_l = beta_area(ceil(i/(Ns/N_cell)))/(4*A_L(i,1)); % global area force coefficient   
%     f_area_l_1 = alpha_area_l.*cross(ksi_2,[b_S(i,:) 0]); % global area force 1 [N]
%     f_area_l_2 = alpha_area_l.*cross(ksi_2,[a_S(i,:) 0]); % global area force 2 [N]    
%     f_area_l_3 = -alpha_area_l.*cross(ksi_2,[c_S(i,:) 0]); % global area force 2 [N]   
%     f_area_loc(a1,:) = f_area_l_3./2; %-f_area_l_1 + f_area_loc(a1,:);% + f_area_l_3./2; % global area force on node 1 [N]
%     f_area_loc(b1,:) = f_area_l_3./2; %-f_area_l_2 + f_area_loc(b1,:);% + f_area_l_3./2; % global area force on node 2 [N]   
% end
%U_area_loc_N = sum(ka_N.*((A_L - A0_N).^2)./(2.*A0_N)); %area induced energy change

F_C_N = f_links +  f_area_g(:,1:2) +  f_area_loc(:,1:2) + f_bending + f_intermediate_N + f_potential_N + f_actin_N; % total conservative force on all nodes [N]
end