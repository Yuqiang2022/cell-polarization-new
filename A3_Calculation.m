clear all 
close all
clc
%% Pre-Simulation Setup:
DATE_START = datestr(now,'mm_dd_yyyy__HH_MM_SS'); % date and time of simulation initiation
global A0 A0t A0_N A0t_N kB ka kb kd kp kp_N kp_inter L0 Lm Lm_N Lp_actin Lp_nucleus m Na Ns Nv Nv_single N_cell...
       T theta0 x0 Lm_N ka_N kd_N Lm_inter Lp_inter Lm_actinfilamet Lp_actinfilamet r_C0
load('CELL_FRAMEWORK','Links','bending_pts','r_C','LEN_0','r_N','LEN_N','Actin','d0_actin','r_MT','d0_MT','d0_inter',...
    'r_C0','r_0','rN_0','r_N0','r_spreading','N_actin','ActinNumber','ind_MT_N','Nv_single','N_cell','Vertex_scale') % triangulation and initial conditions
load('PHYSICAL_PARAMETERS.mat') % parameters

r_1 = r_C;
rN_1 = r_N;
%%
%set the loop steps(and time span)
nt = 10e4;% number of time steps 400,000 steps[unitless]
dt = 5e-3; % length of time step [s] nommalize the time step 2.5e-8 / 250e-9 / 100e-9
%Initialize/Declare Arrays/Vectors/Terms:
K_MM = zeros(nt,1);
Tension = zeros(nt,1);
VELOCITY = zeros(nt,1);
vel = zeros(Nv,2); %initialize current velocity array to zero [m/s]
vel_N = zeros(Nv,2); %initialize current velocity array to zero [m/s]
acc = zeros(Nv,2); %initialize current acceleration array to zero [m/s^2]
acc_N = zeros(Nv,2); %initialize current acceleration array to zero [m/s^2]
vel_p = zeros(Nv,2); %initialize previous velocity array to zero [m/s]
vel_p_N = zeros(Nv,2); %initialize previous velocity array to zero [m/s]
F_C = zeros(Nv,2); %initialize conservative force array
F_C_N = zeros(Nv,2); %initialize conservative force array
F_TOT = zeros(Nv,2); %initialize total force array
F_TOT_N = zeros(Nv,2); %initialize total force array
F_TOT_p = zeros(Nv,2); %initialize previous total force array
F_TOT_p_N = zeros(Nv,2); %initialize previous total force array
F_trac = zeros(Nv,2); %initialize traction force array
%F_C_p = zeros(Nv,2); %initialize previous conservative force array

%Simulation:
U_intermediate = zeros(nt,1); %Intermediate filaments  energy [J]
U_potential = zeros(nt,1); %Interspace potential energy [J]
U_actin = zeros(length(Actin),1); %actin filaments potential energy [J]
U_LINKS = zeros(nt,1); %initialize in-plane potential energy vector [J]
U_BENDING = zeros(nt,1); %initialize bending potential energy vector [J]
U_AREA_loc = zeros(nt,1); %initialize local area potential energy vector [J]
U_AREA_g = zeros(nt,1); %initialize global area potential energy vector [J]
U_INTER =  zeros(nt,1);
U_POTENTIAL = zeros(nt,1);
U_Actin = zeros(nt,1);
U_Microtubule = zeros(nt,1);
U_TOT = zeros(nt,1); %initialize total potential energy vector [J]
V_TOT = zeros(nt,1); %initialize total kinetic energy vector [J]
E_TOT = zeros(nt,1); %initialize total energy vector [J]
U_LINKS_N = zeros(nt,1); %initialize in-plane potential energy vector [J]
U_AREA_loc_N = zeros(nt,1); %initialize local area potential energy vector [J]
U_AREA_g_N = zeros(nt,1); %initialize global area potential energy vector [J]
U_VOLUME_N = zeros(nt,1); %initialize volume potential energy vector [J]
U_BENDING_N = zeros(nt,1); %initialize bending potential energy vector [J]
U_SURFACE = zeros(nt,1); %initialize bending potential energy vector [J]
U_TOT_N = zeros(nt,1); %initialize total potential energy vector [J]
V_TOT_N = zeros(nt,1); %initialize total kinetic energy vector [J]
E_TOT_N = zeros(nt,1); %initialize total energy vector [J]
Indentation_t = zeros(nt,1);%initialize indentation recording
n_leading = zeros(nt,1);%initialize the amount of leading points
F_leading = zeros(nt,1);%initialize the amount of leading forces
r_p0 = r_0; %record initial positions on cytomembrane for potential calculation[m]
rN_p0 = rN_0; %record initial positions on nucleus for potential calculation[m]
A_G_N = A0t_N;
%record the point to fix on the substrate 
r_fix_x = ones(Nv,1);
r_fix_y = ones(Nv,1);
r_fix_x_p = r_fix_x;
r_fix_y_p = r_fix_y;
r_fix_x_N = ones(Nv,1);
r_fix_y_N = ones(Nv,1);
r_fix_x_p_N = r_fix_x_N;
r_fix_y_p_N = r_fix_y_N;
r_C_p = r_C;
r_N_p = r_N;
r_contact_00 = zeros(Nv,2);
r_index = ones(Nv,1);%record the leading edge molecules
tstart = tic; % initialize start time of simulation (wall time)
%% Reaction-diffusion module
Y=2.0.*ones(Nv,1); % inactive Cdc42
X=0.2683312.*ones(Nv,1); % active Cdc42
D_bar = 0.1;
D_bar_ina = 10;
D_bar_finite = D_bar*(mean(LEN_0)*1e6)^2; % the equivalent diffusion rate in finite difference methods
D_bar_ina_finite = D_bar_ina*(mean(LEN_0)*1e6)^2; % the equivalent diffusion rate in finite difference methods
K_plus=0.067; %[s-1] Basal GEF conversion rate
gama=1; %[s-1] GAP hydrolysis rate
K=1; %Saturation parameter
K_minus=1; %Inactivation rate
S=0.15;%0.15; %activation rate of Cdc42
%Y=K_minus.*X./(K_plus+gama.*X.*X./(K*K+X.*X));
vel_head = ones(Nv,1); %leading particle index
X_int=[X(Nv);X(1:(Nv-1))];
X_pre=[X(2:Nv);X(1)];
Y_int=[Y(Nv);Y(1:(Nv-1))];
Y_pre=[Y(2:Nv);Y(1)];
d_int=488.*ones(Nv,1); %[1] Initial integrin density
d_int0=d_int;
d_int_min=488;%[1] bottom limit of integrin number
d_int_max=2500;%[1]up limit of integrin number
d_int_pre = d_int;
d0=4.*ones(Nv,1); %[1] Integrins  added  after  each  reinforcement event
%Fcr=110e-12; %[pN] 
r_substrate=zeros(Nv,1);
r_center=zeros(nt,2);%center of cell
n_lead=ones(Nv,1);%leading particle number
mu_d = 0.003; % 0.003 [Ns/m] damping coefficient
f_area_g = r_C-mean(r_C);
polarity = ones(Nv,1);
F_s = zeros(Nv,1); %traction forces from the elastic ECM
t_polar=zeros(nt,Nv);
t_int=zeros(nt,Nv);
t_X=zeros(nt,Nv);
t_Y=zeros(nt,Nv);
% t_u_load=zeros(Nv,nt);
t_vel=zeros(nt,2);
%t_xy=zeros(nt,2);
%t_leading=zeros(nt,2);

for t=1:1:nt
    if t<=1e4
        Relax=0; %relax time
    else
        Relax=1;
    end
    %Update Numerical Integration Time Span:
    tb = t*dt-dt; %(time before) time at beginning of time step [s]
    ta = t*dt; %(time after) time at end of time step [s]
 %%   
    %Slove by the Multi-scale model
    rr_inter = r_C-mean(r_C);
    de_inter = (rr_inter(:,1).*f_area_g(:,1)+rr_inter(:,2).*f_area_g(:,2));
    de_inter(de_inter>0)=1;
    de_inter(de_inter<0)=-1;    
    vv_inter = f_area_g(:,1:2).*[de_inter de_inter]; %the direction of outward along the area force
    vector_int = vv_inter./sqrt(vv_inter(:,1).*vv_inter(:,1)+vv_inter(:,2).*vv_inter(:,2));
    F_substrate = F_s.*vector_int;
    r_C = r_C +(F_substrate.*[vel_head vel_head] + F_TOT)./mu_d.*dt;% matrix to immobilize attachment on substrate new positions [m]  + F_APP.*vector_int./mu_d
    r_N = r_N + F_TOT_N./mu_d.*dt;% matrix to immobilize attachment on substrate new positions [m] 
    vel_present_C = (r_C-r_C_p)./dt;
    vel_present_N = (r_N-r_N_p)./dt;   
    r_C_p = r_C;
    r_N_p = r_N;
    [F_C,f_links_tensile,f_actin,f_area_g,U_link,U_bending,U_area_loc,U_area_glo,U_inter,U_potential,U_actin,LEN,f_intermediate_C,f_intermediate_N,f_potential_N,f_actin_N] = f_cons_forces_mex(Links,bending_pts,r_C,r_N,Actin,ActinNumber,r_p0,rN_p0,Vertex_scale);
    [F_C_N,U_link_N,U_bending_N,U_area_loc_N,U_area_glo_N,LEN_N,A_G_N] = f_cons_forces_mex_N(Links,bending_pts,r_N,f_intermediate_N,f_potential_N,f_actin_N); % calculate conservative forces and energies of nucleus

%   %substrate soft to stiffness direction
    r_substrate(r_C(:,1)<0.0*r_C0)=-1;    
    r_substrate(r_C(:,1)>=0.0*r_C0)=-1;  
     
    Ksub=10.^r_substrate;    
    F_s = polyval(f_fit,r_substrate).*1e-12; %[N]
%% polarization module
    FF_link = sqrt(f_links_tensile(:,1).*f_links_tensile(:,1)+f_links_tensile(:,2).*f_links_tensile(:,2));
    FF_actin = sqrt(f_actin(:,1).*f_actin(:,1)+f_actin(:,2).*f_actin(:,2));
    FF = FF_link + FF_actin;%  + FF_intermediate;
    polarity = FF./mean(FF);
    if Relax==1
       d_int=d_int+d0.*(polarity-1); % .*(sum(d_int)<Nv*d_int_max) integtin pool is constant .*(sum(d_int)<Nv*d_int_max)
       d_int=round(d_int);
    end
    d_int(d_int>d_int_max,1)=d_int_max; %up limit of integrin number
    d_int(d_int<d_int_min,1)=d_int_min; %up limit of integrin number    

    %finite deference format with variable length
    r_C_int=[r_C(Nv,:);r_C(1:(Nv-1),:)];
    r_C_pre=r_C;
    delt_r_1=sqrt((r_C_int(:,1)-r_C_pre(:,1)).^2+(r_C_int(:,2)-r_C_pre(:,2)).^2).*1e6;%[um]
    delt_r_2=[delt_r_1(2:Nv,1);delt_r_1(1,1)];%[um]
    delt_XY=Y.*(K_plus+gama.*X.*X./(K*K+X.*X))-K_minus.*X+S.*(d_int-d_int_pre).*Y;
    d_int_pre = d_int;  
    X_new=X+dt.*(D_bar_finite*2.*(delt_r_2.*X_int+delt_r_1.*X_pre-(delt_r_1+delt_r_2).*X)./(delt_r_1.^2.*delt_r_2+delt_r_1.*delt_r_2.^2)+delt_XY); %
    Y_new=Y+dt.*(D_bar_ina_finite*2.*(delt_r_2.*Y_int+delt_r_1.*Y_pre-(delt_r_1+delt_r_2).*Y)./(delt_r_1.^2.*delt_r_2+delt_r_1.*delt_r_2.^2)-delt_XY); %
    X = X_new;
    Y = Y_new;   
    X_int=[X(Nv);X(1:(Nv-1))];
    X_pre=[X(2:Nv);X(1)]; 
    Y_int=[Y(Nv);Y(1:(Nv-1))];
    Y_pre=[Y(2:Nv);Y(1)]; 

    %averge X and Y within the adjacent vertexes when negtive values exist   
    X_ID = find(X<=0);
    Y_ID = find(Y<=0);
    for i_X=1:length(X_ID)
        X(Lap_matrix(X_ID(i_X),:)~=0) = mean(X(Lap_matrix(X_ID(i_X),:)~=0));%mean the contentration of X for the regions of X<=0
    end
    for i_Y=1:length(Y_ID)
        Y(Lap_matrix(Y_ID(i_Y),:)~=0) = mean(Y(Lap_matrix(Y_ID(i_Y),:)~=0));%mean the contentration of X for the regions of Y<=0 
    end    
   
    vel_head = (X./mean(X)); % with polarized schem
   
    t_polar(t,:)= polarity;
    t_int(t,:) = d_int;
    t_X(t,:) = X;
    t_Y(t,:) = Y;
    
    t_vel(t,:) = mean(vel_present_C);
    %t_leading(t,:) = r_C(r_C(:,1)==min(r_C),:);
    t

%     d_int_t(t,:)=d_int;
%     delt_XYY(t,:)=delt_XY;
%     D_XY(t,:)=D_bar_finite.*(-Lap_matrix*X);
%     Center_ar(t,:) = mean(polarity.*r_C);
%     Center_real(t,:) = mean(r_C);
    %%
    F_TOT_0 = F_C;% + F_D + Relax.*F_trac ;% + 0.*F_APP + 0.*F_constraint + (-f_MT_v) + 0.*f_doubleCell %new total forces [N]
    F_TOT = F_TOT_0;% + F_friction;
    F_TOT_N = F_C_N;% + F_D_N;%+ f_MT_N %new total forces [N]
    %Solve the potential and kinetic energy
    U_LINKS(t,1) = U_link; %links potential energy [J]
    U_BENDING(t,1) = U_bending; %bending potential energy [J]
    U_AREA_loc(t,1) = U_area_loc; %local area potential energy [J]
    U_AREA_g(t,1) = U_area_glo; % global area potential energy [J]
    U_INTER(t,1) = U_inter; %Intermediate filaments energy[J]
    U_POTENTIAL(t,1)=U_potential;%Potential energy between cytomembrane and nucleus [J]
    U_Actin(t,1) = U_actin;%actin filaments energy [J]
    %    U_Microtubule(t,1) = U_MT; %MT energy [J]
    U_TOT(t,1) = U_LINKS(t,1) + U_BENDING(t,1) + U_AREA_loc(t,1) + U_AREA_g(t,1) + U_INTER(t,1) + U_POTENTIAL(t,1) + U_Actin(t,1) + 1.*U_SURFACE(t,1);% + U_Microtubule(t,1)% total potential energy [J]
    if U_TOT(t,1) > 1e-3 %error if potential energy is greater than threshold value
        n_error = t; %record timestep where error occured
        error('Error: Invalid potential energy on cutomembrane');
    end
   U_LINKS_N(t,1) = U_link_N; %links potential energy [J]
   U_BENDING_N(t,1) = U_bending_N; %bending potential energy [J]
   U_AREA_loc_N(t,1) = U_area_loc_N; %local area potential energy [J]
   U_AREA_g_N(t,1) = U_area_glo; %global area potential energy [J]
   U_TOT_N(t,1) = U_LINKS_N(t,1) + U_AREA_loc_N(t,1) + U_AREA_g_N(t,1) + U_BENDING_N(t,1); % total potential energy [J]
    if U_TOT_N(t,1) > 1e-3 %error if potential energy is greater than threshold value
        n_error_N = t; %#ok<NASGU> %record timestep where error occured
        error('Error: Invalid potential energy on nucleus');
    end
    V_TOT(t,1) = f_ke_mex(vel_present_C,1); %total kinetic energy [J]
    V_TOT_N(t,1) = f_ke_mex(vel_present_N,1); %total kinetic energy [J]
    E_TOT(t,1) = U_TOT(t,1) + V_TOT(t,1); %total energy [J]
    E_TOT_N(t,1) = U_TOT_N(t,1) + V_TOT_N(t,1); %total energy [J]
    
    if max(abs(r_N(:,1))) > max(abs(r_C(:,1)))
        error('Error: Nucleus squeezed out X');
    end
    if max(abs(r_N(:,2))) > max(abs(r_C(:,2)))
        error('Error: Nucleus squeezed out Y');
    end
    %Calculate Simulation Time and Output Progress:
    t_end = toc(tstart); % total elapsed time (real-time) at end of current timestep [s]
    n_error = nt; % record last step of simulation
    if (fix(t/nt*1)==t/nt*1)||t==nt   %save every 1000 loop
        LABEL = datestr(now,'yyyy_mm_dd_HH_MM_SS'); %date and time of completion of parameter derivation
        fname = sprintf('Spread_uniform_%s',LABEL);
        save(fname)
    end
end

BB_Plot_shape