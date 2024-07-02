clear all
close all
clc
% Load Data:
load('CELL_FRAMEWORK.mat') % triangulation and initial node locations

% Given Parameters:
Na = 6.0221415e23; % Avogadro constant [1/mol]
kB = 1.380653e-23; % Boltzmann constant [J/K]
T = 25 + 273.15; % room temperature [K]
kBT = kB*T;
% WLC
x0 = 1/2.2; % L0/Lm ratio [unitless]
x0_inter = 1/2.2;
x0_actinfilament = 1/2.2;
% POW
m = 2; % power law exponent [unitless]
kc = 2.4e-19; % bending rigidity [J]
% bulk mechanical properties
u0 = 6.3e-6; % observed linear elastic shear modulus of cytomembrane [N/m]
up = 2*u0; % actual linear elastic shear modulus [N/m]
Y = 3.92453*up; % Young's modulus [N/m]
K = (Y*up)/(4*up - Y); % Area-compression modulus (2-D Bulk modulus) [N/m]
v = (K - up)/(K + up); % Poisson's ratio (since 2-D, v = 1 gives a incompressible material, instead of v = 0.5) [unitless]
% geometric dimension
theta0 = (Nv-2)*pi/Nv; % 2D equilibrium (initial) angle [rad]
A0 = A0_L; % average area of individual triangle [m^2]
A0t = A0_G; %total cell membrane area
A0_N = A0_L_N; % average area of individual triangle [m^2]
A0t_N = A0_G_N; %total nucleus membrane area
%% calculate the equilibrium length of the filaments
% systematic coarse-graining method based on mean-field theory
L0_fine = 75.5e-9; %[m] %fine model of RBC
Nv_fine = 27344;%[1] %fine model of RBC
Lp_fine = 14.68e-9;%[m] %fine model of RBC
kp_fine = 1.66e-27;%[Nm2] %fine model of RBC
% coarsed cell membrane parameters
%Nv_single_cell=(4*pi*radius_C*radius_C)/0.5/sqrt(3)/(mean(LEN)*mean(LEN))+2;
Nv_single_cell=(Nv_fine-2)*(L0_fine/mean(LEN))*(L0_fine/mean(LEN))+2;
L0_c = L0_fine*sqrt((Nv_fine-2)/(Nv_single_cell-2));
Lm_c = L0_c/x0;
Lp_c = Lp_fine*L0_fine/L0_c;
kp_c = kp_fine*(L0_c/L0_fine)^(m+1);
L0 = L0_c;
Lm = Lm_c;
Lp = Lp_c;
kp = kp_c;
Lp_actin = Lp;
% coarsed nuclear membrane parameters
%Nv_single_nucleus=(4*pi*radius_N*radius_N)/0.5/1.732/(mean(LEN_N)*mean(LEN_N))+2;
Nv_single_nucleus =(Nv_fine-2)*(L0_fine/mean(LEN_N))*(L0_fine/mean(LEN_N))+2;
L0_N = L0_fine*sqrt((Nv_fine-2)/(Nv_single_nucleus-2));
Lm_N = L0_N/x0;
Lp_N = Lp_fine*L0_fine/L0_N;
kp_N = kp_fine*(L0_N/L0_fine)^(m+1);
Lp_nucleus = Lp_N;
% intracelular filaments parameters
Lm_actinfilamet = d0_actin./x0_actinfilament; % maximum length of actin filaments (contour length) [m]
Lp_actinfilamet = Lp_fine.*L0_fine./d0_actin; % persistent length of actin filaments
Lp_actinfilamet = Lp_actinfilamet./Lp_actinfilamet.*Lp;
Lm_inter = d0_inter./x0_inter; % maximum length of intermediate filaments extension (contour length) [m]
Lp_inter = Lp_fine.*L0_fine./d0_inter;  % persistent length of intermediate filaments
Lp_inter = Lp_inter./Lp_inter.*Lp;
%kb = 2*kc/sqrt(3); % bending stiffness [J]
%kb = 1e-14;% bending rigidity [J]
kb = 10000*kc;% the observed bending stiffness [J] 
kb = 1e-14;
kp_inter=5e-16/Vertex_scale; %[J]; % artificial setup of repulsive connection coefficient between nucleus and cell membrane
ka = 4.34e6*A0t;
ka_N = 4.34e6*A0t_N;

%laplace matrix format
Lap_ratio = (LEN.*1e6)./(LEN.*1e6);%LEN_0.*LEN_0./LEN./LEN;         
Lap_matrix = zeros(Nv_single,Nv_single);    
for i=1:Nv_single %use the current length to make Laplace matrix      
    Lap_matrix(Links(i,1),(Links(i,2)))=-Lap_ratio(i,1);
    Lap_matrix(Links(i,2),(Links(i,1)))=-Lap_ratio(i,1);
    Lap_matrix(i,i)=Lap_ratio(Links(:,1)==i,1)+Lap_ratio(Links(:,2)==i,1);
end  

%% motor-cluch fit module
n_myosine = 50; %the myosine number is proportional to length of the filament
K_cluch = 0.8; %set the cluch stiffness as the leading edge stiffness 
K_substrate=logspace(-2,2,20); %define stiffness vector  
n_clutch_max=50;
[v_retrograde,F_retrograde] = F_motor_cluch_FRM(n_myosine,K_cluch,K_substrate,n_clutch_max); %motor cluch to calculate the traction forces
x = log10(K_substrate);
y = -F_retrograde;
f_fit = polyfit(x,y,9);

%load Paper_fig3_400_FRM
xx = log10(logspace(-2,2,50));
figure
plot(xx,polyval(f_fit,xx)) 
%% Save Parameters:
DATE_hRBC_Parameters = datestr(now,'mm_dd_yyyy__HH_MM_SS'); % date and time of completion of parameter derivation
save('PHYSICAL_PARAMETERS.mat','DATE_hRBC_Parameters','kB','kb','ka','kp','L0','Lm','Lp_actin','Lp_nucleus','m','Na','Ns','Nv',...
    'T','theta0','x0','A0','A0t','A0_N','A0t_N','Lm_N','kp_N','ka_N','kp_inter','Lm_inter','Lp_inter','Lm_actinfilamet', ...
    'Lp_actinfilamet','x0_actinfilament','f_fit','Lap_matrix')
