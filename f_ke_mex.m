function KE = f_ke_mex(vel,MASS)

global Nv
TEMP_1 = zeros(Nv,1);
for i = 1:1:Nv
    TEMP_1(i,1) = vel(i,1)*vel(i,1) + vel(i,2)*vel(i,2); %squared velocity of node [m^2/s^2] squared velocity of node [m^2/s^2]
end

KE = 0.5*MASS*sum(TEMP_1); % total kinetic energy [J]

end

