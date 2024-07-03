function [E, nu, rho, alpha, k] = compute_temperature_dependent_material(mat_type, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temperature-dependent material library %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Unit: m - N - kg - s^2 - Pa - K
switch mat_type
    case 0  % Steel
        E = 200e9; nu = 0.3; rho = 8000;
        alpha = 0; k = 0;  % Not in use
    case 1  % Poly methyl methacrylate (PMMA)
        E = 2.5e9; nu = 0.34; rho = 1150;
        alpha = 0; k = 0;  % Not in use
    case 2  % Polystyrene (PS)
        E = 1.9e9; nu = 0.3; rho = 1056.9;
        alpha = 0; k = 0;  % Not in use
    case 3  % Poly[(m-phenylenevinylene)-co-(2,5-dioctoxy-p-phenylenevinylene)] (PmPV)
        T_0 = 300;
        E = (2.1 - 0.0047*(T - T_0))*1e9; nu = 0.34; rho = 1190;
        alpha = 45*(1 + 0.0005*(T - T_0))*1e-6; k = 0;
end
end
