function [E, nu, rho] = compute_material(mat_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Material library %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Unit: m - N - kg - s^2 - Pa
switch mat_type % [E(Pa), nu, rho (kg/m^3)]
    case 1  % Steel
        mat_prop = [200e9, 0.3, 8000];
    case 2  % Alumium
        mat_prop = [70e9, 0.3, 2702];
    case 3  % Titanium
        mat_prop = [106e9, 0.3, 4510];
    case 4  % Copper
        mat_prop = [117e9, 0.3, 8960];
    case 5  % Brass
        mat_prop = [90e9, 0.3, 8730];
    case 6  % Steel 2
        mat_prop = [206.8e9, 0.316, 8730];
    case 7  % Poly methyl methacrylate (PMMA)
        mat_prop = [2.5e9, 0.34, 1150];
%         mat_prop = [1.8e9, 0.34, 1150];
%         mat_prop = [1.9e9, 0.3, 1056.9];
    case 8  % Poly[(m-phenylenevinylene)-co-(2,5-dioctoxy-p-phenylenevinylene)] (PmPV)
%         mat_prop = [2.1e9, 0.34, 1150];
        mat_prop = [2.1e9, 0.34, 1190];
    case 9  % Polystyrene (PS)
        mat_prop = [1.9e9, 0.3, 1056.9];
end
E = mat_prop(1); nu = mat_prop(2); rho = mat_prop(3); 
end
