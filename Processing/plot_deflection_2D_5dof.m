function [CP_deform] = plot_deflection_2D_5dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot deflection of plate %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
p = IGA.NURBS.p; uKnot = IGA.NURBS.uKnot; 
q = IGA.NURBS.q; vKnot = IGA.NURBS.vKnot;
CP = IGA.NURBS.CP;
ndof = IGA.params.ndof; 
U = IGA.result.U;
norm_method = IGA.result.norm_method;

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
q_uniform = Plate.prob.q_uniform;

%% ===== Initial deformation plot =====
CP_deform = CP;

%% ====== Plot deflection ======
size_factor = 1;
count = 0;
for i = 1:size(CP,1)
    for j = 1:size(CP,2) 
        count = count + 1;
        
        % --- Stretching effect deformation function ---
        def = U(ndof*(count-1)+3);
        
        % --- Normalization ---
        switch norm_method
            case 1
                [E_m, nu_m, ~, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
                D_m = E_m*h^(3)/(12*(1-nu_m^2));
                norm = 100*D_m/(q_uniform*L^4);
                def_norm = def*norm;
            case 2
                norm = 1/h;
                def_norm = def*norm;
            case 3 
                [E_m, nu_m, ~, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
                D_m = E_m*h^(3)/(12*(1-nu_m^2));
                norm = 1 / (4.062e-3*q_uniform*L^4/D_m);
                def_norm = def*norm;
        end
        CP_deform(i,j,3) = CP(i,j,3) + size_factor*sign(q_uniform)*def_norm;
    end
end
plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP_deform);
end
