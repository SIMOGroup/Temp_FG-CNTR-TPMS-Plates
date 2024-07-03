%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title: Isogeometric analysis for Temperature-dependent CNTR-FG-TPMS plates with five-variable plate model %%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% ! Please reference to paper: ............................................
% ! This work can be used, modified, and shared under the MIT License
% ! This work can be found in https://github.com/SIMOGroup/Temp_CNTR-FG-TPMS-Plates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================== Initialization =============================
tic
addpath(genpath('./'));

% clc
clear all
% close all
format long

%% ============================ Plate geometry ============================
% === Physical geometric properties ===
Plate.geo.L = 1;
Plate.geo.W = Plate.geo.L;
Plate.geo.h = Plate.geo.L/(10);
% Plate.geo.h = Plate.geo.L/20;

% === Plate theory ===
% [1]: fz = 0, [2]: Reddy (1984), [3]: Shimpi (2002), [4]: Nguyen-Xuan (2013), [5]  Nguyen (2017)
% [6]: Nguyen (2016), [7]: Thai (2014), [8]: Touratier (1991), [9]: Hoang (2023), [10]: Newly developed
Plate.theory.shear_func = 10;

% === Define NURBS functions ===
IGA.NURBS.deg = 3; % Degree of basis functions
IGA.NURBS.ref = 9; % Number of mesh refinement

%% ============================ Material ==================================
% === Base matrix material ===
% [0]: Steel, [1]: PMMA, [2]: PS, [3]: PmPV
Plate.mat.type = 3;
Plate.mat.T_0 = 300;

% === CNT reinforcement ===
Plate.CNT.type = 1;  % Type of CNT
Plate.CNT.fr = 0;  % CNT reinforcement ratio
Plate.CNT.mu = 1;  % Agglomeration ratio
Plate.CNT.eta = 1;  % CNT inside agglomeration/total ratio

% === Porous material ===
% [1]: Primitive, [2]: Gyroid, [3]: IWP, [4]: Closed-cell, [5]: Open-cell (\nu = 0.33), [6]: Mod Open-cell (\nu = 0.3)
Plate.por_mat.type = 3;
% Plate.por_mat.RD = 1;

% === Porosity distribution ===
% [1]: A (asymmetric), [2]: B (symmetric), [3]: C (uniform)
Plate.por_dis.type = 2;
Plate.por_dis.RD_avg = 0.35;
Plate.por_dis.RD_max = 1;
Plate.por_dis.RD_0 = 0.8;  % RD_min = RD_max*(1- RD_0);

% === Temperature distribution ===
% [1]: Uniform (T = T1), [2]: Linear
Plate.temp_dis.type = 1;
Plate.temp_dis.T1 = 400;  % Bot
Plate.temp_dis.T2 = 300;  % Top

%% ========================== Problem type ================================
% [1]: Static, [2]: Vibration
Plate.prob.type = 2;
switch Plate.prob.type
    case 1  % Static
        Plate.prob.q_uniform = 1;  % Static load
%         Plate.prob.q_uniform = 0.1e6;  % Static load
        IGA.result.norm_method = 1;
    case 2  % Vibration
        IGA.result.norm_method = 1;
        IGA.result.nmode = 6;
end

%% ========================= Boundary condition ===========================
% [1]: Fully simply supported (SSSS), [2]: Fully clamped (CCCC)
Plate.bc.bc_case = 1;

%% ======================== Material matrices =============================
[Plate.mat_mat.Db, Plate.mat_mat.Ds, Plate.mat_mat.I, Plate.mat_mat.S_th_0] = cal_Material_Matrices_2D_5dof_CNTR_FG_TPMS(Plate);

%% =============================== IGA mesh ===============================
% === Generate NURBS mesh ===    
IGA.NURBS = Mesh_2D(Plate, IGA.NURBS);
IGA.NURBS = Gen_Ien_Inn_2D(IGA.NURBS);

% === NURBS properties ===
IGA.NURBS.nsd   = 2;                                                             % Number of spatial dimension
IGA.NURBS.nnode = IGA.NURBS.mcp * IGA.NURBS.ncp;                                 % Number of control point
IGA.NURBS.nshl  = (IGA.NURBS.p + 1) * (IGA.NURBS.q + 1);                         % Number of local shape functions (= degree + 1 per element due to k refinement)
IGA.NURBS.nel   = (IGA.NURBS.mcp - IGA.NURBS.p) * (IGA.NURBS.ncp - IGA.NURBS.q); % Number of element

% === IGA properties ===
IGA.params.ndof   = 5;                                                           % Number of dofs of a control point
IGA.params.sdof   = IGA.NURBS.nnode * IGA.params.ndof;                           % Total number of dofs of the structure
IGA.params.nGauss = IGA.NURBS.p + 1;                                             % Number of gauss point in integration

%% ========================= IGA for linear geometry ======================
% === Building global matrices === 
IGA.result.K = cal_Stiffness_Matrices_2D_5dof(IGA,Plate);                    % Stiffness
IGA.result.K_th = cal_Thermal_Stiffness_Matrices_2D_5dof(IGA,Plate);         % Thermal Stiffness
switch Plate.prob.type
    case 1  % Static
        IGA.result.F = cal_Load_Vector_Uniform_2D_5dof(IGA,Plate);           % Extenal load
    case 2  % Free vibration
        IGA.result.M = cal_Mass_Matrices_2D_5dof(IGA,Plate);                 % Mass
end

% === Imposing boundary conditions ===
[IGA.params.bcdof, IGA.params.bcval] = cal_bcdof_2D_5dof(IGA,Plate);
IGA.params.fdof = setdiff((1:IGA.params.sdof)', IGA.params.bcdof');  % Free dofs

%% === Solving weak form ===
switch Plate.prob.type
    case 1  % Static
        bcdof = IGA.params.bcdof; bcval = IGA.params.bcval;
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        IGA.result.U = zeros(sdof, 1); 
        IGA.result.U(bcdof') = bcval';
        IGA.result.F(fdof) = IGA.result.F(fdof) - (IGA.result.K(fdof, bcdof') + IGA.result.K_th(fdof, bcdof'))*bcval';
        IGA.result.U(fdof) = (IGA.result.K(fdof, fdof) + IGA.result.K_th(fdof, fdof)) \ IGA.result.F(fdof);
        
        % --- Normalization ---
        norm_method = IGA.result.norm_method;
        cen_def_norm = cal_central_deflection_2D_5dof(IGA,Plate,norm_method);
        disp("-> Normalized central deflection = " + sprintf('%.4f, ', cen_def_norm))
%         CP_deform = plot_deflection_2D_5dof(IGA,Plate,norm_method);
        
        clear norm_method
        clear sdof bcdof bcval fdof
    case 2  % Free vibration
        bcdof = IGA.params.bcdof;
        KK = full(IGA.result.K); KK_th = full(IGA.result.K_th); MM = full(IGA.result.M);
        [Lambda, ModeShape] = Eigen(KK+KK_th, MM, bcdof);
        [IGA.result.Lambda, sort_index] = sort(Lambda,'ascend');
        IGA.result.ModeShape = ModeShape(:, sort_index);
        IGA.result.Lambda = IGA.result.Lambda(IGA.result.Lambda > 0 & IGA.result.Lambda < inf);
        IGA.result.ModeShape = IGA.result.ModeShape(:, IGA.result.Lambda > 0 & IGA.result.Lambda < inf);
        
        % --- Normalization ---
        norm_method = IGA.result.norm_method; nmode = IGA.result.nmode;
        L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
        [E_m, nu_m, rho_m, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
        I_m = h*rho_m; D_m = E_m*h^(3)/(12*(1-nu_m^2));
        switch norm_method
            case 1
                Lambda_norm = (IGA.result.Lambda(1:nmode)*rho_m*L^4*h/D_m).^0.25;
            case 2
                norm = h*sqrt(rho_m/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 3
                norm = L*sqrt(rho_m*(1-nu_m^(2))/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 4
                norm = L^2/h*sqrt(rho_m/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 5
                norm = W^(2)/pi^(2)*(I_m/D_m)^(1/2);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 6
                norm = W^(2)/h*(rho_m*(1-nu_m^2)/E_m)^(1/2);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
        end
        disp("-> Normalized frequencies = " + sprintf('%.4f, ', Lambda_norm))
        
        clear norm_method nmode E_m nu_m rho_m I_m D_m L W h 
        clear bcdof KK KK_th MM Lambda ModeShape sort_index nmode
end
toc
