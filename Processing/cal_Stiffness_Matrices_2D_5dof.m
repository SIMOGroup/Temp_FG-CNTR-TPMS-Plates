function [K] = cal_Stiffness_Matrices_2D_5dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate stiffness matrices with 5-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Inn = IGA.NURBS.Inn; Ien = IGA.NURBS.Ien; nel = IGA.NURBS.nel;
sdof = IGA.params.sdof; ndof = IGA.params.ndof; nGauss = IGA.params.nGauss; 

%% Used parameters from Plate
Db = Plate.mat_mat.Db;
Ds = Plate.mat_mat.Ds;

%% ===== Initial stiffness matrices =====
K = sparse(sdof,sdof);

%% ===== Initial Gauss integration =====
[Pg, Wg] = gauleg(1,-1,nGauss);
nel_nza = 0; % Number of non-zero-area elements
tol = 1e-8;

%% ===== Stiffness matrices =====
for iel = 1:nel  % Loop over the elements
    % --- Element parameters ---
    sctr = Ien(iel,:);  % Control points indexes
    nn = length(sctr);  % Number of control points in the element
    for idof = 1:ndof  % Dofs of control points
        sctrK(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
    end
    ni = Inn(Ien(iel,1),1);  % Index of the element in parametric domain
    nj = Inn(Ien(iel,1),2);
    
    % --- Gauss integration ---
    if abs(uKnot(ni)-uKnot(ni+1)) > tol && abs(vKnot(nj)-vKnot(nj+1)) > tol  % Check if the current element has nonzero area in the parametric domain
        nel_nza = nel_nza + 1;
        detJ2_xi = (uKnot(ni+1) - uKnot(ni))/2;  % Mapping parametric domain into natural domain of [[-1, 1]; [-1, 1]]
        detJ2_eta = (vKnot(nj+1) - vKnot(nj))/2;

        for iGauss = 1: nGauss  % Loop over the integration points
            for jGauss = 1: nGauss
                % Gauss points & weights
                gpt_xi = Pg(iGauss); gwt_xi = Wg(iGauss);
                gpt_eta = Pg(jGauss); gwt_eta = Wg(jGauss);
                
                % Map the point to parametric domain
                gpt_xi = (uKnot(ni+1)-uKnot(ni))/2*gpt_xi + (uKnot(ni+1)+uKnot(ni))/2; gwt_xi = gwt_xi*detJ2_xi;
                gpt_eta = (vKnot(nj+1)-vKnot(nj))/2*gpt_eta + (vKnot(nj+1)+vKnot(nj))/2; gwt_eta = gwt_eta*detJ2_eta;
                gwt = gwt_xi*gwt_eta;
                
                % Kinematic matrices of NURBS shape function
                [N, dNdxi, dNdxi2, dNdxy, dN2dxy, detJ1] = Kine_Shape_2nd_2D(IGA,ni,nj,gpt_xi,gpt_eta);
                
                % Stiffness matrix
                B0 = zeros(3,ndof*nn);  % Pure membrane
                B0(1,1:ndof:ndof*nn) = dNdxy(:,1)';
                B0(2,2:ndof:ndof*nn) = dNdxy(:,2)';
                B0(3,1:ndof:ndof*nn) = dNdxy(:,2)';
                B0(3,2:ndof:ndof*nn) = dNdxy(:,1)';

                B1 = zeros(3,ndof*nn);  % Pure bending
                B1(1,3:ndof:ndof*nn) = -dN2dxy(:,1)';
                B1(2,3:ndof:ndof*nn) = -dN2dxy(:,2)';
                B1(3,3:ndof:ndof*nn) = -2*dN2dxy(:,3)';
                
                B2 = zeros(3,ndof*nn);  % Shear bending
                B2(1,4:ndof:ndof*nn) = dNdxy(:,1)';
                B2(2,5:ndof:ndof*nn) = dNdxy(:,2)';
                B2(3,4:ndof:ndof*nn) = dNdxy(:,2)';
                B2(3,5:ndof:ndof*nn) = dNdxy(:,1)';
                
                Bs = zeros(2,ndof*nn);  % Pure shear
                Bs(1,5:ndof:ndof*nn) = N';
                Bs(2,4:ndof:ndof*nn) = N';
                
                B{1} = B0; B{2} = B1; B{3} = B2;
                k = 0;
                for i = 1:3
                    for j = 1:3
                        k = k + B{i}'*Db{i,j}*B{j}*gwt*detJ1;
                    end
                end
                k = k + Bs'*Ds*Bs*gwt*detJ1;
                
                K(sctrK,sctrK) = K(sctrK,sctrK) + k;
            end
        end
    end
end
end
