function [E_s, nu_s, rho_s, alpha_s] = compute_CNTR_material(mat_prop, CNT_prop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Material library for CNT reinforcement nanocomposite %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate.mat and Plate.CNT
E_m = mat_prop.E_m; nu_m = mat_prop.nu_m; rho_m = mat_prop.rho_m; alpha_m = mat_prop.alpha_m;
type = CNT_prop.type; f_r = CNT_prop.fr; eta = CNT_prop.eta; mu = CNT_prop.mu;

%% Hill elastic moduli for CNT radius A10
% Unit: m - N - kg - s^2 - Pa - K
if type == 1
    k_r = 30e9; l_r = 10e9; p_r = 1e9; m_r = 1e9; n_r = 450e9; rho_r = 1400;
    alpha_r11 = 3.4584e-6; alpha_r22 = 5.1682e-6;
elseif type == 2
    k_r = 271e9; l_r = 88e9; p_r = 442e9; m_r = 17e9; n_r = 1089e9; rho_r = 1400;
    alpha_r11 = 3.4584e-6; alpha_r22 = 5.1682e-6;
end

%% Check CNT volume fraction ratio 
if f_r == 0
    E_s = E_m; nu_s = nu_m; rho_s = rho_m; alpha_s = alpha_m;
    return
end

%% ===== Calculate nano-composite mechanical properties =====
% --- Properties of matrix material ---
K_m = E_m / (3 - 6*nu_m); G_m = E_m / (2 + 2*nu_m);

% --- Parameters of reinforcement components with agglomerations --- 
alpha_r = (3*K_m + 3*G_m + k_r - l_r) / (3*G_m + 3*k_r);
beta_r = (4*G_m + 2*k_r + l_r) / (15*G_m + 15*k_r) + 4*G_m / (5*G_m + 5*p_r) + ...
         (2*G_m*(3*K_m + G_m) + 2*G_m*(3*K_m + 7*G_m)) / (5*G_m*(3*K_m + G_m) + 5*m_r*(3*K_m + 7*G_m));
delta_r = (3*K_m + 2*G_m - l_r)*(2*k_r + l_r) / (3*G_m + 3*k_r) + 1/3*(n_r + 2*l_r);
eta_r = (8*G_m*p_r)/(5*G_m + 5*p_r) + (8*G_m*m_r*(3*K_m + 4*G_m)) / (15*K_m*(G_m + m_r) + 5*G_m*(7*m_r + G_m)) + ...
        2*(2*G_m + l_r)*(k_r - l_r) / (15*G_m + 15*k_r) + 2/15*(n_r - l_r);

% --- Agglomeration and non-agglomeration phases properties ---
if mu ~= 0  % Not Zero agglomeration
    K_in = K_m + f_r*eta*(delta_r - 3*K_m*alpha_r)/3/(mu - f_r*eta + f_r*eta*alpha_r);
    G_in = G_m + f_r*eta*(eta_r - 2*G_m*beta_r)/2/(mu - f_r*eta + f_r*eta*beta_r);
end
if mu ~= 1  % Not Full agglomeration
    K_out = K_m + f_r*(1-eta)*(delta_r - 3*K_m*alpha_r)/3/(1 - mu - f_r*(1-eta) + f_r*(1-eta)*alpha_r);
    G_out = G_m + f_r*(1-eta)*(eta_r - 2*G_m*beta_r)/2/(1 - mu - f_r*(1-eta) + f_r*(1-eta)*beta_r);
    nu_out = (3*K_out - 2*G_out) / (6*K_out + 2*G_out);
end

% --- Nano-composite properties --- 
if mu == 0  % Zero agglomeration
%     K_in = 0; G_in = 0;
    K_s = K_out; G_s = G_out; 
elseif mu == 1  % Full agglomeration
%     K_out = 0; G_out = 0;
    K_s = K_in; G_s = G_in;
else
    alpha = (1 + nu_out) / (3 - 3*nu_out);
    beta = (8 - 10*nu_out) / (15 - 15*nu_out);

    K_s = K_out * (1 + mu*(K_in/K_out-1) / (1 + alpha*(1-mu)*(K_in/K_out-1)));
    G_s = G_out * (1 + mu*(G_in/G_out-1) / (1 + beta*(1-mu)*(G_in/G_out-1)));
end

E_s = 9*K_s*G_s / (3*K_s + G_s);
nu_s = (3*K_s - 2*G_s) / (6*K_s + 2*G_s);
rho_s = f_r*rho_r + (1-f_r)*rho_m;

%% ===== Calculate nano-composite thermal properties =====
% --- Properties of random orientation CNT phase ---
K_rRO = (4*k_r + 4*l_r + n_r) / 9;
alpha_rRO = ((n_r + 2*l_r)*alpha_r11 + 2*(l_r + 2*k_r)*alpha_r22) / (4*k_r + 4*l_r + n_r);

% --- Agglomeration and non-agglomeration phases properties ---
if mu ~= 0  % Not Zero agglomeration
    alpha_in = alpha_m + f_r*eta*(alpha_rRO - alpha_m) / (mu + (mu - f_r*eta)*(K_m/K_rRO-1)*(4*G_m/(4*G_m+3*K_m)));
end
if mu ~= 1  % Not Full agglomeration
    alpha_out = alpha_m + f_r*(1-eta)*(alpha_rRO - alpha_m) / (1 - mu + (1 - mu - f_r*(1-eta))*(K_m/K_rRO-1)*(4*G_m/(4*G_m+3*K_m)));
end

% --- Nano-composite properties --- 
if mu == 0  % Zero agglomeration
%     alpha_in = 0;
    alpha_s = alpha_out;
elseif mu == 1  % Full agglomeration
%     alpha_out = 0;
    alpha_s = alpha_in;
else
    alpha_s = alpha_out + K_in/K_out*mu*(alpha_in - alpha_out) / (1 + (K_in/K_out-1)*((1 - mu)*(1+nu_out)/3/(1-nu_out) + mu));
end

end
