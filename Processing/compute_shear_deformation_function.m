function [fz, dfz, ddfz] = compute_shear_deformation_function(z,h,shear_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute shear deformation function at point z %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% === Shear deformation function ===
switch shear_func
    case 1 % Reissner plate model
        fz = 0; % 3DOF with no higher-order term
        dfz = 0;
        ddfz = 0;
    case 2 % Reddy et al. (1984)
        fz = z-4*z^3/(3*h^2);
        dfz = 1-4*z^2/h^2;
        ddfz = -4*2*z/h^2;
    case 3 % Shimpi et al. (2002)
        fz  = z-5*z^3/(3*h^2)+z/4;
        dfz = 1-5*z^2/h^2+1/4;
        ddfz = -5*2*z/h^2;
    case 4 % Nguyen-Xuan et al. (2013)
        fz = z-1/8*z-2/h^2*z^3+2/h^4*z^5;
        dfz = 1-1/8-6/h^2*z^2+10/h^4*z^4;
        ddfz = -6/h^2*2*z+10/h^4*4*z^3;
    case 5 % Nguyen et al. (2017)
        fz = z-9*z+(10/h^2)*z^3+6/(5*h^4)*z^5+8/(7*h^6)*z^7;
        dfz = 1-9+(10/h^2)*3*z^2+6/(5*h^4)*5*z^4+8/(7*h^6)*7*z^6;
        ddfz = (10/h^2)*3*2*z+6/(5*h^4)*5*4*z^3+8/(7*h^6)*7*6*z^5;
    case 6 % Nguyen et al. (2016)
        fz = z - (17*z^3)/(10*h^2) + (22*z^5)/(25*h^4);
        dfz = 1 - (17*3*z^2)/(10*h^2) + (22*5*z^4)/(25*h^4);
        ddfz = (17*3*2*z)/(10*h^2) + (22*5*4*z^3)/(25*h^4);
    case 7 % Thai et al. (2014)
        fz = atan(sin(pi/h*z));
        dfz = pi/h*cos(pi/h*z)/(1+(sin(pi/h*z))^2);
        ddfz = - (pi^2*sin((pi*z)/h))/(h^2*(sin((pi*z)/h)^2 + 1)) - (2*pi^2*cos((pi*z)/h)^2*sin((pi*z)/h))/(h^2*(sin((z*pi)/h)^2 + 1)^2);
    case 8 % Touratier et al. (1991)
        fz = (h/pi)*sin(pi*z/h);
        dfz = cos(pi*z/h);
        ddfz = -(pi/h)*sin(pi*z/h);
    case 9 % Hoang and Pham (2023)
        fz = h/3 * (atan(2*z/h) + 1/(2*pi)*sin(2*pi*z/h));
        dfz = 1/3 * (2/(1+4*z^2/h^2) + cos(2*pi*z/h));
        ddfz = 1/3 * (-(16*z)/(h^2*((4*z^2)/h^2 + 1)^2) - 2*pi/h*sin(2*pi*z/h));
    case 10 % Newly developed
        fz = z + h/(2*pi)*sin(2*pi*z/h) - h*atan(sin(pi/h*z));
        dfz = 1 + cos(2*pi*z/h) - pi*cos(pi/h*z)/(1+(sin(pi/h*z))^2);
        ddfz = - 2*pi/h*sin(2*pi*z/h) + h/pi*(pi^2*sin((pi*z)/h))/(h^2*(sin((pi*z)/h)^2 + 1)) + (2*pi^2*cos((pi*z)/h)^2*sin((pi*z)/h))/(h^2*(sin((z*pi)/h)^2 + 1)^2);
    otherwise
        disp('Error')
        pause
end

end
