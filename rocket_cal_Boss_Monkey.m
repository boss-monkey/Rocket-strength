%% rocket ansys
% A Boss Monkey code 2024.04.01
% Baiht0201@163.com
% 仅供教学，禁止商用
%% main data
m_t0 = 48856.2;
delta = deg2rad(3.5);
alpha = atan(40/1100);
T_e = 1160000;
T_x = T_e*cos(delta);
T_y = T_e*sin(delta);
rho = 0.0084634;
a = 308.2997;
velocity = (1100^2+40^2)^0.5;
M = velocity/a;

%% Geometric data
x0 = 0; x1 = 2.2; x2 = 5.7; x3 = 7.4; x4 = 8.06;
x5 = 8.41; x6 = 9.4; x7 = 10.9; x8 = 13.3; x9 = 14;
x10 = 17.3; x11=19.5; x12 = 27.8; x13 = 29; x14=30.9;
scale = 1000;
N_1 = x14*scale+1;

x = linspace(x0,x14,N_1);
x_frame = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14];
r1 = 1.6/2; r2 = 2.4/2; S1 = pi*r1^2; S2 = pi*r2^2;
beta_1 = atan(r1/x1); beta_2 = atan((r2-r1)/(x8-x7)) ;
xi_t0 = 21.556;
I_t0 = 0.222*10^7;

%% Mass
% Distribution mass
q_m = (x.*60/x1) .*(x<x1) + 50 *(x>=x1 & x<x2)+ 85.*(x>=x2 & x<x7)+...
        (50+(x-x7)*(115-50)/(x8-x7)).*(x>=x7 & x<x8)+100.*(x>=x8 & x<x11)+...
        130.*(x>=x11 & x<x12)+150.*(x>=x12 & x<=x14);

l_r1 = 1181/(890*pi*r1^2); l_o1 = 1020/(1450*pi*r1^2);
l_r2 = 6321.4/(790*pi*r2^2); l_o2 = 26462.8/(1144*pi*r2^2);    
q_fuel = 0*(x>=x1 & x<x3)+(890*pi*r1^2)*(x>=x3 & x<x4)+(1450*pi*r1^2)*(x>=x4 & x<x5)+...
         0*(x>=x5 & x<(x10-l_r2))+(790*pi*r2^2)*(x>=(x10-l_r2) & x<x10)+0*(x>=x10 & x<(x12-l_o2))+...
         (1144*pi*r2^2)*(x>=(x12-l_o2) & x<x12)+0*(x>=x12 & x<x14);
       
% Concentrated mass
m_concentrated_0 = [1289, 262, 80, 100, 80, 120, 0, 690, 100, 100, 100, 100, 1148];
m_fuel_0 = [0,0,0,1181,(1020+824),0,0,0,0,(6321.4+2383),0,(26462.8+3450),0];
m_fule_bilge = [0,0,0,0,824,0,0,0,0,2383,0,3450,0];  %fuel mass in the bilge of the tank

xi_center_mc = [(x1+1.8),(x2+0.85),x3,x4,x5,(x6+1.0),x7,(x8-0.8),x9,x10,x11,x12,(x13+1)]; %mass center 修改了发动机距离数据
xi_center_mb = [x1,x2,x3,x4,(x5+0.175),x6,x7,x8,x9,(x10+0.375),x11,(x12+0.375),x13];

I_fule_bilge = [0,0,0,0,(12+3.04),0,0,0,0,(141+8.55),0,(205+8.55),0]; 


%% aerodynamic
q_dynamic = 1/2*rho*velocity^2;
% x_axis
q_ax_p = (-2*pi*q_dynamic*(alpha^2+2*beta_1^2)* x*r1/x1 *tan(beta_1)).*(x<x1) + 0 *(x>=x1 & x<x7)+...
         (-2*pi*q_dynamic*(alpha^2+2*beta_2^2)* (r1+(x-x7)*(r2-r1)/(x8-x7))*tan(beta_2)).*(x>=x7 & x<=x8)+...
         0*(x>x8 & x<=x14);
   
X_p = q_dynamic*(alpha^2+2*beta_1^2)*(pi*r1^2)+q_dynamic*(alpha^2+2*beta_2^2)*(pi*(r2^2-r1^2));   
X_f = 0.3*X_p;
F = (r1*x1)/2+r1*(x7-x1)+(r1+r2)*(x8-x7)/2+r2*(x14-x8);
q_ax_f = (-X_f/F * x * r1/x1).*(x<x1) + (-X_f/F * r1) *(x>=x1 & x<x7)+...
         (-X_f/F * (r1+(x-x7)*(r2-r1)/(x8-x7))).*(x>=x7 & x<=x8)+...
         (-X_f/F * r2)*(x>x8 & x<=x14);
  
q_ax = q_ax_p+q_ax_f;
X_tail = 0.2*(X_p+X_f);
X_a = X_p+X_f+X_tail;

% y_axis
Y_nose1 = (4*pi*q_dynamic*alpha*beta_1*r1)*x1/2;
Y_body1 = 1.5*alpha^2*(x7-x1)/(2*r1)*q_dynamic*S1;
Y_nose2 = (4*pi*q_dynamic*alpha*beta_2)*(r2+r1)*(x8-x7)/2;
Y_body2 = 1.5*alpha^2*(x14-x8)/(2*r2)*q_dynamic*S2;
Y_a = Y_nose1+Y_body1+Y_nose2+Y_body2;
q_ay = (4*pi*q_dynamic*alpha*beta_1* x*r1/x1).*(x<x1) + Y_body1/(x7-x1) *(x>=x1 & x<x7)+...
       (4*pi*q_dynamic*alpha*beta_2* (r1+(x-x7)*(r2-r1)/(x8-x7))).*(x>=x7 & x<=x8)+...
         Y_body2/(x14-x8) *(x>x8 & x<=x14);
 
xi_nose1 = 2/3*x1;
xi_body1 = x1+(x7-x1)/2;
xi_nose2 = x7+(x8-x7)*(r1+2*r2)/(r1+r2)/3;
xi_body_2 = x8+(x14-x8)/2;

xi_a = 1/Y_a*(Y_nose1*xi_nose1+Y_body1*xi_body1+Y_nose2*xi_nose2+Y_body2*xi_body_2);
%% acceleration
g = 9.81;
n_x = (T_x-X_a)/(m_t0*g);
n_y = (T_y+Y_a)/(m_t0*g);
epsilon = (Y_a*(xi_t0-xi_a)+T_y*(xi_t0-xi_center_mc(13)))/I_t0;
%%
m_concentrated = zeros(1,N_1);   j=1;
for i = 2:N_1   
      m_concentrated(i) =m_concentrated(i-1);
    if j<=13
        if i/scale == x_frame(j)
            m_concentrated(i) = m_concentrated(i)-m_concentrated_0(j);  
            j=j+1;
        end
    end
end

m_distributed= zeros(1,N_1); 
for i = 2:N_1   
      m_distributed(i) =m_distributed(i-1)-q_m(i)/scale;
end

m_structure = m_distributed+m_concentrated;

m_fuel = zeros(1,N_1);   j=1;
for i = 2:N_1   
      m_fuel(i) =m_fuel(i-1);
    if j<=13
        if i/scale == x_frame(j)
            m_fuel(i) = m_fuel(i)-m_fuel_0(j);  
            j=j+1;
        end
    end
end

%% force to x aixs
p_r1 = 0.24*10^6*pi*r1^2; p_o1 = 0.28*10^6*pi*r1^2;
p_r2 = 0.20*10^6*pi*r2^2; p_o2 = 0.24*10^6*pi*r2^2;
f_fuel_p0 = [0,0,p_r1,(-p_r1+p_o1),-p_o1,0,0,0,p_r2,-p_r2,p_o1,-p_o1,0];
f_fuel_p = zeros(1,N_1);   j=1;  % Pressure at both ends of fuel tank
for i = 2:N_1   
      f_fuel_p(i) =f_fuel_p(i-1);
    if j<=13
        if i/scale == x_frame(j)
            f_fuel_p(i) = f_fuel_p(i)+f_fuel_p0(j);  
            j=j+1;
        end
    end
end

f_ax = zeros(1,N_1);
for i = 2:N_1   
      f_ax(i) =f_ax(i-1)+q_ax(i)/scale;
      if i/scale == x14
            f_ax(i) = f_ax(i)-X_tail;  
      end
end

f_tx = zeros(1,N_1);
for i = 2:N_1   
      f_tx(i) =f_tx(i-1);
      if i/scale == x_frame(13)
            f_tx(i) = f_tx(i)+T_x;  
      end
end

f_distributed_x = m_distributed * n_x *g;
f_concentrated_x = m_concentrated * n_x *g;
f_fuel_rho_x = m_fuel* n_x *g;
N_x = f_distributed_x+f_concentrated_x+f_fuel_rho_x+f_fuel_p+f_ax+f_tx;
figure(1)
plot(x, N_x)
%% force to y aixs
f_distributed_y= zeros(1,N_1);
for i = 2:N_1   
      f_distributed_y(i) =f_distributed_y(i-1)-q_m(i)/scale*(n_y+epsilon/g*(xi_t0-x(i)))*g;
end

f_fuel_rho_y = zeros(1,N_1);
for i = 2:N_1   
      f_fuel_rho_y(i) =f_fuel_rho_y(i-1)-q_fuel(i)/scale*(n_y+epsilon/g*(xi_t0-x(i)))*g;
end

f_concentrated_y = zeros(1,N_1); j = 1;
for i = 2:N_1   
      f_concentrated_y(i) = f_concentrated_y(i-1);
      if j<=13
          if i/scale == x_frame(j)
                f_concentrated_y(i) = f_concentrated_y(i)- m_concentrated_0(j)*(n_y+epsilon/g*(xi_t0-xi_center_mc(j)))*g;  
                j = j+1;
          end
      end
end

f_fule_bilge_y = zeros(1,N_1); j = 1;
for i = 2:N_1   
      f_fule_bilge_y(i) = f_fule_bilge_y(i-1);
      if j<=13
          if i/scale == x_frame(j)
                f_fule_bilge_y(i) = f_fule_bilge_y(i)- m_fule_bilge(j)*(n_y+epsilon/g*(xi_t0-xi_center_mb(j)))*g;  
                j = j+1;
          end
      end
end

f_ay = zeros(1,N_1);
for i = 2:N_1   
      f_ay(i) =f_ay(i-1)+q_ay(i)/scale;
end

f_ty = zeros(1,N_1);
for i = 2:N_1   
      f_ty(i) =f_ty(i-1);
      if i/scale == x_frame(13)
            f_ty(i) = f_ty(i)+T_y;  
      end
end
Q_y = f_distributed_y+f_concentrated_y+f_fuel_rho_y+f_fule_bilge_y+f_ay+f_ty;
figure(2)
plot(x,Q_y)
%% moment
M_Qy= zeros(1,N_1);
for i = 2:N_1   
      M_Qy(i) =M_Qy(i-1)-Q_y(i)/scale;
end

M_concentrated_y = zeros(1,N_1); j = 1;
for i = 2:N_1   
      M_concentrated_y(i) = M_concentrated_y(i-1);
      if j<=13
          if i/scale == x_frame(j)
              M_concentrated_y(i) = M_concentrated_y(i)- ...
                                    m_concentrated_0(j)*(n_y+epsilon/g*(xi_t0-xi_center_mc(j)))*g*(x_frame(j)-xi_center_mc(j));  
              j = j+1;
          end
      end
end

M_fule_bilge_y = zeros(1,N_1); j = 1;
for i = 2:N_1   
      M_fule_bilge_y(i) = M_fule_bilge_y(i-1);
      if j<=13
          if i/scale == x_frame(j)
              M_fule_bilge_y(i) = M_fule_bilge_y(i)- ...
                                  m_fule_bilge(j)*(n_y+epsilon/g*(xi_t0-xi_center_mb(j)))*g*(x_frame(j)-xi_center_mb(j))+...
                                  -epsilon*I_fule_bilge(j);  
              j = j+1;
          end
      end
end

M_ty = zeros(1,N_1);
for i = 2:N_1   
      M_ty(i) =M_ty(i-1);
      if i/scale == x_frame(13)
            M_ty(i) = M_ty(i)-T_y*(x_frame(13)-xi_center_mc(13));  
      end
end

M_z = M_Qy+M_concentrated_y+M_fule_bilge_y+M_ty;
figure(3)
plot(x,M_z)
