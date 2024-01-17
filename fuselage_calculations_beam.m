clc
clear
close


%mass data
m_m = 3.1;                                  % mass of motor     [kg]
m_w = 27.42;                                % mass of wings     [kg]
m_p = 25;                                   % mass of payload   [kg]
m_b = 20.8;                                 % mass of battery   [kg]
m_t = 7.5;                                    % mass of tail      [kg]
m_g = 10.8;                                 % mass of gear      [kg]

%dimensions
D = 0.25;                                    % outer diameter fuselage [m]
d = D - 0.001;                              % inner diameter fuselage [m]
L = 3;                                      % length fuselage         [m]
t = 0.001;                                  % thickness fuselage wall [m]
A_str = 0.0001;                             % area stringer           [m^2]
density = 2780;                             % density of aluminium    [kg/m^2]

%mass calculations
V_skin = (D^2-d^2)/4*pi*L + (d/2)^2*pi*t;   % volume of skin fuselage [m^2]
m_skin = V_skin * density;                  % mass of skin fuselage   [kg]
V_stringer = 0.0001 * L * 8;                % volume stringers        [m^2]
m_stringer = V_stringer * density;          % mass stringers          [kg]
m_f = m_skin + m_stringer;                  % total mass fuselage     [kg]

mass = [m_m m_w m_g m_f m_p m_b m_t];       % all masses in array     [kg]
mass_total = sum(mass);                     % total mass aeroplane    [kg]

%weight loads
W = mass *9.81;                             % gravitational forces    [N]
W_total = sum(W);                           % total grav. force       [N]

%data loads
L_w = W_total;                              % wing Lift               [N]
L_t = L_w - W_total;                        % tail lift               [N]

%moment of inertia
I = ((D^4 - d^4)*pi)/64;                    % inertia tube            [m^4]

%centre of gravity
x_m = 0.037;
x_w = 1.2;
x_g = x_w;
x_b = 0.6;
x_t = L;
x_f = L/2;
X_cg = ((m_f * x_f) + (m_m * L) + (m_w *x_w) + (m_b*x_b) + ... 
    + (m_t*x_t) + (m_g*x_g))/(m_m+m_w+m_b+m_f+m_t + m_g);
x_p = X_cg;

cg = ((m_f * x_f) + (m_m * L) + (m_w *x_w) + (m_b*x_b) + ... 
    + (m_t*x_t) + (m_g*x_g) + (m_p*x_p))/(m_m+m_w+m_b+m_f+m_t + m_g+m_p);

M_ac = W(4)*(x_f - x_w) + W(5) * (x_p-x_w) + W(6) * (x_b - x_w) + (W(6)+L_t) * (x_t-x_w) - W(1) * x_w;

%distance to clear up the upcoming formulas
W_lift = L_w-W(2)-W(3);
a = x_m - x_w;
b = x_b - x_w;
c = x_p - x_w;
d = x_f - x_w;
e = x_t - x_w;

E = 72.4e9;                             % E-modulus Aluminium

%plot
% formulas Moment and displacement for every part of the fuselage
x_03 = linspace(-x_w, x_m-x_w,10);
f_03 = ((x_03).^3/6);
g_03 = ((x_03).^2/2);
M_03 = -M_ac - (W_lift)*x_03  + W(6)*(x_03-(b)) + W(1)*(x_03 - (a));
V_03 = (1/(E*I))*(-M_ac*g_03 - (W_lift*f_03) + W(6)*(f_03-b*g_03) + W(1)*(f_03-a*g_03) ...
    + (W(6)*(b^2/2)+W(1)*(a^2/2))*x_03 - W(6)*(b^3/6)- W(1)*(a^3/6)); 

x_02 = linspace(x_m-x_w, x_b-x_w,10);
f_02 = ((x_02).^3/6);
g_02 = ((x_02).^2/2);
M_02 = -M_ac - (W_lift)*x_02  + W(6)*(x_02-(x_b-x_w));
V_02 = (1/(E*I))*(-M_ac*g_02 - (W_lift*f_02) + W(6)*(f_02-b*g_02) ...
    + (W(6)*(b^2/2))*x_02 - W(6)*(b^3/6));

x_01 = linspace(-x_w + x_b, 0, 10);
f_01 = ((x_01).^3/6);
g_01 = ((x_01).^2/2);
M_01 = -M_ac - (W_lift)*x_01;
V_01 = (1/(E*I))*(-M_ac*g_01 - (W_lift*f_01));

x_1 = linspace(0,x_p-x_w, 10);
f_1 = ((x_1).^3/6);
g_1 = ((x_1).^2/2);
M_1 = -M_ac + (W_lift)*x_1;
V_1 = (1/(E*I))*(-M_ac*g_1 + (W_lift*f_1));

x_2 = linspace(x_p-x_w, x_f-x_w,10);
f_2 = ((x_2).^3/6);
g_2 = ((x_2).^2/2);
M_2 = -M_ac + (W_lift)*x_2  - W(5)*(x_2-(x_p-x_w));
V_2 = (1/(E*I))*(-M_ac*g_2 + (W_lift*f_2) - W(5)*(f_2-c*g_2) ...
     - (W(5)*(c^2/2))*x_2 + W(5)*(c^3/6));

x_3 = linspace(x_f-x_w, x_t-x_w,10);
f_3 = ((x_3).^3/6);
g_3 = ((x_3).^2/2);
M_3 = -M_ac + (W_lift)*x_3  - W(5) * (x_3 - (x_p-x_w)) ...
      - W(4) * (x_3 - (x_f - x_w));
V_3 = (1/(E*I))*(-M_ac*g_3 + (W_lift*f_3) - W(5)*(f_3-c*g_3) - W(4)*(f_3-d*g_3) ...
    - (W(5)*(c^2/2)+W(4)*(d^2/2))*x_3 + W(5)*(c^3/6)+ W(4)*(d^3/6));

x_4 = linspace(x_t-x_w, L,10);
f_4 = ((x_4).^3/6);
g_4 = ((x_4).^2/2);
M_4 = -M_ac + (W_lift)*x_4  - W(5) * (x_4 - (x_p-x_w)) ...
      - W(4) * (x_4 - (x_f - x_w)) - W(7) * (x_4 - (x_t - x_w));
V_4 = (1/(E*I))*(-M_ac*g_4 + (W_lift*f_4) - W(5)*(f_4-c*g_4) - W(4)*(f_4-d*g_4) - W(7)*(f_4-e*g_4) ...
    - (W(5)*(c^2/2)+W(4)*(d^2/2)+W(7)*(e^2/2))*x_4 + W(5)*(c^3/6)+ W(4)*(d^3/6) + W(7)*(e^3/6));

% making plots of the displacement and momement
V = [V_03 V_02 V_01 V_1 V_2 V_3 V_4];
M = [M_03 M_02 M_01 M_1 M_2 M_3 M_4];
x = [x_03 x_02 x_01 x_1 x_2 x_3 x_4];


figure(1)
plot(x,M)

figure(2)
plot(x, V)


%determining maximum stress at wing
displacement = max(V);
sigma = (-M_ac/I)*(-displacement);





