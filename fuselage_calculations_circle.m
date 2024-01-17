clc
clear
close

% changeable data
weight = 104;                       % weight airplane               [kg]
n_s = linspace(1, 8, 8);            % number of stringers
d = 0.25;                            % diameter of fuselage          [m]
t = 0.001;                          % thickness plates              [m]
A_str = 0.0001;                  % area stringers                [m^2]

% calculated values from data
angle = 2 * pi /length(n_s);        % angle between stringers       [rad]
force = weight*9.81;                % downward force due to weight  [N]
r = d/2;                            % radius                        [m]
moment = force*r;                   % moment created by force       [Nm]
circ = d * pi;                      % circumfrence fuselage         [m]
b = circ / length(n_s);             % distance between stringers    [m]
heights = cos((n_s-1)*angle)*r;     % heights from middelpoint      [m]

%zero values for the loop
I = 0;                              % ero value for inertia         [m^4]  
b_all = zeros(1, length(n_s));      % zero array for boom           [m^2]
stress_z = zeros(1, length(n_s));   % zero array for bendingstress  [N/m^2]
q_s = zeros(1, length(n_s));        % zero array for shear stress   [N/m

%calculate boom values
for n = 1:(length(n_s)/2)
    if n == 1
        b_1 = A_str + ((t*b)/6)*(2 + heights(n+1)/heights(n)) * 2;
        b_all(1) = b_1;
        b_all((length(n_s)/2)+1) = b_1;
    else 
        b_n = A_str + ((t*b)/6)*(2 + heights(n+1)/heights(n)) + ((t*b)/6)*(2 + heights(n-1)/heights(n));
        b_all(n) = b_n;
        b_all(length(n_s)-(n-2)) = b_n;
    end
    b_all((length(n_s)/4)+1) = 0;
    b_all((length(n_s)-(length(n_s)/4)+1)) = 0;
end

%calculate inertia
for n = 1:(length(n_s)/4)
    if n == 1
        I_1 = 2 * b_all(n) * heights(n)^2;
        I = I + I_1;
    else
        I_n = 4 * b_all(n) * heights(n)^2;
        I = I + I_n;
    end
end

%calculate bending stress
for n = 1:(length(n_s)/2)
    if n == 1
        stress_z1 = (moment*heights(n))/I;
        stress_z(n) = stress_z1;
        stress_z((length(n_s)/2)+1) = stress_z1;
    else
        stress_zn = (moment*heights(n))/I;
        stress_z(n) = stress_zn;
        stress_z(length(n_s)-(n-2)) = stress_zn;
    end
    stress_z((length(n_s)/4)+1) = 0;
    stress_z((length(n_s)-(length(n_s)/4)+1)) = 0;
end

%calculate shear stress
for n = 1:(length(n_s))
    if n == 1
        q_s1 = (-force/I)*(b_all(n))*heights(n);
        q_s(n) =q_s1;
    else
        q_sn = (-force/I)*(b_all(n))*heights(n);
        q_s(n) = q_s(n-1) + q_sn;
    end
end

Stess_max = max(stress_z)

