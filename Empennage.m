% Martin Lewis Corona 2024
% Module 10 Aircraft Engineering
% Empennage sizing 
clc
clear
close all

% Given data
W = 100*9.81;                                              % Total weight estimate (N)
Clmax = 2.23;                                               % Maximum lift coefficient
rho = 1.225;                                               % Air density at sea level (kg/m^3)
V = 41.667;                                                % Velocity (m/s)
c_r = 0.85;                                                % Root chord (m)
c_t = c_r*0.3;                                               % Tip chord (m)
lambda = c_t/c_r;                                          % Taper ratio
Span = 2.26*2;                                               % Main Wing span (m)
AreaW = c_t*Span+(c_r-c_t)*Span*0.5;                             % Main wing area (m^2)
AR = (Span^2)/AreaW;                                       % Main wing aspect ratio
MAC = c_r*(2/3)*((1+lambda+lambda^2)/(1+lambda));          % Mean aerodynamic cord

% Assumed data for empennage sizing
%WoF = 0.5;                                                 % Ratio between main wing and fuselage length (assumed as 0.5)
%FLength = Span*WoF;                                        % Fuselage length(m) 
FLength = 3;
LeverLength = 0.6*FLength;                                 % Distance between main wing and tail (0.6 is the value for a tractor propeller engine)
CVertVolume = 0.04;                                        % Vertical tail Volume ratio (assumed for homebuilt props)
CHorVolume = 0.5;                                          % Horizontal tail volume ratio (assumed for homebuilt props)

Vstall = sqrt(W/(0.5*rho*Clmax*AreaW));                    % Stall speed at sea level


% Empennage sizing
AreaVert = (AreaW * Span * CVertVolume)/LeverLength ;      % Surface of vertical stabilizer
AreaHor =  (AreaW * MAC * CHorVolume)/LeverLength ;        % Surface of horizontal stabilizer

AR2= 4;                                                    % Aspect ratio of horizontal stabilizer (usually half of main wing to allow higher aoa without stall)
AR3 = 2;                                                   % Aspect ratio of vertical stabilizer (half of the horizontal since it only has one part, while the horizontal stabilizer hase two parts)
VSpan = sqrt(AreaVert*AR3);                                % Vertical stabilizer span (height)
HSpan = sqrt(AreaHor*AR2);                                 % Horizontal stabilizer span

ctH_r = (2*AreaHor)/((lambda+1)*HSpan);                    % Horizontal stabilizer cord root
ctH_t=lambda*ctH_r;                                        % Horizontal stabilizer cord tip

ctV_r = (2*AreaVert)/((lambda+1)*VSpan);                   % Vertical stabilizer cord root
ctV_t=lambda*ctV_r;                                        % Vertical stabilizer cord tip


% Elevator sizing 
AreaEL = 0.25*AreaHor;                                      % Area elevator
cEL = AreaEL/HSpan;                                        % Cord elevator (the elevator will be along the whole stabilizer length)

%Rudder sizing
AreaRud = 0.4*AreaVert;                                    % Area rudder
cRud = AreaRud/VSpan;                                      % Cord rudder (the rudder will be along the whole stabilizer length)


a = [0 1 2 3 4 5 6 7 8 9 10 11];                                                     % Angle of attack in degrees
alphaw = a*(pi / 180) ;                                    % Angle of attack of the main wing in RADIANS (angle in degrees *(pi/180) 
it = 2*(pi/180);                                           % Angle of attack of the tail in RADIANS (angle in degrees *(pi/180) FREE CHOICE
at = 6.3;                                                  % Lift curve of the tail slope
aW = 6.334;                                                % Lift curve of the main wing
Cmacw = 0.0152;                                            % Moment coefficient around the aerodynamic center of the main wing
h = 0.12;                                                   % Distance between main wing and center of gravity FREE CHOICE
hac = 0.25*c_r;                                            % Aerodynamic center
for i = 1: 12

Cl = [0.0925 0.0128 0.3136 0.4238 0.534 0.6424 0.7480... 
    0.8624 0.9938 1.1293 1.2608 1.3337];                   % Values of lift coefficient for angles of attack going from 0 to 11 degrees)
deda(i) = (2*Cl(a(i)+1))/(pi*AR);                                % Derivative of the downwash angle in respect to the angle of attack 
epsilon0 = (2*0.0925)/(pi*AR);                             % Downwash for a 0 angle of attack
Cm0cg = Cmacw+CHorVolume*at*(it+epsilon0);
Cmcg(i) = Cmacw+aW*alphaw(i)*(h-hac-CHorVolume*(at/aW)* ...
    (1-deda(i)))+CHorVolume*at*(it+epsilon0);                 % Moment coefficient around the center of gravity
dcdalfa = aW*(h-hac-CHorVolume*(at/aW)*(1-deda(i)));          % Derivative of the moment coefficient in respect to the angle of attack
hn = hac+CHorVolume*(at/aW)*(1-deda(i));                      % Neutral point location
deltaH = hn-h;                                             % Static margin

end
plot(a,Cmcg,'Color','r','LineWidth',1)
xlabel('AoA','FontSize',16)
ylabel('Cmcg','FontSize',16)
title('Coefficient of moment around center of gravity for different AoA','FontSize',12)
grid on
% Longitudinal stability analysis
a = 4;                                                     % Angle of attack in degrees
alphaw = a*(pi / 180) ;                                    % Angle of attack of the main wing in RADIANS (angle in degrees *(pi/180) 
it = 2*(pi/180);                                           % Angle of attack of the tail in RADIANS (angle in degrees *(pi/180) FREE CHOICE
at = 6.3;                                                  % Lift curve of the tail slope
aW = 6.334;                                                % Lift curve of the main wing
Cmacw = 0.0152;                                            % Moment coefficient around the aerodynamic center of the main wing
h = 0.12;                                                  % Distance between main wing and center of gravity FREE CHOICE
hac = 0.25*c_r;                                            % Aerodynamic center


Cl = [0.0925 0.0128 0.3136 0.4238 0.534 0.6424 0.7480... 
    0.8624 0.9938 1.1293 1.2608 1.3337];                   % Values of lift coefficient for angles of attack going from 0 to 11 degrees)
deda = (2*Cl(a+1))/(pi*AR);                                % Derivative of the downwash angle in respect to the angle of attack 
epsilon0 = (2*0.0925)/(pi*AR);                             % Downwash for a 0 angle of attack
Cm0cg = Cmacw+CHorVolume*at*(it+epsilon0);
Cmcg = Cmacw+aW*alphaw*(h-hac-CHorVolume*(at/aW)* ...
    (1-deda))+CHorVolume*at*(it+epsilon0);                 % Moment coefficient around the center of gravity
dcdalfa = aW*(h-hac-CHorVolume*(at/aW)*(1-deda));          % Derivative of the moment coefficient in respect to the angle of attack
hn = hac+CHorVolume*(at/aW)*(1-deda);                      % Neutral point location
deltaH = hn-h;                                             % Static margin


% Drag calculation
e = 0.8;                                                   % Oswald efficiency
Cltail = it*at;                                            % Lift coefficient of tail
Cd0 = 0.01;                                                % Drag coefficient for 0 angle of attack
Cdi = (Cltail^2)/(pi*AR2*e);                               % Lift induced drag coefficient
DragHor = 0.5 * (Cdi+Cd0) * rho * V^2 * AreaHor;           % Drag produced by the Horizontal stabilizer
DragVert = 0.5 * (Cd0) * rho * V^2 * AreaVert;             % Drag produced by the vertical stabilizer
Drag = DragVert + DragHor;                                 % Total drag produced by the tail 
Lifttail = 0.5 * Cltail * rho * V^2 * AreaHor;             % Lift generated by the tail