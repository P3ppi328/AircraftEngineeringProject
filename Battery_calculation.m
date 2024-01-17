clc
clear
close all

%% Power consumption script

n = 100;
rho = [1.225 * 1, 1.225 * 0.6689];
alt = [0, 4000];

Q = @(vel,y) 1/2 * interp1(alt, rho, y) * vel^2;

W = 100 * 9.81; % N
S = 2.5; % m^2

e = 0.8;
AR = 8.16;
ClMax = 2.23;
vStall = sqrt( (2 * W) ./ (S * rho(1) * ClMax) ); % stall velocity
vLO = 1.2 * vStall; % lift-off velocity
vTD = 1.3 * vStall; % touch-down velocity
headWind = 12.3;
vMax = 150/3.6 + headWind;
vCruise = 150/3.6;


flying = true;
x = 0;
y = 0;
v = vLO;
dt = 1;
Pr = 0;
Pt = 0;
runwayLength = 150;
cruiseAlt = 1000;
dropAlt = 50; % 50 metres altitude
dropPos = 50e3; % 50km away

loitAlt = 250;
stage = 'takeOff';
while flying
    % v(end) = v(end) + headWind;
    q = Q(v(end),y(end)); % Local Dynamic Pressure
    Cl = W / (q * S);
    Cd = Cl / 15;
    C_D0 = 0.00857;
    vStall = sqrt( (2 * W) / (S * interp1(alt, rho, y(end)) * ClMax) ); % stall velocity
    Tr = q * S * Cd; % N Thrust required
    propEff = 0.83; % Estimation
    Pa = 7.5e3 * propEff; % W Power available

    switch stage
        case 'takeOff'
            if x(end) < runwayLength
                v(end + 1) = vLO; % Special case IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pt(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                %Pr(end+1) = Pa;
                RC = (Pa - Pt(end)) / W;
                y(end + 1) = y(end);
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
                Pr(end+1) = Pa;
            else
                stage = 'initialClimb';
            end
        case 'initialClimb'
            if y(end) < cruiseAlt
                v(end + 1) = vStall + headWind;
               % Pr(end + 1) = Tr * v(end); % W Power required
                Pt(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                %Pr(end+1) = Pa;
                RC = (Pa - Pt(end)) / W;
                y(end + 1) = y(end) + RC * dt;
                x(end + 1) = x(end) + (v(end) - headWind) * dt;
                Pr(end+1) = Pa;
            else
                stage = 'initialCruise';
            end

        case 'initialCruise'
            if x(end) < 47.466e3
                v(end + 1) = vCruise + headWind; % IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end);
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
            else
                stage = 'initialDescent';
            end

        case 'initialDescent'
            %descentRate = sqrt(vMax^2 - vStall^2);
            if y(end) > dropAlt
                v(end + 1) = vStall + headWind; % IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                %descentRate = -(Pa - Pr(end)) / W;
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end) - RC*dt;%;descentRate * dt;
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
            else
                stage = 'dropOff';
            end
        case 'dropOff'
            if x(end) < 60e3
                v(end + 1) = vStall + headWind; % IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end);
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
            else
                W = 75 * 9.81; 
                stage = 'secondClimb';
            end

        case 'secondClimb'
            if y(end) < cruiseAlt
                v(end + 1) = vStall + headWind; % IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pt(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                % Pr(end+1) = Pa;
                RC = (Pa - Pt(end)) / W;
                y(end + 1) = y(end) + RC * dt;
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
                Pr(end+1) = Pa;
            else
                stage = 'secondCruise';
            end

        case 'secondCruise'
            if x(end) < 110e3
                v(end + 1) = vCruise + headWind; % IAS
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end);
                x(end + 1) = x(end) + (v(end) - headWind) * dt; % GS
            else
                stage = 'secondDescent';
            end

        case 'secondDescent'
            %descentRate = sqrt(vMax^2 - vStall^2);
            if y(end) > loitAlt
                v(end + 1) = vStall + headWind;
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                %descentRate = -(Pa - Pr(end)) / W;
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end) - RC*dt;%descentRate * dt;
                x(end + 1) = x(end) + (v(end) - headWind) * dt;
            else
                stage = 'loiter';
            end

        case 'loiter'
            vLoiter = vStall;%0.9 * vCruise + headWind;
            loiterTime = 15*60; % loiter time in seconds

            for i = 1:loiterTime
                v(end + 1) = vLoiter;
                %Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end);
                x(end + 1) = x(end) + (v(end) - headWind) * dt;
            end
            stage = 'finalDescent';

        case 'finalDescent'
            %descentRate = sqrt(vMax^2 - vStall^2);
            if y(end) > 0
                v(end + 1) = vStall + headWind;
               % Pr(end + 1) = Tr * v(end); % W Power required
                Pr(end+1) = (0.5 * interp1(alt, rho, y(end)) * (v(end).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 *interp1(alt, rho, y(end))* v(end) * S)); 
               %descentRate = -(Pa - Pr(end)) / W;
                RC = (Pa - Pr(end)) / W;
                y(end + 1) = y(end) -RC*dt;%descentRate * dt;
                x(end + 1) = x(end) + (v(end) - headWind) * dt;
            else
                stage = 'landed';
                Pr(end + 1) = 0;
                v(end + 1) = 0;
                x(end + 1) = x(end) + runwayLength;
                y(end + 1) = 0;
            end

        otherwise
            flying = false;
    end
end
flightTime = length(x) * dt;
wattHour = sum(Pr) / 3600;
battDensity = 200; % Wh/kg
battMass = wattHour / battDensity;
battVolume = 330; % Wh/L
battVolume = wattHour / battVolume;
% Motor https://www.elmofo.com.au/qs138-90h-v3-electric-motor.html
disp(['The plane has covered a distance of ' num2str(x(end) / 1e3, 3) ' km in ' num2str(flightTime /60,3) ' minutes and used ' num2str(wattHour / 1e3,2) ' kWh'])
%v = linspace(10, 60, 61);
%Pr = []; % Initialize Pr array
%for i = 1:61
%    Pr(end+1) = (0.5 * 1.225 * (v(i).^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 * 1.225 * v(i) * S));
%end
%plot(v, Pr)
%xlabel('Velocity (m/s)')
%ylabel('Power Required')
%title('Power Required vs Velocity')
%grid on
%KK = (0.5 * 1.225 * (54.^3) * S * C_D0) + ((W^2) ./ (pi * e * AR * 0.5 * 1.225 * 54 * S))
%% plotting
figure
subplot(2,1,1)
%yyaxis left
plot(x,y,'LineWidth',2)
ylabel('Altitude [m]')
hold on

%yyaxis right
%plot(x,v,'-','LineWidth',1.3,'DisplayName','Velocity')
%ylabel('Velocity [m/s]')
xlabel('Distance [m]')
%legend('Flight path','Velocity','Location','southwest')
legend('Flight path','Location','southwest')


subplot(2,1,2)
%yyaxis right
%plot(v,'-','LineWidth',1.3)
hold on
%ylabel('Velocity [m/s]')
xlabel('Time [s]')

%yyaxis left
plot(Pr,'-','LineWidth',1.3)
ylabel('Power [W]')
xlabel('Time [s]')
legend('Power','Velocity','Location','southwest')
