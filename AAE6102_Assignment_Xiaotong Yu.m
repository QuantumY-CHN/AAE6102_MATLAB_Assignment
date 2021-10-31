%% AAE6102 Satellite Communication and Navigation
%  Assignment
%  Author: Xiaotong Yu 21065068R
%  Date: 30/10/2021
%  ****************************************************

%%  Load the data
rcvr = importdata('rcvr.dat');
eph = importdata('eph.dat');
[eph_row, eph_column] = size(eph);
[rcvr_row, rcvr_column] = size(rcvr);

%  Reorder the matrix
rcvr_reorder = ones(eph_row, rcvr_column);
for i = 1:eph_row
    for g = 1:rcvr_row
        if rcvr(g,2) == eph(i,2)
            rcvr_reorder(i,:) = rcvr(g,:);
        end
    end
end

%%  Definitions and initializations
c = 299792458;
mu = 3.986005*10^(14);
omega_e = 7.2921151467*10^(-5);
Time_elapsed = eph(:,1)-eph(:,4)-rcvr_reorder(:,3)./c; 
Clock_correction_tpara = eph(:,1) - eph(:,3);
Mean_anomaly = eph(:, 12);
e = eph(:, 9);
sqrta = eph(:,10);
motion_correction = eph(:,11);
Odot = eph(:,16);
idot = eph(:,17);
Cus = eph(:,18);
Cuc = eph(:,19);
Cis = eph(:,20);
Cic = eph(:,21);
Crs = eph(:,22);
Crc = eph(:,23);
omega0 = eph(:,13);
Omega0 = eph(:,14);
i0 = eph(:,15);
Ecc_anomaly = ones(1,eph_row);
True_anomaly = ones(1,eph_row);
Perigee = ones(1,eph_row);
Radial_distance = ones(1,eph_row);
Incilination = ones(1,eph_row);
Right_ascension = ones(1,eph_row);
xp = ones(1,eph_row);
yp = ones(1,eph_row);
x = ones(1,eph_row);
y = ones(1,eph_row);
z = ones(1,eph_row);
Satellite_position = ones(eph_row,4);
Clock_correction = ones(eph_row,1);
Max_iteration = 100;

%%  Satellite position
for i = 1:eph_row
    temp(1) = Mean_anomaly(1)+(motion_correction(1)+sqrt(mu./(sqrta(1).^6)))*Time_elapsed(1); 
    for j = 2:Max_iteration
        temp(j) = Mean_anomaly(i)+(motion_correction(i)+sqrt(mu./(sqrta(i).^6)))*Time_elapsed(i)...
            + e(i) * sin(temp(j-1));
        if (abs(temp(j) - temp(j-1))) <= 0.00000001
            break
        end
    end
    %  Eccentric anomaly
    Ecc_anomaly(i) = temp(j);
    %  True anomaly
    True_anomaly(i) = atan(sqrt(1-e(i).^2).*sin(Ecc_anomaly(i))./(cos(Ecc_anomaly(i))-e(i)));
    if cos(Ecc_anomaly(i)) < 0
        if sin(Ecc_anomaly(i)) >= 0
            True_anomaly(i) = True_anomaly(i)+pi;
        else
            True_anomaly(i) = True_anomaly(i)-pi;
        end
    end
    %  Perigee
    Perigee(i) = omega0(i)+True_anomaly(i)+Cuc(i).*cos(2.*(omega0(i)+True_anomaly(i)))+...
        Cus(i).*sin(2.*(omega0(i)+True_anomaly(i)));
    %  Radial distance
    Radial_distance(i) = sqrta(i)^2.*(1-e(i).*cos(Ecc_anomaly(i)))+Crc(i).*cos(2.*(omega0(i)...
        +True_anomaly(i)))+Crs(i).*sin(2.*(omega0(i)+True_anomaly(i)));
    %  Incilination
    Incilination(i) = i0(i)+idot(i).*(Time_elapsed(i))+Cic(i).*cos(2.*(omega0(i)+True_anomaly(i)))...
        +Cis(i).*sin(2.*(omega0(i)+True_anomaly(i)));
    %  Right ascension
    Right_ascension(i) = Omega0(i)+ (Odot(i)-omega_e).*Time_elapsed(i)-omega_e.*eph(i,4);
    xp(i) = Radial_distance(i).*cos(Perigee(i));
    yp(i) = Radial_distance(i).*sin(Perigee(i));
    %  [X,Y,Z] in WGS 84
    x(i) = xp(i)*cos(Right_ascension(i))-yp(i)*cos(Incilination(i))*sin(Right_ascension(i));
    y(i) = xp(i)*sin(Right_ascension(i))+yp(i)*cos(Incilination(i))*cos(Right_ascension(i));
    z(i) = yp(i)*sin(Incilination(i));
    Satellite_position(i,:) = [eph(i,2) x(i) y(i) z(i)];
    %%  Clock correction
    Clock_correction(i,1) = eph(i,5)+eph(i,6)*Clock_correction_tpara(i)+eph(i,7)*...
        Clock_correction_tpara(i).^2;
end
%%  User position
User_position = [-2694685.473 -4293642.366 3857878.924];
R = ones(eph_row,1);
err = 100;
iteration_num = 0;
Clock_bias = 0;
while (err > 10^(-4) && iteration_num < Max_iteration)
    for k = 1:eph_row
        R(k,1) = sqrt((Satellite_position(k,2:4)-User_position)*(Satellite_position(k,2:4)-...
            User_position)') + Clock_bias*c;
    end
    %  Position offset
    Delta_rho = R - rcvr_reorder(:,3) - c.*Clock_correction;
    Hp = ones(eph_row,3);
    for h = 1:eph_row
        Hp(h,:) = (Satellite_position(h,2:4) - User_position)./R(h,1)';
    end
    %  H matrix 
    H = [Hp ones(eph_row,1)];
    %  Iteration function
    DeltaX = inv(H'*H) * H' * Delta_rho;
    TempDeltaX = DeltaX(1:3,:);
    %  User position error term update
    err = max(abs(TempDeltaX));
    %  User position update
    User_position = User_position + DeltaX(1:3,:)';
    %  Clock bias update
    Clock_bias = DeltaX(4,1)/(-c);
    iteration_num = iteration_num + 1;
end
