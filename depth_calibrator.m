clc; clear; close all;

%% User factors
start_set = 1;
end_set = 3;        % Number of scene candidates

min_depth = 1.5;    % Rough target distance
max_depth = 2.1;

sz_img_x = 512;     % Depth image size
sz_img_y = 640;

lb = [0.02 0.03 -0.02 88 -2 88]';       % lower bound of extrinsic camera parameters
Eparam_0 = [0.05 0.06 0 90 0 90]';      % initial guess
ub = [0.08 0.09 0.02 92 2 92]';         % upper bound

%% Optimization of Extrinsic Parameters b/w thermal and radar

% define optimization options
oldoptions = optimoptions(@lsqnonlin,...
    'Algorithm','Levenberg-Marquardt',... % choose LM over Gauss-Newton
    'Diagnostics','on',...                % print diagnostic info
    'Jacobian','off',...                  % jacobian function defined by user
    'MaxFunEvals',1e10000, ...            % let it iterate indefinitely
    'MaxIter',100, ...                    % if initialized correctly, should coverge very quickly (was5e8)
    'Display','iter', ...                 % level of display
    'TolFun', 1e-6, ...
    'FiniteDifferenceStepSize',1e-3);                         

Eparam_estm = lsqnonlin(@(Eparam) exrinsic_loss(Eparam, start_set ,end_set, sz_img_x, sz_img_y, min_depth, max_depth, lb, ub), Eparam_0, lb, ub, oldoptions);

Eparam_estm

%% loss function
function loss = exrinsic_loss(Eparam, start_set ,end_set, sz_img_x, sz_img_y, min_depth, max_depth, lb, ub)

    loss = [];
    for k = start_set:end_set
    %     Eparam = Eparam_0;
        file_name = ['scene_' num2str(k)];
        load(file_name);
        [row_depthimage, col_depthimage, v_depthimage] = find(depth_resized);

        translation = [Eparam(1); Eparam(2); Eparam(3)];      % translation: x, y, z (unit: m)
        rotation = [Eparam(4) Eparam(5) Eparam(6)];           % rotation: yaw (z-axis), pitch (y-axis), roll (x-axis)

        a = deg2rad(rotation(1));
        b = deg2rad(rotation(2));
        c = deg2rad(rotation(3));

        Rz = [cos(a) -sin(a) 0 ;sin(a) cos(a) 0; 0 0 1];
        Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
        Rx = [1 0 0; 0 cos(c) -sin(c); 0 sin(c) cos(c)];
        Rrc = Rx*Ry*Rz;
        Trc = translation;
        E = [Rrc Trc; 0 0 0 1];

        pt_projected = K_therm * E * radarDat4Proj;
        pt_projected = pt_projected(:,:) ./ pt_projected(3,:);

        numPt = size(pt_projected); % number of radar points

        % find depth-points projected in the thermal image with the right distance range
        points = [];
        for i = 1:numPt(2)
            % check if the projected depth is within thermal images
            if (1 < pt_projected(1,i)) && (pt_projected(1,i) < sz_img_x-1) && (1 < pt_projected(2,i)) && (pt_projected(2,i) < sz_img_y-1)
                depth = 1*radarDat4Proj(1,i);
                % check if the radar depth is in the right depth range
                if (depth > min_depth && depth < max_depth)
                    points = [points ; round(pt_projected(1,i)) round(pt_projected(2,i)) depth];
                end
            end
        end

        loss_temp = zeros(6,1);
        loss_temp(1) = (median(row_depthimage) - median(points(:,2,:)));
        loss_temp(2) = (median(col_depthimage) - median(points(:,1,:)));
        loss_temp(3) = min(row_depthimage) - min(points(:,2,:));
        loss_temp(4) = max(row_depthimage) - max(points(:,2,:));
        loss_temp(5) = min(col_depthimage) - min(points(:,1,:));
        loss_temp(6) = max(col_depthimage) - max(points(:,1,:));

        loss = [loss ; loss_temp];
%         img_proj = depth_resized;
%         figure;
%         imshow(img_proj);
%         hold on;
%         plot(points(:,1), points(:,2), 'r*', 'LineWidth', 2, 'MarkerSize', 5);
%         text(points(:,1), points(:,2),num2str(points(:,3)),'Color', 'w')
%         drawnow;
    end
end