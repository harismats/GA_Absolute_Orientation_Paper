%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script registers a 3D point cloud against its transformed version 
% using three methods:
%    1) SVD-based rigid registration
%    2) A Geometric Algebra (Characteristic Multivector, CM) method
%    3) MATLAB's in-built Procrustes function
%
% The script runs a number of test cases (each with a random rotation and
% translation) and computes error metrics (RMSE, Hausdorff distance, angular
% convergence, etc.). It then plots the results and saves the figures to an 
% output folder.
%
% IMPORTANT: This code uses the Clifford Multivector Toolbox. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Workspace and Setup
clear; close all; clc;
clifford_signature(4,0);

% Create output folder for saving figures
output_folder = 'output';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Set random seed for reproducibility
rng(100);

%% Read and Display the Model Point Cloud (example: apple.ply)
ptCloud = pcread('apple.ply');
f = figure;
colormap(f, 'copper');
pcshow(ptCloud, 'BackgroundColor', [0.9 0.9 0.9]);
xlabel('X'); ylabel('Y'); zlabel('Z');
set(gca, 'FontSize', 16);
title('Model Point Cloud');

%% Extract the Coordinates
% Convert the Nx3 point cloud into a 3xN matrix
xyz_coordinates = ptCloud.Location;  % Nx3 matrix
model_Data = xyz_coordinates';         % 3xN matrix
size_model_Data = size(model_Data);

%% Define Number of Test Cases and Preallocate Variables
no_test_cases = 200;  % Number of iterations

% Preallocate variables 
measured_Data       = zeros(3, size_model_Data(2), no_test_cases);
x                   = zeros(1,1,no_test_cases);
y                   = zeros(1,1,no_test_cases);
z                   = zeros(1,1,no_test_cases);
t                   = zeros(1,3,no_test_cases);
Rotation_GT         = zeros(3,3,no_test_cases);
Translation_GT      = zeros(1,3,no_test_cases);
P2                  = zeros(3, size_model_Data(2), no_test_cases);
R_SVD               = zeros(3,3,no_test_cases);
R_CM                = zeros(3,3,no_test_cases);
R_MATLAB            = zeros(3,3,no_test_cases);
t_SVD               = zeros(1,3,no_test_cases);
t_CM                = zeros(1,3,no_test_cases);
t_MATLAB            = zeros(1,3,no_test_cases);

RMSE_SVD            = zeros(1,no_test_cases);
RMSE_CM             = zeros(1,no_test_cases);
RMSE_MATLAB         = zeros(1,no_test_cases);
haus_dist_SVD       = zeros(1,no_test_cases);
haus_dist_CM        = zeros(1,no_test_cases);
haus_dist_MATLAB    = zeros(1,no_test_cases);
mean_ae_SVD         = zeros(1,no_test_cases);
mean_ae_CM          = zeros(1,no_test_cases);
mean_ae_MATLAB      = zeros(1,no_test_cases);
RMSLE_SVD           = zeros(1,no_test_cases);
RMSLE_CM            = zeros(1,no_test_cases);
RMSLE_MATLAB        = zeros(1,no_test_cases);

convergence_error_SVD         = zeros(1,no_test_cases);
convergence_error_CM          = zeros(1,no_test_cases);
convergence_error_MATLAB      = zeros(1,no_test_cases);
convergence_error_SVD_chordal = zeros(1,no_test_cases);
convergence_error_CM_chordal  = zeros(1,no_test_cases);
convergence_error_MATLAB_chordal = zeros(1,no_test_cases);

time_finish_SVD    = zeros(1,no_test_cases);
time_finish_CM     = zeros(1,no_test_cases);
time_finish_MATLAB = zeros(1,no_test_cases);

%% Main Simulation Loop: Apply Random Transformations and Register
for k = 1:no_test_cases
    fprintf('Test Case: %d\n', k);
    
    %% Generate Random Rotation Angles (in degrees) and Translation Vector
    x(:,:,k) = randi([0, 180]);
    y(:,:,k) = randi([0, 180]);
    z(:,:,k) = randi([0, 180]);
    t(:,:,k) = [500 * rand(1,1), 500 * rand(1,1), 500 * rand(1,1)];
    
    %% Transform the Model Data and Obtain Ground Truth Transformation
    [measured_Data(:,:,k), ~, Rotation_GT(:,:,k), Translation_GT(:,:,k)] = ...
        transformation(model_Data, x(:,:,k), y(:,:,k), z(:,:,k), t(:,:,k));
    
    %% Prepare Data for Registration (Convert to Nx3 matrices)
    P1 = model_Data';            % Model: Nx3
    P2(:,:,k) = measured_Data(:,:,k);  % Measured: still 3xN (will be used as is)
    
    %% Registration Using Different Methods
    % SVD-based method (method == 1)
    t_start = tic;
    [R_SVD(:,:,k), t_SVD(:,:,k), RMSE_SVD(k), haus_dist_SVD(k), ...
        mean_ae_SVD(k)] = rigidtform(P1, P2(:,:,k), 1, 1);
    time_finish_SVD(k) = toc(t_start);
    
    % CM method (method == 2)
    t_start = tic;
    [R_CM(:,:,k), t_CM(:,:,k), RMSE_CM(k), haus_dist_CM(k), ...
        mean_ae_CM(k)] = rigidtform(P1, P2(:,:,k), 1, 2);
    time_finish_CM(k) = toc(t_start);
    
    % MATLAB's Procrustes (method == 3)
    t_start = tic;
    [R_MATLAB(:,:,k), ~, RMSE_MATLAB(k), haus_dist_MATLAB(k), ...
        mean_ae_MATLAB(k)] = rigidtform(P1, P2(:,:,k), 1, 3);
    time_finish_MATLAB(k) = toc(t_start);
    
    %% Compute Convergence Errors (Angular and Chordal Distances)
    convergence_error_SVD(k) = 1/sqrt(2) * norm(logm(Rotation_GT(:,:,k) * R_SVD(:,:,k)'), 'fro');
    convergence_error_CM(k)  = 1/sqrt(2) * norm(logm(Rotation_GT(:,:,k) * R_CM(:,:,k)'), 'fro');
    convergence_error_MATLAB(k) = 1/sqrt(2) * norm(logm(Rotation_GT(:,:,k) * R_MATLAB(:,:,k)'), 'fro');
    
    convergence_error_SVD_chordal(k) = norm(Rotation_GT(:,:,k) - R_SVD(:,:,k), 'fro');
    convergence_error_CM_chordal(k)  = norm(Rotation_GT(:,:,k) - R_CM(:,:,k), 'fro');
    convergence_error_MATLAB_chordal(k) = norm(Rotation_GT(:,:,k) - R_MATLAB(:,:,k), 'fro');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting and Saving Figures
% (For demonstration we show plots for the "apple" point cloud case. In your 
% full study you likely have similar variables for helix, tree, etc. Here we
% assume the apple case variables are the ones computed.)

% For convenience, assign these variables:
RMSE_SVD_apple = RMSE_SVD;
RMSE_CM_apple  = RMSE_CM;
RMSE_MATLAB_apple = RMSE_MATLAB;
haus_dist_SVD_apple = haus_dist_SVD;
haus_dist_CM_apple  = haus_dist_CM;
haus_dist_MATLAB_apple = haus_dist_MATLAB;
convergence_error_SVD_apple = convergence_error_SVD;
convergence_error_CM_apple  = convergence_error_CM;
convergence_error_MATLAB_apple = convergence_error_MATLAB;
time_finish_SVD_apple = time_finish_SVD;
time_finish_CM_apple  = time_finish_CM;
time_finish_MATLAB_apple = time_finish_MATLAB;

%% 1. RMSE Comparison: SVD vs. CM (Apple)
fig1 = figure;
plot(RMSE_SVD_apple, 'ko', 'LineWidth',1, 'Marker','o', 'MarkerSize',9);
hold on;
plot(RMSE_CM_apple, 'r+', 'LineWidth',1, 'Marker','+', 'MarkerSize',9);
set(gca, 'FontSize',15, 'XTick', []);
legend(sprintf('SVD (μ: %.3e)', mean(RMSE_SVD_apple)), sprintf('CM (μ: %.3e)', mean(RMSE_CM_apple)), 'FontSize',14, 'Location','best');
grid on;
% Title on two lines:
title({'RMSE: SVD vs. CM', '(Apple)'});
han = axes(fig1,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 18);
ylabel(han, 'RMSE', 'FontSize', 18);
saveas(fig1, fullfile(output_folder, 'RMSE_SVD_vs_CM_Apple.png'));

%% 2. RMSE Comparison: MATLAB vs. CM (Apple)
fig2 = figure;
plot(RMSE_MATLAB_apple, 'ko', 'LineWidth',1, 'Marker','o', 'MarkerSize',9);
hold on;
plot(RMSE_CM_apple, 'r+', 'LineWidth',1, 'Marker','+', 'MarkerSize',9);
set(gca, 'FontSize',15, 'XTick', []);
legend(sprintf('MATLAB (μ: %.3e)', mean(RMSE_MATLAB_apple)), sprintf('CM (μ: %.3e)', mean(RMSE_CM_apple)), 'FontSize',14, 'Location','best');
grid on;
title({'RMSE: MATLAB vs. CM', '(Apple)'});
han = axes(fig2,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 18);
ylabel(han, 'RMSE', 'FontSize', 18);
saveas(fig2, fullfile(output_folder, 'RMSE_MATLAB_vs_CM_Apple.png'));

%% 3. Hausdorff Distance: SVD vs. CM (Apple)
fig3 = figure;
plot(haus_dist_SVD_apple, 'ko', 'LineWidth',1, 'Marker','o', 'MarkerSize',9);
hold on;
plot(haus_dist_CM_apple, 'rx', 'LineWidth',1, 'Marker','x', 'MarkerSize',9);
set(gca, 'FontSize',15, 'XTick', []);
legend(sprintf('SVD (μ: %.3e)', mean(haus_dist_SVD_apple)), sprintf('CM (μ: %.3e)', mean(haus_dist_CM_apple)), 'FontSize',14, 'Location','best');
grid on;
title({'Hausdorff Distance: SVD vs. CM', '(Apple)'});
han = axes(fig3,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 18);
ylabel(han, 'Hausdorff Distance', 'FontSize', 18);
saveas(fig3, fullfile(output_folder, 'Hausdorff_SVD_vs_CM_Apple.png'));

%% 4. Hausdorff Distance: MATLAB vs. CM (Apple)
fig4 = figure;
plot(haus_dist_MATLAB_apple, 'ko', 'LineWidth',1, 'Marker','o', 'MarkerSize',9);
hold on;
plot(haus_dist_CM_apple, 'rx', 'LineWidth',1, 'Marker','x', 'MarkerSize',9);
set(gca, 'FontSize',15, 'XTick', []);
legend(sprintf('MATLAB (μ: %.3e)', mean(haus_dist_MATLAB_apple)), sprintf('CM (μ: %.3e)', mean(haus_dist_CM_apple)), 'FontSize',14, 'Location','best');
grid on;
title({'Hausdorff Distance: MATLAB vs. CM', '(Apple)'});
han = axes(fig4,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 18);
ylabel(han, 'Hausdorff Distance', 'FontSize', 18);
saveas(fig4, fullfile(output_folder, 'Hausdorff_MATLAB_vs_CM_Apple.png'));

%% 5. Angular (Riemannian) Distance Comparison (Apple)
fig5 = figure;
plot(convergence_error_SVD_apple, 'k-', 'LineWidth',1.8);
hold on;
plot(convergence_error_CM_apple, 'r--', 'LineWidth',2.5);
hold on;
plot(convergence_error_MATLAB_apple, ':', 'LineWidth',2.5, 'Color', [0.3020 0.7373 0.8980]);  % Color: #4DBEEE
set(gca, 'FontSize',15, 'XTick', []);
legend(sprintf('SVD (μ: %.1e)', mean(convergence_error_SVD_apple)), ...
       sprintf('CM (μ: %.1e)', mean(convergence_error_CM_apple)), ...
       sprintf('MATLAB (μ: %.1e)', mean(convergence_error_MATLAB_apple)), 'FontSize',14, 'Location','best');
grid on;
title({'Angular Distance (radians)', 'Comparison (Apple)'});
han = axes(fig5,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 18);
ylabel(han, 'Angular Distance (radians)', 'FontSize', 18);
saveas(fig5, fullfile(output_folder, 'AngularDistance_Apple.png'));

%% 6. Running Times Comparison (Apple)
fig6 = figure;
plot(time_finish_SVD_apple, 'k-', 'LineWidth',2.5);
hold on;
plot(time_finish_CM_apple, 'r:', 'LineWidth',2.5);
hold on;
plot(time_finish_MATLAB_apple, '-.', 'LineWidth',2.5, 'Color', [0.6350 0.0780 0.1840]); % Color: #A2142F
set(gca, 'FontSize',14, 'XTick', []);
legend('SVD', 'CM', 'MATLAB', 'FontSize',14, 'Location','best');
grid on;
title({'Running Times (s)', 'Comparison (Apple)'});
han = axes(fig6,'visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Test Case', 'FontSize', 20);
ylabel(han, 'Time (s)', 'FontSize', 20);
saveas(fig6, fullfile(output_folder, 'RunningTimes_Apple.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of Main Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, t, RMSE, haus_dist, mean_ae] = rigidtform(Point_Cloud_1, Point_Cloud_2, dim, method)
% rigidtform Compute the optimal rigid body transformation between two point clouds.
%
%   [R, t, RMSE, haus_dist, mean_ae, RMSLE] = rigidtform(P1, P2, dim, method)
%
%   Input:
%     - Point_Cloud_1: Nx3 matrix representing the model point cloud.
%     - Point_Cloud_2: Nx3 matrix representing the measured point cloud.
%     - dim         : Dimension flag (1 for 3D; only 3D is supported here).
%     - method      : 
%                     1 = SVD-based method,
%                     2 = Characteristic Multivector (CM) method,
%                     3 = MATLAB's in-built Procrustes.
%
%   Output:
%     - R         : 3x3 rotation matrix (returned as the transpose for
%                   consistency with ground truth).
%     - t         : 3x1 translation vector.
%     - RMSE      : Root Mean Squared Error.
%     - haus_dist : Hausdorff distance.
%     - mean_ae   : Mean Absolute Error.
%     - RMSLE     : Root Mean Squared Log Error.
%
%   The function first converts the inputs to single precision, centers the
%   point clouds, and computes the cross-covariance matrix H.
%   Depending on the selected method, it computes the optimal rotation R via:
%     (1) SVD, (2) CM method using F/G matrices, or (3) MATLAB's Procrustes.
%

Point_Cloud_1 = single(Point_Cloud_1);
Point_Cloud_2 = single(Point_Cloud_2);
s = size(Point_Cloud_1);
if nargin < 3
    dim = 1;
else
    if ~isscalar(dim)
        error('rigidtform:NonscalarDimension','Dimension must be a scalar.');
    end
    if dim ~= 1 && dim ~= 2
        error('rigidtform:InvalidDimension','Dimension must be 1 or 2.');
    end
    if dim == 2
        Point_Cloud_1 = Point_Cloud_1';
        Point_Cloud_2 = Point_Cloud_2';
        s = s([2 1]);
    end
end
if s(2) > s(1) && s(2) > 2^15
    warning('rigidtform:LargeDimension','Number of dimensions (%d) is very large. Datasets may be transposed.', s(2));
end
s2 = size(Point_Cloud_2);
if s(2) ~= 3
    Point_Cloud_1 = Point_Cloud_1';
end
if s2(2) ~= 3
    Point_Cloud_2 = Point_Cloud_2';
end

% Compute centroids
c1 = mean(Point_Cloud_1,1);
c2 = mean(Point_Cloud_2,1);

% Center data and form cross-covariance matrix
H = bsxfun(@minus,Point_Cloud_1,c1)' * bsxfun(@minus,Point_Cloud_2,c2);

if method == 1
    % SVD-based method
    [U,~,V] = svd(H);
    R = U * V';
    if s(2) == 3 && det(R) < 0
        V(:,end) = -V(:,end);
        R = U * V';
    end
elseif method == 2
    % Characteristic Multivector (CM) method using F/G matrices
    % Cross-Correlation Matrix (we can use H as calculated above)
    F_matrix = matrix_F(Point_Cloud_1, Point_Cloud_2);
    % Auto-Correlation Matrix from measured data
    G_matrix = matrix_G(Point_Cloud_2);
    [~, ~, R] = CM_FG_matrices(F_matrix, G_matrix);
    R = single(R);
elseif method == 3
    % MATLAB's in-built Procrustes method
    [~, transformed_data_MATLAB, transformation_MATLAB] = procrustes(Point_Cloud_2, Point_Cloud_1);
    R = transformation_MATLAB.T;
else
    error('rigidtform:UnknownMethod','Method must be 1, 2, or 3.');
end

% Calculate optimal translation vector
t = c2-c1*R;

% Registration Result
registration_result = Point_Cloud_1 * R + t;

% Calculate root mean squared distance error for each point
if method == 1 || method ==2 
    % Calculate squared error between actual and prediction
    score = (Point_Cloud_2(:) - (registration_result(:))).^2;
    % Calculate Mean Squared Error
    MSE = mean(score);
    % Calculate RMSE
    RMSE = sqrt(MSE);
else
    % MATLAB
    score = (Point_Cloud_2(:) - (transformed_data_MATLAB(:))).^2;
    % Calculate Mean Squared Error
    MSE = mean(score);
    % Calculate RMSE
    RMSE = sqrt(MSE);
end

% Hausdorff Distance
if method == 1 || method ==2 
    haus_dist = hausdorff_distance(Point_Cloud_2,registration_result);
else
    haus_dist = hausdorff_distance(Point_Cloud_2,transformed_data_MATLAB);
end

% Mean Absolute Error (MAE)
if method == 1 || method ==2 
    mean_ae = mean_absolute_error(Point_Cloud_2,registration_result);
else
    mean_ae = mean_absolute_error(Point_Cloud_2,transformed_data_MATLAB);
end

R = R';  % Transpose for consistency with ground truth convention
end

function [transformed_data, transform_matrix, Rotation_GT, Translation_GT] = transformation(data, x, y, z, translation)
% transformation Applies a 3D rigid transformation (rotation + translation) to the data.
%
%   [transformed_data, transform_matrix, Rotation_GT, Translation_GT] =
%       transformation(data, x, y, z, translation)
%
%   Input:
%     - data        : 3xN matrix (if not, transposed accordingly)
%     - x, y, z     : Rotation angles in degrees about X, Y, and Z axes.
%     - translation : 1x3 translation vector.
%
%   Output:
%     - transformed_data : 3xN matrix of transformed points.
%     - transform_matrix : 4x4 homogeneous transformation matrix.
%     - Rotation_GT      : Ground truth 3x3 rotation matrix.
%     - Translation_GT   : Ground truth 3x1 translation vector.
%
%   Note: The angles are converted to radians for computation.

if size(data,1) ~= 3
    data = data';
end

x = deg2rad(x);
y = deg2rad(y);
z = deg2rad(z);

% Rotation matrices about each axis
Rx = [1      0       0;
      0 cos(x) -sin(x);
      0 sin(x)  cos(x)];
Ry = [cos(y)  0 sin(y);
      0       1  0;
     -sin(y)  0 cos(y)];
Rz = [cos(z) -sin(z) 0;
      sin(z)  cos(z) 0;
      0       0      1];

% Combined rotation matrix
Rotation_GT = Rx * Ry * Rz;

% Form the homogeneous transformation matrix
transform_matrix = [Rotation_GT, translation(:);
                    0 0 0 1];

Translation_GT = transform_matrix(1:3,4);

% Apply transformation
num_points = size(data, 2);
homog_data = [data; ones(1, num_points)];
transformed_homog = transform_matrix * homog_data;
transformed_data = transformed_homog(1:3, :);
end

function distance = hausdorff_distance(A, B)
% HAUSDORFF_DISTANCE Computes the Hausdorff distance between two point sets.
%
%   distance = hausdorff_distance(A, B)
%
%   INPUT:
%       A, B - Two sets of points (each Nx3 for 3D data).
%
%   OUTPUT:
%       distance - The Hausdorff distance between A and B.
%
%   The function computes the distance from every point in A to its closest point in B 
%   (and vice versa) and returns the maximum of these directed distances.
%

    D = pdist2(A, B);
    dist_ab = max(min(D, [], 2)); % Directed distance: from each point in A to the closest in B
    dist_ba = max(min(D, [], 1)); % Directed distance: from each point in B to the closest in A
    distance = max([dist_ab, dist_ba]);
end

function result = mean_absolute_error(X, Y)
% MEAN_ABSOLUTE_ERROR Computes the mean absolute error between two datasets.
%
%   result = mean_absolute_error(X, Y)
%
%   INPUT:
%       X, Y - Two matrices of the same size representing corresponding points.
%
%   OUTPUT:
%       result - The mean absolute error (MAE) computed as the average of the absolute differences.
%

    result = sum(sum(abs(X - Y))) / numel(X);
end

function F = matrix_F(u, v)
% Computes the F matrix used in the Characteristic
% Multivector (CM) method.
%
%   F = matrix_F(u, v)
%
%   Input:
%     u, v - 3xN matrices representing two point clouds.
%
%   Output:
%     F - A 3x3 matrix whose (j,k) entry equals the sum over i=1:N of the
%         product of the scalar components (obtained via the Clifford product)
%         corresponding to e_j (from u) and e_k (from v).
%

    %% Checking size of data
    [rows_u, ~] = size(u);
    [rows_v, ~] = size(v);
    if rows_u ~= 3
        u = u';
    end
    if rows_v ~= 3
        v = v';
    end

    %% Convert to double
    u = double(u);
    v = double(v);

    %% Subtract the mean from each point set as this is what GA method does
    u = u - mean(u, 2);
    v = v - mean(v, 2);

    %% Vector containing e1,e2,e3
    vec = [e1, e2, e3];

    %% Multiplication of u and v (COLUMN) - Convert them to Clifford multivectors
    % This multiplication uses the toolbox operator such that vec * u produces
    % a 1xN vector of Clifford multivectors.
    clif_u = vec * u;
    clif_v = vec * v;

    %% Forming Matrix F --> F_{jk} = sum_{i=1}^{N} (u_i . e_j)(v_i . e_k)
    %% 1st Column calculations
    F11a = bsxfun(@scalar_product, clif_u(1,:), e1);
    F11a_coef = cell2mat(coefficients(F11a));
    F11b = bsxfun(@scalar_product, clif_v(1,:), e1);
    F11b_coef = cell2mat(coefficients(F11b));
    if isempty(F11a_coef) || isempty(F11b_coef)
        F11 = 0;
    else 
        F11 = sum(F11a_coef .* F11b_coef);
    end

    F21a = bsxfun(@scalar_product, clif_u(1,:), e2);
    F21a_coef = cell2mat(coefficients(F21a));
    F21b = bsxfun(@scalar_product, clif_v(1,:), e1);
    F21b_coef = cell2mat(coefficients(F21b));
    if isempty(F21a_coef) || isempty(F21b_coef)
        F21 = 0;
    else 
        F21 = sum(F21a_coef .* F21b_coef);
    end

    F31a = bsxfun(@scalar_product, clif_u(1,:), e3);
    F31a_coef = cell2mat(coefficients(F31a));
    F31b = bsxfun(@scalar_product, clif_v(1,:), e1);
    F31b_coef = cell2mat(coefficients(F31b));
    if isempty(F31a_coef) || isempty(F31b_coef)
        F31 = 0;
    else 
        F31 = sum(F31a_coef .* F31b_coef);
    end

    %% 2nd Column calculations
    F12a = bsxfun(@scalar_product, clif_u(1,:), e1);
    F12a_coef = cell2mat(coefficients(F12a));
    F12b = bsxfun(@scalar_product, clif_v(1,:), e2);
    F12b_coef = cell2mat(coefficients(F12b));
    if isempty(F12a_coef) || isempty(F12b_coef)
        F12 = 0;
    else 
        F12 = sum(F12a_coef .* F12b_coef);
    end

    F22a = bsxfun(@scalar_product, clif_u(1,:), e2);
    F22a_coef = cell2mat(coefficients(F22a));
    F22b = bsxfun(@scalar_product, clif_v(1,:), e2);
    F22b_coef = cell2mat(coefficients(F22b));
    if isempty(F22a_coef) || isempty(F22b_coef)
        F22 = 0;
    else 
        F22 = sum(F22a_coef .* F22b_coef);
    end

    F32a = bsxfun(@scalar_product, clif_u(1,:), e3);
    F32a_coef = cell2mat(coefficients(F32a));
    F32b = bsxfun(@scalar_product, clif_v(1,:), e2);
    F32b_coef = cell2mat(coefficients(F32b));
    if isempty(F32a_coef) || isempty(F32b_coef)
        F32 = 0;
    else 
        F32 = sum(F32a_coef .* F32b_coef);
    end

    %% 3rd Column calculations
    F13a = bsxfun(@scalar_product, clif_u(1,:), e1);
    F13a_coef = cell2mat(coefficients(F13a));
    F13b = bsxfun(@scalar_product, clif_v(1,:), e3);
    F13b_coef = cell2mat(coefficients(F13b));
    if isempty(F13a_coef) || isempty(F13b_coef)
        F13 = 0;
    else 
        F13 = sum(F13a_coef .* F13b_coef);
    end

    F23a = bsxfun(@scalar_product, clif_u(1,:), e2);
    F23a_coef = cell2mat(coefficients(F23a));
    F23b = bsxfun(@scalar_product, clif_v(1,:), e3);
    F23b_coef = cell2mat(coefficients(F23b));
    if isempty(F23a_coef) || isempty(F23b_coef)
        F23 = 0;
    else 
        F23 = sum(F23a_coef .* F23b_coef);
    end

    F33a = bsxfun(@scalar_product, clif_u(1,:), e3);
    F33a_coef = cell2mat(coefficients(F33a));
    F33b = bsxfun(@scalar_product, clif_v(1,:), e3);
    F33b_coef = cell2mat(coefficients(F33b));
    if isempty(F33a_coef) || isempty(F33b_coef)
        F33 = 0;
    else 
        F33 = sum(F33a_coef .* F33b_coef);
    end

    % Assemble the F matrix by columns
    F = [F11, F12, F13; F21, F22, F23; F31, F32, F33];
end

function G = matrix_G(v)
% Computes the G matrix used in the Characteristic
% Multivector (CM) method using explicit loops.
%
%   G = matrix_G(v)
%
%   Input:
%     v - A 3xN matrix representing a point cloud.
%
%   Output:
%     G - A 3x3 matrix where each entry is computed as:
%         G_{jk} = sum_{i=1}^{N} (v_i . e_j)(v_i . e_k)
%

    %% Checking size of data
    [rows_v, ~] = size(v);
    if rows_v ~= 3
        v = v';
    end

    %% Convert to double precision
    v = double(v);

    %% Subtract the mean from v (as in the GA method)
    v = v - mean(v, 2);

    %% Vector containing e1,e2,e3
    vec = [e1, e2, e3];

    %% Multiply v by the basis to convert it to a Clifford multivector
    clif_v = vec * v;

    %% Forming Matrix G --> G_{jk} = sum_{i=1}^{N} (v_i . e_j)(v_i . e_k)
    % 1st Column calculations
    G11a = bsxfun(@scalar_product, clif_v(1,:), e1);
    G11b = bsxfun(@scalar_product, clif_v(1,:), e1);
    G11 = sum(G11a .* G11b);

    G21a = bsxfun(@scalar_product, clif_v(1,:), e2);
    G21b = bsxfun(@scalar_product, clif_v(1,:), e1);
    G21 = sum(G21a .* G21b);

    G31a = bsxfun(@scalar_product, clif_v(1,:), e3);
    G31b = bsxfun(@scalar_product, clif_v(1,:), e1);
    G31 = sum(G31a .* G31b);

    %% 2nd Column calculations
    G12a = bsxfun(@scalar_product, clif_v(1,:), e1);
    G12b = bsxfun(@scalar_product, clif_v(1,:), e2);
    G12 = sum(G12a .* G12b);

    G22a = bsxfun(@scalar_product, clif_v(1,:), e2);
    G22b = bsxfun(@scalar_product, clif_v(1,:), e2);
    G22 = sum(G22a .* G22b);

    G32a = bsxfun(@scalar_product, clif_v(1,:), e3);
    G32b = bsxfun(@scalar_product, clif_v(1,:), e2);
    G32 = sum(G32a .* G32b);

    %% 3rd Column calculations
    G13a = bsxfun(@scalar_product, clif_v(1,:), e1);
    G13b = bsxfun(@scalar_product, clif_v(1,:), e3);
    G13 = sum(G13a .* G13b);

    G23a = bsxfun(@scalar_product, clif_v(1,:), e2);
    G23b = bsxfun(@scalar_product, clif_v(1,:), e3);
    G23 = sum(G23a .* G23b);

    G33a = bsxfun(@scalar_product, clif_v(1,:), e3);
    G33b = bsxfun(@scalar_product, clif_v(1,:), e3);
    G33 = sum(G33a .* G33b);

    % Assemble the G matrix by columns and convert the result to a numeric array
    G = coefficients([G11, G12, G13; G21, G22, G23; G31, G32, G33]);
    G = cell2mat(G);
end

function [New_Rotor, New_Rotor_Reverse, Rnew] = CM_FG_matrices(F_matrix, G_matrix)
% CM_FG_MATRICES Computes the rotor using the Characteristic Multivector (CM) method.
%
%   [New_Rotor, New_Rotor_Reverse, Rnew] = CM_FG_matrices(F_matrix, G_matrix)
%
%   This function recovers the rotation operator (rotor) from the two 3x3 matrices
%   F_matrix and G_matrix. These matrices are derived from the cross-correlation and
%   auto-correlation of the point clouds, respectively. The function uses Clifford
%   Algebra operations (e.g. wedge, unit, reverse, coefficients) to compute three invariants 
%   (Inv1, Inv2, Inv3) based on the reciprocal frame of F_matrix.
%
%   INPUT:
%       F_matrix - 3x3 matrix (must be of type double). Typically derived from the model.
%       G_matrix - 3x3 matrix (must be of type double). Typically derived from the measured data.
%
%   OUTPUT:
%       New_Rotor         - The computed rotor (normalized) as a Clifford multivector.
%       New_Rotor_Reverse - The reverse of the computed rotor.
%       Rnew              - The recovered rotation matrix extracted from the rotor.
%
%   NOTES:
%   - This function preserves all original calculations (e.g. wedge products and cross‐products).
%   - The commented sections (e.g. normalization of F_matrix and G_matrix) are retained for reference.
%   - In rare cases when near-zero values occur, zeros in F_matrix are replaced by a very small number.
%

% Convert input matrices to double
F_matrix = double(F_matrix);
G_matrix = double(G_matrix);

%% Ensure no zero entries (to avoid division problems)
if any(F_matrix == 0)
    F_matrix(F_matrix == 0) = 10e-30;
end

%% Convert F_matrix and G_matrix into Clifford vectors
% Matrix F
vector1_Fclif = F_matrix(1,1) * e1 + F_matrix(2,1) * e2 + F_matrix(3,1) * e3;
vector2_Fclif = F_matrix(1,2) * e1 + F_matrix(2,2) * e2 + F_matrix(3,2) * e3;
vector3_Fclif = F_matrix(1,3) * e1 + F_matrix(2,3) * e2 + F_matrix(3,3) * e3;

% Matrix G
vector1_Gclif = G_matrix(1,1) * e1 + G_matrix(2,1) * e2 + G_matrix(3,1) * e3;
vector2_Gclif = G_matrix(1,2) * e1 + G_matrix(2,2) * e2 + G_matrix(3,2) * e3;
vector3_Gclif = G_matrix(1,3) * e1 + G_matrix(2,3) * e2 + G_matrix(3,3) * e3;

%% Calculate Reciprocals for F (the reciprocal frame)
% These calculations recover the reciprocal frame vectors for F_matrix.
vector1_Fclif_up = - (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * ...
                    wedge(vector2_Fclif, vector3_Fclif) * wedge(e1, e2, e3);
vector2_Fclif_up =   (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * ...
                    wedge(vector1_Fclif, vector3_Fclif) * wedge(e1, e2, e3);
vector3_Fclif_up = - (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * ...
                    wedge(vector1_Fclif, vector2_Fclif) * wedge(e1, e2, e3);


%% Handle potential Inf values in the reciprocals for vector1_Fclif_up
coef_vector1_FC_up = cell2mat(coefficients(vector1_Fclif_up));
if any(coef_vector1_FC_up == Inf) || any(coef_vector1_FC_up == -Inf)
    coef_vector1_FC_up(coef_vector1_FC_up == Inf) = 10e30;
    coef_vector1_FC_up(coef_vector1_FC_up == -Inf) = -10e30;
    coef_vector1_FC_up = [coef_vector1_FC_up, 10e30];
    vector1_Fclif_up = coef_vector1_FC_up(1) * e1 + coef_vector1_FC_up(2) * e2 + ...
                       coef_vector1_FC_up(3) * e3;
end

coef_vector2_FC_up = cell2mat(coefficients(vector2_Fclif_up));
if any(coef_vector2_FC_up == Inf) || any(coef_vector2_FC_up == -Inf)
    coef_vector2_FC_up(coef_vector2_FC_up == Inf) = 10e30;
    coef_vector2_FC_up(coef_vector2_FC_up == -Inf) = -10e30;
    coef_vector2_FC_up = [coef_vector2_FC_up, -10e30];
    vector2_Fclif_up = coef_vector2_FC_up(1) * e1 + coef_vector2_FC_up(2) * e2 + ...
                       coef_vector2_FC_up(3) * e3;
end

coef_vector3_FC_up = cell2mat(coefficients(vector3_Fclif_up));
if any(coef_vector3_FC_up == Inf) || any(coef_vector3_FC_up == -Inf)
    coef_vector3_FC_up(coef_vector3_FC_up == Inf) = 10e30;
    coef_vector3_FC_up(coef_vector3_FC_up == -Inf) = -10e30;
    coef_vector3_FC_up = [coef_vector3_FC_up, -10e30];
    vector3_Fclif_up = coef_vector3_FC_up(1) * e1 + coef_vector3_FC_up(2) * e2 + ...
                       coef_vector3_FC_up(3) * e3;
end

%% Calculate the Invariants based on the reciprocal frame of F and G
Inv1 = vector1_Fclif_up * vector1_Gclif + vector2_Fclif_up * vector2_Gclif + vector3_Fclif_up * vector3_Gclif;
Inv2 = wedge(vector2_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector2_Gclif) + ...
       wedge(vector3_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector3_Gclif) + ...
       wedge(vector3_Fclif_up, vector2_Fclif_up) * wedge(vector2_Gclif, vector3_Gclif);
Inv3 = wedge(vector3_Fclif_up, vector2_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector2_Gclif, vector3_Gclif);
Inv_all = Inv1 + Inv2 + Inv3;

%% Compute the Rotor from the invariants
Rotor_reverse = 1 + Inv_all;              % Rotor reverse (unnormalized)
Rotor_reverse_normalized = unit(Rotor_reverse);  % Normalize the rotor reverse
New_Rotor = reverse(Rotor_reverse_normalized);   % The rotor (normalized)
New_Rotor_Reverse = Rotor_reverse_normalized;      % Its reverse

% Check (should be identity)
R_Rtilde_3D = New_Rotor * New_Rotor_Reverse; %#ok<NASGU>

%% Recover the Rotation Matrix from the Rotor
f1_F = New_Rotor * e1 * New_Rotor_Reverse;
f2_F = New_Rotor * e2 * New_Rotor_Reverse;
f3_F = New_Rotor * e3 * New_Rotor_Reverse;

coefficients_f1_F = cell2mat(coefficients(f1_F));
coefficients_f2_F = cell2mat(coefficients(f2_F));
coefficients_f3_F = cell2mat(coefficients(f3_F));

% Depending on the number of coefficients returned, build the rotation matrix
if length(coefficients_f1_F) && length(coefficients_f2_F) && length(coefficients_f3_F) == 1
    try
        F1 = coefficients_f1_F(1:3);
        F2 = coefficients_f2_F(1:3);
        F3 = coefficients_f3_F(1:3);
    catch
        F1 = coefficients_f1_F(1);
        F2 = coefficients_f2_F(1);
        F3 = coefficients_f3_F(1);
        F = [F1 0 0; 0 F2 0; 0 0 F3];
    end
elseif length(coefficients_f1_F) == 1
    F1 = coefficients_f1_F(1);
    F2 = coefficients_f2_F(1:2);
    F3 = coefficients_f3_F(1:2);
    F = [F1 0 0; 0 F2; 0 F3];
elseif length(coefficients_f1_F) == 2
    F1 = coefficients_f1_F(1:2);
    F2 = coefficients_f2_F(1:3);
    F3 = coefficients_f3_F(1:3);
    F = [F1 0; F2; F3];
elseif length(coefficients_f2_F) == 1
    F2 = coefficients_f2_F(1);
    F1 = coefficients_f1_F(1:2);
    F3 = coefficients_f3_F(1:2);
    F = [F1 0; 0 F2 0; 0 F3];
else
    F1 = coefficients_f1_F(1:3);
    F2 = coefficients_f2_F(1:3);
    F3 = coefficients_f3_F(1:3);
    F = [F1; F2; F3];  % By row
end

Rnew = F;
end
