%% Finite-size scaling collapse of level-spacing data for the 7NN lattice
%
% Input:
%   7NN/data/paper/level_spacing_7nn_paper.log
%
% Scaling variable:
%   (W - Wc) L^(1/nu)
%
% Level-spacing parameters reported for the 7NN lattice:
%   Wc = 18.2
%   nu = 1.49

clear;
clc;
close all;


%% =====================================================
% Critical parameters
% ======================================================

Wc = 18.2;
nu = 1.49;


%% =====================================================
% Locate the original paper data
% ======================================================

script_directory = fileparts(mfilename('fullpath'));

data_file = fullfile( ...
    script_directory, ...
    '..', ...
    'data', ...
    'paper', ...
    'level_spacing_7nn_paper.log' ...
);

if ~isfile(data_file)
    error('Data file not found:\n%s', data_file);
end

text_data = fileread(data_file);


%% =====================================================
% Identify lattice-size blocks
% ======================================================

blocks = regexp( ...
    text_data, ...
    'Running\s+for\s+L\s*=\s*(\d+)\s*(.*?)(?=Running\s+for\s+L\s*=|$)', ...
    'tokens', ...
    'dotall' ...
);

% Alternative block format
if isempty(blocks)
    blocks = regexp( ...
        text_data, ...
        '==========\s*L\s*=\s*(\d+)\s*==========\s*(.*?)(?==========|$)', ...
        'tokens', ...
        'dotall' ...
    );
end

if isempty(blocks)
    error('No lattice-size blocks were found in the log file.');
end


%% =====================================================
% Extract W and mean level-spacing ratio
% ======================================================

number_of_blocks = length(blocks);

lattice_sizes = zeros(1, number_of_blocks);
all_W = cell(1, number_of_blocks);
all_mean_r = cell(1, number_of_blocks);

% Accepts either:
%   L=12, W=5.00, r=0.5313
% or
%   L = 12, W = 5.00, <r> = 0.5313

data_pattern = [ ...
    'L\s*=\s*\d+\s*,\s*' ...
    'W\s*=\s*([\d.+\-eE]+)\s*,\s*' ...
    '(?:<r>|r)\s*=\s*([\d.+\-eE]+)' ...
];

for block_index = 1:number_of_blocks

    lattice_sizes(block_index) = str2double( ...
        blocks{block_index}{1} ...
    );

    current_block = blocks{block_index}{2};

    extracted_data = regexp( ...
        current_block, ...
        data_pattern, ...
        'tokens' ...
    );

    if isempty(extracted_data)
        warning( ...
            'No level-spacing data found for L = %d.', ...
            lattice_sizes(block_index) ...
        );

        all_W{block_index} = [];
        all_mean_r{block_index} = [];

        continue;
    end

    number_of_points = length(extracted_data);

    W = zeros(number_of_points, 1);
    mean_r = zeros(number_of_points, 1);

    for point_index = 1:number_of_points

        W(point_index) = str2double( ...
            extracted_data{point_index}{1} ...
        );

        mean_r(point_index) = str2double( ...
            extracted_data{point_index}{2} ...
        );
    end

    all_W{block_index} = W;
    all_mean_r{block_index} = mean_r;
end


%% =====================================================
% Plot the finite-size scaling collapse
% ======================================================

figure;
hold on;
box on;

plot_colors = lines(number_of_blocks);
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+'};

for block_index = 1:number_of_blocks

    if isempty(all_W{block_index})
        continue;
    end

    L = lattice_sizes(block_index);
    W = all_W{block_index};
    mean_r = all_mean_r{block_index};

    scaled_disorder = (W - Wc) .* L.^(1.0 / nu);

    marker = markers{ ...
        mod(block_index - 1, length(markers)) + 1 ...
    };

    plot( ...
        scaled_disorder, ...
        mean_r, ...
        'LineWidth', 1.5, ...
        'Marker', marker, ...
        'MarkerSize', 6, ...
        'Color', plot_colors(block_index, :), ...
        'DisplayName', sprintf('L = %d', L) ...
    );
end


%% =====================================================
% Figure formatting
% ======================================================

xlabel( ...
    '$(W-W_c)L^{1/\nu}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14 ...
);

ylabel( ...
    '$\langle r \rangle$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14 ...
);

title( ...
    sprintf( ...
        '7NN level-spacing collapse: W_c = %.1f, \\nu = %.2f', ...
        Wc, nu ...
    ), ...
    'Interpreter', 'latex', ...
    'FontSize', 14 ...
);

legend('show', 'Location', 'best');

set(gca, ...
    'FontWeight', 'bold', ...
    'LineWidth', 1.5, ...
    'FontSize', 12, ...
    'TickDir', 'in', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on' ...
);

xlim([-100, 175]);

hold off;