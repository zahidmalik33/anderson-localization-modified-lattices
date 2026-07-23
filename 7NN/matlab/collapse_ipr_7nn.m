%% Finite-size scaling collapse of the typical IPR for the 7NN lattice
%
% Input:
%   7NN/data/paper/ipr_7nn_paper.log
%
% Scaling form:
%   IPR_typ(W,L) = L^(-tau) F[(W-Wc)L^(1/nu)]
%
% Parameters used in the manuscript:
%   Wc  = 18.4
%   nu  = 1.49
%   tau = 1.18

clear;
clc;
close all;


%% =====================================================
% Critical parameters
% ======================================================

Wc = 18.4;
nu = 1.49;
tau = 1.18;


%% =====================================================
% Locate the original paper data
% ======================================================

script_directory = fileparts(mfilename('fullpath'));

data_file = fullfile( ...
    script_directory, ...
    '..', ...
    'data', ...
    'paper', ...
    'ipr_7nn_paper.log' ...
);

if ~isfile(data_file)
    error('Data file not found:\n%s', data_file);
end

text_data = fileread(data_file);


%% =====================================================
% Separate the data into lattice-size blocks
% ======================================================

blocks = regexp( ...
    text_data, ...
    '==========\s*L\s*=\s*(\d+)\s*==========\s*(.*?)(?==========|$)', ...
    'tokens', ...
    'dotall' ...
);

if isempty(blocks)
    error('No lattice-size blocks were found in the log file.');
end


%% =====================================================
% Extract W and typical IPR
% ======================================================

number_of_blocks = length(blocks);

lattice_sizes = zeros(1, number_of_blocks);
all_W = cell(1, number_of_blocks);
all_typical_ipr = cell(1, number_of_blocks);

data_pattern = [ ...
    'W\s*=\s*([\d.+\-eE]+),\s*' ...
    'Mean\s+ln\(IPR\)\s*=\s*([\d.+\-eE]+),\s*' ...
    'IPR_typ\s*=\s*([\d.+\-eE]+)' ...
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
            'No IPR data found for L = %d.', ...
            lattice_sizes(block_index) ...
        );

        all_W{block_index} = [];
        all_typical_ipr{block_index} = [];

        continue;
    end

    number_of_points = length(extracted_data);

    W = zeros(number_of_points, 1);
    typical_ipr = zeros(number_of_points, 1);

    for point_index = 1:number_of_points

        W(point_index) = str2double( ...
            extracted_data{point_index}{1} ...
        );

        typical_ipr(point_index) = str2double( ...
            extracted_data{point_index}{3} ...
        );
    end

    all_W{block_index} = W;
    all_typical_ipr{block_index} = typical_ipr;
end


%% =====================================================
% Construct and plot the scaled variables
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
    typical_ipr = all_typical_ipr{block_index};

    scaled_disorder = (W - Wc) .* L.^(1.0 / nu);
    scaled_ipr = typical_ipr .* L.^tau;

    marker = markers{ ...
        mod(block_index - 1, length(markers)) + 1 ...
    };

    plot( ...
        scaled_disorder, ...
        scaled_ipr, ...
        'LineStyle', '-', ...
        'LineWidth', 1.5, ...
        'Marker', marker, ...
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
    '$\mathrm{IPR}_{\mathrm{typ}}L^{\tau}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14 ...
);

title( ...
    sprintf( ...
        '7NN IPR collapse: W_c = %.1f, \\nu = %.2f, \\tau = %.2f', ...
        Wc, nu, tau ...
    ), ...
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

xlim([-100, 100]);
ylim([0, 2]);

hold off;