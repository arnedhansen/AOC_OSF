%% Visualization of microsaccades for AOC Sternberg

% Visualizations:
%   Boxplots of microsaccade rate

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load gaze data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);
    load('ms_rate_sternberg.mat');
    ms2(subj) = l2msrate;
    ms4(subj) = l4msrate;
    ms6(subj) = l6msrate;
end

%% Plot Euclidean deviations BOXPLOTS
dataMSrate = [ms2', ms4', ms6'];
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'k', 'r'};
hold on;

% Boxplots
boxplot(dataMSrate, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataMSrate(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)), dataMSrate(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Mean Microsaccade Rate [Microsaccades / Second]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
set(gca, 'YLim', [min(dataMSrate(:)) * 0.85 max(dataMSrate(:)) * 1.15]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Mean Microsaccade Rate', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/microsaccades/AOC_gaze_microsaccades_boxplot_sternberg.png');

% Stats
means = mean(dataMSrate, 'omitnan');
Stds = std(dataMSrate, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);
