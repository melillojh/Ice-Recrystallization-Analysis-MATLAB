%% Ice Recrystallization Analysis
%
% MATLAB workflow for quantitative analysis of ice recrystallization
% from microscopy image sequences.
%
% The script:
%   1. Selects the last image as reference.
%   2. Allows the user to select a region of interest (ROI).
%   3. Segments ice crystals using Cellpose.
%   4. Measures area, perimeter, radius, and total ice area.
%   5. Estimates ice area fraction as an approximation of volume fraction.
%   6. Tracks crystals between frames.
%   7. Exports tables and annotated images.
%
% Author: Jorge H. Melillo

clear; clc; close all;

%% =================== User parameters ===================

% Camera and microscope calibration
cameraPixelSize_um = 4.4;      % Camera pixel size [µm]
objectiveMagnification = 10;   % Objective magnification

um_per_pixel = cameraPixelSize_um / objectiveMagnification;

% Cellpose parameters
averageCrystalDiameter_px = 18;
cellposeModel = "cyto";

% Scale bar
scaleBar_um = 10;

%% =================== Select reference image ===================

[referenceFile, inputFolder] = uigetfile({'*.png;*.jpg;*.tif'}, ...
    'Select the LAST image (latest in time)');

if isequal(referenceFile, 0)
    error('No image selected.');
end

outputFolder = uigetdir(inputFolder, ...
    'Select folder to save results');

if isequal(outputFolder, 0)
    error('No output folder selected.');
end

referenceImage = imread(fullfile(inputFolder, referenceFile));

%% =================== Select ROI ===================

spatialRef = imref2d(size(referenceImage), um_per_pixel, um_per_pixel);

fig = figure('Name', 'Select ROI', 'NumberTitle', 'off');
ax = axes('Parent', fig);

imshow(referenceImage, spatialRef, 'Parent', ax, ...
    'InitialMagnification', 100);

axis(ax, 'image');

title(ax, sprintf(['Select the region of interest.\n', ...
    'Press Enter to confirm. Scale: %.3f µm/pixel'], um_per_pixel));

roi = drawrectangle(ax);

set(fig, 'CurrentCharacter', '@');
waitfor(fig, 'CurrentCharacter', char(13));

roiPosition_um = roi.Position;
close(fig);

roiPosition_px = [ ...
    roiPosition_um(1:2) ./ um_per_pixel, ...
    roiPosition_um(3:4) ./ um_per_pixel];

roiPosition_px = round(roiPosition_px);
roiPosition_px = clampROI(roiPosition_px, size(referenceImage));

save(fullfile(outputFolder, 'roi_and_calibration.mat'), ...
    'roiPosition_px', ...
    'cameraPixelSize_um', ...
    'objectiveMagnification', ...
    'um_per_pixel');

%% =================== Initialize Cellpose ===================

cp = cellpose(Model = cellposeModel);

%% =================== Load image sequence ===================

[~, ~, imageExtension] = fileparts(referenceFile);

imageFiles = dir(fullfile(inputFolder, ['*' imageExtension]));

[~, sortingIndex] = sort({imageFiles.name});
imageFiles = imageFiles(sortingIndex);

referenceIndex = find(strcmp({imageFiles.name}, referenceFile), 1);

if isempty(referenceIndex)
    error('Reference image not found in image list.');
end

fprintf('Reference image: %s\n', referenceFile);
fprintf('Reference index: %d of %d\n', referenceIndex, numel(imageFiles));

%% =================== Initialize results ===================

allCrystalResults = table();
summaryRows = [];

previousCentroids = [];
previousIDs = [];
nextID = 1;

%% =================== Main analysis loop ===================

for frameIndex = referenceIndex:-1:1

    currentFile = imageFiles(frameIndex).name;
    currentPath = fullfile(inputFolder, currentFile);

    fprintf('\nProcessing %s, frame %d...\n', currentFile, frameIndex);

    rawImage = imread(currentPath);

    currentROI = clampROI(roiPosition_px, size(rawImage));
    croppedImage = imcrop(rawImage, currentROI);

    nRows = size(croppedImage, 1);
    nCols = size(croppedImage, 2);

    xWorld_um = [0, (nCols - 1) * um_per_pixel];
    yWorld_um = [0, (nRows - 1) * um_per_pixel];

    roiWidth_um = nCols * um_per_pixel;
    roiHeight_um = nRows * um_per_pixel;
    roiArea_um2 = roiWidth_um * roiHeight_um;

    grayImage = im2gray(croppedImage);
    invertedGrayImage = imcomplement(grayImage);

    labels = segmentCells2D(cp, invertedGrayImage, ...
        ImageCellDiameter = averageCrystalDiameter_px);

    stats = regionprops(labels, ...
        'Area', ...
        'Perimeter', ...
        'Centroid', ...
        'EquivDiameter');

    %% =================== No crystals detected ===================

    if isempty(stats)

        summaryRows = [summaryRows; ...
            frameIndex, ...
            0, ...
            roiArea_um2, ...
            0, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN]; %#ok<AGROW>

        saveDetectionImage( ...
            croppedImage, labels, xWorld_um, yWorld_um, ...
            outputFolder, frameIndex, um_per_pixel, scaleBar_um);

        continue;
    end

    %% =================== Crystal properties ===================

    area_px2 = [stats.Area]';
    area_um2 = area_px2 * um_per_pixel^2;

    perimeter_px = [stats.Perimeter]';
    perimeter_um = perimeter_px * um_per_pixel;

    diameter_px = [stats.EquivDiameter]';
    diameter_um = diameter_px * um_per_pixel;

    radius_px = diameter_px / 2;
    radius_um = diameter_um / 2;

    centroids_px = vertcat(stats.Centroid);
    nCrystals = size(centroids_px, 1);

    totalArea_um2 = sum(area_um2);

    meanArea_um2 = mean(area_um2);
    medianArea_um2 = median(area_um2);
    maxArea_um2 = max(area_um2);

    meanRadius_um = mean(radius_um);
    stdRadius_um = std(radius_um);
    stdArea_um2 = std(area_um2);

    % IceAreaFraction:
    % Fraction of the analyzed ROI occupied by detected ice crystals.
    % Assuming approximately constant sample thickness, this can be
    % interpreted as an approximation of the ice volume fraction.
    iceAreaFraction = totalArea_um2 / roiArea_um2;

    %% =================== Simple crystal tracking ===================

    trackedIDs = zeros(nCrystals, 1);
    localIDs = (1:nCrystals)';

    if isempty(previousCentroids)

        trackedIDs = (nextID : nextID + nCrystals - 1)';
        nextID = nextID + nCrystals;

    else

        searchRadius_px = 2 * mean(radius_px);
        nPrevious = size(previousCentroids, 1);

        for j = 1:nPrevious

            distances = vecnorm(centroids_px - previousCentroids(j, :), 2, 2);

            matchingIndices = find(distances <= searchRadius_px & trackedIDs == 0);

            trackedIDs(matchingIndices) = previousIDs(j);
        end

        unmatched = find(trackedIDs == 0);

        if ~isempty(unmatched)
            trackedIDs(unmatched) = ...
                (nextID : nextID + numel(unmatched) - 1)';
            nextID = nextID + numel(unmatched);
        end
    end

    previousCentroids = centroids_px;
    previousIDs = trackedIDs;

    %% =================== Save per-crystal results ===================

    fileColumn = repmat(string(currentFile), nCrystals, 1);
    frameColumn = repmat(frameIndex, nCrystals, 1);

    frameResults = table( ...
        fileColumn, ...
        frameColumn, ...
        trackedIDs, ...
        localIDs, ...
        centroids_px(:,1), ...
        centroids_px(:,2), ...
        area_px2, ...
        area_um2, ...
        perimeter_px, ...
        perimeter_um, ...
        radius_px, ...
        radius_um, ...
        'VariableNames', { ...
        'FileName', ...
        'FrameIndex', ...
        'TrackedID', ...
        'LocalID', ...
        'CentroidX_px', ...
        'CentroidY_px', ...
        'Area_px2', ...
        'Area_um2', ...
        'Perimeter_px', ...
        'Perimeter_um', ...
        'Radius_px', ...
        'Radius_um'});

    allCrystalResults = [allCrystalResults; frameResults]; %#ok<AGROW>

    summaryRows = [summaryRows; ...
        frameIndex, ...
        totalArea_um2, ...
        roiArea_um2, ...
        iceAreaFraction, ...
        meanArea_um2, ...
        medianArea_um2, ...
        maxArea_um2, ...
        meanRadius_um, ...
        stdRadius_um, ...
        stdArea_um2]; %#ok<AGROW>

    writetable(frameResults, ...
        fullfile(outputFolder, sprintf('results_frame_%d.csv', frameIndex)));

    %% =================== Save detection image ===================

    saveDetectionImage( ...
        croppedImage, labels, xWorld_um, yWorld_um, ...
        outputFolder, frameIndex, um_per_pixel, scaleBar_um);

end

%% =================== Save summary tables ===================

summaryRows = sortrows(summaryRows, 1);

summaryTable = array2table(summaryRows, ...
    'VariableNames', { ...
    'FrameIndex', ...
    'TotalArea_um2', ...
    'ROI_Area_um2', ...
    'IceAreaFraction', ...
    'MeanArea_um2', ...
    'MedianArea_um2', ...
    'MaxArea_um2', ...
    'MeanRadius_um', ...
    'StdRadius_um', ...
    'StdArea_um2'});

writetable(summaryTable, ...
    fullfile(outputFolder, 'summary_per_frame.csv'));

allCrystalResults = sortrows(allCrystalResults, ...
    {'FrameIndex', 'TrackedID'});

writetable(allCrystalResults, ...
    fullfile(outputFolder, 'all_crystals_all_frames.csv'));

fprintf('\nProcessing finished successfully.\n');

%% =================== Local functions ===================

function roi_px = clampROI(roi_px, imageSize)

    roi_px(1) = max(1, roi_px(1));
    roi_px(2) = max(1, roi_px(2));

    roi_px(3) = max(1, min(roi_px(3), imageSize(2) - roi_px(1) + 1));
    roi_px(4) = max(1, min(roi_px(4), imageSize(1) - roi_px(2) + 1));
end

function saveDetectionImage( ...
    croppedImage, labels, xWorld_um, yWorld_um, ...
    outputFolder, frameIndex, um_per_pixel, scaleBar_um)

    overlayImage = labeloverlay(croppedImage, labels);

    fig = figure('Visible', 'off');
    ax = axes(fig);

    imshow(overlayImage, ...
        'XData', xWorld_um, ...
        'YData', yWorld_um, ...
        'Parent', ax);

    axis(ax, 'image');
    hold(ax, 'on');

    title(ax, sprintf('Ice crystal detection - frame %d', frameIndex), ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    drawScaleBar_um(ax, size(overlayImage), um_per_pixel, scaleBar_um);

    exportgraphics(fig, ...
        fullfile(outputFolder, sprintf('crystals_frame_%d.png', frameIndex)), ...
        'Resolution', 300);

    close(fig);
end

function drawScaleBar_um(ax, imageSize, um_per_pixel, scaleBar_um)

    rows = imageSize(1);
    cols = imageSize(2);

    imageWidth_um = cols * um_per_pixel;
    imageHeight_um = rows * um_per_pixel;

    inset_um = 0.08 * min(imageWidth_um, imageHeight_um);
    barThickness_um = 0.01 * min(imageWidth_um, imageHeight_um);

    x0_um = imageWidth_um - inset_um - scaleBar_um;
    y0_um = imageHeight_um - inset_um - barThickness_um;

    rectangle(ax, ...
        'Position', [x0_um, y0_um, scaleBar_um, barThickness_um], ...
        'FaceColor', 'w', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1);

    text(ax, ...
        x0_um + scaleBar_um/2, ...
        y0_um - 0.02 * imageHeight_um, ...
        sprintf('%g \\mum', scaleBar_um), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold', ...
        'Color', 'w', ...
        'BackgroundColor', 'k', ...
        'Margin', 1, ...
        'Clipping', 'on');
end