clear;
clc;
image = imread('einstein.tif');
image = double(image);
[row, column] = size(image);

%% Main
for D = [3, 7]
    main(image, D, row, column);
end

function [] = main(image, D, row, column)
    % Assign
    downsampled = downsample(image, D);
    upsampled = upsampleZero(downsampled, D, row, column);
    upsampledNear = near(upsampled, D, row, column);
    upsampledBili = bili(upsampled, D, row, column);
    % Figures
    figure('Name', ['D = ', num2str(D)]);
    subplot(2, 2, 1);
    imshow(uint8(image));
    title('Original');
    subplot(2, 2, 2);
    imshow(uint8(upsampled));
    title('Upsampled');
    subplot(2, 2, 3);
    imshow(uint8(upsampledNear));
    title('Nearest Neighbor');
    subplot(2, 2, 4);
    imshow(uint8(upsampledBili));
    title('Bilinear');
end

%% Functions
% Nearest Neighbor
function upsampledNear = near(upsampled, D, row, column)
    kernel = nearGen(D);
    upsampledNear = upsampled;
    for i = 1:row
        upsampledNear(i, :) = conv(upsampledNear(i, :), kernel, 'same');
    end
    kernel = kernel';
    for i = 1:column
        upsampledNear(:, i) = conv(upsampledNear(:, i), kernel, 'same');
    end
end

% Bilinear Interpolation
function upsampledBili = bili(upsampled, D, row, column)
    kernel = biliGen(D);
    upsampledBili = upsampled;
    for i = 1:row
        upsampledBili(i, :) = conv(upsampledBili(i, :), kernel, 'same');
    end
    kernel = kernel';
    for i = 1:column
        upsampledBili(:, i) = conv(upsampledBili(:, i), kernel, 'same');
    end
end

% Downsample
function output = downsample(input, D)
    [row, column] = size(input);
    sizeOutput = floor((row-(D+1)/2)/D)+1;
    output = zeros(sizeOutput, sizeOutput);
    i = (D + 1) / 2;
    while (i < row)
        j = (D + 1) / 2;
        while (j < column)
            output(ceil(i/3), ceil(j/3)) = input(i, j);
            j = j + D;
        end
        i = i + D;
    end
end

% Upsample
function output = upsampleZero(input, D, oriRow, oriCol)
    output = zeros(oriRow, oriCol);
    i = (D + 1) / 2;
    while (i < oriRow)
        j = (D + 1) / 2;
        while (j < oriCol)
            output(i, j) = input(ceil(i/3), ceil(j/3));
            j = j + D;
        end
        i = i + D;
    end
end

% Near kernel generator
function output = nearGen(D)
    output = ones(1, D);
end

% Bili kernel generator
function output = biliGen(D)
    output = ones(1, 2 * D - 1);
    for i = 1:2 * D - 1
        output(1, i) = output(1, i) - abs(i - D) / D;
    end
end