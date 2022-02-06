clear;
clc;
image = imread('einstein.tif');
image = double(image);

D = 3;
downsampled = downsample(image, D);
upsampled = upsample(downsampled, D);
imshow(uint8(upsampled));

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

function output = upsample(input, D)
    [row, column] = size(input);
    row = D * row;
    column = D * column;
    output = zeros(row, column);
    i = (D + 1) / 2;
    while (i < row)
        j = (D + 1) / 2;
        while (j < column)
            output(i, j) = input(ceil(i/3), ceil(j/3));
            j = j + D;
        end
        i = i + D;
    end
end