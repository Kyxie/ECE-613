clear;
clc;
book = imread('books.tif');
bit = 8;

Gamma = 0.5;
gamma = gammaCorrection(book, bit, Gamma);

contrast =  fullScaleContrast(book, bit);

histogram = histogramEqual(book, bit);

figure('Name', 'Images');
subplot(2,2,1);
imshow(book);
title('Origional Image');
subplot(2,2,2);
imshow(gamma);
title('Gamma-mapped Image');
subplot(2,2,3);
imshow(contrast);
title('Full-scale Contast Stretched Image');
subplot(2,2,4);
imshow(histogram);
title('Histogram Equalized Image');


figure('Name', 'Histograms');
subplot(2,2,1);
imhist(book);
title('Origional Image');
subplot(2,2,2);
imhist(gamma);
title('Gamma-mapped Image');
subplot(2,2,3);
imhist(contrast);
title('Full-scale Contast Stretched Image');
subplot(2,2,4);
imhist(histogram);
title('Histogram Equalized Image');

% Gamma Correction
function output = gammaCorrection(input, bit, Gamma)
    input = double(input);
    output = input / 2^bit;
    output = imadjust(output,[],[],Gamma);
    output = uint8(output * 2^bit);
end

% Full Scale Contrast
function output = fullScaleContrast(input, bit)
    input = double(input);
    reshapeInput = input(:);
    minNum = min(reshapeInput);
    maxNum = max(reshapeInput);
    output = round((2 ^ bit - 1) * (input - minNum) / (maxNum - minNum));
    output = uint8(output);
end

% Histogram Equalization
function output = histogramEqual(input, bit)
    input = double(input);
    [row, column] = size(input);
    reshapeInput = input(:);
    reshapeInput = sort(reshapeInput);
    % Difference is 0 means same number
    difference = diff([reshapeInput; max(reshapeInput) + 1]);

    % count is pdf
    count = diff(find([1;difference]));

    % sumCount is cdf
    sumCount = count;
    component = 0;
    for i = 1:length(count)
        component = component + count(i, 1);
        sumCount(i, 1) = component;
    end

    % Replace form
    form = [reshapeInput(find(difference)) sumCount];

    %  Replace
    [rowForm, ~] = size(form);
    newReshapeInput = input(:);
    equalized = input(:);
    for i = 1:length(newReshapeInput)
        for j = 1:rowForm
            if (newReshapeInput(i) == form(j,1))
                equalized(i) = form(j,2);
            end
        end
    end
    equalized = reshape(equalized,[row, column]);

    % Full-scale contrast stretch
    output = fullScaleContrast(equalized, bit);
    output = uint8(output);
end