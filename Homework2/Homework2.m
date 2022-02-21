%%
% * Name: Kunyang Xie
% * Student ID: 20958936
% * Email: k47xie@uwaterloo.ca

%% Question 1

% Read the image
clear;
clc;
slope = imread('slope.tif');

% Process
bitPlane1 = double(bitget(slope, 8));
bitPlane2 = double(bitget(slope, 7));
bitPlane3 = double(bitget(slope, 6));
bitPlane4 = double(bitget(slope, 5));
bitPlane5 = double(bitget(slope, 4));
bitPlane6 = double(bitget(slope, 3));
bitPlane7 = double(bitget(slope, 2));
bitPlane8 = double(bitget(slope, 1));

version1 = bitshift(bitPlane1, 7);
version2 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6);
version3 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5);
version4 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5) + bitshift(bitPlane4, 4);
version5 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5) + bitshift(bitPlane4, 4) + bitshift(bitPlane5, 3);
version6 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5) + bitshift(bitPlane4, 4) + bitshift(bitPlane5, 3) + bitshift(bitPlane6, 2);
version7 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5) + bitshift(bitPlane4, 4) + bitshift(bitPlane5, 3) + bitshift(bitPlane6, 2) + bitshift(bitPlane7, 1);
version8 = bitshift(bitPlane1, 7) + bitshift(bitPlane2, 6) + bitshift(bitPlane3, 5) + bitshift(bitPlane4, 4) + bitshift(bitPlane5, 3) + bitshift(bitPlane6, 2) + bitshift(bitPlane7, 1) + bitPlane8;

% Plot the figures
figure('Name', 'Upper 4 bitplanes');
subplot(2, 2, 1);
imshow(bitPlane1);
title('1st Bitplane');
subplot(2, 2, 2);
imshow(bitPlane2);
title('2nd Bitplane');
subplot(2, 2, 3);
imshow(bitPlane3);
title('3rd Bitplane');
subplot(2, 2, 4);
imshow(bitPlane4);
title('4th Bitplane');

figure('Name', 'Lower 4 bitplanes');
subplot(2, 2, 1);
imshow(bitPlane5);
title('5th Bitplane');
subplot(2, 2, 2);
imshow(bitPlane6);
title('6th Bitplane');
subplot(2, 2, 3);
imshow(bitPlane7);
title('7th Bitplane');
subplot(2, 2, 4);
imshow(bitPlane8);
title('8th Bitplane');

figure('Name', 'First 4 images');
subplot(2, 2, 1);
imshow(uint8(version1));
title('First 1 image');
subplot(2, 2, 2);
imshow(uint8(version2));
title('First 2 image');
subplot(2, 2, 3);
imshow(uint8(version3));
title('First 3 image');
subplot(2, 2, 4);
imshow(uint8(version4));
title('First 4 image');

figure('Name', 'Remaining 4 images');
subplot(2, 2, 1);
imshow(uint8(version5));
title('First 5 image');
subplot(2, 2, 2);
imshow(uint8(version6));
title('First 6 image');
subplot(2, 2, 3);
imshow(uint8(version7));
title('First 7 image');
subplot(2, 2, 4);
imshow(uint8(version8));
title('First 8 image');

%%
% 
% * As for the image of bit-plane, it changes more and more frequently 
% from front to end, this is because the more significant bits change 
% slower than less significant bits.
% * As for the images of reconstructed images, the more bits the 
% reconstructed image contains, the richer the information and the closer 
% it is to the original image.
% 

%% Question 2

% Read the image
clear;
clc;
book = imread('books.tif');
bit = 8;

% Process
Gamma = 0.5;
gamma = gammaCorrection(book, bit, Gamma);
contrast =  fullScaleContrast(book, bit);
histogram = histogramEqual(book, bit);

% Plot the figures
figure('Name', 'Images');
subplot(2,2,1);
imshow(book);
title('Origional Image');
subplot(2,2,2);
imshow(gamma);
title('Gamma-mapped Image');
subplot(2,2,3);
imshow(contrast);
title('Full-scale Contrast Stretched Image');
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
title('Full-scale Contrast Stretched Image');
subplot(2,2,4);
imhist(histogram);
title('Histogram Equalized Image');

%%
% 
% * The original image is little bit dark, after gamma mapping, it becomes 
% brighter, then after the Full-scale contrast stretch, its contrast 
% enhanced, finally after the histogram equalization, its contrast 
% enhanced again, the white goes whiter and the black goes blacker.
% * As for the histograms, after the gamma map the data shifted to the 
% middle of the chart, while the width does not change, then after the 
% full-scale contrast stretch, the figure is wider than the second image, 
% however there is still blank between 20 ¨C 50, then at the last image, 
% the data is stretched even further. 
% 

%% Question 3

% Read the image
clear;
clc;
bridge = imread('bridge.tif');
bridge = double(bridge);
[side, ~] = size(bridge);

% Process
% fft
Bridge = fft2(bridge);
Bridge = fftshift(Bridge);

% Frequencies
lowFreq = 1/8;
bandFreq = [1/8, 1/2];
highFreq = 1/2;

% Generate filters
lowFilter = lowPass(side, lowFreq);
bandFilter = bandPass(side, bandFreq);
highFilter = highPass(side, highFreq);

% Filtering
lowFiltered = Bridge .* lowFilter;
bandFiltered = Bridge .* bandFilter;
highFiltered = Bridge .* highFilter;

% Shift
lowFilteredShifted = fftshift(lowFiltered);
bandFilteredShifted = fftshift(bandFiltered);
highFilteredShifted = fftshift(highFiltered);

% ifft
lowFilteredSpatial = real(ifft2(lowFilteredShifted));
bandFilteredSpatial = real(ifft2(bandFilteredShifted)) + 128;
highFilteredSpatial = real(ifft2(highFilteredShifted)) + 128;

% Plot the figures
% Spatial Domain Figure
figure('Name', 'Spatial Domain');
subplot(2, 2, 1);
imshow(uint8(bridge));
title('Original Image');
subplot(2, 2, 2);
imshow(uint8(lowFilteredSpatial));
title('Low-pass Filtered Image');
subplot(2, 2, 3);
imshow(uint8(bandFilteredSpatial));
title('Band-pass Filtered Image');
subplot(2, 2, 4);
imshow(uint8(highFilteredSpatial));
title('High-pass Filtered Image');

% Frequency Domain Figure
origin = log(1 + Bridge);
maxOrigin = max(max(origin));
origin = origin / maxOrigin;

lowFiltered = log(1 + lowFiltered);
maxLowFiltered = max(max(lowFiltered));
lowFiltered = lowFiltered / maxLowFiltered;

bandFiltered = log(1 + bandFiltered);
maxBandFiltered = max(max(bandFiltered));
bandFiltered = bandFiltered / maxBandFiltered;

highFiltered = log(1 + highFiltered);
maxHighFiltered = max(max(highFiltered));
highFiltered = highFiltered / maxHighFiltered;

figure('Name', 'Frequency Domain');
subplot(2, 2, 1);
imshow(real(origin));
title('Original Image');
subplot(2, 2, 2);
imshow(real(lowFiltered));
title('Low-pass Filtered Image');
subplot(2, 2, 3);
imshow(real(bandFiltered));
title('Band-pass Filtered Image');
subplot(2, 2, 4);
imshow(real(highFiltered));
title('High-pass Filtered Image');

%%
% 
% * The center of the image is brightest, which means it contains more
% information in the center, in this case, conpare to the band-pass and
% high-pass filter, the low-pass filtered image is more similar to the
% original image.
% 

%% Question 4

% Read the image
clear;
clc;
text = imread('text.tif');
text = double(text);
TextMag = process(text);

% Process
side = 21;
sigma = 1;
filter = gaussianFilter(side, sigma);
FilterMag = process(filter);
filtered = conv2(text, filter, 'same');
FilteredMag = process(filtered);

% Frequency filter
Filter = fft2(filter);
% Reciprocal frequency filter
reciFilter = 1 ./ Filter;
% Reciprocal spatial filter
recifilter = real(ifft2(reciFilter));
RecifilterMag = process(recifilter);

% Deblur
deblurred = conv2(filtered, recifilter, 'same');
DeblurredMag = process(deblurred);

% Plot the figures
figure('Name', 'Images');
subplot(2, 3, 1);
imshow(uint8(text));
title('Original Image');
subplot(2, 3, 2);
imshow(uint8(filtered));
title('Blurred Image');
subplot(2, 3, 3);
imshow(uint8(deblurred));
title('Deblurred Image');
subplot(2, 3, 4);
imshow(TextMag);
title('Original Magnitude Spectra');
subplot(2, 3, 5);
imshow(FilteredMag);
title('Filtered Magnitude Spectra');
subplot(2, 3, 6);
imshow(DeblurredMag);
title('Deblurred Magnitude Spectra');

figure('Name', 'Filters');
subplot(2, 2, 1);
mesh(filter);
title('Spatial Gaussian Filter');
subplot(2, 2, 2);
mesh(FilterMag);
title('Frequency Gaussian Filter Magnitude');
subplot(2, 2, 3);
mesh(recifilter);
title('Spatial Inverse Gaussian Filter');
subplot(2, 2, 4);
mesh(RecifilterMag);
title('Frequency Inverse Gaussian Filter Magnitude');

%%
% 
% * We can see that the filtered magnitude spectra is the original spectra 
% times the gaussian filter in frequency, while the deblur filter is 
% opposite to the gaussian filter in frequency.
% 

%% Question 5

% Read the image
clear;
clc;
image = imread('einstein.tif');
image = double(image);
[row, column] = size(image);

% Process
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

%%
% 
% * The performance of bilinear method is better than nearest neighbor 
% method, the larger D is, the more information is lost, and therefore the 
% more blurred the interpolated image is.
% 

%% Question 2 Functions

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

%% Question 3 Functions
function output = lowPass(side, frequency)
    frequency = frequency * side / 2;
    output = zeros(side, side);
    center = side / 2;
    
    for i = 1:side
        for j = 1:side
            if ((i - center)^2 + (j - center)^2 <= frequency^2)
                output(i, j) = 1;
            end
        end
    end
end

function output = bandPass(side, frequency)
    frequency = frequency * side / 2;
    output = zeros(side, side);
    center = side / 2;
    
    for i = 1:side
        for j = 1:side
            if (((i - center)^2 + (j - center)^2 <= frequency(1,2)^2) && ((i - center)^2 + (j - center)^2 >= frequency(1,1)^2))
                output(i, j) = 1;
            end
        end
    end
end

function output = highPass(side, frequency)
    frequency = frequency * side / 2;
    output = zeros(side, side);
    center = side / 2;
    
    for i = 1:side
        for j = 1:side
            if ((i - center)^2 + (j - center)^2 >= frequency^2)
                output(i, j) = 1;
            end
        end
    end
end

%% Question 4 Functions
function output = gaussianFilter(side, sigma)
    output = ones(side, side);
    center = ceil(side / 2);
    for i = 1:side
        for j = 1:side
            x = i - center;
            y = -(j - center);
            output(i, j) = 1 / (2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2));
        end
    end
    output = output / sum(sum(output));
end

function output = process(image)
    Image = fft2(image);
    ImageMag = abs(Image);
    ImageMag = log(1 + ImageMag);
    ImageMagNormal = ImageMag / max(max(ImageMag));
    output = fftshift(ImageMagNormal);
end

%% Question 5 Functions

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