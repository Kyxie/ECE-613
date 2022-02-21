clear;
clc;
bridge = imread('bridge.tif');
bridge = double(bridge);
[side, ~] = size(bridge);

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