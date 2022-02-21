clear;
clc;
text = imread('text.tif');
text = double(text);
TextMag = process(text);

%% Assign
side = 21;
sigma = 1;
filter = gaussianFilter(side, sigma);
FilterMag = fftshift(abs(fft2(filter)));
filtered = conv2(text, filter, 'same');
FilteredMag = process(filtered);

% Frequency filter
Filter = fft2(filter);
% Reciprocal frequency filter
reciFilter = 1 ./ Filter;
% Reciprocal spatial filter
recifilter = real(ifft2(reciFilter));
RecifilterMag = fftshift(abs(reciFilter));

% Deblur
deblurred = conv2(filtered, recifilter, 'same');
DeblurredMag = process(deblurred);

%% Figures
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

%% Functions
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