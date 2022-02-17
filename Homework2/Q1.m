clear;
clc;
slope = imread('slope.tif');
bit = 8;

%% Assignment
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

%% Figures
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