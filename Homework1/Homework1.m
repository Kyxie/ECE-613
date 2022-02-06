%% MAE MSE PSNR
clear;
clc;
r = [7 2 7 2;
    2 7 2 7;
    7 2 7 2;
    2 7 2 7];
s = [3 4 3 3;
    4 4 5 3;
    3 5 4 4;
    3 3 4 3];
% s = [7 7 2 2;
%     7 7 2 2;
%     2 2 7 7;
%     2 2 7 7];
MAE = mae(r, s);
MSE = mse(r, s);
PSNR = psnr(r, s);

%% Full-Scale Contrast Stretch
clear;
clc;
r = [2 14 11 14;
    11 2 5 11;
    14 5 5 11;
    16 7 7 16];
B = 4;
minNum = min(r);
minNum = min(minNum);
maxNum = max(r);
maxNum = max(maxNum);
s = round((2^B - 1) * (r - minNum) / (maxNum - minNum));

%% Matrix Convolution
clear;
clc;
coe = 1/6;
h = [0 1 0;
    1 2 1;
    0 1 0];
r = [1 8 6 6;
    6 3 11 8;
    8 8 9 10;
    9 10 10 7];
s = conv2(r, h, 'same');
s = s * coe;
s = round(s);

%% 2D DFT
clear;
clc;
x = [20 -2-2i 0 -2+2i;
    -2-2i 2i 0 2;
    0 0 0 0;
    -2+2i 2 0 -2i];
X = fft2(x);
1 ./ X;

%% Inverse 2D DFT
clear;
clc;
X = [1 -2.5+2.5i 0 -2.5-2.5i;
    -2.5+2.5i 0 0 0;
    0 0 0 0;
    -2.5-2.5i 0 0 0];
x = ifft2(X);

%% Filter
clear;
clc;
% H = [1 -0.2-0.2i 0 -0.2+0.2i;
%     -0.2-0.2i 0.05i 0 0.05;
%     0 0 0 0;
%     -0.2+0.2i 0.05 0 -0.05i];
h = [1 1 1 1;
    1 2 2 1;
    1 2 2 1;
    1 1 1 1];
H = fft2(h);
Inv = 1 ./ H;

sigmawsq = 100;
sigmaxsq = 400;
K = sigmawsq / sigmaxsq;
Wiener = conj(H) ./ (abs(H) .* abs(H) + K);
h = ifft2(H);

%% Median Min & Max Filter
clear;
clc;
r = [5 6 7 8;
    0 6 7 8;
    5 6 15 8;
    5 6 7 8];
mask = [1 1 1;
    1 1 1;
    1 1 1];
% 第二个参数的个数与mask的1的个数有关
% 1就是最小
% middle of 1 and mask中1的个数就是median
% mask中1的个数就是最大
Minr = ordfilt2(r, 1, mask, 'symmetric');
Medr = ordfilt2(r, 5, mask, 'symmetric');
Maxr = ordfilt2(r, 5, mask, 'symmetric');

%% Edge
clear;
clc;
Robertg1 = [1 0;
    0 -1];
Robertg2 = [0 1;
    -1 0];
Prewittg1 = [-1 0 1;
    -1 0 1;
    -1 0 1];
Prewittg2 = [-1 -1 -1;
    0 0 0;
    1 1 1];
Sobelg1 = [-1 0 1;
    -2 0 2;
    -1 0 1];
Sobelg2 = [-1 -2 -1;
    0 0 0;
    1 2 1];
r = [3 3 1 3 3 3 4;
    0 3 3 3 3 3 3;
    3 3 3 2 3 3 12;
    12 3 3 3 3 12 12;
    10 12 2 3 3 12 12;
    12 14 12 12 12 12 11;
    11 12 12 12 10 12 12];
midg1 = -conv2(r, Prewittg1, 'valid');
midg2 = -conv2(r, Prewittg2, 'valid');
s = abs(midg1) + abs(midg2);
threshold = 22;
s(s < threshold) = 0;
s(s > 0) = 1;