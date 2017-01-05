function [ n_cap, n_inf, n_corr, r ] = getFD(img)
% GETFD      Compute fractal dimensions.
%
% Input:     img     a gray scale image (2D matrix, 0..255)
%
% Output:    FDcap   Capacity Fractal Dimension (box counting)
%            FDinf   Information Fractal Dimension
%            FDcor   Correlation Fractal Dimension
%
%
% Author:    Martin Reuter
% Copyright: Martin Reuter
% Date:      May , 2016
% License:   BSD 3-clause (see LICENSE)
%


% if color image is passed reduce to gray
if (size(img,3) > 1) 
    img = sum(img(:,:,1:3),3)/3.0;
end

% if range is 0..1, scale to 0 ..255:
if (max(img(:)) <= 1)
    img = 255 * img;
end
img = uint8(img);


width = max(size(img));    % largest size of the box
p = log(width)/log(2);   % nbre of generations

% remap the array if the sizes are not all equal,
% or if they are not power of two
% (this slows down the computation!)
if p~=round(p) || any(size(img)~=width)
    % get the highest power of 2
    p = ceil(p);
    % compute the new width
    width = 2^p;
    % generate the new image
    mz = uint8(zeros(width, width));
    % put current image in the very corner of the new image
    mz(1:size(img,1), 1:size(img,2)) = img;
    % reassign to c, c_inf and c_cor
    img = mz;
end

[h, w] = size(img);
lmax = 2^p;

wn = floor(w / lmax) * lmax;
hn = floor(h / lmax) * lmax;


%fprintf(1,'Calculating over h=%d and w=%d\n', hn, wn);
L = zeros(p,1);
NL = zeros(p,1);
InfL = zeros(p,1);
SqrFreqL = zeros(p,1);
count = 0;
hnwn = hn*wn;
S = log(hnwn); % number of pixels, S in Reuters thesis
for boxsize = pow2(1:floor(log2(lmax)))
    Nsum = 0;
    i = 0;
    Inf = 0;
    SqrFreq = 0;
    bs2 = boxsize*boxsize;
    count = count+1;
    for (k = 1:boxsize:hn-boxsize+1)
      for (l = 1:boxsize:wn-boxsize+1)
        ibox = img(k:k+boxsize-1,l:l+boxsize-1);
        maxi = single(max(ibox(:)));
        mini = single(min(ibox(:)));
        %fprintf('max %d  min %d   bs %d  add %d\n',maxi,mini,boxsize,floor((maxi - mini) / boxsize)+1);
        N = floor((maxi - mini) / boxsize)+1;  %number of boxes (round up)
        
        % accumulate the number of boxes
        Nsum = Nsum+ N; 
        
        
        Inf = Inf + bs2 * (log(bs2 / N) - S );
        
        
        SqrFreq = SqrFreq + bs2 *bs2 / N;
        
        
        i = i+1;    % count squares
      end
    end
    %fprintf('Quadrate: %5d ', i);
    %fprintf('N=%7.1f and -ln(N)=%10.6f for boxsize %3d\n', Nsum, -log(single(Nsum)), boxsize);
    L(count) = boxsize;
    NL(count) = Nsum;
    InfL(count) = Inf / hnwn;
    SqrFreqL(count) = SqrFreq / hnwn^2;
end

r = L;
n_cap = -log(NL);
n_inf = InfL;
n_corr = log(SqrFreqL);

% istart = 2;
% iend = length(NL);
% X = [ ones(iend-istart+1,1) log(L(istart:iend))];
% 
% logNL = -log(NL);
% Y = logNL(istart:iend);
% [B,BINT,R,RINT, STATS] = regress(Y,X);
% FDcap = B(2);
% 
% Y = InfL(istart:iend);
% [B,BINT,R,RINT, STATS] = regress(Y,X);
% FDinf = B(2);
% 
% logsqr = log(SqrFreqL);
% Y = logsqr(istart:iend);
% [B,BINT,R,RINT, STATS] = regress(Y,X);
% FDcor = B(2);
% 
% fprintf('\nCapacity Fractal Dimension    : %0.5g\nInformation Fractal Dimension : %0.5g\nCorrelation Fractal Dimension : %0.5g\n',FDcap,FDinf, FDcor);
