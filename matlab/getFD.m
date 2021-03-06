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

% get height and width
[h, w] = size(img);
lmax = 2^p;

% initialize the output variables
rs = zeros(p,1);
Nr = zeros(p,1);
Infr = zeros(p,1);
SqrFreqr = zeros(p,1);
% % ALTERNATIVE
%Inf_alternative = zeros(p, 1);
%SqrFreqr_alternative = zeros(p, 1);

% initialize a counter to iterate per each of the scales
current_r_id = 0;

% precompute some values to accelerate the computation
hw = h*w; % size of the image
loghw = log(hw); % log of the size of the image

% for each grid size
for r = pow2(1:floor(log2(lmax)))
    
    % initialize accumulators
    Nsum = 0;
    Inf_ = 0;
    SqrFreq = 0;
    % precompute r^2
    r2 = r*r;
    % initialize box id
    box_id = 0;
    % increment current r id
    current_r_id = current_r_id+1;
    
    % % ALTERNATIVE
    %Inf_sub_alternative = 0;
    %SqrFreqr_sub_alternative = 0;
    
    % for each of the cells in the grid
    for k = 1 : r : h-r+1
        for l = 1 : r : w-r+1

            % get current grid
            ibox = img(k:k+r-1,l:l+r-1);
            % retrieve max and minimum values
            maxi = single(max(ibox(:)));
            mini = single(min(ibox(:)));

            % compute N
            N = floor((maxi - mini) / r)+1;  %number of boxes (round up)

            % sum for the capacity measure
            Nsum = Nsum + N; 
            % sum for the information measure
            Inf_ = Inf_ + r2 * (log(r2 / N) - loghw );
            % sum for the correlation measure
            SqrFreq = SqrFreq + r2 *r2 / N;
            
            % % ALTERNATIVE
            %P = N / (hw);
            %Inf_sub_alternative = Inf_sub_alternative - P * log(P);
            %SqrFreqr_sub_alternative = SqrFreqr_sub_alternative + P^2;
            
            % increment current box id
            box_id = box_id+1;    
            
        end
    end

    
    rs(current_r_id) = r;
    Nr(current_r_id) = Nsum;
    Infr(current_r_id) = Inf_ / hw;
    SqrFreqr(current_r_id) = SqrFreq / hw^2;
    
    % % ALTERNATIVE
    %Inf_alternative(current_r_id) = - Inf_sub_alternative;
    %SqrFreqr_alternative(current_r_id) = SqrFreqr_sub_alternative;
    
    
end

r = rs;
n_cap = Nr;
n_inf = Infr;
n_corr = SqrFreqr;

% r = rs;
% n_cap = Nr;
% n_inf = Inf_alternative;
% n_corr = SqrFreqr_alternative;
