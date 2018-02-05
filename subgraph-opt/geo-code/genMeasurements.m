function [yy, S] = genMeasurements(pts, K, l_0, l_1, numSamples, thinRatio)
% Given pts uniform in a [-1,1]^d hypercube, generate an ellipsoid anomaly of ~K
% points centered at the origin. thinRatio gives the ratio of max semi axis
% length to min semi axis length. Non-anomalous are Poisson(l_0), anomalous
% are Poisson(l_1) with numSamples IID copies.

[n, d] = size(pts);

if nargin < 6
    thinRatio = 1;
end

% radii = zeros(1, d);
numLargerDims = 1;
while numLargerDims < d
%     r = (K/n / thinRatio^numLargerDims)^(1/d);
    r = 2/sqrt(pi) * (K/n / thinRatio^numLargerDims * gamma(d/2 + 1))^(1/d);
    if r*thinRatio <= 1
        radii = r * [ones(1, d - numLargerDims) thinRatio*ones(1, numLargerDims)];
        break
    end
    numLargerDims = numLargerDims + 1;
end

dists = sum((pts./radii(ones(1,n),:)).^2, 2);

S_ind = dists <= 1;
S = find(S_ind);
k = length(S);

yy = zeros(n, numSamples);
yy(~S_ind, :) = poissrnd(l_0, n-k, numSamples);
yy(S_ind, :) = poissrnd(l_1, k, numSamples);
