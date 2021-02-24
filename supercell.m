clear variables;
nMax = 20;
mMax = 20;
a = 1;
cutMaxA = 10* a;

theta1 = zeros(nMax, mMax);
theta2 = zeros(nMax, mMax);
length1 = zeros(nMax, mMax);
length2 = zeros(nMax, mMax);

for numIdxN = 1: nMax
    for numIdxM = 1: mMax
        theta1(numIdxN, numIdxM) = 2*atand((2*numIdxN - 1)/(sqrt(3)* (2*numIdxM - 1)));
        theta2(numIdxN, numIdxM) = 2*atand(numIdxN/(sqrt(3)* numIdxM));
        length1(numIdxN, numIdxM) = a/2 * sqrt((2*numIdxN - 1)^2 + 3*(2*numIdxM - 1)^2);
        length2(numIdxN, numIdxM) = a* sqrt(numIdxN^2 + 3* numIdxM^2);
    end
end

% [sortedTheta1, sortedIndex1] = sort(theta1(:));
% [sortedn1, sortedm1] = ind2sub(size(theta1), sortedIndex1);
% [sortedTheta2, sortedIndex2] = sort(theta2(:));
% [sortedn2, sortedm2] = ind2sub(size(theta2), sortedIndex2);
% 
% figure();
% hold on;
% scatter(sortedTheta1, length1(sortedIndex1), 5, 'b');
% scatter(sortedTheta2, length2(sortedIndex2), 5, 'r');
% hold off;

thetaMix = [theta1, theta2];
lengthMix = [length1, length2];
[sortedThetaMix, sortedIndexMix] = sort(thetaMix(:));
[sortednMix, sortedmMix] = ind2sub(size(thetaMix), sortedIndexMix);
sortedmMix(sortedmMix > 20) = sortedmMix(sortedmMix > 20) - 20;
sortedMixType = nan(size(thetaMix));
sortedMixType = 2*(sortedIndexMix(:) > (numel(thetaMix)/2)) + 1*(sortedIndexMix(:) <= (numel(thetaMix)/2));
scatter(sortedThetaMix, lengthMix(sortedIndexMix), 5);

[thetaMixUnique, sortedIndexMixUnique] = uniquetol(thetaMix, 'ByRows', true);
lengthMixUnique = lengthMix(sortedIndexMixUnique);
figure;
scatter(thetaMixUnique, lengthMixUnique, 5);