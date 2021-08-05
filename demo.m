clc; close all; clear;
Start = tic;

%% Inputs
numSig = 5;
t = linspace(0,1,1000); % Time vector
nbins = 5;

methods = {'Euclidean', "CityBlock", "Minkowski", "Chebyshev", "Sorensen", "Gower", ...
           "Soergel", "Kulczynski", "Canberra", "Lorentzian", 'Intersection', 'WaveHedges', ...
           'Czekanowski', 'Motyka', 'Kulczynski', 'Ruzicka', 'Tanimoto', "InnerProduct", "HarmonicMean", ...
           "Cosine", "Kumar-Hassebrook", "Jaccard", "Dice", 'Fidelity', 'Bhattacharya', 'Hellinger', ...
           'Matusita', 'Squared-chord', 'SquaredEuclidean', 'Pearson', 'Neyman', 'Squared', 'ProbabilisticSymmetrix', ...
           'Divergence', 'Clark', 'AdditiveSymmetric', 'Kullback-Leibler', 'Jeffreys', 'KDivergence', 'Topsoe', ...
           'Jensen-Shannon', 'JensenDifference', 'Taneja', 'Kumar-Johnson', 'Average'};
       
%% Signals 1 and 2
phi = rand(1,numSig)*2*pi; % Random phase
A = rand(1,numSig); % Random ampltiude

x = zeros(size(t));
for i = 1:numSig
    x = x + A(i)*sin(rand(1)*randi(500).*t+phi(i));
end

phi = rand(1,numSig)*2*pi; % Random phase
A = rand(1,numSig); % Random ampltiude

y = zeros(size(t));
for i = 1:numSig
    y = y + A(i)*cos(rand(1)*randi(500)*t+phi(i));
end

%% Histogram counts
[Nx,edgesx] = histcounts(x, nbins);
[Ny,edgesy] = histcounts(y, nbins);

Phist = Nx/sum(Nx);
Qhist = Ny/sum(Ny);


%% Distance/ similarity
for i = 1:length(methods)
    method = methods{i};
    [dhist, shist, R] = pdfsDistanceSimilarity(Phist, Qhist , method, 2);
    Output(i).method = method;
    Output(i).distance = dhist;
    Output(i).similarity = shist;
    Output(i).corr = R;

end

%% Plot the P and Q
figure;
subplot(2,2,1)
plot(t,x)
xlabel('Time (t)'); ylabel('Amplitude'); title('Random signal 1')

subplot(2,2,2)
histogram(Phist,nbins)
xlabel('Values'); ylabel('Count'); title('Probability Density Function 1 (Histogram)')

subplot(2,2,3)
plot(t,y)
xlabel('Time (t)'); ylabel('Amplitude'); title('Random signal 2')

subplot(2,2,4)
histogram(Qhist,nbins)
xlabel('Values'); ylabel('Count'); title('Probability Density Function 2 (Histogram)')

%% End parameters
%--------------------------------------------------------------------------
Runtime = toc(Start);

