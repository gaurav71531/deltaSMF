% The code to reproduce the Figure 2 of the paper titled:
% Approximate submodular functions and performance guarantees"
% Authors: Gaurav Gupta, Sergio Pequito and Paul Bogdan

% The mathematical details are provided in the Section 4.1, and the
% notations are same unless stated otherwise

clear
n = 10; % dimensions of the sensor vector
N = 30; % total number of sensors

betaVec = 0.05:0.1:20; %beta as appearing in the Section-
perB1 = zeros(1,length(betaVec));
perB2 = zeros(1,length(betaVec));

for betaInd  = 1:length(betaVec)
    beta = betaVec(betaInd);
    X = randn(n,N); %randomly generate the entries of X matrix and normalize
    for i = 1:N
        X(:,i) = X(:,i)/norm(X(:,i));
    end

    W = @(S) eye(n)*beta^2 + X(:,S)*X(:,S)'; %define Gramian

    f1 = @(S) log(det(W(S))); % define surrogate submodular function

    delta_u = 1/beta^2;

    eigTemp = eig(W([1:N]));
    delta_l = 1/max(eigTemp);

    valMin = Inf;
    val = zeros(1,N);
    for i = 1:N
        setTemp = 1:N;
        setTemp(i) = [];
        val(i) = (f1(1:N) - f1(setTemp))/f1(i);
    end
    curvF1 = 1- min(val);  %curvature of the surrogate submodular function

    eigTemp = eig(X*X');
    eigUse = max(eigTemp);
    
    %lower- and upper-bound of submodularity ratio and curvature as derived
    %in paper titled "Guarantees for Greedy Maximization of
    %Non-submodular Functions with Applications"
    %link: https://arxiv.org/abs/1703.02100
    
    gammaB = beta^2/eigUse/(beta^2 + eigUse); 
    alphaB = 1 - gammaB;

%     curvature and submodularity ratio equivalence for the proposed
%     delta-submodularity using Table-1 of the paper
    alphaUse = 1 - delta_l/delta_u*(1-curvF1);
    gammaUse = delta_l/delta_u;
    
    perB1(betaInd) = 1/alphaUse*(1-exp(-alphaUse*gammaUse));
    perB2(betaInd) = 1/alphaB*(1 - exp(-alphaB*gammaB));
    
end
figure;

plot(betaVec, perB1, betaVec, perB2);grid;
hold on;
plot(betaVec, ones(1,length(betaVec))*(1-1/exp(1)), 'k');
xlabel('\beta');
ylabel('Performance bound');
legend({'delta-approximate', 'parameters bound'})