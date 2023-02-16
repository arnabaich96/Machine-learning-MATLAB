clear;
load('DataFile1.mat');
D1 = D1;

%Initialization
[m,n] = size(D1);
T = 100;
sweep_num = 20;
sweep_list = [5,10,15,20];
sigma1 = 100;
sigma2 = 30;

%Viewing image
figure(1);
clf;
subplot(3,2,1);
imagesc(D1);
colormap(gray);
title('Initial image');
hold on;

%Modifying matrix for easy calculation
D = zeros(m+2,n+2);
D(2:m+1,2:n+1) = D1;
I = D;

%Loop
k = 1;
for s=1:sweep_num
    for i=2:m+1
        for j=2:n+1
            neighbor_val = [I(i-1,j),I(i+1,j),I(i,j-1),I(i,j+1)];
            mu = mean(neighbor_val);
            I_gibbs(1) = normrnd(mu,sigma1);
            for t=1:(T-1)
                I_prop = normrnd(mu,sigma1);
                frac = (f(I_prop,mu,D(i,j),sigma1,sigma2)*q(I_gibbs(t),mu,sigma1))/(f(I_gibbs(t),mu,D(i,j),sigma1,sigma2)*q(I_prop,mu,sigma1));
                rho = min(frac,1);
                if binornd(1,rho)==1
                    I_gibbs(t+1) = I_prop;
                else
                    I_gibbs(t+1) = I_gibbs(t);
                end
            end
            I(i,j) = I_gibbs(end);
        end
    end
    
    if ismember(s,sweep_list)==1
        k = k+1;
        subplot(3,2,k);
        imagesc(I);
        colormap(gray);
        title(sprintf('Sweep %d',s));
    end
end

% Target density
function y=f(x,m,d,sigma1,sigma2)
y = exp(-(1/2)*((x-m)/sigma1)^2) * exp(-(1/2)*((x-d)/sigma2)^2);
end

% Proposal density
function z=q(x,m,sigma1)
z = exp(-(1/2)*((x-m)/sigma1)^2);
end
