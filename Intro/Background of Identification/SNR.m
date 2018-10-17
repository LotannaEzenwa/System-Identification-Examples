%% Signal-to-Noise Ratio
% The signal to noise ratio is defined as the variance of the signal over
% the variance of the noise

%%
% 
% $$SNR = \frac{Variance of Signal}{Variance of Noise} $$
% 

%% Example
% This example shows the relationsship between SNR and identification
% (regression).
%
% Let 
%
% $$ y[k] = b_{1}u[k] + b_{0} $$
%
% $$ y_{m}[k] = y[k] + v[k] $$
%
% where $b_{1} = 3$ and $b_{0} = 4$, $y_{m}[k]$ is the measured output, and
% $v[k]$ is the sensor noise

%%
% Create Random Input-Output Data
clf
u = rand(200,1);
y = 5*u + 2;

% Add noise

e = randn(200,1);
ym1 = y + e*sqrt(var(y)/100); % SNR 100
ym2 = y + e*sqrt(var(y)/10); % SNR 10

% Plot Measurements

subplot(1,2,1); 
plot(u, ym1, 'x', 'linewidth', 1);
subplot(1,2,2);
plot(u, ym2, 'x', 'linewidth', 2);

polyfit(u, ym1, 1)
polyfit(u, ym2, 1)

% Error Estimates

Phi = [u ones(200,1)];
thetah1 = Phi \ ym1;
err1 = ym1 - Phi*thetah1;
errth1 = sqrt(diag(inv(Phi'*Phi)*sum(err1.^2)/198))


thetah2 = Phi \ ym2;
err2 = ym2 - Phi*thetah2;
errth2 = sqrt(diag(inv(Phi'*Phi)*sum(err2.^2)/198))


