% Parameters
freq = 10; 
strengthScaler = 0.3; 
window = 20; 
k = 2; % # of std signal appears

% Step 1: Simulate factor returns 
rng(42); 
dates = datetime(2023,1,1) + caldays(0:299); 
FactorRtns = randn(1, 300); 
FactorRtns = FactorRtns + strengthScaler * sin(linspace(0, freq, 300)); 

% Step 2: Calculate rolling mean and std
rolling_mean = movmean(FactorRtns, window); 
rolling_std = movstd(FactorRtns, window); 
upper_threshold = rolling_mean + k * rolling_std; 
lower_threshold = rolling_mean - k * rolling_std; 

% Step 3: Generate signals
signals = zeros(1, 300); 
for i = window:300 
    if FactorRtns(i) > upper_threshold(i)
        signals(i) = -1; 
    elseif FactorRtns(i) < lower_threshold(i)
        signals(i) = 1; 
    end
end

% Step 4: Calculate strategy returns
StratRtns = zeros(1, 300); 
StratRtns(2:end) = signals(1:end-1) .* diff(FactorRtns); 

% Step 5: Plot Factor Returns and Signals
figure;
plot(dates, FactorRtns, 'DisplayName', 'Factor Returns'); hold on;
plot(dates, rolling_mean, '--', 'DisplayName', 'Rolling Mean');
plot(dates, upper_threshold, ':', 'DisplayName', 'Upper Threshold');
plot(dates, lower_threshold, ':', 'DisplayName', 'Lower Threshold');
title('Mean Reversion Strategy Using Factor Returns');
xlabel('Date');
ylabel('Factor Returns');
legend('Location', 'best');
grid on;

% Step 6: Plot Cumulative Strategy Returns
CumStratRtns = cumprod(1 + StratRtns); % Cumulative returns
figure;
plot(dates, CumStratRtns, 'g', 'DisplayName', 'Cumulative Strategy Returns');
title('Cumulative Returns of Factor-Based Mean Reversion Strategy');
xlabel('Date');
ylabel('Cumulative Returns');
legend('Location', 'best');
grid on;
