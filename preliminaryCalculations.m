%Flicker photometry reliability measures

% Results are from dce's self test on 7/23/19. 
% Settings: flickerFrequency = 20Hz, adjustCone = M, steadyConeContrast =
% 0.09

% Final selected M contrast values for each of the 10 trials in order
finalAdjustments = [0.0249 0.0580 0.0166 0.0332 0.0497 0.0332 0.0332 0.0414 0.0414 0.0414]; 

% Mean and standard deviation: 0.0373 +- 0.0119
contrastMean = mean(finalAdjustments)
contrastStd = std(finalAdjustments)

% Mean and standard deviation of last 7 values: 0.0391 +- 0.0062
shortArr = finalAdjustments(4:end);
contrastMeanLast7 = mean(shortArr)
contrastStdLast7 = std(shortArr)