function y = zscore_per_channel(X)
mu = mean(X,1); sd = std(X,[],1); y = (X - mu) ./ (sd + eps);
end
