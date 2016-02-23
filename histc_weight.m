function [N,BIN] = histc_weight(X,W,edges)

[N,BIN] = histc(X,edges);

for i = 1:length(N)
    N(i) = sum(W(BIN==i));
end
