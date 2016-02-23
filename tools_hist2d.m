function [mHist,x_mid,y_mid] = tools_hist2d (x, y, xEdge, yEdge)

% Usages:
%  mHist = tools_hist2d (x, y, xEdge, yEdge) computes a 2d
%  histogram, which counts number of points in the bins defined by "xEdge" and
%  "yEdge".
%
% Inputs:
%  x, y are column vectors
%
% Outputs:
%  mHist = size [length(yEdge)-1, length(xEdge)-1]
%
% Tian Tsong Ng, July 2005

x_len = length(xEdge)-1;
y_len = length(yEdge)-1;
mHist = zeros(y_len,x_len);


x_mid = (xEdge(1:end-1) + xEdge(2:end))/2;
y_mid = (yEdge(1:end-1) + yEdge(2:end))/2;

if isempty(x)
    mHist = zeros(y_mid,x_mid);
    return
end


[NY,binY] = histc(y,yEdge);

for yC = 1:y_len       
    x_selected = x(binY==yC);
    if ~isempty(x_selected)
        [NX,binX] = histc(x_selected,xEdge);
        NX = NX(:)';
        mHist(yC,:) = NX(1:end-1);
    end    
end