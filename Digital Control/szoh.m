function [tzoh, xzoh] = szoh(td, xd)
%% Simplified ZOH
% xd(i) is the i'th sampled data at time td(i).
% td: time points
% xd: discrete x

% tzoh
% xzoh: zoh'd x

N = numel(xd);
xzoh = zeros(1,2*N);
xzoh(1:2:2*N-1) = xd';
xzoh(2:2:2*N) = xd';
xzoh(end) = [];

tzoh = zeros(1,2*N);
tzoh(1:2:2*N-1) = td;
tzoh(2:2:2*N) = td;
tzoh(1) = [];

end

