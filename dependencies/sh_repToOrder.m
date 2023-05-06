function out = sh_repToOrder(in)
% replicate order weights from size (n-1) to size (n+1)^2 along first
% dimension
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

n = size(in,1) - 1;
l = (n+1)^2;
out = zeros(l, size(in,2), 'like', in);
%nm2acn = @(n_,m_) n_.^2 + m_ + n_ + 1;

for nn = 0:n
    for mm = -nn:nn
        %out(nm2acn(nn,mm),:) = in(nn+1,:);
        out(nn.^2+mm+nn+1, :) = in(nn+1,:);
    end
end

end
