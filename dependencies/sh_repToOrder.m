function out = sh_repToOrder(in)
% replicate order weights from size (n-1) to size (n+1)^2 along first
% dimension
n = size(in,1) - 1;
l = (n+1)^2;
out = zeros(l, size(in,2));
%nm2acn = @(n_,m_) n_.^2 + m_ + n_ + 1;

for nn = 0:n
    for mm = -nn:nn
        %out(nm2acn(nn,mm),:) = in(nn+1,:);
        out(nn.^2+mm+nn+1, :) = in(nn+1,:);
    end
end

end