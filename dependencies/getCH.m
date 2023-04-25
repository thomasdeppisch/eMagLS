function Y = getCH(order, gridAziRad, shDefinition)
% yield circular harmonics basis functions
% 
% Hannes Helmholz, 2022

if strcmpi(shDefinition, 'complex')
    m = ch_stackOrder(order);
    Y = exp(-1i .* m .* gridAziRad);
elseif strcmpi(shDefinition, 'real')
    Y = ones(2*order+1, length(gridAziRad));
    for m = 1 : order
        Y(2*m, :) = sqrt(2) .* sin(m .* gridAziRad);
        Y(2*m+1, :) = sqrt(2) .* cos(m .* gridAziRad);
    end
else
    error('Unknown shDefinition "%s".', shDefinition);
end

end
