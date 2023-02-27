function ids = ch_fromShIds(order)
    ids = 1;
    for n = 1 : order
        ids  = [ids, n^2+1, n^2+n+n+1]; %#ok<AGROW>
    end
end
