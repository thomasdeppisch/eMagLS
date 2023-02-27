function ids = eh_fromShIds(order)
    ids = 1;
    for n = 1 : order
        for m = -n : 2 : n % each (n+m)==even
            ids = [ids, n^2+n+m+1]; %#ok<AGROW> 
        end
    end
end
