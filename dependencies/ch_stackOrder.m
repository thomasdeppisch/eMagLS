function ms = ch_stackOrder(order)
    ms = 0;
    for m = 1 : order
        ms = [ms, -m, m]; %#ok<*AGROW> 
    end
end
