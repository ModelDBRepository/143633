function [H] = shannon_entropy(X,bins)

    XH = histc(X,bins);

    p = XH/sum(XH);
    h = zeros(size(p));
    ind = find(p~=0);
    h(ind) = -(p(ind)).*log2(p(ind));
    H = sum(h);
end
