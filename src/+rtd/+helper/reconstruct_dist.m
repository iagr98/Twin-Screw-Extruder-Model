function [n, p, f] = reconstruct_dist(mue_0,mue_1,mue_2)
% reconstructs the molecular weight distribution from moments 0 to 2 using
% a modified gamma distribution (ref: Yu2014)
    N = 10e4; % maximal chain length
    n = 1:1:N; % chain length vector
    
    a = mue_1(end)/mue_0(end);
    b = a^2/(mue_2(end)/mue_0(end)-a^2);
    
    z = b/a*n;
    
    p = 1/factorial(round(b-1))*z.^(b-1).*exp(-z);
    f = p/sum(p)*mue_0(end); % distribution

end