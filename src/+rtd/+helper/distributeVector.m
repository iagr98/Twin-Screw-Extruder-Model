function T = distributeVector(N, n, Temps)
    if length(Temps) ~= n
        error('The length of vector Temps must coincide with n.');
    end    
    baseSize = floor(N / n);    
    extraElements = mod(N, n);
    T = zeros(1, N);
    startIdx = 1;
    
    for i = 1:n        
        currentPartSize = baseSize + (i <= extraElements);       
        endIdx = startIdx + currentPartSize - 1;            
        T(startIdx:endIdx) = Temps(i);        
        startIdx = endIdx + 1;
    end
end
