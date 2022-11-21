%% rdmOrder
% Use subject ID for assignment of pseudorandom task Order (Sternberg & N-back)
if mod(subject.ID,2) == 0       
    SternbergNBack = 0;
elseif mod(subject.ID,2) == 1
    SternbergNBack = 1;
end
