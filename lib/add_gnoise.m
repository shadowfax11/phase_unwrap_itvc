function noisy_map = add_gnoise(map,noise_lvl)
%ADD_GNOSIE Summary of this function goes here
%   Detailed explanation goes here
noisy_map = map + noise_lvl*randn(size(map));
end

