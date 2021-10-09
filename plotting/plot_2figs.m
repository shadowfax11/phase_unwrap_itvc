function plot_2figs(map1,map2,colorbar_flag,title_str)
%PLOT_2FIGS Summary of this function goes here
%   Detailed explanation goes here
figure;
if nargin>3
    if colorbar_flag
        subplot(1,2,1);
        imagesc(map1); title(title_str{1}); colorbar;
        subplot(1,2,2); 
        imagesc(map2); title(title_str{2}); colorbar;
    else
        subplot(1,2,1);
        imagesc(map1); title(title_str{1});
        subplot(1,2,2); 
        imagesc(map2); title(title_str{2});
    end
else
    if colorbar_flag
        subplot(1,2,1);
        imagesc(map1); colorbar;
        subplot(1,2,2); 
        imagesc(map2); colorbar;
    else
        subplot(1,2,1);
        imagesc(map1);
        subplot(1,2,2); 
        imagesc(map2);
    end
end
end

