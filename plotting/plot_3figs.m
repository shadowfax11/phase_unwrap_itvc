function plot_3figs(map1,map2,map3,colorbar_flag,title_str)
%PLOT_2FIGS Summary of this function goes here
%   Detailed explanation goes here
figure;
if nargin>4
    if colorbar_flag
        subplot(1,3,1);
        imagesc(map1); title(title_str{1}); colorbar;
        subplot(1,3,2); 
        imagesc(map2); title(title_str{2}); colorbar;
        subplot(1,3,3); 
        imagesc(map3); title(title_str{3}); colorbar;
    else
        subplot(1,3,1);
        imagesc(map1); title(title_str{1});
        subplot(1,3,2); 
        imagesc(map2); title(title_str{2});
        subplot(1,3,3);
        imagesc(map3); title(title_str{3});
    end
else
    if colorbar_flag
        subplot(1,3,1);
        imagesc(map1); colorbar;
        subplot(1,3,2); 
        imagesc(map2); colorbar;
        subplot(1,3,3); 
        imagesc(map3); colorbar;
    else
        subplot(1,3,1);
        imagesc(map1);
        subplot(1,3,2); 
        imagesc(map2);
        subplot(1,3,3);
        imagesc(map3);
    end
end
end

