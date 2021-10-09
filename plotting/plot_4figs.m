function plot_4figs(map1,map2,map3,map4,colorbar_flag,title_str)
%PLOT_2FIGS Summary of this function goes here
%   Detailed explanation goes here
figure;
if nargin>5
    if colorbar_flag
        subplot(2,2,1);
        imagesc(map1); title(title_str{1}); colorbar;
        subplot(2,2,2); 
        imagesc(map2); title(title_str{2}); colorbar;
        subplot(2,2,3);
        imagesc(map3); title(title_str{3}); colorbar;
        subplot(2,2,4);
        imagesc(map4); title(title_str{4}); colorbar;
    else
        subplot(2,2,1);
        imagesc(map1); title(title_str{1});
        subplot(2,2,2); 
        imagesc(map2); title(title_str{2});
        subplot(2,2,3);
        imagesc(map3); title(title_str{3});
        subplot(2,2,4);
        imagesc(map4); title(title_str{4});
    end
else
    if colorbar_flag
        subplot(2,2,1);
        imagesc(map1); colorbar;
        subplot(2,2,2); 
        imagesc(map2); colorbar;
        subplot(2,2,3);
        imagesc(map3); colorbar;
        subplot(2,2,4);
        imagesc(map4); colorbar;
    else
        subplot(2,2,1);
        imagesc(map1);
        subplot(2,2,2); 
        imagesc(map2);
        subplot(2,2,3);
        imagesc(map3);
        subplot(2,2,4);
        imagesc(map4);
    end
end
end

