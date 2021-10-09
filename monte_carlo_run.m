NUM_MAPS = 50;

addpath ./lib
for n=1:NUM_MAPS
    [~,~,~,P,~,~] = generate_terrain(7,256,0,25,10);
    for noise=[pi/6,pi/4]
        Pn = P + noise*randn(size(P));
        if noise==pi/6
            noise_str = 'pi6';
        else
            noise_str = 'pi4';
        end
        save_str = [noise_str,'-map-',num2str(n),'-fine'];
        wrapper_unwrap_itv3(P,Pn,'IHTV',save_str);
        wrapper_unwrap_itv3(P,Pn,'ITVC',save_str);
        wrapper_unwrap_itv(P,Pn,'ITV',save_str);
        save_str = [noise_str,'-map-',num2str(n),'-coarse'];
        wrapper_unwrap_itv3(P(1:4:end,1:4:end),Pn(1:4:end,1:4:end),'IHTV',save_str);
        wrapper_unwrap_itv3(P(1:4:end,1:4:end),Pn(1:4:end,1:4:end),'ITVC',save_str);
        wrapper_unwrap_itv(P(1:4:end,1:4:end),Pn(1:4:end,1:4:end),'ITV',save_str);
    end
end