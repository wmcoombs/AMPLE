function plot_rects(mesh,mpData,crit)
    nmp = length(mpData);
    pos = [mpData.mpC];
    lp = [mpData.lp];
    pos_x = pos(1:2:end);
    pos_y = pos(2:2:end);
    lp_x = lp(1:2:end)*2;
    lp_y = lp(2:2:end)*2;
    data_split_l = split_critera(mesh,mpData,crit);
    data_split_l(data_split_l==1) = 1;
    data_split_l(data_split_l==2) = 1;
    data_split_l(data_split_l==3) = 1;
    data_split = data_split_l;
    colours = zeros(nmp,3);
    positions = [(pos_x-lp_x*0.5)', (pos_y-lp_y*0.5)', lp_x',lp_y'];
    colours(:,1) = data_split;
    colours(:,2) = data_split==0;
    cla;
    for i=1:nmp
    rectangle('Position', positions(i,:), 'FaceColor', colours(i,:));
    end
    colormap(colours);
    xlim([0,mesh.h(1)]);
    caxis([0,1]);
    xlim([0,32]);
    ylim([0,16]);
    drawnow;
    pause(0.5)
end

