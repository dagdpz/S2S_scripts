function plot_convex_target_shape(convexities,convex_sides,filenames)
% plot_convex_target_shape({-0.5,-0.2,0.2,0.5,-0.5,-0.5},{'LR','LR','LR','LR','L','R'},{'tar1','tar2','tar3','tar4','tar5','tar6'})
for k=1:numel(convexities)
convexity=convexities{k};
convex_side=convex_sides{k};
    
figure
pointList=CalculateConvexPointList([0,0],1,convexity,convex_side);
patch(pointList(:,1),pointList(:,2),'r');
set(gca,'xlim',[-2,2],'ylim',[-2,2]);
export_fig(filenames{k}, '-pdf','-transparent') % append to existing pdf
close(gcf)
end
