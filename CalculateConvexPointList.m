function pointList=CalculateConvexPointList(center,a_ellipse,convexity,convex_sides)
convexity_sign=sign(convexity);
b_ellipse=abs(convexity)*a_ellipse;
b_ellipse_reference=0.3*a_ellipse;

steps_for_interpolation=100;
euclidian_distance=(-a_ellipse):2*a_ellipse/steps_for_interpolation:(a_ellipse);
%corner_positions=[-1,-1;1,-1;1,1;-1,1].*a_ellipse;

Half_circle_Area=a_ellipse^2*pi/2;
rect_b=(Half_circle_Area-a_ellipse*b_ellipse*pi*convexity_sign/2)/(2*a_ellipse);
rect_b_reference=(Half_circle_Area-a_ellipse*b_ellipse_reference*pi*convexity_sign/2)/(2*a_ellipse);
%rect_ratio=Half_circle_Area/(2*a_ellipse^2);

tmp_bow_parameter=sqrt((b_ellipse)^2.*(1-euclidian_distance'.^2/a_ellipse^2));
tmp_bow_reference=sqrt((b_ellipse_reference)^2.*(1-euclidian_distance'.^2/a_ellipse^2));

bow_vector=[euclidian_distance',(tmp_bow_parameter.*convexity_sign.*-1-rect_b)];
bow_vector_reference=[euclidian_distance',(tmp_bow_reference.*convexity_sign.*-1-rect_b_reference)];
         
if strcmp(convex_sides,'LR')|| strcmp(convex_sides,'R') || strcmp(convex_sides,'L')
   bow_vector=[bow_vector(:,2)*-1,bow_vector(:,1)];
   bow_vector_reference=[bow_vector_reference(:,2)*-1,bow_vector_reference(:,1)];
   %corner_positions=[corner_positions(:,1)*rect_ratio,corner_positions(:,2)];
else
   %corner_positions=[corner_positions(:,1),corner_positions(:,2)*rect_ratio];
end

switch convex_sides
    case 'T'       
        pointList=[bow_vector;bow_vector_reference*-1];
    case 'B'       
        pointList=[bow_vector_reference;bow_vector*-1];
    case 'TB'
        pointList=[bow_vector;bow_vector*-1];
    case 'R'       
        pointList=[bow_vector;bow_vector_reference*-1];
    case 'L'       
        pointList=[bow_vector_reference;bow_vector*-1];   
    case 'LR'
        pointList=[bow_vector;bow_vector*-1];
end
%Position_correction=[(max(pointList(:,1))+min(pointList(:,1)))/2,(max(pointList(:,2))+min(pointList(:,2)))/2];
%pointList=pointList+repmat(center-Position_correction,size(pointList,1),1);
pointList=pointList+repmat(center,size(pointList,1),1);
