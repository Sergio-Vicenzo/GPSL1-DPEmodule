 cnslxyz = llh2xyz([22.328337/180 * pi, 114.171328/180 * pi, 7]);

for i=1:length(navSolutions.latitude)
    
pos_llh(i,:) = [navSolutions.latitude(i) navSolutions.longitude(i) navSolutions.height(i)];
pos_xyz(i,:) = [navSolutions.X(i), navSolutions.Y(i), navSolutions.Z(i)];
pos_enu(i,:) = xyz2enu(pos_xyz(i,:),cnslxyz);
end