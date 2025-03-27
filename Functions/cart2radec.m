function [mag, RA, DEC] = cart2radec(r)
% Converts x,y,z Cartesian coordinates to spherical RA and DEC coordinates

mag = norm(r);
RA = atan2(r(2), r(1));
DEC = asin(r(3)/mag);

end