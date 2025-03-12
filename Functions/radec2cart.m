function [r] = radec2cart(mag, RA, DEC)

x = mag*cos(RA)*cos(DEC);
y = mag*sin(RA)*cos(DEC);
z = mag*sin(DEC);
r = [x,y,z];

end