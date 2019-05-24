function val = comp(tode,xode)
% load param
% phi   = xode(4)
% theta = xode(5)
% psi   = xode(6)
% p     = xode(1)
% q     = xode(2)
% r     = xode(3)
% u     = xode(12)
% v     = xode(10)
% w     = xode(8)
% x     = xode(11)
% y     = xode(9)
% z     = xode(7)

% val = zeros(ncgl+1,12);
% global itr
% itr = itr+1;
% if(itr == ncgl+2)
%     itr = 1;
% end
global m g Ix Iy Iz fwx fwy fwz
global U ptspan
global k
% u_ = spline(ptspan, U, tode);
u_ = -k*(xode);
u_(1) = u_(1) + m*g;
val = zeros(1,12);
if(tode< 0.367)
    f_pert_x = 0.5;
else
    f_pert_x = 0.0;
end

% for i = 1:41
    val(4)  = xode(1) + xode(3)*cos(xode(4))*tan(xode(5)) + xode(2)*sin(xode(4))*tan(xode(5));
    val(5)  = xode(2)*cos(xode(4)) - xode(3)*sin(xode(4));
    val(6)  = xode(3)*(cos(xode(4))/cos(xode(5))) + xode(2)*(sin(xode(4))/cos(xode(5)));
    val(1)  = ((Iy - Iz)*xode(3)*xode(2))/Ix + u_(2)/Ix;  %(tow_x + tow_wx)/Ix; 
    val(2)  = ((Iz - Ix)*xode(1)*xode(3))/Iy + u_(3)/Iy; %(tow_y + tow_wy)/Iy;
    val(3)  = ((Ix - Iy)*xode(1)*xode(2))/Iz + u_(4)/Iz; %(tow_z + tow_wz)/Iz;
    val(12)  = xode(3)*xode(10) - xode(2)*xode(8) - g*sin(xode(5)) + fwx/m;
    val(10)  = xode(1)*xode(8) - xode(3)*xode(12) + g*sin(xode(4))*cos(xode(5)) + fwy/m;
    val(8)  = xode(2)*xode(12) - xode(1)*xode(10) + g*cos(xode(5))*cos(xode(4)) + (fwz - u_(1))/m - f_pert_x;  %-ft)/m;
    val(11) = xode(8)*(sin(xode(4))*sin(xode(6)) + cos(xode(4))*cos(xode(6))*sin(xode(5))) - xode(10)*(cos(xode(4))*sin(xode(6)) - cos(xode(6))*sin(xode(4))*sin(xode(5))) + xode(12)*cos(xode(6))*cos(xode(5));
    val(9) = xode(10)*(cos(xode(4))*cos(xode(6)) + sin(xode(4))*sin(xode(6))*sin(xode(5))) - xode(8)*(cos(xode(6))*sin(xode(5)) - cos(xode(4))*sin(xode(6))*sin(xode(5))) + xode(12)*cos(xode(5))*sin(xode(6));
    val(7) = xode(8)*cos(xode(4))*cos(xode(5)) - xode(12)*sin(xode(5)) + xode(10)*cos(xode(5))*sin(xode(4));
% end
% val = reshape(val',[12*(ncgl+1),1]);
val = val';
val;