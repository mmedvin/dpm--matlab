% smoothstep(t) returns:
%     0   if t < 0
%     1   if t > 1
%    p(t) if 0 < t < 1, where p(t) is polynomial that smoothly goes from 0 to 1
%
function [y,dy] = smoothstep(t)
    s = min((t>0).*t,1);
    y = (126+(-420+(540+(-315+70*s).*s).*s).*s).*s.^5;
    dy = (630 + (-2520 + (3780 + (630*s - 2520).*s).*s ).*s).*s.^4;
end