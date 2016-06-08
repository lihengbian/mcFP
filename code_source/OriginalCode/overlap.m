function y = overlap(NA, H,LEDp)
%%function y = overlap(NA, H,LEDp)
%%calculate overlap in pupil
a=NA;
b=sin(atan(LEDp/H/2));
c=abs(sqrt(a^2-b^2));
Stri=c*b;
Sfan=acos(b/a)*a^2;
y=2*(Sfan-Stri)/pi/a^2;
end
