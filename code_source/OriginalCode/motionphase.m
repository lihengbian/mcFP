function addphase = motionphase(x0, y0, M, N)

temp1 = [0:M-1]';
temp1 = repmat(temp1, [1, N]);
temp2 = [0:N-1];
temp2 = repmat(temp2, [M, 1]);

addphase = exp(1j*2*pi*(temp1*x0/M + temp2*y0/N));

end