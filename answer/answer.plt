series(x,y,n) = sum[k=1:n] u0*2.0*(1-(-1)**k)*sin(k*pi*x/L)*(exp(k*pi*y/L)-exp(-k*pi*y/L))/(k*pi*(exp(k*pi)-exp(-k*pi)))

u0 = 10.0
L = 1.0

set xrange [0:L]                                                                                                                                                                                                                               
set yrange [0:L]

splot series(x,y,100)
