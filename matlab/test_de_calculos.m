clc;
B1=0.7899944821689744;
B2=0.7084696621689741;
C1=(-0.4193241427264186-0.09040962826817195j);
C2=(0.22787597862192258-0.3186513314157572j) ;

r_max_in =  ((B1 - sqrt((B1^2)-(4*(abs(C1)^2))))/(2*(abs(C1))))
abs(r_max_in)
%r_max_in =  ((B1 + sqrt((B1^2)-(4*(abs(C1)^2))))/(2*(abs(C1))))

[thetha,rho] = cart2pol(real(r_max_in),imag(r_max_in));
rho
thetha = rad2deg(thetha)

r_max_out =  ((B2 - sqrt((B2^2)-(4*(abs(C2)^2))))/(2*(abs(C2))))
%r_max_out =  ((B2 + sqrt((B2^2)-(4*(abs(C2)^2))))/(2*(abs(C2))))

[thetha,rho] = cart2pol(real(r_max_out),imag(r_max_out));
rho
thetha = rad2deg(thetha)