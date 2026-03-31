function [rho,theta2,N1,N2] = refl(n1,k1,n2,k2,theta1)

N1 = sqrt(0.5*(sqrt((n1.^2-k1.^2).^2+((2.*n1.*k1)./cos(theta1)).^2)+(n1.^2-k1.^2)));
N2 = sqrt(0.5*(sqrt((n2.^2-k2.^2-N1.^2.*sin(theta1).^2).^2+4.*n2.^2.*k2.^2)+(n2.^2-k2.^2+N1.^2.*sin(theta1).^2)));
K1 = (n1.*k1)./(N1.*cos(theta1));

theta2 = asin(N1.*sin(theta1)./N2);

K2 = (n2.*k2)./(N2.*cos(theta2));

r_perp  = ((N1.*cos(theta1)-1i.*K1)-(N2.*cos(theta2)-1i.*K2))./((N1.*cos(theta1)-1i.*K1)+(N2.*cos(theta2)-1i.*K2));
r_par  = ((n2-1i.*k2).^2.*(N1.*cos(theta1)-1i.*K1)-(n1-1i.*k1).^2.*(N2.*cos(theta2)-1i.*K2))./((n2-1i.*k2).^2.*(N1.*cos(theta1)-1i.*K1)+(n1-1i.*k1).^2.*(N2.*cos(theta2)-1i.*K2));

rho_perp = r_perp.*conj(r_perp);
rho_par = r_par.*conj(r_par);

rho = 0.5.*(abs(r_perp).^2 + abs(r_par).^2);


end

