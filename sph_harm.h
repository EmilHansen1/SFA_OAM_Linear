#include <cmath>
#include <complex>

using dcmplx = std::complex<double>;
constexpr dcmplx I(0., 1.);
constexpr double Pi = 3.14159265358979323846;

dcmplx sph_harm(dcmplx x, dcmplx y, dcmplx z, dcmplx r, int l, int m){
if(std::abs(m) > l || l < 0){return 0.;}

if(l == 0){
return 1./(2.*std::sqrt(Pi));
}

else if (l == 1){
if(m == -1){
return (std::sqrt(3./(2.*Pi))*(x-I*y)/r)/(2.*1.);
}
else if (m == 0){
return (std::sqrt(3./Pi)*z/r)/2.;
}
else if (m == 1){
return -0.5*(1.*std::sqrt(3./(2.*Pi))*(x+I*y)/r);
}
else{return 0.;}
}

else if (l == 2){
if(m == -2){
return (std::sqrt(15./(2.*Pi))*std::pow((x-I*y)/r,2.))/(4.*1.);
}
else if (m == -1){
return (std::sqrt(15./(2.*Pi))*z/r*(x-I*y)/r)/(2.*1.);
}
else if (m == 0){
return (std::sqrt(5./Pi)*(-1. + 3.*std::pow(z/r,2.)))/4.;
}
else if (m == 1){
return -0.5*(1.*std::sqrt(15./(2.*Pi))*z/r*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(15./(2.*Pi))*std::pow((x+I*y)/r,2.))/4.;
}
else{return 0.;}
}

else if (l == 3){
if(m == -3){
return (std::sqrt(35./Pi)*std::pow((x-I*y)/r,3.))/(8.*1.);
}
else if (m == -2){
return (std::sqrt(105./(2.*Pi))*z/r*std::pow((x-I*y)/r,2.))/(4.*1.);
}
else if (m == -1){
return (std::sqrt(21./Pi)*(-1. + 5.*std::pow(z/r,2.))*(x-I*y)/r)/(8.*1.);
}
else if (m == 0){
return (std::sqrt(7./Pi)*(-3.*z/r + 5.*std::pow(z/r,3.)))/4.;
}
else if (m == 1){
return -0.125*(1.*std::sqrt(21./Pi)*(-1. + 5.*std::pow(z/r,2.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(105./(2.*Pi))*z/r*std::pow((x+I*y)/r,2.))/4.;
}
else if (m == 3){
return -0.125*(1.*std::sqrt(35./Pi)*std::pow((x+I*y)/r,3.));
}
else{return 0.;}
}

else if (l == 4){
if(m == -4){
return (3.*std::sqrt(35./(2.*Pi))*std::pow((x-I*y)/r,4.))/(16.*1.);
}
else if (m == -3){
return (3.*std::sqrt(35./Pi)*z/r*std::pow((x-I*y)/r,3.))/(8.*1.);
}
else if (m == -2){
return (3.*std::sqrt(5./(2.*Pi))*(-1. + 7.*std::pow(z/r,2.))*std::pow((x-I*y)/r,2.))/(8.*1.);
}
else if (m == -1){
return (3.*std::sqrt(5./Pi)*z/r*(-3. + 7.*std::pow(z/r,2.))*(x-I*y)/r)/(8.*1.);
}
else if (m == 0){
return (3.*(3. - 30.*std::pow(z/r,2.) + 35.*std::pow(z/r,4.)))/(16.*std::sqrt(Pi));
}
else if (m == 1){
return (-3.*1.*std::sqrt(5./Pi)*z/r*(-3. + 7.*std::pow(z/r,2.))*(x+I*y)/r)/8.;
}
else if (m == 2){
return (3.*1.*std::sqrt(5./(2.*Pi))*(-1. + 7.*std::pow(z/r,2.))*std::pow((x+I*y)/r,2.))/8.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(35./Pi)*z/r*std::pow((x+I*y)/r,3.))/8.;
}
else if (m == 4){
return (3.*1.*std::sqrt(35./(2.*Pi))*std::pow((x+I*y)/r,4.))/16.;
}
else{return 0.;}
}

else if (l == 5){
if(m == -5){
return (3.*std::sqrt(77./Pi)*std::pow((x-I*y)/r,5.))/(32.*1.);
}
else if (m == -4){
return (3.*std::sqrt(385./(2.*Pi))*z/r*std::pow((x-I*y)/r,4.))/(16.*1.);
}
else if (m == -3){
return (std::sqrt(385./Pi)*(-1. + 9.*std::pow(z/r,2.))*std::pow((x-I*y)/r,3.))/(32.*1.);
}
else if (m == -2){
return (std::sqrt(1155./(2.*Pi))*z/r*(-1. + 3.*std::pow(z/r,2.))*std::pow((x-I*y)/r,2.))/(8.*1.);
}
else if (m == -1){
return (std::sqrt(165./(2.*Pi))*(1. - 14.*std::pow(z/r,2.) + 21.*std::pow(z/r,4.))*(x-I*y)/r)/(16.*1.);
}
else if (m == 0){
return (std::sqrt(11./Pi)*(15.*z/r - 70.*std::pow(z/r,3.) + 63.*std::pow(z/r,5.)))/16.;
}
else if (m == 1){
return -0.0625*(1.*std::sqrt(165./(2.*Pi))*(1. - 14.*std::pow(z/r,2.) + 21.*std::pow(z/r,4.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(1155./(2.*Pi))*z/r*(-1. + 3.*std::pow(z/r,2.))*std::pow((x+I*y)/r,2.))/8.;
}
else if (m == 3){
return -0.03125*(1.*std::sqrt(385./Pi)*(-1. + 9.*std::pow(z/r,2.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(385./(2.*Pi))*z/r*std::pow((x+I*y)/r,4.))/16.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(77./Pi)*std::pow((x+I*y)/r,5.))/32.;
}
else{return 0.;}
}

else if (l == 6){
if(m == -6){
return (std::sqrt(3003./Pi)*std::pow((x-I*y)/r,6.))/(64.*1.);
}
else if (m == -5){
return (3.*std::sqrt(1001./Pi)*z/r*std::pow((x-I*y)/r,5.))/(32.*1.);
}
else if (m == -4){
return (3.*std::sqrt(91./(2.*Pi))*(-1. + 11.*std::pow(z/r,2.))*std::pow((x-I*y)/r,4.))/(32.*1.);
}
else if (m == -3){
return (std::sqrt(1365./Pi)*z/r*(-3. + 11.*std::pow(z/r,2.))*std::pow((x-I*y)/r,3.))/(32.*1.);
}
else if (m == -2){
return (std::sqrt(1365./Pi)*(1. - 18.*std::pow(z/r,2.) + 33.*std::pow(z/r,4.))*std::pow((x-I*y)/r,2.))/(64.*1.);
}
else if (m == -1){
return (std::sqrt(273./(2.*Pi))*z/r*(5. - 30.*std::pow(z/r,2.) + 33.*std::pow(z/r,4.))*(x-I*y)/r)/(16.*1.);
}
else if (m == 0){
return (std::sqrt(13./Pi)*(-5. + 105.*std::pow(z/r,2.) - 315.*std::pow(z/r,4.) + 231.*std::pow(z/r,6.)))/32.;
}
else if (m == 1){
return -0.0625*(1.*std::sqrt(273./(2.*Pi))*z/r*(5. - 30.*std::pow(z/r,2.) + 33.*std::pow(z/r,4.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(1365./Pi)*(1. - 18.*std::pow(z/r,2.) + 33.*std::pow(z/r,4.))*std::pow((x+I*y)/r,2.))/64.;
}
else if (m == 3){
return -0.03125*(1.*std::sqrt(1365./Pi)*z/r*(-3. + 11.*std::pow(z/r,2.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(91./(2.*Pi))*(-1. + 11.*std::pow(z/r,2.))*std::pow((x+I*y)/r,4.))/32.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(1001./Pi)*z/r*std::pow((x+I*y)/r,5.))/32.;
}
else if (m == 6){
return (1.*std::sqrt(3003./Pi)*std::pow((x+I*y)/r,6.))/64.;
}
else{return 0.;}
}

else if (l == 7){
if(m == -7){
return (3.*std::sqrt(715./(2.*Pi))*std::pow((x-I*y)/r,7.))/(64.*1.);
}
else if (m == -6){
return (3.*std::sqrt(5005./Pi)*z/r*std::pow((x-I*y)/r,6.))/(64.*1.);
}
else if (m == -5){
return (3.*std::sqrt(385./(2.*Pi))*(-1. + 13.*std::pow(z/r,2.))*std::pow((x-I*y)/r,5.))/(64.*1.);
}
else if (m == -4){
return (3.*std::sqrt(385./(2.*Pi))*z/r*(-3. + 13.*std::pow(z/r,2.))*std::pow((x-I*y)/r,4.))/(32.*1.);
}
else if (m == -3){
return (3.*std::sqrt(35./(2.*Pi))*(3. - 66.*std::pow(z/r,2.) + 143.*std::pow(z/r,4.))*std::pow((x-I*y)/r,3.))/(64.*1.);
}
else if (m == -2){
return (3.*std::sqrt(35./Pi)*z/r*(15. - 110.*std::pow(z/r,2.) + 143.*std::pow(z/r,4.))*std::pow((x-I*y)/r,2.))/(64.*1.);
}
else if (m == -1){
return (std::sqrt(105./(2.*Pi))*(-5. + 135.*std::pow(z/r,2.) - 495.*std::pow(z/r,4.) + 429.*std::pow(z/r,6.))*(x-I*y)/r)/(64.*1.);
}
else if (m == 0){
return (std::sqrt(15./Pi)*(-35.*z/r + 315.*std::pow(z/r,3.) - 693.*std::pow(z/r,5.) + 429.*std::pow(z/r,7.)))/32.;
}
else if (m == 1){
return -0.015625*(1.*std::sqrt(105./(2.*Pi))*(-5. + 135.*std::pow(z/r,2.) - 495.*std::pow(z/r,4.) + 429.*std::pow(z/r,6.))*(x+I*y)/r);
}
else if (m == 2){
return (3.*1.*std::sqrt(35./Pi)*z/r*(15. - 110.*std::pow(z/r,2.) + 143.*std::pow(z/r,4.))*std::pow((x+I*y)/r,2.))/64.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(35./(2.*Pi))*(3. - 66.*std::pow(z/r,2.) + 143.*std::pow(z/r,4.))*std::pow((x+I*y)/r,3.))/64.;
}
else if (m == 4){
return (3.*1.*std::sqrt(385./(2.*Pi))*z/r*(-3. + 13.*std::pow(z/r,2.))*std::pow((x+I*y)/r,4.))/32.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(385./(2.*Pi))*(-1. + 13.*std::pow(z/r,2.))*std::pow((x+I*y)/r,5.))/64.;
}
else if (m == 6){
return (3.*1.*std::sqrt(5005./Pi)*z/r*std::pow((x+I*y)/r,6.))/64.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(715./(2.*Pi))*std::pow((x+I*y)/r,7.))/64.;
}
else{return 0.;}
}

else if (l == 8){
if(m == -8){
return (3.*std::sqrt(12155./(2.*Pi))*std::pow((x-I*y)/r,8.))/(256.*1.);
}
else if (m == -7){
return (3.*std::sqrt(12155./(2.*Pi))*z/r*std::pow((x-I*y)/r,7.))/(64.*1.);
}
else if (m == -6){
return (std::sqrt(7293./Pi)*(-1. + 15.*std::pow(z/r,2.))*std::pow((x-I*y)/r,6.))/(128.*1.);
}
else if (m == -5){
return (3.*std::sqrt(17017./(2.*Pi))*z/r*(-1. + 5.*std::pow(z/r,2.))*std::pow((x-I*y)/r,5.))/(64.*1.);
}
else if (m == -4){
return (3.*std::sqrt(1309./(2.*Pi))*(1. - 26.*std::pow(z/r,2.) + 65.*std::pow(z/r,4.))*std::pow((x-I*y)/r,4.))/(128.*1.);
}
else if (m == -3){
return (std::sqrt(19635./(2.*Pi))*z/r*(3. - 26.*std::pow(z/r,2.) + 39.*std::pow(z/r,4.))*std::pow((x-I*y)/r,3.))/(64.*1.);
}
else if (m == -2){
return (3.*std::sqrt(595./Pi)*(-1. + 33.*std::pow(z/r,2.) - 143.*std::pow(z/r,4.) + 143.*std::pow(z/r,6.))*std::pow((x-I*y)/r,2.))/(128.*1.);
}
else if (m == -1){
return (3.*std::sqrt(17./(2.*Pi))*z/r*(-35. + 385.*std::pow(z/r,2.) - 1001.*std::pow(z/r,4.) + 715.*std::pow(z/r,6.))*(x-I*y)/r)/(64.*1.);
}
else if (m == 0){
return (std::sqrt(17./Pi)*(35. - 1260.*std::pow(z/r,2.) + 6930.*std::pow(z/r,4.) - 12012.*std::pow(z/r,6.) + 6435.*std::pow(z/r,8.)))/256.;
}
else if (m == 1){
return (-3.*1.*std::sqrt(17./(2.*Pi))*z/r*(-35. + 385.*std::pow(z/r,2.) - 1001.*std::pow(z/r,4.) + 715.*std::pow(z/r,6.))*(x+I*y)/r)/64.;
}
else if (m == 2){
return (3.*1.*std::sqrt(595./Pi)*(-1. + 33.*std::pow(z/r,2.) - 143.*std::pow(z/r,4.) + 143.*std::pow(z/r,6.))*std::pow((x+I*y)/r,2.))/128.;
}
else if (m == 3){
return -0.015625*(1.*std::sqrt(19635./(2.*Pi))*z/r*(3. - 26.*std::pow(z/r,2.) + 39.*std::pow(z/r,4.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(1309./(2.*Pi))*(1. - 26.*std::pow(z/r,2.) + 65.*std::pow(z/r,4.))*std::pow((x+I*y)/r,4.))/128.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(17017./(2.*Pi))*z/r*(-1. + 5.*std::pow(z/r,2.))*std::pow((x+I*y)/r,5.))/64.;
}
else if (m == 6){
return (1.*std::sqrt(7293./Pi)*(-1. + 15.*std::pow(z/r,2.))*std::pow((x+I*y)/r,6.))/128.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(12155./(2.*Pi))*z/r*std::pow((x+I*y)/r,7.))/64.;
}
else if (m == 8){
return (3.*1.*std::sqrt(12155./(2.*Pi))*std::pow((x+I*y)/r,8.))/256.;
}
else{return 0.;}
}

else if (l == 9){
if(m == -9){
return (std::sqrt(230945./Pi)*std::pow((x-I*y)/r,9.))/(512.*1.);
}
else if (m == -8){
return (3.*std::sqrt(230945./(2.*Pi))*z/r*std::pow((x-I*y)/r,8.))/(256.*1.);
}
else if (m == -7){
return (3.*std::sqrt(13585./Pi)*(-1. + 17.*std::pow(z/r,2.))*std::pow((x-I*y)/r,7.))/(512.*1.);
}
else if (m == -6){
return (std::sqrt(40755./Pi)*z/r*(-3. + 17.*std::pow(z/r,2.))*std::pow((x-I*y)/r,6.))/(128.*1.);
}
else if (m == -5){
return (3.*std::sqrt(2717./Pi)*(1. - 30.*std::pow(z/r,2.) + 85.*std::pow(z/r,4.))*std::pow((x-I*y)/r,5.))/(256.*1.);
}
else if (m == -4){
return (3.*std::sqrt(95095./(2.*Pi))*z/r*(1. - 10.*std::pow(z/r,2.) + 17.*std::pow(z/r,4.))*std::pow((x-I*y)/r,4.))/(128.*1.);
}
else if (m == -3){
return (std::sqrt(21945./Pi)*(-1. + 39.*std::pow(z/r,2.) - 195.*std::pow(z/r,4.) + 221.*std::pow(z/r,6.))*std::pow((x-I*y)/r,3.))/(256.*1.);
}
else if (m == -2){
return (3.*std::sqrt(1045./Pi)*z/r*(-7. + 91.*std::pow(z/r,2.) - 273.*std::pow(z/r,4.) + 221.*std::pow(z/r,6.))*std::pow((x-I*y)/r,2.))/(128.*1.);
}
else if (m == -1){
return (3.*std::sqrt(95./(2.*Pi))*(7. - 308.*std::pow(z/r,2.) + 2002.*std::pow(z/r,4.) - 4004.*std::pow(z/r,6.) + 2431.*std::pow(z/r,8.))*(x-I*y)/r)/(256.*1.);
}
else if (m == 0){
return (std::sqrt(19./Pi)*(315.*z/r - 4620.*std::pow(z/r,3.) + 18018.*std::pow(z/r,5.) - 25740.*std::pow(z/r,7.) + 12155.*std::pow(z/r,9.)))/256.;
}
else if (m == 1){
return (-3.*1.*std::sqrt(95./(2.*Pi))*(7. - 308.*std::pow(z/r,2.) + 2002.*std::pow(z/r,4.) - 4004.*std::pow(z/r,6.) + 2431.*std::pow(z/r,8.))*(x+I*y)/r)/256.;
}
else if (m == 2){
return (3.*1.*std::sqrt(1045./Pi)*z/r*(-7. + 91.*std::pow(z/r,2.) - 273.*std::pow(z/r,4.) + 221.*std::pow(z/r,6.))*std::pow((x+I*y)/r,2.))/128.;
}
else if (m == 3){
return -0.00390625*(1.*std::sqrt(21945./Pi)*(-1. + 39.*std::pow(z/r,2.) - 195.*std::pow(z/r,4.) + 221.*std::pow(z/r,6.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(95095./(2.*Pi))*z/r*(1. - 10.*std::pow(z/r,2.) + 17.*std::pow(z/r,4.))*std::pow((x+I*y)/r,4.))/128.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(2717./Pi)*(1. - 30.*std::pow(z/r,2.) + 85.*std::pow(z/r,4.))*std::pow((x+I*y)/r,5.))/256.;
}
else if (m == 6){
return (1.*std::sqrt(40755./Pi)*z/r*(-3. + 17.*std::pow(z/r,2.))*std::pow((x+I*y)/r,6.))/128.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(13585./Pi)*(-1. + 17.*std::pow(z/r,2.))*std::pow((x+I*y)/r,7.))/512.;
}
else if (m == 8){
return (3.*1.*std::sqrt(230945./(2.*Pi))*z/r*std::pow((x+I*y)/r,8.))/256.;
}
else if (m == 9){
return -0.001953125*(1.*std::sqrt(230945./Pi)*std::pow((x+I*y)/r,9.));
}
else{return 0.;}
}

else if (l == 10){
if(m == -10){
return (std::sqrt(969969./Pi)*std::pow((x-I*y)/r,10.))/(1024.*1.);
}
else if (m == -9){
return (std::sqrt(4849845./Pi)*z/r*std::pow((x-I*y)/r,9.))/(512.*1.);
}
else if (m == -8){
return (std::sqrt(255255./(2.*Pi))*(-1. + 19.*std::pow(z/r,2.))*std::pow((x-I*y)/r,8.))/(512.*1.);
}
else if (m == -7){
return (3.*std::sqrt(85085./Pi)*z/r*(-3. + 19.*std::pow(z/r,2.))*std::pow((x-I*y)/r,7.))/(512.*1.);
}
else if (m == -6){
return (3.*std::sqrt(5005./Pi)*(3. - 102.*std::pow(z/r,2.) + 323.*std::pow(z/r,4.))*std::pow((x-I*y)/r,6.))/(1024.*1.);
}
else if (m == -5){
return (3.*std::sqrt(1001./Pi)*z/r*(15. - 170.*std::pow(z/r,2.) + 323.*std::pow(z/r,4.))*std::pow((x-I*y)/r,5.))/(256.*1.);
}
else if (m == -4){
return (3.*std::sqrt(5005./(2.*Pi))*(-1. + 45.*std::pow(z/r,2.) - 255.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x-I*y)/r,4.))/(256.*1.);
}
else if (m == -3){
return (3.*std::sqrt(5005./Pi)*z/r*(-7. + 105.*std::pow(z/r,2.) - 357.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x-I*y)/r,3.))/(256.*1.);
}
else if (m == -2){
return (3.*std::sqrt(385./(2.*Pi))*(7. - 364.*std::pow(z/r,2.) + 2730.*std::pow(z/r,4.) - 6188.*std::pow(z/r,6.) + 4199.*std::pow(z/r,8.))*std::pow((x-I*y)/r,2.))/(512.*1.);
}
else if (m == -1){
return (std::sqrt(1155./(2.*Pi))*z/r*(63. - 1092.*std::pow(z/r,2.) + 4914.*std::pow(z/r,4.) - 7956.*std::pow(z/r,6.) + 4199.*std::pow(z/r,8.))*(x-I*y)/r)/(256.*1.);
}
else if (m == 0){
return (std::sqrt(21./Pi)*(-63. + 3465.*std::pow(z/r,2.) - 30030.*std::pow(z/r,4.) + 90090.*std::pow(z/r,6.) - 109395.*std::pow(z/r,8.) + 46189.*std::pow(z/r,10.)))/512.;
}
else if (m == 1){
return -0.00390625*(1.*std::sqrt(1155./(2.*Pi))*z/r*(63. - 1092.*std::pow(z/r,2.) + 4914.*std::pow(z/r,4.) - 7956.*std::pow(z/r,6.) + 4199.*std::pow(z/r,8.))*(x+I*y)/r);
}
else if (m == 2){
return (3.*1.*std::sqrt(385./(2.*Pi))*(7. - 364.*std::pow(z/r,2.) + 2730.*std::pow(z/r,4.) - 6188.*std::pow(z/r,6.) + 4199.*std::pow(z/r,8.))*std::pow((x+I*y)/r,2.))/512.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(5005./Pi)*z/r*(-7. + 105.*std::pow(z/r,2.) - 357.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x+I*y)/r,3.))/256.;
}
else if (m == 4){
return (3.*1.*std::sqrt(5005./(2.*Pi))*(-1. + 45.*std::pow(z/r,2.) - 255.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x+I*y)/r,4.))/256.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(1001./Pi)*z/r*(15. - 170.*std::pow(z/r,2.) + 323.*std::pow(z/r,4.))*std::pow((x+I*y)/r,5.))/256.;
}
else if (m == 6){
return (3.*1.*std::sqrt(5005./Pi)*(3. - 102.*std::pow(z/r,2.) + 323.*std::pow(z/r,4.))*std::pow((x+I*y)/r,6.))/1024.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(85085./Pi)*z/r*(-3. + 19.*std::pow(z/r,2.))*std::pow((x+I*y)/r,7.))/512.;
}
else if (m == 8){
return (1.*std::sqrt(255255./(2.*Pi))*(-1. + 19.*std::pow(z/r,2.))*std::pow((x+I*y)/r,8.))/512.;
}
else if (m == 9){
return -0.001953125*(1.*std::sqrt(4849845./Pi)*z/r*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (1.*std::sqrt(969969./Pi)*std::pow((x+I*y)/r,10.))/1024.;
}
else{return 0.;}
}

else if (l == 11){
if(m == -11){
return (std::sqrt(2028117./(2.*Pi))*std::pow((x-I*y)/r,11.))/(1024.*1.);
}
else if (m == -10){
return (std::sqrt(22309287./Pi)*z/r*std::pow((x-I*y)/r,10.))/(1024.*1.);
}
else if (m == -9){
return (std::sqrt(1062347./(2.*Pi))*(-1. + 21.*std::pow(z/r,2.))*std::pow((x-I*y)/r,9.))/(1024.*1.);
}
else if (m == -8){
return (std::sqrt(15935205./(2.*Pi))*z/r*(-1. + 7.*std::pow(z/r,2.))*std::pow((x-I*y)/r,8.))/(512.*1.);
}
else if (m == -7){
return (std::sqrt(838695./(2.*Pi))*(1. - 38.*std::pow(z/r,2.) + 133.*std::pow(z/r,4.))*std::pow((x-I*y)/r,7.))/(1024.*1.);
}
else if (m == -6){
return (std::sqrt(167739./Pi)*z/r*(15. - 190.*std::pow(z/r,2.) + 399.*std::pow(z/r,4.))*std::pow((x-I*y)/r,6.))/(1024.*1.);
}
else if (m == -5){
return (3.*std::sqrt(3289./(2.*Pi))*(-5. + 255.*std::pow(z/r,2.) - 1615.*std::pow(z/r,4.) + 2261.*std::pow(z/r,6.))*std::pow((x-I*y)/r,5.))/(1024.*1.);
}
else if (m == -4){
return (3.*std::sqrt(23023./(2.*Pi))*z/r*(-5. + 85.*std::pow(z/r,2.) - 323.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x-I*y)/r,4.))/(256.*1.);
}
else if (m == -3){
return (std::sqrt(345345./Pi)*(1. - 60.*std::pow(z/r,2.) + 510.*std::pow(z/r,4.) - 1292.*std::pow(z/r,6.) + 969.*std::pow(z/r,8.))*std::pow((x-I*y)/r,3.))/(1024.*1.);
}
else if (m == -2){
return (std::sqrt(49335./(2.*Pi))*z/r*(21. - 420.*std::pow(z/r,2.) + 2142.*std::pow(z/r,4.) - 3876.*std::pow(z/r,6.) + 2261.*std::pow(z/r,8.))*std::pow((x-I*y)/r,2.))/(512.*1.);
}
else if (m == -1){
return (std::sqrt(759./Pi)*(-21. + 1365.*std::pow(z/r,2.) - 13650.*std::pow(z/r,4.) + 46410.*std::pow(z/r,6.) - 62985.*std::pow(z/r,8.) + 29393.*std::pow(z/r,10.))*(x-I*y)/r)/(1024.*1.);
}
else if (m == 0){
return (std::sqrt(23./Pi)*(-693.*z/r + 15015.*std::pow(z/r,3.) - 90090.*std::pow(z/r,5.) + 218790.*std::pow(z/r,7.) - 230945.*std::pow(z/r,9.) + 88179.*std::pow(z/r,11.)))/512.;
}
else if (m == 1){
return -0.0009765625*(1.*std::sqrt(759./Pi)*(-21. + 1365.*std::pow(z/r,2.) - 13650.*std::pow(z/r,4.) + 46410.*std::pow(z/r,6.) - 62985.*std::pow(z/r,8.) + 29393.*std::pow(z/r,10.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(49335./(2.*Pi))*z/r*(21. - 420.*std::pow(z/r,2.) + 2142.*std::pow(z/r,4.) - 3876.*std::pow(z/r,6.) + 2261.*std::pow(z/r,8.))*std::pow((x+I*y)/r,2.))/512.;
}
else if (m == 3){
return -0.0009765625*(1.*std::sqrt(345345./Pi)*(1. - 60.*std::pow(z/r,2.) + 510.*std::pow(z/r,4.) - 1292.*std::pow(z/r,6.) + 969.*std::pow(z/r,8.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(23023./(2.*Pi))*z/r*(-5. + 85.*std::pow(z/r,2.) - 323.*std::pow(z/r,4.) + 323.*std::pow(z/r,6.))*std::pow((x+I*y)/r,4.))/256.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(3289./(2.*Pi))*(-5. + 255.*std::pow(z/r,2.) - 1615.*std::pow(z/r,4.) + 2261.*std::pow(z/r,6.))*std::pow((x+I*y)/r,5.))/1024.;
}
else if (m == 6){
return (1.*std::sqrt(167739./Pi)*z/r*(15. - 190.*std::pow(z/r,2.) + 399.*std::pow(z/r,4.))*std::pow((x+I*y)/r,6.))/1024.;
}
else if (m == 7){
return -0.0009765625*(1.*std::sqrt(838695./(2.*Pi))*(1. - 38.*std::pow(z/r,2.) + 133.*std::pow(z/r,4.))*std::pow((x+I*y)/r,7.));
}
else if (m == 8){
return (1.*std::sqrt(15935205./(2.*Pi))*z/r*(-1. + 7.*std::pow(z/r,2.))*std::pow((x+I*y)/r,8.))/512.;
}
else if (m == 9){
return -0.0009765625*(1.*std::sqrt(1062347./(2.*Pi))*(-1. + 21.*std::pow(z/r,2.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (1.*std::sqrt(22309287./Pi)*z/r*std::pow((x+I*y)/r,10.))/1024.;
}
else if (m == 11){
return -0.0009765625*(1.*std::sqrt(2028117./(2.*Pi))*std::pow((x+I*y)/r,11.));
}
else{return 0.;}
}

else if (l == 12){
if(m == -12){
return (5.*std::sqrt(676039./Pi)*std::pow((x-I*y)/r,12.))/(4096.*1.);
}
else if (m == -11){
return (5.*std::sqrt(2028117./(2.*Pi))*z/r*std::pow((x-I*y)/r,11.))/(1024.*1.);
}
else if (m == -10){
return (5.*std::sqrt(88179./Pi)*(-1. + 23.*std::pow(z/r,2.))*std::pow((x-I*y)/r,10.))/(2048.*1.);
}
else if (m == -9){
return (5.*std::sqrt(323323./(2.*Pi))*z/r*(-3. + 23.*std::pow(z/r,2.))*std::pow((x-I*y)/r,9.))/(1024.*1.);
}
else if (m == -8){
return (5.*std::sqrt(138567./(2.*Pi))*(1. - 42.*std::pow(z/r,2.) + 161.*std::pow(z/r,4.))*std::pow((x-I*y)/r,8.))/(2048.*1.);
}
else if (m == -7){
return (5.*std::sqrt(138567./(2.*Pi))*z/r*(5. - 70.*std::pow(z/r,2.) + 161.*std::pow(z/r,4.))*std::pow((x-I*y)/r,7.))/(1024.*1.);
}
else if (m == -6){
return (5.*std::sqrt(2431./Pi)*(-5. + 285.*std::pow(z/r,2.) - 1995.*std::pow(z/r,4.) + 3059.*std::pow(z/r,6.))*std::pow((x-I*y)/r,6.))/(2048.*1.);
}
else if (m == -5){
return (15.*std::sqrt(17017./(2.*Pi))*z/r*(-5. + 95.*std::pow(z/r,2.) - 399.*std::pow(z/r,4.) + 437.*std::pow(z/r,6.))*std::pow((x-I*y)/r,5.))/(1024.*1.);
}
else if (m == -4){
return (15.*std::sqrt(1001./Pi)*(5. - 340.*std::pow(z/r,2.) + 3230.*std::pow(z/r,4.) - 9044.*std::pow(z/r,6.) + 7429.*std::pow(z/r,8.))*std::pow((x-I*y)/r,4.))/(4096.*1.);
}
else if (m == -3){
return (5.*std::sqrt(1001./Pi)*z/r*(45. - 1020.*std::pow(z/r,2.) + 5814.*std::pow(z/r,4.) - 11628.*std::pow(z/r,6.) + 7429.*std::pow(z/r,8.))*std::pow((x-I*y)/r,3.))/(1024.*1.);
}
else if (m == -2){
return (5.*std::sqrt(3003./(2.*Pi))*(-3. + 225.*std::pow(z/r,2.) - 2550.*std::pow(z/r,4.) + 9690.*std::pow(z/r,6.) - 14535.*std::pow(z/r,8.) + 7429.*std::pow(z/r,10.))*std::pow((x-I*y)/r,2.))/(1024.*1.);
}
else if (m == -1){
return (5.*std::sqrt(39./Pi)*z/r*(-231. + 5775.*std::pow(z/r,2.) - 39270.*std::pow(z/r,4.) + 106590.*std::pow(z/r,6.) - 124355.*std::pow(z/r,8.) + 52003.*std::pow(z/r,10.))*(x-I*y)/r)/(1024.*1.);
}
else if (m == 0){
return (5.*(231. - 18018.*std::pow(z/r,2.) + 225225.*std::pow(z/r,4.) - 1021020.*std::pow(z/r,6.) + 2078505.*std::pow(z/r,8.) - 1939938.*std::pow(z/r,10.) + 676039.*std::pow(z/r,12.)))/(2048.*std::sqrt(Pi));
}
else if (m == 1){
return (-5.*1.*std::sqrt(39./Pi)*z/r*(-231. + 5775.*std::pow(z/r,2.) - 39270.*std::pow(z/r,4.) + 106590.*std::pow(z/r,6.) - 124355.*std::pow(z/r,8.) + 52003.*std::pow(z/r,10.))*(x+I*y)/r)/1024.;
}
else if (m == 2){
return (5.*1.*std::sqrt(3003./(2.*Pi))*(-3. + 225.*std::pow(z/r,2.) - 2550.*std::pow(z/r,4.) + 9690.*std::pow(z/r,6.) - 14535.*std::pow(z/r,8.) + 7429.*std::pow(z/r,10.))*std::pow((x+I*y)/r,2.))/1024.;
}
else if (m == 3){
return (-5.*1.*std::sqrt(1001./Pi)*z/r*(45. - 1020.*std::pow(z/r,2.) + 5814.*std::pow(z/r,4.) - 11628.*std::pow(z/r,6.) + 7429.*std::pow(z/r,8.))*std::pow((x+I*y)/r,3.))/1024.;
}
else if (m == 4){
return (15.*1.*std::sqrt(1001./Pi)*(5. - 340.*std::pow(z/r,2.) + 3230.*std::pow(z/r,4.) - 9044.*std::pow(z/r,6.) + 7429.*std::pow(z/r,8.))*std::pow((x+I*y)/r,4.))/4096.;
}
else if (m == 5){
return (-15.*1.*std::sqrt(17017./(2.*Pi))*z/r*(-5. + 95.*std::pow(z/r,2.) - 399.*std::pow(z/r,4.) + 437.*std::pow(z/r,6.))*std::pow((x+I*y)/r,5.))/1024.;
}
else if (m == 6){
return (5.*1.*std::sqrt(2431./Pi)*(-5. + 285.*std::pow(z/r,2.) - 1995.*std::pow(z/r,4.) + 3059.*std::pow(z/r,6.))*std::pow((x+I*y)/r,6.))/2048.;
}
else if (m == 7){
return (-5.*1.*std::sqrt(138567./(2.*Pi))*z/r*(5. - 70.*std::pow(z/r,2.) + 161.*std::pow(z/r,4.))*std::pow((x+I*y)/r,7.))/1024.;
}
else if (m == 8){
return (5.*1.*std::sqrt(138567./(2.*Pi))*(1. - 42.*std::pow(z/r,2.) + 161.*std::pow(z/r,4.))*std::pow((x+I*y)/r,8.))/2048.;
}
else if (m == 9){
return (-5.*1.*std::sqrt(323323./(2.*Pi))*z/r*(-3. + 23.*std::pow(z/r,2.))*std::pow((x+I*y)/r,9.))/1024.;
}
else if (m == 10){
return (5.*1.*std::sqrt(88179./Pi)*(-1. + 23.*std::pow(z/r,2.))*std::pow((x+I*y)/r,10.))/2048.;
}
else if (m == 11){
return (-5.*1.*std::sqrt(2028117./(2.*Pi))*z/r*std::pow((x+I*y)/r,11.))/1024.;
}
else if (m == 12){
return (5.*1.*std::sqrt(676039./Pi)*std::pow((x+I*y)/r,12.))/4096.;
}
else{return 0.;}
}

else if (l == 13){
if(m == -13){
return (15.*std::sqrt(156009./(2.*Pi))*std::pow((x-I*y)/r,13.))/(4096.*1.);
}
else if (m == -12){
return (15.*std::sqrt(2028117./Pi)*z/r*std::pow((x-I*y)/r,12.))/(4096.*1.);
}
else if (m == -11){
return (3.*std::sqrt(2028117./(2.*Pi))*(-1. + 25.*std::pow(z/r,2.))*std::pow((x-I*y)/r,11.))/(4096.*1.);
}
else if (m == -10){
return (3.*std::sqrt(2028117./Pi)*z/r*(-3. + 25.*std::pow(z/r,2.))*std::pow((x-I*y)/r,10.))/(2048.*1.);
}
else if (m == -9){
return (3.*std::sqrt(88179./Pi)*(3. - 138.*std::pow(z/r,2.) + 575.*std::pow(z/r,4.))*std::pow((x-I*y)/r,9.))/(4096.*1.);
}
else if (m == -8){
return (3.*std::sqrt(4849845./(2.*Pi))*z/r*(3. - 46.*std::pow(z/r,2.) + 115.*std::pow(z/r,4.))*std::pow((x-I*y)/r,8.))/(2048.*1.);
}
else if (m == -7){
return (3.*std::sqrt(692835./Pi)*(-1. + 63.*std::pow(z/r,2.) - 483.*std::pow(z/r,4.) + 805.*std::pow(z/r,6.))*std::pow((x-I*y)/r,7.))/(4096.*1.);
}
else if (m == -6){
return (3.*std::sqrt(969969./Pi)*z/r*(-5. + 105.*std::pow(z/r,2.) - 483.*std::pow(z/r,4.) + 575.*std::pow(z/r,6.))*std::pow((x-I*y)/r,6.))/(2048.*1.);
}
else if (m == -5){
return (3.*std::sqrt(51051./(2.*Pi))*(5. - 380.*std::pow(z/r,2.) + 3990.*std::pow(z/r,4.) - 12236.*std::pow(z/r,6.) + 10925.*std::pow(z/r,8.))*std::pow((x-I*y)/r,5.))/(4096.*1.);
}
else if (m == -4){
return (3.*std::sqrt(51051./Pi)*z/r*(45. - 1140.*std::pow(z/r,2.) + 7182.*std::pow(z/r,4.) - 15732.*std::pow(z/r,6.) + 10925.*std::pow(z/r,8.))*std::pow((x-I*y)/r,4.))/(4096.*1.);
}
else if (m == -3){
return (3.*std::sqrt(15015./(2.*Pi))*(-9. + 765.*std::pow(z/r,2.) - 9690.*std::pow(z/r,4.) + 40698.*std::pow(z/r,6.) - 66861.*std::pow(z/r,8.) + 37145.*std::pow(z/r,10.))*std::pow((x-I*y)/r,3.))/(4096.*1.);
}
else if (m == -2){
return (3.*std::sqrt(1365./(2.*Pi))*z/r*(-99. + 2805.*std::pow(z/r,2.) - 21318.*std::pow(z/r,4.) + 63954.*std::pow(z/r,6.) - 81719.*std::pow(z/r,8.) + 37145.*std::pow(z/r,10.))*std::pow((x-I*y)/r,2.))/(1024.*1.);
}
else if (m == -1){
return (3.*std::sqrt(273./(2.*Pi))*(33. - 2970.*std::pow(z/r,2.) + 42075.*std::pow(z/r,4.) - 213180.*std::pow(z/r,6.) + 479655.*std::pow(z/r,8.) - 490314.*std::pow(z/r,10.) + 185725.*std::pow(z/r,12.))*(x-I*y)/r)/(2048.*1.);
}
else if (m == 0){
return (3.*std::sqrt(3./Pi)*(3003.*z/r - 90090.*std::pow(z/r,3.) + 765765.*std::pow(z/r,5.) - 2771340.*std::pow(z/r,7.) + 4849845.*std::pow(z/r,9.) - 4056234.*std::pow(z/r,11.) + 1300075.*std::pow(z/r,13.)))/2048.;
}
else if (m == 1){
return (-3.*1.*std::sqrt(273./(2.*Pi))*(33. - 2970.*std::pow(z/r,2.) + 42075.*std::pow(z/r,4.) - 213180.*std::pow(z/r,6.) + 479655.*std::pow(z/r,8.) - 490314.*std::pow(z/r,10.) + 185725.*std::pow(z/r,12.))*(x+I*y)/r)/2048.;
}
else if (m == 2){
return (3.*1.*std::sqrt(1365./(2.*Pi))*z/r*(-99. + 2805.*std::pow(z/r,2.) - 21318.*std::pow(z/r,4.) + 63954.*std::pow(z/r,6.) - 81719.*std::pow(z/r,8.) + 37145.*std::pow(z/r,10.))*std::pow((x+I*y)/r,2.))/1024.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(15015./(2.*Pi))*(-9. + 765.*std::pow(z/r,2.) - 9690.*std::pow(z/r,4.) + 40698.*std::pow(z/r,6.) - 66861.*std::pow(z/r,8.) + 37145.*std::pow(z/r,10.))*std::pow((x+I*y)/r,3.))/4096.;
}
else if (m == 4){
return (3.*1.*std::sqrt(51051./Pi)*z/r*(45. - 1140.*std::pow(z/r,2.) + 7182.*std::pow(z/r,4.) - 15732.*std::pow(z/r,6.) + 10925.*std::pow(z/r,8.))*std::pow((x+I*y)/r,4.))/4096.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(51051./(2.*Pi))*(5. - 380.*std::pow(z/r,2.) + 3990.*std::pow(z/r,4.) - 12236.*std::pow(z/r,6.) + 10925.*std::pow(z/r,8.))*std::pow((x+I*y)/r,5.))/4096.;
}
else if (m == 6){
return (3.*1.*std::sqrt(969969./Pi)*z/r*(-5. + 105.*std::pow(z/r,2.) - 483.*std::pow(z/r,4.) + 575.*std::pow(z/r,6.))*std::pow((x+I*y)/r,6.))/2048.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(692835./Pi)*(-1. + 63.*std::pow(z/r,2.) - 483.*std::pow(z/r,4.) + 805.*std::pow(z/r,6.))*std::pow((x+I*y)/r,7.))/4096.;
}
else if (m == 8){
return (3.*1.*std::sqrt(4849845./(2.*Pi))*z/r*(3. - 46.*std::pow(z/r,2.) + 115.*std::pow(z/r,4.))*std::pow((x+I*y)/r,8.))/2048.;
}
else if (m == 9){
return (-3.*1.*std::sqrt(88179./Pi)*(3. - 138.*std::pow(z/r,2.) + 575.*std::pow(z/r,4.))*std::pow((x+I*y)/r,9.))/4096.;
}
else if (m == 10){
return (3.*1.*std::sqrt(2028117./Pi)*z/r*(-3. + 25.*std::pow(z/r,2.))*std::pow((x+I*y)/r,10.))/2048.;
}
else if (m == 11){
return (-3.*1.*std::sqrt(2028117./(2.*Pi))*(-1. + 25.*std::pow(z/r,2.))*std::pow((x+I*y)/r,11.))/4096.;
}
else if (m == 12){
return (15.*1.*std::sqrt(2028117./Pi)*z/r*std::pow((x+I*y)/r,12.))/4096.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(156009./(2.*Pi))*std::pow((x+I*y)/r,13.))/4096.;
}
else{return 0.;}
}

else if (l == 14){
if(m == -14){
return (15.*std::sqrt(646323./(2.*Pi))*std::pow((x-I*y)/r,14.))/(8192.*1.);
}
else if (m == -13){
return (15.*std::sqrt(4524261./(2.*Pi))*z/r*std::pow((x-I*y)/r,13.))/(4096.*1.);
}
else if (m == -12){
return (5.*std::sqrt(1508087./Pi)*(-1. + 27.*std::pow(z/r,2.))*std::pow((x-I*y)/r,12.))/(8192.*1.);
}
else if (m == -11){
return (5.*std::sqrt(58815393./(2.*Pi))*z/r*(-1. + 9.*std::pow(z/r,2.))*std::pow((x-I*y)/r,11.))/(4096.*1.);
}
else if (m == -10){
return (std::sqrt(58815393./(2.*Pi))*(1. - 50.*std::pow(z/r,2.) + 225.*std::pow(z/r,4.))*std::pow((x-I*y)/r,10.))/(8192.*1.);
}
else if (m == -9){
return (std::sqrt(98025655./Pi)*z/r*(3. - 50.*std::pow(z/r,2.) + 135.*std::pow(z/r,4.))*std::pow((x-I*y)/r,9.))/(4096.*1.);
}
else if (m == -8){
return (std::sqrt(12785955./(2.*Pi))*(-1. + 69.*std::pow(z/r,2.) - 575.*std::pow(z/r,4.) + 1035.*std::pow(z/r,6.))*std::pow((x-I*y)/r,8.))/(4096.*1.);
}
else if (m == -7){
return (std::sqrt(20092215./Pi)*z/r*(-7. + 161.*std::pow(z/r,2.) - 805.*std::pow(z/r,4.) + 1035.*std::pow(z/r,6.))*std::pow((x-I*y)/r,7.))/(4096.*1.);
}
else if (m == -6){
return (std::sqrt(46881835./(2.*Pi))*(1. - 84.*std::pow(z/r,2.) + 966.*std::pow(z/r,4.) - 3220.*std::pow(z/r,6.) + 3105.*std::pow(z/r,8.))*std::pow((x-I*y)/r,6.))/(8192.*1.);
}
else if (m == -5){
return (3.*std::sqrt(9376367./(2.*Pi))*z/r*(5. - 140.*std::pow(z/r,2.) + 966.*std::pow(z/r,4.) - 2300.*std::pow(z/r,6.) + 1725.*std::pow(z/r,8.))*std::pow((x-I*y)/r,5.))/(4096.*1.);
}
else if (m == -4){
return (3.*std::sqrt(2467465./Pi)*(-1. + 95.*std::pow(z/r,2.) - 1330.*std::pow(z/r,4.) + 6118.*std::pow(z/r,6.) - 10925.*std::pow(z/r,8.) + 6555.*std::pow(z/r,10.))*std::pow((x-I*y)/r,4.))/(8192.*1.);
}
else if (m == -3){
return (std::sqrt(224315./(2.*Pi))*z/r*(-99. + 3135.*std::pow(z/r,2.) - 26334.*std::pow(z/r,4.) + 86526.*std::pow(z/r,6.) - 120175.*std::pow(z/r,8.) + 58995.*std::pow(z/r,10.))*std::pow((x-I*y)/r,3.))/(4096.*1.);
}
else if (m == -2){
return (std::sqrt(39585./(2.*Pi))*(33. - 3366.*std::pow(z/r,2.) + 53295.*std::pow(z/r,4.) - 298452.*std::pow(z/r,6.) + 735471.*std::pow(z/r,8.) - 817190.*std::pow(z/r,10.) + 334305.*std::pow(z/r,12.))*std::pow((x-I*y)/r,2.))/(8192.*1.);
}
else if (m == -1){
return (std::sqrt(3045./(2.*Pi))*z/r*(429. - 14586.*std::pow(z/r,2.) + 138567.*std::pow(z/r,4.) - 554268.*std::pow(z/r,6.) + 1062347.*std::pow(z/r,8.) - 965770.*std::pow(z/r,10.) + 334305.*std::pow(z/r,12.))*(x-I*y)/r)/(2048.*1.);
}
else if (m == 0){
return (std::sqrt(29./Pi)*(-429. + 45045.*std::pow(z/r,2.) - 765765.*std::pow(z/r,4.) + 4849845.*std::pow(z/r,6.) - 14549535.*std::pow(z/r,8.) + 22309287.*std::pow(z/r,10.) - 16900975.*std::pow(z/r,12.) + 5014575.*std::pow(z/r,14.)))/4096.;
}
else if (m == 1){
return -0.00048828125*(1.*std::sqrt(3045./(2.*Pi))*z/r*(429. - 14586.*std::pow(z/r,2.) + 138567.*std::pow(z/r,4.) - 554268.*std::pow(z/r,6.) + 1062347.*std::pow(z/r,8.) - 965770.*std::pow(z/r,10.) + 334305.*std::pow(z/r,12.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(39585./(2.*Pi))*(33. - 3366.*std::pow(z/r,2.) + 53295.*std::pow(z/r,4.) - 298452.*std::pow(z/r,6.) + 735471.*std::pow(z/r,8.) - 817190.*std::pow(z/r,10.) + 334305.*std::pow(z/r,12.))*std::pow((x+I*y)/r,2.))/8192.;
}
else if (m == 3){
return -0.000244140625*(1.*std::sqrt(224315./(2.*Pi))*z/r*(-99. + 3135.*std::pow(z/r,2.) - 26334.*std::pow(z/r,4.) + 86526.*std::pow(z/r,6.) - 120175.*std::pow(z/r,8.) + 58995.*std::pow(z/r,10.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(2467465./Pi)*(-1. + 95.*std::pow(z/r,2.) - 1330.*std::pow(z/r,4.) + 6118.*std::pow(z/r,6.) - 10925.*std::pow(z/r,8.) + 6555.*std::pow(z/r,10.))*std::pow((x+I*y)/r,4.))/8192.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(9376367./(2.*Pi))*z/r*(5. - 140.*std::pow(z/r,2.) + 966.*std::pow(z/r,4.) - 2300.*std::pow(z/r,6.) + 1725.*std::pow(z/r,8.))*std::pow((x+I*y)/r,5.))/4096.;
}
else if (m == 6){
return (1.*std::sqrt(46881835./(2.*Pi))*(1. - 84.*std::pow(z/r,2.) + 966.*std::pow(z/r,4.) - 3220.*std::pow(z/r,6.) + 3105.*std::pow(z/r,8.))*std::pow((x+I*y)/r,6.))/8192.;
}
else if (m == 7){
return -0.000244140625*(1.*std::sqrt(20092215./Pi)*z/r*(-7. + 161.*std::pow(z/r,2.) - 805.*std::pow(z/r,4.) + 1035.*std::pow(z/r,6.))*std::pow((x+I*y)/r,7.));
}
else if (m == 8){
return (1.*std::sqrt(12785955./(2.*Pi))*(-1. + 69.*std::pow(z/r,2.) - 575.*std::pow(z/r,4.) + 1035.*std::pow(z/r,6.))*std::pow((x+I*y)/r,8.))/4096.;
}
else if (m == 9){
return -0.000244140625*(1.*std::sqrt(98025655./Pi)*z/r*(3. - 50.*std::pow(z/r,2.) + 135.*std::pow(z/r,4.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (1.*std::sqrt(58815393./(2.*Pi))*(1. - 50.*std::pow(z/r,2.) + 225.*std::pow(z/r,4.))*std::pow((x+I*y)/r,10.))/8192.;
}
else if (m == 11){
return (-5.*1.*std::sqrt(58815393./(2.*Pi))*z/r*(-1. + 9.*std::pow(z/r,2.))*std::pow((x+I*y)/r,11.))/4096.;
}
else if (m == 12){
return (5.*1.*std::sqrt(1508087./Pi)*(-1. + 27.*std::pow(z/r,2.))*std::pow((x+I*y)/r,12.))/8192.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(4524261./(2.*Pi))*z/r*std::pow((x+I*y)/r,13.))/4096.;
}
else if (m == 14){
return (15.*1.*std::sqrt(646323./(2.*Pi))*std::pow((x+I*y)/r,14.))/8192.;
}
else{return 0.;}
}

else if (l == 15){
if(m == -15){
return (3.*std::sqrt(33393355./Pi)*std::pow((x-I*y)/r,15.))/(16384.*1.);
}
else if (m == -14){
return (15.*std::sqrt(20036013./(2.*Pi))*z/r*std::pow((x-I*y)/r,14.))/(8192.*1.);
}
else if (m == -13){
return (15.*std::sqrt(690897./Pi)*(-1. + 29.*std::pow(z/r,2.))*std::pow((x-I*y)/r,13.))/(16384.*1.);
}
else if (m == -12){
return (15.*std::sqrt(1612093./Pi)*z/r*(-3. + 29.*std::pow(z/r,2.))*std::pow((x-I*y)/r,12.))/(8192.*1.);
}
else if (m == -11){
return (5.*std::sqrt(4836279./Pi)*(1. - 54.*std::pow(z/r,2.) + 261.*std::pow(z/r,4.))*std::pow((x-I*y)/r,11.))/(16384.*1.);
}
else if (m == -10){
return (std::sqrt(314358135./(2.*Pi))*z/r*(5. - 90.*std::pow(z/r,2.) + 261.*std::pow(z/r,4.))*std::pow((x-I*y)/r,10.))/(8192.*1.);
}
else if (m == -9){
return (std::sqrt(104786045./Pi)*(-1. + 75.*std::pow(z/r,2.) - 675.*std::pow(z/r,4.) + 1305.*std::pow(z/r,6.))*std::pow((x-I*y)/r,9.))/(16384.*1.);
}
else if (m == -8){
return (std::sqrt(44908305./(2.*Pi))*z/r*(-7. + 175.*std::pow(z/r,2.) - 945.*std::pow(z/r,4.) + 1305.*std::pow(z/r,6.))*std::pow((x-I*y)/r,8.))/(4096.*1.);
}
else if (m == -7){
return (std::sqrt(1952535./Pi)*(7. - 644.*std::pow(z/r,2.) + 8050.*std::pow(z/r,4.) - 28980.*std::pow(z/r,6.) + 30015.*std::pow(z/r,8.))*std::pow((x-I*y)/r,7.))/(16384.*1.);
}
else if (m == -6){
return (std::sqrt(21477885./(2.*Pi))*z/r*(21. - 644.*std::pow(z/r,2.) + 4830.*std::pow(z/r,4.) - 12420.*std::pow(z/r,6.) + 10005.*std::pow(z/r,8.))*std::pow((x-I*y)/r,6.))/(8192.*1.);
}
else if (m == -5){
return (3.*std::sqrt(10023013./Pi)*(-1. + 105.*std::pow(z/r,2.) - 1610.*std::pow(z/r,4.) + 8050.*std::pow(z/r,6.) - 15525.*std::pow(z/r,8.) + 10005.*std::pow(z/r,10.))*std::pow((x-I*y)/r,5.))/(16384.*1.);
}
else if (m == -4){
return (3.*std::sqrt(4555915./Pi)*z/r*(-11. + 385.*std::pow(z/r,2.) - 3542.*std::pow(z/r,4.) + 12650.*std::pow(z/r,6.) - 18975.*std::pow(z/r,8.) + 10005.*std::pow(z/r,10.))*std::pow((x-I*y)/r,4.))/(8192.*1.);
}
else if (m == -3){
return (std::sqrt(719355./Pi)*(11. - 1254.*std::pow(z/r,2.) + 21945.*std::pow(z/r,4.) - 134596.*std::pow(z/r,6.) + 360525.*std::pow(z/r,8.) - 432630.*std::pow(z/r,10.) + 190095.*std::pow(z/r,12.))*std::pow((x-I*y)/r,3.))/(16384.*1.);
}
else if (m == -2){
return (std::sqrt(55335./(2.*Pi))*z/r*(429. - 16302.*std::pow(z/r,2.) + 171171.*std::pow(z/r,4.) - 749892.*std::pow(z/r,6.) + 1562275.*std::pow(z/r,8.) - 1533870.*std::pow(z/r,10.) + 570285.*std::pow(z/r,12.))*std::pow((x-I*y)/r,2.))/(8192.*1.);
}
else if (m == -1){
return (std::sqrt(465./Pi)*(-429. + 51051.*std::pow(z/r,2.) - 969969.*std::pow(z/r,4.) + 6789783.*std::pow(z/r,6.) - 22309287.*std::pow(z/r,8.) + 37182145.*std::pow(z/r,10.) - 30421755.*std::pow(z/r,12.) + 9694845.*std::pow(z/r,14.))*(x-I*y)/r)/(16384.*1.);
}
else if (m == 0){
return (std::sqrt(31./Pi)*(-6435.*z/r + 255255.*std::pow(z/r,3.) - 2909907.*std::pow(z/r,5.) + 14549535.*std::pow(z/r,7.) - 37182145.*std::pow(z/r,9.) + 50702925.*std::pow(z/r,11.) - 35102025.*std::pow(z/r,13.) + 9694845.*std::pow(z/r,15.)))/4096.;
}
else if (m == 1){
return -0.00006103515625*(1.*std::sqrt(465./Pi)*(-429. + 51051.*std::pow(z/r,2.) - 969969.*std::pow(z/r,4.) + 6789783.*std::pow(z/r,6.) - 22309287.*std::pow(z/r,8.) + 37182145.*std::pow(z/r,10.) - 30421755.*std::pow(z/r,12.) + 9694845.*std::pow(z/r,14.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(55335./(2.*Pi))*z/r*(429. - 16302.*std::pow(z/r,2.) + 171171.*std::pow(z/r,4.) - 749892.*std::pow(z/r,6.) + 1562275.*std::pow(z/r,8.) - 1533870.*std::pow(z/r,10.) + 570285.*std::pow(z/r,12.))*std::pow((x+I*y)/r,2.))/8192.;
}
else if (m == 3){
return -0.00006103515625*(1.*std::sqrt(719355./Pi)*(11. - 1254.*std::pow(z/r,2.) + 21945.*std::pow(z/r,4.) - 134596.*std::pow(z/r,6.) + 360525.*std::pow(z/r,8.) - 432630.*std::pow(z/r,10.) + 190095.*std::pow(z/r,12.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(4555915./Pi)*z/r*(-11. + 385.*std::pow(z/r,2.) - 3542.*std::pow(z/r,4.) + 12650.*std::pow(z/r,6.) - 18975.*std::pow(z/r,8.) + 10005.*std::pow(z/r,10.))*std::pow((x+I*y)/r,4.))/8192.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(10023013./Pi)*(-1. + 105.*std::pow(z/r,2.) - 1610.*std::pow(z/r,4.) + 8050.*std::pow(z/r,6.) - 15525.*std::pow(z/r,8.) + 10005.*std::pow(z/r,10.))*std::pow((x+I*y)/r,5.))/16384.;
}
else if (m == 6){
return (1.*std::sqrt(21477885./(2.*Pi))*z/r*(21. - 644.*std::pow(z/r,2.) + 4830.*std::pow(z/r,4.) - 12420.*std::pow(z/r,6.) + 10005.*std::pow(z/r,8.))*std::pow((x+I*y)/r,6.))/8192.;
}
else if (m == 7){
return -0.00006103515625*(1.*std::sqrt(1952535./Pi)*(7. - 644.*std::pow(z/r,2.) + 8050.*std::pow(z/r,4.) - 28980.*std::pow(z/r,6.) + 30015.*std::pow(z/r,8.))*std::pow((x+I*y)/r,7.));
}
else if (m == 8){
return (1.*std::sqrt(44908305./(2.*Pi))*z/r*(-7. + 175.*std::pow(z/r,2.) - 945.*std::pow(z/r,4.) + 1305.*std::pow(z/r,6.))*std::pow((x+I*y)/r,8.))/4096.;
}
else if (m == 9){
return -0.00006103515625*(1.*std::sqrt(104786045./Pi)*(-1. + 75.*std::pow(z/r,2.) - 675.*std::pow(z/r,4.) + 1305.*std::pow(z/r,6.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (1.*std::sqrt(314358135./(2.*Pi))*z/r*(5. - 90.*std::pow(z/r,2.) + 261.*std::pow(z/r,4.))*std::pow((x+I*y)/r,10.))/8192.;
}
else if (m == 11){
return (-5.*1.*std::sqrt(4836279./Pi)*(1. - 54.*std::pow(z/r,2.) + 261.*std::pow(z/r,4.))*std::pow((x+I*y)/r,11.))/16384.;
}
else if (m == 12){
return (15.*1.*std::sqrt(1612093./Pi)*z/r*(-3. + 29.*std::pow(z/r,2.))*std::pow((x+I*y)/r,12.))/8192.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(690897./Pi)*(-1. + 29.*std::pow(z/r,2.))*std::pow((x+I*y)/r,13.))/16384.;
}
else if (m == 14){
return (15.*1.*std::sqrt(20036013./(2.*Pi))*z/r*std::pow((x+I*y)/r,14.))/8192.;
}
else if (m == 15){
return (-3.*1.*std::sqrt(33393355./Pi)*std::pow((x+I*y)/r,15.))/16384.;
}
else{return 0.;}
}

else if (l == 16){
if(m == -16){
return (3.*std::sqrt(1101980715./(2.*Pi))*std::pow((x-I*y)/r,16.))/(65536.*1.);
}
else if (m == -15){
return (3.*std::sqrt(1101980715./Pi)*z/r*std::pow((x-I*y)/r,15.))/(16384.*1.);
}
else if (m == -14){
return (3.*std::sqrt(35547765./(2.*Pi))*(-1. + 31.*std::pow(z/r,2.))*std::pow((x-I*y)/r,14.))/(16384.*1.);
}
else if (m == -13){
return (15.*std::sqrt(7109553./Pi)*z/r*(-3. + 31.*std::pow(z/r,2.))*std::pow((x-I*y)/r,13.))/(16384.*1.);
}
else if (m == -12){
return (15.*std::sqrt(245157./Pi)*(3. - 174.*std::pow(z/r,2.) + 899.*std::pow(z/r,4.))*std::pow((x-I*y)/r,12.))/(32768.*1.);
}
else if (m == -11){
return (3.*std::sqrt(8580495./Pi)*z/r*(15. - 290.*std::pow(z/r,2.) + 899.*std::pow(z/r,4.))*std::pow((x-I*y)/r,11.))/(16384.*1.);
}
else if (m == -10){
return (std::sqrt(8580495./(2.*Pi))*(-5. + 405.*std::pow(z/r,2.) - 3915.*std::pow(z/r,4.) + 8091.*std::pow(z/r,6.))*std::pow((x-I*y)/r,10.))/(16384.*1.);
}
else if (m == -9){
return (std::sqrt(15935205./Pi)*z/r*(-35. + 945.*std::pow(z/r,2.) - 5481.*std::pow(z/r,4.) + 8091.*std::pow(z/r,6.))*std::pow((x-I*y)/r,9.))/(16384.*1.);
}
else if (m == -8){
return (std::sqrt(15935205./(2.*Pi))*(7. - 700.*std::pow(z/r,2.) + 9450.*std::pow(z/r,4.) - 36540.*std::pow(z/r,6.) + 40455.*std::pow(z/r,8.))*std::pow((x-I*y)/r,8.))/(32768.*1.);
}
else if (m == -7){
return (3.*std::sqrt(5311735./Pi)*z/r*(21. - 700.*std::pow(z/r,2.) + 5670.*std::pow(z/r,4.) - 15660.*std::pow(z/r,6.) + 13485.*std::pow(z/r,8.))*std::pow((x-I*y)/r,7.))/(16384.*1.);
}
else if (m == -6){
return (3.*std::sqrt(46189./(2.*Pi))*(-21. + 2415.*std::pow(z/r,2.) - 40250.*std::pow(z/r,4.) + 217350.*std::pow(z/r,6.) - 450225.*std::pow(z/r,8.) + 310155.*std::pow(z/r,10.))*std::pow((x-I*y)/r,6.))/(16384.*1.);
}
else if (m == -5){
return (3.*std::sqrt(46189./Pi)*z/r*(-231. + 8855.*std::pow(z/r,2.) - 88550.*std::pow(z/r,4.) + 341550.*std::pow(z/r,6.) - 550275.*std::pow(z/r,8.) + 310155.*std::pow(z/r,10.))*std::pow((x-I*y)/r,5.))/(16384.*1.);
}
else if (m == -4){
return (3.*std::sqrt(323323./Pi)*(11. - 1386.*std::pow(z/r,2.) + 26565.*std::pow(z/r,4.) - 177100.*std::pow(z/r,6.) + 512325.*std::pow(z/r,8.) - 660330.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x-I*y)/r,4.))/(32768.*1.);
}
else if (m == -3){
return (3.*std::sqrt(124355./Pi)*z/r*(143. - 6006.*std::pow(z/r,2.) + 69069.*std::pow(z/r,4.) - 328900.*std::pow(z/r,6.) + 740025.*std::pow(z/r,8.) - 780390.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x-I*y)/r,3.))/(16384.*1.);
}
else if (m == -2){
return (3.*std::sqrt(935./(2.*Pi))*(-143. + 19019.*std::pow(z/r,2.) - 399399.*std::pow(z/r,4.) + 3062059.*std::pow(z/r,6.) - 10935925.*std::pow(z/r,8.) + 19684665.*std::pow(z/r,10.) - 17298645.*std::pow(z/r,12.) + 5892945.*std::pow(z/r,14.))*std::pow((x-I*y)/r,2.))/(16384.*1.);
}
else if (m == -1){
return (std::sqrt(561./Pi)*z/r*(-6435. + 285285.*std::pow(z/r,2.) - 3594591.*std::pow(z/r,4.) + 19684665.*std::pow(z/r,6.) - 54679625.*std::pow(z/r,8.) + 80528175.*std::pow(z/r,10.) - 59879925.*std::pow(z/r,12.) + 17678835.*std::pow(z/r,14.))*(x-I*y)/r)/(16384.*1.);
}
else if (m == 0){
return (std::sqrt(33./Pi)*(6435. - 875160.*std::pow(z/r,2.) + 19399380.*std::pow(z/r,4.) - 162954792.*std::pow(z/r,6.) + 669278610.*std::pow(z/r,8.) - 1487285800.*std::pow(z/r,10.) + 1825305300.*std::pow(z/r,12.) - 1163381400.*std::pow(z/r,14.) + 300540195.*std::pow(z/r,16.)))/65536.;
}
else if (m == 1){
return -0.00006103515625*(1.*std::sqrt(561./Pi)*z/r*(-6435. + 285285.*std::pow(z/r,2.) - 3594591.*std::pow(z/r,4.) + 19684665.*std::pow(z/r,6.) - 54679625.*std::pow(z/r,8.) + 80528175.*std::pow(z/r,10.) - 59879925.*std::pow(z/r,12.) + 17678835.*std::pow(z/r,14.))*(x+I*y)/r);
}
else if (m == 2){
return (3.*1.*std::sqrt(935./(2.*Pi))*(-143. + 19019.*std::pow(z/r,2.) - 399399.*std::pow(z/r,4.) + 3062059.*std::pow(z/r,6.) - 10935925.*std::pow(z/r,8.) + 19684665.*std::pow(z/r,10.) - 17298645.*std::pow(z/r,12.) + 5892945.*std::pow(z/r,14.))*std::pow((x+I*y)/r,2.))/16384.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(124355./Pi)*z/r*(143. - 6006.*std::pow(z/r,2.) + 69069.*std::pow(z/r,4.) - 328900.*std::pow(z/r,6.) + 740025.*std::pow(z/r,8.) - 780390.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x+I*y)/r,3.))/16384.;
}
else if (m == 4){
return (3.*1.*std::sqrt(323323./Pi)*(11. - 1386.*std::pow(z/r,2.) + 26565.*std::pow(z/r,4.) - 177100.*std::pow(z/r,6.) + 512325.*std::pow(z/r,8.) - 660330.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x+I*y)/r,4.))/32768.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(46189./Pi)*z/r*(-231. + 8855.*std::pow(z/r,2.) - 88550.*std::pow(z/r,4.) + 341550.*std::pow(z/r,6.) - 550275.*std::pow(z/r,8.) + 310155.*std::pow(z/r,10.))*std::pow((x+I*y)/r,5.))/16384.;
}
else if (m == 6){
return (3.*1.*std::sqrt(46189./(2.*Pi))*(-21. + 2415.*std::pow(z/r,2.) - 40250.*std::pow(z/r,4.) + 217350.*std::pow(z/r,6.) - 450225.*std::pow(z/r,8.) + 310155.*std::pow(z/r,10.))*std::pow((x+I*y)/r,6.))/16384.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(5311735./Pi)*z/r*(21. - 700.*std::pow(z/r,2.) + 5670.*std::pow(z/r,4.) - 15660.*std::pow(z/r,6.) + 13485.*std::pow(z/r,8.))*std::pow((x+I*y)/r,7.))/16384.;
}
else if (m == 8){
return (1.*std::sqrt(15935205./(2.*Pi))*(7. - 700.*std::pow(z/r,2.) + 9450.*std::pow(z/r,4.) - 36540.*std::pow(z/r,6.) + 40455.*std::pow(z/r,8.))*std::pow((x+I*y)/r,8.))/32768.;
}
else if (m == 9){
return -0.00006103515625*(1.*std::sqrt(15935205./Pi)*z/r*(-35. + 945.*std::pow(z/r,2.) - 5481.*std::pow(z/r,4.) + 8091.*std::pow(z/r,6.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (1.*std::sqrt(8580495./(2.*Pi))*(-5. + 405.*std::pow(z/r,2.) - 3915.*std::pow(z/r,4.) + 8091.*std::pow(z/r,6.))*std::pow((x+I*y)/r,10.))/16384.;
}
else if (m == 11){
return (-3.*1.*std::sqrt(8580495./Pi)*z/r*(15. - 290.*std::pow(z/r,2.) + 899.*std::pow(z/r,4.))*std::pow((x+I*y)/r,11.))/16384.;
}
else if (m == 12){
return (15.*1.*std::sqrt(245157./Pi)*(3. - 174.*std::pow(z/r,2.) + 899.*std::pow(z/r,4.))*std::pow((x+I*y)/r,12.))/32768.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(7109553./Pi)*z/r*(-3. + 31.*std::pow(z/r,2.))*std::pow((x+I*y)/r,13.))/16384.;
}
else if (m == 14){
return (3.*1.*std::sqrt(35547765./(2.*Pi))*(-1. + 31.*std::pow(z/r,2.))*std::pow((x+I*y)/r,14.))/16384.;
}
else if (m == 15){
return (-3.*1.*std::sqrt(1101980715./Pi)*z/r*std::pow((x+I*y)/r,15.))/16384.;
}
else if (m == 16){
return (3.*1.*std::sqrt(1101980715./(2.*Pi))*std::pow((x+I*y)/r,16.))/65536.;
}
else{return 0.;}
}

else if (l == 17){
if(m == -17){
return (15.*std::sqrt(90751353./Pi)*std::pow((x-I*y)/r,17.))/(131072.*1.);
}
else if (m == -16){
return (15.*std::sqrt(1542773001./(2.*Pi))*z/r*std::pow((x-I*y)/r,16.))/(65536.*1.);
}
else if (m == -15){
return (15.*std::sqrt(46750697./Pi)*(-1. + 33.*std::pow(z/r,2.))*std::pow((x-I*y)/r,15.))/(131072.*1.);
}
else if (m == -14){
return (15.*std::sqrt(140252091./(2.*Pi))*z/r*(-1. + 11.*std::pow(z/r,2.))*std::pow((x-I*y)/r,14.))/(16384.*1.);
}
else if (m == -13){
return (15.*std::sqrt(4524261./(2.*Pi))*(1. - 62.*std::pow(z/r,2.) + 341.*std::pow(z/r,4.))*std::pow((x-I*y)/r,13.))/(32768.*1.);
}
else if (m == -12){
return (15.*std::sqrt(1508087./Pi)*z/r*(15. - 310.*std::pow(z/r,2.) + 1023.*std::pow(z/r,4.))*std::pow((x-I*y)/r,12.))/(32768.*1.);
}
else if (m == -11){
return (15.*std::sqrt(156009./(2.*Pi))*(-5. + 435.*std::pow(z/r,2.) - 4495.*std::pow(z/r,4.) + 9889.*std::pow(z/r,6.))*std::pow((x-I*y)/r,11.))/(32768.*1.);
}
else if (m == -10){
return (15.*std::sqrt(156009./(2.*Pi))*z/r*(-35. + 1015.*std::pow(z/r,2.) - 6293.*std::pow(z/r,4.) + 9889.*std::pow(z/r,6.))*std::pow((x-I*y)/r,10.))/(16384.*1.);
}
else if (m == -9){
return (5.*std::sqrt(52003./Pi)*(35. - 3780.*std::pow(z/r,2.) + 54810.*std::pow(z/r,4.) - 226548.*std::pow(z/r,6.) + 267003.*std::pow(z/r,8.))*std::pow((x-I*y)/r,9.))/(65536.*1.);
}
else if (m == -8){
return (15.*std::sqrt(676039./(2.*Pi))*z/r*(35. - 1260.*std::pow(z/r,2.) + 10962.*std::pow(z/r,4.) - 32364.*std::pow(z/r,6.) + 29667.*std::pow(z/r,8.))*std::pow((x-I*y)/r,8.))/(32768.*1.);
}
else if (m == -7){
return (3.*std::sqrt(3380195./Pi)*(-7. + 875.*std::pow(z/r,2.) - 15750.*std::pow(z/r,4.) + 91350.*std::pow(z/r,6.) - 202275.*std::pow(z/r,8.) + 148335.*std::pow(z/r,10.))*std::pow((x-I*y)/r,7.))/(65536.*1.);
}
else if (m == -6){
return (std::sqrt(111546435./(2.*Pi))*z/r*(-21. + 875.*std::pow(z/r,2.) - 9450.*std::pow(z/r,4.) + 39150.*std::pow(z/r,6.) - 67425.*std::pow(z/r,8.) + 40455.*std::pow(z/r,10.))*std::pow((x-I*y)/r,6.))/(16384.*1.);
}
else if (m == -5){
return (3.*std::sqrt(1616615./(2.*Pi))*(7. - 966.*std::pow(z/r,2.) + 20125.*std::pow(z/r,4.) - 144900.*std::pow(z/r,6.) + 450225.*std::pow(z/r,8.) - 620310.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x-I*y)/r,5.))/(32768.*1.);
}
else if (m == -4){
return (3.*std::sqrt(11305./Pi)*z/r*(1001. - 46046.*std::pow(z/r,2.) + 575575.*std::pow(z/r,4.) - 2960100.*std::pow(z/r,6.) + 7153575.*std::pow(z/r,8.) - 8064030.*std::pow(z/r,10.) + 3411705.*std::pow(z/r,12.))*std::pow((x-I*y)/r,4.))/(32768.*1.);
}
else if (m == -3){
return (std::sqrt(33915./(2.*Pi))*(-143. + 21021.*std::pow(z/r,2.) - 483483.*std::pow(z/r,4.) + 4029025.*std::pow(z/r,6.) - 15540525.*std::pow(z/r,8.) + 30045015.*std::pow(z/r,10.) - 28224105.*std::pow(z/r,12.) + 10235115.*std::pow(z/r,14.))*std::pow((x-I*y)/r,3.))/(32768.*1.);
}
else if (m == -2){
return (3.*std::sqrt(11305./(2.*Pi))*z/r*(-715. + 35035.*std::pow(z/r,2.) - 483483.*std::pow(z/r,4.) + 2877875.*std::pow(z/r,6.) - 8633625.*std::pow(z/r,8.) + 13656825.*std::pow(z/r,10.) - 10855425.*std::pow(z/r,12.) + 3411705.*std::pow(z/r,14.))*std::pow((x-I*y)/r,2.))/(16384.*1.);
}
else if (m == -1){
return (3.*std::sqrt(595./(2.*Pi))*(715. - 108680.*std::pow(z/r,2.) + 2662660.*std::pow(z/r,4.) - 24496472.*std::pow(z/r,6.) + 109359250.*std::pow(z/r,8.) - 262462200.*std::pow(z/r,10.) + 345972900.*std::pow(z/r,12.) - 235717800.*std::pow(z/r,14.) + 64822395.*std::pow(z/r,16.))*(x-I*y)/r)/(65536.*1.);
}
else if (m == 0){
return (std::sqrt(35./Pi)*(109395.*z/r - 5542680.*std::pow(z/r,3.) + 81477396.*std::pow(z/r,5.) - 535422888.*std::pow(z/r,7.) + 1859107250.*std::pow(z/r,9.) - 3650610600.*std::pow(z/r,11.) + 4071834900.*std::pow(z/r,13.) - 2404321560.*std::pow(z/r,15.) + 583401555.*std::pow(z/r,17.)))/65536.;
}
else if (m == 1){
return (-3.*1.*std::sqrt(595./(2.*Pi))*(715. - 108680.*std::pow(z/r,2.) + 2662660.*std::pow(z/r,4.) - 24496472.*std::pow(z/r,6.) + 109359250.*std::pow(z/r,8.) - 262462200.*std::pow(z/r,10.) + 345972900.*std::pow(z/r,12.) - 235717800.*std::pow(z/r,14.) + 64822395.*std::pow(z/r,16.))*(x+I*y)/r)/65536.;
}
else if (m == 2){
return (3.*1.*std::sqrt(11305./(2.*Pi))*z/r*(-715. + 35035.*std::pow(z/r,2.) - 483483.*std::pow(z/r,4.) + 2877875.*std::pow(z/r,6.) - 8633625.*std::pow(z/r,8.) + 13656825.*std::pow(z/r,10.) - 10855425.*std::pow(z/r,12.) + 3411705.*std::pow(z/r,14.))*std::pow((x+I*y)/r,2.))/16384.;
}
else if (m == 3){
return -0.000030517578125*(1.*std::sqrt(33915./(2.*Pi))*(-143. + 21021.*std::pow(z/r,2.) - 483483.*std::pow(z/r,4.) + 4029025.*std::pow(z/r,6.) - 15540525.*std::pow(z/r,8.) + 30045015.*std::pow(z/r,10.) - 28224105.*std::pow(z/r,12.) + 10235115.*std::pow(z/r,14.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(11305./Pi)*z/r*(1001. - 46046.*std::pow(z/r,2.) + 575575.*std::pow(z/r,4.) - 2960100.*std::pow(z/r,6.) + 7153575.*std::pow(z/r,8.) - 8064030.*std::pow(z/r,10.) + 3411705.*std::pow(z/r,12.))*std::pow((x+I*y)/r,4.))/32768.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(1616615./(2.*Pi))*(7. - 966.*std::pow(z/r,2.) + 20125.*std::pow(z/r,4.) - 144900.*std::pow(z/r,6.) + 450225.*std::pow(z/r,8.) - 620310.*std::pow(z/r,10.) + 310155.*std::pow(z/r,12.))*std::pow((x+I*y)/r,5.))/32768.;
}
else if (m == 6){
return (1.*std::sqrt(111546435./(2.*Pi))*z/r*(-21. + 875.*std::pow(z/r,2.) - 9450.*std::pow(z/r,4.) + 39150.*std::pow(z/r,6.) - 67425.*std::pow(z/r,8.) + 40455.*std::pow(z/r,10.))*std::pow((x+I*y)/r,6.))/16384.;
}
else if (m == 7){
return (-3.*1.*std::sqrt(3380195./Pi)*(-7. + 875.*std::pow(z/r,2.) - 15750.*std::pow(z/r,4.) + 91350.*std::pow(z/r,6.) - 202275.*std::pow(z/r,8.) + 148335.*std::pow(z/r,10.))*std::pow((x+I*y)/r,7.))/65536.;
}
else if (m == 8){
return (15.*1.*std::sqrt(676039./(2.*Pi))*z/r*(35. - 1260.*std::pow(z/r,2.) + 10962.*std::pow(z/r,4.) - 32364.*std::pow(z/r,6.) + 29667.*std::pow(z/r,8.))*std::pow((x+I*y)/r,8.))/32768.;
}
else if (m == 9){
return (-5.*1.*std::sqrt(52003./Pi)*(35. - 3780.*std::pow(z/r,2.) + 54810.*std::pow(z/r,4.) - 226548.*std::pow(z/r,6.) + 267003.*std::pow(z/r,8.))*std::pow((x+I*y)/r,9.))/65536.;
}
else if (m == 10){
return (15.*1.*std::sqrt(156009./(2.*Pi))*z/r*(-35. + 1015.*std::pow(z/r,2.) - 6293.*std::pow(z/r,4.) + 9889.*std::pow(z/r,6.))*std::pow((x+I*y)/r,10.))/16384.;
}
else if (m == 11){
return (-15.*1.*std::sqrt(156009./(2.*Pi))*(-5. + 435.*std::pow(z/r,2.) - 4495.*std::pow(z/r,4.) + 9889.*std::pow(z/r,6.))*std::pow((x+I*y)/r,11.))/32768.;
}
else if (m == 12){
return (15.*1.*std::sqrt(1508087./Pi)*z/r*(15. - 310.*std::pow(z/r,2.) + 1023.*std::pow(z/r,4.))*std::pow((x+I*y)/r,12.))/32768.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(4524261./(2.*Pi))*(1. - 62.*std::pow(z/r,2.) + 341.*std::pow(z/r,4.))*std::pow((x+I*y)/r,13.))/32768.;
}
else if (m == 14){
return (15.*1.*std::sqrt(140252091./(2.*Pi))*z/r*(-1. + 11.*std::pow(z/r,2.))*std::pow((x+I*y)/r,14.))/16384.;
}
else if (m == 15){
return (-15.*1.*std::sqrt(46750697./Pi)*(-1. + 33.*std::pow(z/r,2.))*std::pow((x+I*y)/r,15.))/131072.;
}
else if (m == 16){
return (15.*1.*std::sqrt(1542773001./(2.*Pi))*z/r*std::pow((x+I*y)/r,16.))/65536.;
}
else if (m == 17){
return (-15.*1.*std::sqrt(90751353./Pi)*std::pow((x+I*y)/r,17.))/131072.;
}
else{return 0.;}
}

else if (l == 18){
if(m == -18){
return (5.*std::sqrt(3357800061./Pi)*std::pow((x-I*y)/r,18.))/(262144.*1.);
}
else if (m == -17){
return (15.*std::sqrt(3357800061./Pi)*z/r*std::pow((x-I*y)/r,17.))/(131072.*1.);
}
else if (m == -16){
return (3.*std::sqrt(2398428615./(2.*Pi))*(-1. + 35.*std::pow(z/r,2.))*std::pow((x-I*y)/r,16.))/(131072.*1.);
}
else if (m == -15){
return (3.*std::sqrt(13591095485./Pi)*z/r*(-3. + 35.*std::pow(z/r,2.))*std::pow((x-I*y)/r,15.))/(131072.*1.);
}
else if (m == -14){
return (3.*std::sqrt(3706662405./Pi)*(1. - 66.*std::pow(z/r,2.) + 385.*std::pow(z/r,4.))*std::pow((x-I*y)/r,14.))/(262144.*1.);
}
else if (m == -13){
return (15.*std::sqrt(741332481./(2.*Pi))*z/r*(1. - 22.*std::pow(z/r,2.) + 77.*std::pow(z/r,4.))*std::pow((x-I*y)/r,13.))/(32768.*1.);
}
else if (m == -12){
return (15.*std::sqrt(7971317./Pi)*(-1. + 93.*std::pow(z/r,2.) - 1023.*std::pow(z/r,4.) + 2387.*std::pow(z/r,6.))*std::pow((x-I*y)/r,12.))/(65536.*1.);
}
else if (m == -11){
return (3.*std::sqrt(836988285./(2.*Pi))*z/r*(-5. + 155.*std::pow(z/r,2.) - 1023.*std::pow(z/r,4.) + 1705.*std::pow(z/r,6.))*std::pow((x-I*y)/r,11.))/(32768.*1.);
}
else if (m == -10){
return (3.*std::sqrt(28861665./Pi)*(5. - 580.*std::pow(z/r,2.) + 8990.*std::pow(z/r,4.) - 39556.*std::pow(z/r,6.) + 49445.*std::pow(z/r,8.))*std::pow((x-I*y)/r,10.))/(131072.*1.);
}
else if (m == -9){
return (std::sqrt(4123095./Pi)*z/r*(315. - 12180.*std::pow(z/r,2.) + 113274.*std::pow(z/r,4.) - 356004.*std::pow(z/r,6.) + 346115.*std::pow(z/r,8.))*std::pow((x-I*y)/r,9.))/(65536.*1.);
}
else if (m == -8){
return (15.*std::sqrt(274873./(2.*Pi))*(-7. + 945.*std::pow(z/r,2.) - 18270.*std::pow(z/r,4.) + 113274.*std::pow(z/r,6.) - 267003.*std::pow(z/r,8.) + 207669.*std::pow(z/r,10.))*std::pow((x-I*y)/r,8.))/(65536.*1.);
}
else if (m == -7){
return (15.*std::sqrt(39306839./Pi)*z/r*(-7. + 315.*std::pow(z/r,2.) - 3654.*std::pow(z/r,4.) + 16182.*std::pow(z/r,6.) - 29667.*std::pow(z/r,8.) + 18879.*std::pow(z/r,10.))*std::pow((x-I*y)/r,7.))/(65536.*1.);
}
else if (m == -6){
return (std::sqrt(117920517./Pi)*(7. - 1050.*std::pow(z/r,2.) + 23625.*std::pow(z/r,4.) - 182700.*std::pow(z/r,6.) + 606825.*std::pow(z/r,8.) - 890010.*std::pow(z/r,10.) + 471975.*std::pow(z/r,12.))*std::pow((x-I*y)/r,6.))/(131072.*1.);
}
else if (m == -5){
return (3.*std::sqrt(3023603./(2.*Pi))*z/r*(91. - 4550.*std::pow(z/r,2.) + 61425.*std::pow(z/r,4.) - 339300.*std::pow(z/r,6.) + 876525.*std::pow(z/r,8.) - 1051830.*std::pow(z/r,10.) + 471975.*std::pow(z/r,12.))*std::pow((x-I*y)/r,5.))/(32768.*1.);
}
else if (m == -4){
return (3.*std::sqrt(920227./Pi)*(-13. + 2093.*std::pow(z/r,2.) - 52325.*std::pow(z/r,4.) + 470925.*std::pow(z/r,6.) - 1950975.*std::pow(z/r,8.) + 4032015.*std::pow(z/r,10.) - 4032015.*std::pow(z/r,12.) + 1550775.*std::pow(z/r,14.))*std::pow((x-I*y)/r,4.))/(65536.*1.);
}
else if (m == -3){
return (std::sqrt(1254855./(2.*Pi))*z/r*(-429. + 23023.*std::pow(z/r,2.) - 345345.*std::pow(z/r,4.) + 2220075.*std::pow(z/r,6.) - 7153575.*std::pow(z/r,8.) + 12096045.*std::pow(z/r,10.) - 10235115.*std::pow(z/r,12.) + 3411705.*std::pow(z/r,14.))*std::pow((x-I*y)/r,3.))/(32768.*1.);
}
else if (m == -2){
return (3.*std::sqrt(59755./(2.*Pi))*(143. - 24024.*std::pow(z/r,2.) + 644644.*std::pow(z/r,4.) - 6446440.*std::pow(z/r,6.) + 31081050.*std::pow(z/r,8.) - 80120040.*std::pow(z/r,10.) + 112896420.*std::pow(z/r,12.) - 81880920.*std::pow(z/r,14.) + 23881935.*std::pow(z/r,16.))*std::pow((x-I*y)/r,2.))/(131072.*1.);
}
else if (m == -1){
return (3.*std::sqrt(703./(2.*Pi))*z/r*(12155. - 680680.*std::pow(z/r,2.) + 10958948.*std::pow(z/r,4.) - 78278200.*std::pow(z/r,6.) + 293543250.*std::pow(z/r,8.) - 619109400.*std::pow(z/r,10.) + 738168900.*std::pow(z/r,12.) - 463991880.*std::pow(z/r,14.) + 119409675.*std::pow(z/r,16.))*(x-I*y)/r)/(65536.*1.);
}
else if (m == 0){
return (std::sqrt(37./Pi)*(-12155. + 2078505.*std::pow(z/r,2.) - 58198140.*std::pow(z/r,4.) + 624660036.*std::pow(z/r,6.) - 3346393050.*std::pow(z/r,8.) + 10039179150.*std::pow(z/r,10.) - 17644617900.*std::pow(z/r,12.) + 18032411700.*std::pow(z/r,14.) - 9917826435.*std::pow(z/r,16.) + 2268783825.*std::pow(z/r,18.)))/131072.;
}
else if (m == 1){
return (-3.*1.*std::sqrt(703./(2.*Pi))*z/r*(12155. - 680680.*std::pow(z/r,2.) + 10958948.*std::pow(z/r,4.) - 78278200.*std::pow(z/r,6.) + 293543250.*std::pow(z/r,8.) - 619109400.*std::pow(z/r,10.) + 738168900.*std::pow(z/r,12.) - 463991880.*std::pow(z/r,14.) + 119409675.*std::pow(z/r,16.))*(x+I*y)/r)/65536.;
}
else if (m == 2){
return (3.*1.*std::sqrt(59755./(2.*Pi))*(143. - 24024.*std::pow(z/r,2.) + 644644.*std::pow(z/r,4.) - 6446440.*std::pow(z/r,6.) + 31081050.*std::pow(z/r,8.) - 80120040.*std::pow(z/r,10.) + 112896420.*std::pow(z/r,12.) - 81880920.*std::pow(z/r,14.) + 23881935.*std::pow(z/r,16.))*std::pow((x+I*y)/r,2.))/131072.;
}
else if (m == 3){
return -0.000030517578125*(1.*std::sqrt(1254855./(2.*Pi))*z/r*(-429. + 23023.*std::pow(z/r,2.) - 345345.*std::pow(z/r,4.) + 2220075.*std::pow(z/r,6.) - 7153575.*std::pow(z/r,8.) + 12096045.*std::pow(z/r,10.) - 10235115.*std::pow(z/r,12.) + 3411705.*std::pow(z/r,14.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(920227./Pi)*(-13. + 2093.*std::pow(z/r,2.) - 52325.*std::pow(z/r,4.) + 470925.*std::pow(z/r,6.) - 1950975.*std::pow(z/r,8.) + 4032015.*std::pow(z/r,10.) - 4032015.*std::pow(z/r,12.) + 1550775.*std::pow(z/r,14.))*std::pow((x+I*y)/r,4.))/65536.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(3023603./(2.*Pi))*z/r*(91. - 4550.*std::pow(z/r,2.) + 61425.*std::pow(z/r,4.) - 339300.*std::pow(z/r,6.) + 876525.*std::pow(z/r,8.) - 1051830.*std::pow(z/r,10.) + 471975.*std::pow(z/r,12.))*std::pow((x+I*y)/r,5.))/32768.;
}
else if (m == 6){
return (1.*std::sqrt(117920517./Pi)*(7. - 1050.*std::pow(z/r,2.) + 23625.*std::pow(z/r,4.) - 182700.*std::pow(z/r,6.) + 606825.*std::pow(z/r,8.) - 890010.*std::pow(z/r,10.) + 471975.*std::pow(z/r,12.))*std::pow((x+I*y)/r,6.))/131072.;
}
else if (m == 7){
return (-15.*1.*std::sqrt(39306839./Pi)*z/r*(-7. + 315.*std::pow(z/r,2.) - 3654.*std::pow(z/r,4.) + 16182.*std::pow(z/r,6.) - 29667.*std::pow(z/r,8.) + 18879.*std::pow(z/r,10.))*std::pow((x+I*y)/r,7.))/65536.;
}
else if (m == 8){
return (15.*1.*std::sqrt(274873./(2.*Pi))*(-7. + 945.*std::pow(z/r,2.) - 18270.*std::pow(z/r,4.) + 113274.*std::pow(z/r,6.) - 267003.*std::pow(z/r,8.) + 207669.*std::pow(z/r,10.))*std::pow((x+I*y)/r,8.))/65536.;
}
else if (m == 9){
return -0.0000152587890625*(1.*std::sqrt(4123095./Pi)*z/r*(315. - 12180.*std::pow(z/r,2.) + 113274.*std::pow(z/r,4.) - 356004.*std::pow(z/r,6.) + 346115.*std::pow(z/r,8.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (3.*1.*std::sqrt(28861665./Pi)*(5. - 580.*std::pow(z/r,2.) + 8990.*std::pow(z/r,4.) - 39556.*std::pow(z/r,6.) + 49445.*std::pow(z/r,8.))*std::pow((x+I*y)/r,10.))/131072.;
}
else if (m == 11){
return (-3.*1.*std::sqrt(836988285./(2.*Pi))*z/r*(-5. + 155.*std::pow(z/r,2.) - 1023.*std::pow(z/r,4.) + 1705.*std::pow(z/r,6.))*std::pow((x+I*y)/r,11.))/32768.;
}
else if (m == 12){
return (15.*1.*std::sqrt(7971317./Pi)*(-1. + 93.*std::pow(z/r,2.) - 1023.*std::pow(z/r,4.) + 2387.*std::pow(z/r,6.))*std::pow((x+I*y)/r,12.))/65536.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(741332481./(2.*Pi))*z/r*(1. - 22.*std::pow(z/r,2.) + 77.*std::pow(z/r,4.))*std::pow((x+I*y)/r,13.))/32768.;
}
else if (m == 14){
return (3.*1.*std::sqrt(3706662405./Pi)*(1. - 66.*std::pow(z/r,2.) + 385.*std::pow(z/r,4.))*std::pow((x+I*y)/r,14.))/262144.;
}
else if (m == 15){
return (-3.*1.*std::sqrt(13591095485./Pi)*z/r*(-3. + 35.*std::pow(z/r,2.))*std::pow((x+I*y)/r,15.))/131072.;
}
else if (m == 16){
return (3.*1.*std::sqrt(2398428615./(2.*Pi))*(-1. + 35.*std::pow(z/r,2.))*std::pow((x+I*y)/r,16.))/131072.;
}
else if (m == 17){
return (-15.*1.*std::sqrt(3357800061./Pi)*z/r*std::pow((x+I*y)/r,17.))/131072.;
}
else if (m == 18){
return (5.*1.*std::sqrt(3357800061./Pi)*std::pow((x+I*y)/r,18.))/262144.;
}
else{return 0.;}
}

else if (l == 19){
if(m == -19){
return (15.*std::sqrt(765814049./(2.*Pi))*std::pow((x-I*y)/r,19.))/(262144.*1.);
}
else if (m == -18){
return (15.*std::sqrt(14550466931./Pi)*z/r*std::pow((x-I*y)/r,18.))/(262144.*1.);
}
else if (m == -17){
return (15.*std::sqrt(393255863./(2.*Pi))*(-1. + 37.*std::pow(z/r,2.))*std::pow((x-I*y)/r,17.))/(262144.*1.);
}
else if (m == -16){
return (15.*std::sqrt(1179767589./(2.*Pi))*z/r*(-3. + 37.*std::pow(z/r,2.))*std::pow((x-I*y)/r,16.))/(131072.*1.);
}
else if (m == -15){
return (3.*std::sqrt(842691135./(2.*Pi))*(3. - 210.*std::pow(z/r,2.) + 1295.*std::pow(z/r,4.))*std::pow((x-I*y)/r,15.))/(262144.*1.);
}
else if (m == -14){
return (15.*std::sqrt(2865149859./Pi)*z/r*(3. - 70.*std::pow(z/r,2.) + 259.*std::pow(z/r,4.))*std::pow((x-I*y)/r,14.))/(262144.*1.);
}
else if (m == -13){
return (15.*std::sqrt(260468169./(2.*Pi))*(-1. + 99.*std::pow(z/r,2.) - 1155.*std::pow(z/r,4.) + 2849.*std::pow(z/r,6.))*std::pow((x-I*y)/r,13.))/(262144.*1.);
}
else if (m == -12){
return (15.*std::sqrt(1823277183./Pi)*z/r*(-1. + 33.*std::pow(z/r,2.) - 231.*std::pow(z/r,4.) + 407.*std::pow(z/r,6.))*std::pow((x-I*y)/r,12.))/(65536.*1.);
}
else if (m == -11){
return (15.*std::sqrt(58815393./(2.*Pi))*(1. - 124.*std::pow(z/r,2.) + 2046.*std::pow(z/r,4.) - 9548.*std::pow(z/r,6.) + 12617.*std::pow(z/r,8.))*std::pow((x-I*y)/r,11.))/(131072.*1.);
}
else if (m == -10){
return (3.*std::sqrt(98025655./Pi)*z/r*(45. - 1860.*std::pow(z/r,2.) + 18414.*std::pow(z/r,4.) - 61380.*std::pow(z/r,6.) + 63085.*std::pow(z/r,8.))*std::pow((x-I*y)/r,10.))/(131072.*1.);
}
else if (m == -9){
return (15.*std::sqrt(676039./(2.*Pi))*(-9. + 1305.*std::pow(z/r,2.) - 26970.*std::pow(z/r,4.) + 178002.*std::pow(z/r,6.) - 445005.*std::pow(z/r,8.) + 365893.*std::pow(z/r,10.))*std::pow((x-I*y)/r,9.))/(131072.*1.);
}
else if (m == -8){
return (15.*std::sqrt(1062347./(2.*Pi))*z/r*(-63. + 3045.*std::pow(z/r,2.) - 37758.*std::pow(z/r,4.) + 178002.*std::pow(z/r,6.) - 346115.*std::pow(z/r,8.) + 232841.*std::pow(z/r,10.))*std::pow((x-I*y)/r,8.))/(65536.*1.);
}
else if (m == -7){
return (15.*std::sqrt(1062347./(2.*Pi))*(7. - 1134.*std::pow(z/r,2.) + 27405.*std::pow(z/r,4.) - 226548.*std::pow(z/r,6.) + 801009.*std::pow(z/r,8.) - 1246014.*std::pow(z/r,10.) + 698523.*std::pow(z/r,12.))*std::pow((x-I*y)/r,7.))/(131072.*1.);
}
else if (m == -6){
return (15.*std::sqrt(1062347./Pi)*z/r*(91. - 4914.*std::pow(z/r,2.) + 71253.*std::pow(z/r,4.) - 420732.*std::pow(z/r,6.) + 1157013.*std::pow(z/r,8.) - 1472562.*std::pow(z/r,10.) + 698523.*std::pow(z/r,12.))*std::pow((x-I*y)/r,6.))/(131072.*1.);
}
else if (m == -5){
return (3.*std::sqrt(7436429./(2.*Pi))*(-13. + 2275.*std::pow(z/r,2.) - 61425.*std::pow(z/r,4.) + 593775.*std::pow(z/r,6.) - 2629575.*std::pow(z/r,8.) + 5785065.*std::pow(z/r,10.) - 6135675.*std::pow(z/r,12.) + 2494725.*std::pow(z/r,14.))*std::pow((x-I*y)/r,5.))/(131072.*1.);
}
else if (m == -4){
return (3.*std::sqrt(37182145./Pi)*z/r*(-39. + 2275.*std::pow(z/r,2.) - 36855.*std::pow(z/r,4.) + 254475.*std::pow(z/r,6.) - 876525.*std::pow(z/r,8.) + 1577745.*std::pow(z/r,10.) - 1415925.*std::pow(z/r,12.) + 498945.*std::pow(z/r,14.))*std::pow((x-I*y)/r,4.))/(65536.*1.);
}
else if (m == -3){
return (3.*std::sqrt(1616615./Pi)*(39. - 7176.*std::pow(z/r,2.) + 209300.*std::pow(z/r,4.) - 2260440.*std::pow(z/r,6.) + 11705850.*std::pow(z/r,8.) - 32256120.*std::pow(z/r,10.) + 48384180.*std::pow(z/r,12.) - 37218600.*std::pow(z/r,14.) + 11475735.*std::pow(z/r,16.))*std::pow((x-I*y)/r,3.))/(262144.*1.);
}
else if (m == -2){
return (3.*std::sqrt(8645./(2.*Pi))*z/r*(7293. - 447304.*std::pow(z/r,2.) + 7827820.*std::pow(z/r,4.) - 60386040.*std::pow(z/r,6.) + 243221550.*std::pow(z/r,8.) - 548354040.*std::pow(z/r,10.) + 695987820.*std::pow(z/r,12.) - 463991880.*std::pow(z/r,14.) + 126233085.*std::pow(z/r,16.))*std::pow((x-I*y)/r,2.))/(131072.*1.);
}
else if (m == -1){
return (std::sqrt(3705./Pi)*(-2431. + 459459.*std::pow(z/r,2.) - 14090076.*std::pow(z/r,4.) + 164384220.*std::pow(z/r,6.) - 951080130.*std::pow(z/r,8.) + 3064591530.*std::pow(z/r,10.) - 5757717420.*std::pow(z/r,12.) + 6263890380.*std::pow(z/r,14.) - 3653936055.*std::pow(z/r,16.) + 883631595.*std::pow(z/r,18.))*(x-I*y)/r)/(262144.*1.);
}
else if (m == 0){
return (std::sqrt(39./Pi)*(-230945.*z/r + 14549535.*std::pow(z/r,3.) - 267711444.*std::pow(z/r,5.) + 2230928700.*std::pow(z/r,7.) - 10039179150.*std::pow(z/r,9.) + 26466926850.*std::pow(z/r,11.) - 42075627300.*std::pow(z/r,13.) + 39671305740.*std::pow(z/r,15.) - 20419054425.*std::pow(z/r,17.) + 4418157975.*std::pow(z/r,19.)))/131072.;
}
else if (m == 1){
return -3.814697265625e-6*(1.*std::sqrt(3705./Pi)*(-2431. + 459459.*std::pow(z/r,2.) - 14090076.*std::pow(z/r,4.) + 164384220.*std::pow(z/r,6.) - 951080130.*std::pow(z/r,8.) + 3064591530.*std::pow(z/r,10.) - 5757717420.*std::pow(z/r,12.) + 6263890380.*std::pow(z/r,14.) - 3653936055.*std::pow(z/r,16.) + 883631595.*std::pow(z/r,18.))*(x+I*y)/r);
}
else if (m == 2){
return (3.*1.*std::sqrt(8645./(2.*Pi))*z/r*(7293. - 447304.*std::pow(z/r,2.) + 7827820.*std::pow(z/r,4.) - 60386040.*std::pow(z/r,6.) + 243221550.*std::pow(z/r,8.) - 548354040.*std::pow(z/r,10.) + 695987820.*std::pow(z/r,12.) - 463991880.*std::pow(z/r,14.) + 126233085.*std::pow(z/r,16.))*std::pow((x+I*y)/r,2.))/131072.;
}
else if (m == 3){
return (-3.*1.*std::sqrt(1616615./Pi)*(39. - 7176.*std::pow(z/r,2.) + 209300.*std::pow(z/r,4.) - 2260440.*std::pow(z/r,6.) + 11705850.*std::pow(z/r,8.) - 32256120.*std::pow(z/r,10.) + 48384180.*std::pow(z/r,12.) - 37218600.*std::pow(z/r,14.) + 11475735.*std::pow(z/r,16.))*std::pow((x+I*y)/r,3.))/262144.;
}
else if (m == 4){
return (3.*1.*std::sqrt(37182145./Pi)*z/r*(-39. + 2275.*std::pow(z/r,2.) - 36855.*std::pow(z/r,4.) + 254475.*std::pow(z/r,6.) - 876525.*std::pow(z/r,8.) + 1577745.*std::pow(z/r,10.) - 1415925.*std::pow(z/r,12.) + 498945.*std::pow(z/r,14.))*std::pow((x+I*y)/r,4.))/65536.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(7436429./(2.*Pi))*(-13. + 2275.*std::pow(z/r,2.) - 61425.*std::pow(z/r,4.) + 593775.*std::pow(z/r,6.) - 2629575.*std::pow(z/r,8.) + 5785065.*std::pow(z/r,10.) - 6135675.*std::pow(z/r,12.) + 2494725.*std::pow(z/r,14.))*std::pow((x+I*y)/r,5.))/131072.;
}
else if (m == 6){
return (15.*1.*std::sqrt(1062347./Pi)*z/r*(91. - 4914.*std::pow(z/r,2.) + 71253.*std::pow(z/r,4.) - 420732.*std::pow(z/r,6.) + 1157013.*std::pow(z/r,8.) - 1472562.*std::pow(z/r,10.) + 698523.*std::pow(z/r,12.))*std::pow((x+I*y)/r,6.))/131072.;
}
else if (m == 7){
return (-15.*1.*std::sqrt(1062347./(2.*Pi))*(7. - 1134.*std::pow(z/r,2.) + 27405.*std::pow(z/r,4.) - 226548.*std::pow(z/r,6.) + 801009.*std::pow(z/r,8.) - 1246014.*std::pow(z/r,10.) + 698523.*std::pow(z/r,12.))*std::pow((x+I*y)/r,7.))/131072.;
}
else if (m == 8){
return (15.*1.*std::sqrt(1062347./(2.*Pi))*z/r*(-63. + 3045.*std::pow(z/r,2.) - 37758.*std::pow(z/r,4.) + 178002.*std::pow(z/r,6.) - 346115.*std::pow(z/r,8.) + 232841.*std::pow(z/r,10.))*std::pow((x+I*y)/r,8.))/65536.;
}
else if (m == 9){
return (-15.*1.*std::sqrt(676039./(2.*Pi))*(-9. + 1305.*std::pow(z/r,2.) - 26970.*std::pow(z/r,4.) + 178002.*std::pow(z/r,6.) - 445005.*std::pow(z/r,8.) + 365893.*std::pow(z/r,10.))*std::pow((x+I*y)/r,9.))/131072.;
}
else if (m == 10){
return (3.*1.*std::sqrt(98025655./Pi)*z/r*(45. - 1860.*std::pow(z/r,2.) + 18414.*std::pow(z/r,4.) - 61380.*std::pow(z/r,6.) + 63085.*std::pow(z/r,8.))*std::pow((x+I*y)/r,10.))/131072.;
}
else if (m == 11){
return (-15.*1.*std::sqrt(58815393./(2.*Pi))*(1. - 124.*std::pow(z/r,2.) + 2046.*std::pow(z/r,4.) - 9548.*std::pow(z/r,6.) + 12617.*std::pow(z/r,8.))*std::pow((x+I*y)/r,11.))/131072.;
}
else if (m == 12){
return (15.*1.*std::sqrt(1823277183./Pi)*z/r*(-1. + 33.*std::pow(z/r,2.) - 231.*std::pow(z/r,4.) + 407.*std::pow(z/r,6.))*std::pow((x+I*y)/r,12.))/65536.;
}
else if (m == 13){
return (-15.*1.*std::sqrt(260468169./(2.*Pi))*(-1. + 99.*std::pow(z/r,2.) - 1155.*std::pow(z/r,4.) + 2849.*std::pow(z/r,6.))*std::pow((x+I*y)/r,13.))/262144.;
}
else if (m == 14){
return (15.*1.*std::sqrt(2865149859./Pi)*z/r*(3. - 70.*std::pow(z/r,2.) + 259.*std::pow(z/r,4.))*std::pow((x+I*y)/r,14.))/262144.;
}
else if (m == 15){
return (-3.*1.*std::sqrt(842691135./(2.*Pi))*(3. - 210.*std::pow(z/r,2.) + 1295.*std::pow(z/r,4.))*std::pow((x+I*y)/r,15.))/262144.;
}
else if (m == 16){
return (15.*1.*std::sqrt(1179767589./(2.*Pi))*z/r*(-3. + 37.*std::pow(z/r,2.))*std::pow((x+I*y)/r,16.))/131072.;
}
else if (m == 17){
return (-15.*1.*std::sqrt(393255863./(2.*Pi))*(-1. + 37.*std::pow(z/r,2.))*std::pow((x+I*y)/r,17.))/262144.;
}
else if (m == 18){
return (15.*1.*std::sqrt(14550466931./Pi)*z/r*std::pow((x+I*y)/r,18.))/262144.;
}
else if (m == 19){
return (-15.*1.*std::sqrt(765814049./(2.*Pi))*std::pow((x+I*y)/r,19.))/262144.;
}
else{return 0.;}
}

else if (l == 20){
if(m == -20){
return (3.*std::sqrt(156991880045./Pi)*std::pow((x-I*y)/r,20.))/(1.048576e6*1.);
}
else if (m == -19){
return (15.*std::sqrt(31398376009./(2.*Pi))*z/r*std::pow((x-I*y)/r,19.))/(262144.*1.);
}
else if (m == -18){
return (5.*std::sqrt(7245779079./Pi)*(-1. + 39.*std::pow(z/r,2.))*std::pow((x-I*y)/r,18.))/(524288.*1.);
}
else if (m == -17){
return (15.*std::sqrt(45889934167./(2.*Pi))*z/r*(-1. + 13.*std::pow(z/r,2.))*std::pow((x-I*y)/r,17.))/(262144.*1.);
}
else if (m == -16){
return (15.*std::sqrt(1240268491./(2.*Pi))*(1. - 74.*std::pow(z/r,2.) + 481.*std::pow(z/r,4.))*std::pow((x-I*y)/r,16.))/(524288.*1.);
}
else if (m == -15){
return (3.*std::sqrt(6201342455./(2.*Pi))*z/r*(15. - 370.*std::pow(z/r,2.) + 1443.*std::pow(z/r,4.))*std::pow((x-I*y)/r,15.))/(262144.*1.);
}
else if (m == -14){
return (15.*std::sqrt(531543639./Pi)*(-1. + 105.*std::pow(z/r,2.) - 1295.*std::pow(z/r,4.) + 3367.*std::pow(z/r,6.))*std::pow((x-I*y)/r,14.))/(524288.*1.);
}
else if (m == -13){
return (15.*std::sqrt(63253693041./(2.*Pi))*z/r*(-1. + 35.*std::pow(z/r,2.) - 259.*std::pow(z/r,4.) + 481.*std::pow(z/r,6.))*std::pow((x-I*y)/r,13.))/(262144.*1.);
}
else if (m == -12){
return (15.*std::sqrt(1916778577./Pi)*(1. - 132.*std::pow(z/r,2.) + 2310.*std::pow(z/r,4.) - 11396.*std::pow(z/r,6.) + 15873.*std::pow(z/r,8.))*std::pow((x-I*y)/r,12.))/(1.048576e6*1.);
}
else if (m == -11){
return (15.*std::sqrt(1916778577./(2.*Pi))*z/r*(3. - 132.*std::pow(z/r,2.) + 1386.*std::pow(z/r,4.) - 4884.*std::pow(z/r,6.) + 5291.*std::pow(z/r,8.))*std::pow((x-I*y)/r,11.))/(131072.*1.);
}
else if (m == -10){
return (3.*std::sqrt(309157835./Pi)*(-3. + 465.*std::pow(z/r,2.) - 10230.*std::pow(z/r,4.) + 71610.*std::pow(z/r,6.) - 189255.*std::pow(z/r,8.) + 164021.*std::pow(z/r,10.))*std::pow((x-I*y)/r,10.))/(262144.*1.);
}
else if (m == -9){
return (5.*std::sqrt(2040441711./(2.*Pi))*z/r*(-9. + 465.*std::pow(z/r,2.) - 6138.*std::pow(z/r,4.) + 30690.*std::pow(z/r,6.) - 63085.*std::pow(z/r,8.) + 44733.*std::pow(z/r,10.))*std::pow((x-I*y)/r,9.))/(131072.*1.);
}
else if (m == -8){
return (15.*std::sqrt(23453353./(2.*Pi))*(3. - 522.*std::pow(z/r,2.) + 13485.*std::pow(z/r,4.) - 118668.*std::pow(z/r,6.) + 445005.*std::pow(z/r,8.) - 731786.*std::pow(z/r,10.) + 432419.*std::pow(z/r,12.))*std::pow((x-I*y)/r,8.))/(262144.*1.);
}
else if (m == -7){
return (15.*std::sqrt(43556227./(2.*Pi))*z/r*(21. - 1218.*std::pow(z/r,2.) + 18879.*std::pow(z/r,4.) - 118668.*std::pow(z/r,6.) + 346115.*std::pow(z/r,8.) - 465682.*std::pow(z/r,10.) + 232841.*std::pow(z/r,12.))*std::pow((x-I*y)/r,7.))/(131072.*1.);
}
else if (m == -6){
return (5.*std::sqrt(914680767./Pi)*(-1. + 189.*std::pow(z/r,2.) - 5481.*std::pow(z/r,4.) + 56637.*std::pow(z/r,6.) - 267003.*std::pow(z/r,8.) + 623007.*std::pow(z/r,10.) - 698523.*std::pow(z/r,12.) + 299367.*std::pow(z/r,14.))*std::pow((x-I*y)/r,6.))/(262144.*1.);
}
else if (m == -5){
return (3.*std::sqrt(117266765./(2.*Pi))*z/r*(-65. + 4095.*std::pow(z/r,2.) - 71253.*std::pow(z/r,4.) + 525915.*std::pow(z/r,6.) - 1928355.*std::pow(z/r,8.) + 3681405.*std::pow(z/r,10.) - 3492615.*std::pow(z/r,12.) + 1297257.*std::pow(z/r,14.))*std::pow((x-I*y)/r,5.))/(131072.*1.);
}
else if (m == -4){
return (3.*std::sqrt(117266765./(2.*Pi))*(13. - 2600.*std::pow(z/r,2.) + 81900.*std::pow(z/r,4.) - 950040.*std::pow(z/r,6.) + 5259150.*std::pow(z/r,8.) - 15426840.*std::pow(z/r,10.) + 24542700.*std::pow(z/r,12.) - 19957800.*std::pow(z/r,14.) + 6486285.*std::pow(z/r,16.))*std::pow((x-I*y)/r,4.))/(524288.*1.);
}
else if (m == -3){
return (std::sqrt(20694135./Pi)*z/r*(663. - 44200.*std::pow(z/r,2.) + 835380.*std::pow(z/r,4.) - 6921720.*std::pow(z/r,6.) + 29801850.*std::pow(z/r,8.) - 71524440.*std::pow(z/r,10.) + 96282900.*std::pow(z/r,12.) - 67856520.*std::pow(z/r,14.) + 19458855.*std::pow(z/r,16.))*std::pow((x-I*y)/r,3.))/(262144.*1.);
}
else if (m == -2){
return (std::sqrt(899745./(2.*Pi))*(-221. + 45747.*std::pow(z/r,2.) - 1524900.*std::pow(z/r,4.) + 19213740.*std::pow(z/r,6.) - 119399670.*std::pow(z/r,8.) + 411265530.*std::pow(z/r,10.) - 822531060.*std::pow(z/r,12.) + 949074300.*std::pow(z/r,14.) - 585262485.*std::pow(z/r,16.) + 149184555.*std::pow(z/r,18.))*std::pow((x-I*y)/r,2.))/(262144.*1.);
}
else if (m == -1){
return (std::sqrt(4305./Pi)*z/r*(-46189. + 3187041.*std::pow(z/r,2.) - 63740820.*std::pow(z/r,4.) + 573667380.*std::pow(z/r,6.) - 2772725670.*std::pow(z/r,8.) + 7814045070.*std::pow(z/r,10.) - 13223768580.*std::pow(z/r,12.) + 13223768580.*std::pow(z/r,14.) - 7195285845.*std::pow(z/r,16.) + 1641030105.*std::pow(z/r,18.))*(x-I*y)/r)/(262144.*1.);
}
else if (m == 0){
return (std::sqrt(41./Pi)*(46189. - 9699690.*std::pow(z/r,2.) + 334639305.*std::pow(z/r,4.) - 4461857400.*std::pow(z/r,6.) + 30117537450.*std::pow(z/r,8.) - 116454478140.*std::pow(z/r,10.) + 273491577450.*std::pow(z/r,12.) - 396713057400.*std::pow(z/r,14.) + 347123925225.*std::pow(z/r,16.) - 167890003050.*std::pow(z/r,18.) + 34461632205.*std::pow(z/r,20.)))/524288.;
}
else if (m == 1){
return -3.814697265625e-6*(1.*std::sqrt(4305./Pi)*z/r*(-46189. + 3187041.*std::pow(z/r,2.) - 63740820.*std::pow(z/r,4.) + 573667380.*std::pow(z/r,6.) - 2772725670.*std::pow(z/r,8.) + 7814045070.*std::pow(z/r,10.) - 13223768580.*std::pow(z/r,12.) + 13223768580.*std::pow(z/r,14.) - 7195285845.*std::pow(z/r,16.) + 1641030105.*std::pow(z/r,18.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(899745./(2.*Pi))*(-221. + 45747.*std::pow(z/r,2.) - 1524900.*std::pow(z/r,4.) + 19213740.*std::pow(z/r,6.) - 119399670.*std::pow(z/r,8.) + 411265530.*std::pow(z/r,10.) - 822531060.*std::pow(z/r,12.) + 949074300.*std::pow(z/r,14.) - 585262485.*std::pow(z/r,16.) + 149184555.*std::pow(z/r,18.))*std::pow((x+I*y)/r,2.))/262144.;
}
else if (m == 3){
return -3.814697265625e-6*(1.*std::sqrt(20694135./Pi)*z/r*(663. - 44200.*std::pow(z/r,2.) + 835380.*std::pow(z/r,4.) - 6921720.*std::pow(z/r,6.) + 29801850.*std::pow(z/r,8.) - 71524440.*std::pow(z/r,10.) + 96282900.*std::pow(z/r,12.) - 67856520.*std::pow(z/r,14.) + 19458855.*std::pow(z/r,16.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(117266765./(2.*Pi))*(13. - 2600.*std::pow(z/r,2.) + 81900.*std::pow(z/r,4.) - 950040.*std::pow(z/r,6.) + 5259150.*std::pow(z/r,8.) - 15426840.*std::pow(z/r,10.) + 24542700.*std::pow(z/r,12.) - 19957800.*std::pow(z/r,14.) + 6486285.*std::pow(z/r,16.))*std::pow((x+I*y)/r,4.))/524288.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(117266765./(2.*Pi))*z/r*(-65. + 4095.*std::pow(z/r,2.) - 71253.*std::pow(z/r,4.) + 525915.*std::pow(z/r,6.) - 1928355.*std::pow(z/r,8.) + 3681405.*std::pow(z/r,10.) - 3492615.*std::pow(z/r,12.) + 1297257.*std::pow(z/r,14.))*std::pow((x+I*y)/r,5.))/131072.;
}
else if (m == 6){
return (5.*1.*std::sqrt(914680767./Pi)*(-1. + 189.*std::pow(z/r,2.) - 5481.*std::pow(z/r,4.) + 56637.*std::pow(z/r,6.) - 267003.*std::pow(z/r,8.) + 623007.*std::pow(z/r,10.) - 698523.*std::pow(z/r,12.) + 299367.*std::pow(z/r,14.))*std::pow((x+I*y)/r,6.))/262144.;
}
else if (m == 7){
return (-15.*1.*std::sqrt(43556227./(2.*Pi))*z/r*(21. - 1218.*std::pow(z/r,2.) + 18879.*std::pow(z/r,4.) - 118668.*std::pow(z/r,6.) + 346115.*std::pow(z/r,8.) - 465682.*std::pow(z/r,10.) + 232841.*std::pow(z/r,12.))*std::pow((x+I*y)/r,7.))/131072.;
}
else if (m == 8){
return (15.*1.*std::sqrt(23453353./(2.*Pi))*(3. - 522.*std::pow(z/r,2.) + 13485.*std::pow(z/r,4.) - 118668.*std::pow(z/r,6.) + 445005.*std::pow(z/r,8.) - 731786.*std::pow(z/r,10.) + 432419.*std::pow(z/r,12.))*std::pow((x+I*y)/r,8.))/262144.;
}
else if (m == 9){
return (-5.*1.*std::sqrt(2040441711./(2.*Pi))*z/r*(-9. + 465.*std::pow(z/r,2.) - 6138.*std::pow(z/r,4.) + 30690.*std::pow(z/r,6.) - 63085.*std::pow(z/r,8.) + 44733.*std::pow(z/r,10.))*std::pow((x+I*y)/r,9.))/131072.;
}
else if (m == 10){
return (3.*1.*std::sqrt(309157835./Pi)*(-3. + 465.*std::pow(z/r,2.) - 10230.*std::pow(z/r,4.) + 71610.*std::pow(z/r,6.) - 189255.*std::pow(z/r,8.) + 164021.*std::pow(z/r,10.))*std::pow((x+I*y)/r,10.))/262144.;
}
else if (m == 11){
return (-15.*1.*std::sqrt(1916778577./(2.*Pi))*z/r*(3. - 132.*std::pow(z/r,2.) + 1386.*std::pow(z/r,4.) - 4884.*std::pow(z/r,6.) + 5291.*std::pow(z/r,8.))*std::pow((x+I*y)/r,11.))/131072.;
}
else if (m == 12){
return (15.*1.*std::sqrt(1916778577./Pi)*(1. - 132.*std::pow(z/r,2.) + 2310.*std::pow(z/r,4.) - 11396.*std::pow(z/r,6.) + 15873.*std::pow(z/r,8.))*std::pow((x+I*y)/r,12.))/1.048576e6;
}
else if (m == 13){
return (-15.*1.*std::sqrt(63253693041./(2.*Pi))*z/r*(-1. + 35.*std::pow(z/r,2.) - 259.*std::pow(z/r,4.) + 481.*std::pow(z/r,6.))*std::pow((x+I*y)/r,13.))/262144.;
}
else if (m == 14){
return (15.*1.*std::sqrt(531543639./Pi)*(-1. + 105.*std::pow(z/r,2.) - 1295.*std::pow(z/r,4.) + 3367.*std::pow(z/r,6.))*std::pow((x+I*y)/r,14.))/524288.;
}
else if (m == 15){
return (-3.*1.*std::sqrt(6201342455./(2.*Pi))*z/r*(15. - 370.*std::pow(z/r,2.) + 1443.*std::pow(z/r,4.))*std::pow((x+I*y)/r,15.))/262144.;
}
else if (m == 16){
return (15.*1.*std::sqrt(1240268491./(2.*Pi))*(1. - 74.*std::pow(z/r,2.) + 481.*std::pow(z/r,4.))*std::pow((x+I*y)/r,16.))/524288.;
}
else if (m == 17){
return (-15.*1.*std::sqrt(45889934167./(2.*Pi))*z/r*(-1. + 13.*std::pow(z/r,2.))*std::pow((x+I*y)/r,17.))/262144.;
}
else if (m == 18){
return (5.*1.*std::sqrt(7245779079./Pi)*(-1. + 39.*std::pow(z/r,2.))*std::pow((x+I*y)/r,18.))/524288.;
}
else if (m == 19){
return (-15.*1.*std::sqrt(31398376009./(2.*Pi))*z/r*std::pow((x+I*y)/r,19.))/262144.;
}
else if (m == 20){
return (3.*1.*std::sqrt(156991880045./Pi)*std::pow((x+I*y)/r,20.))/1.048576e6;
}
else{return 0.;}
}

else if (l == 21){
if(m == -21){
return (std::sqrt(2893136075115./(2.*Pi))*std::pow((x-I*y)/r,21.))/(1.048576e6*1.);
}
else if (m == -20){
return (3.*std::sqrt(6750650841935./Pi)*z/r*std::pow((x-I*y)/r,20.))/(1.048576e6*1.);
}
else if (m == -19){
return (3.*std::sqrt(164650020535./(2.*Pi))*(-1. + 41.*std::pow(z/r,2.))*std::pow((x-I*y)/r,19.))/(1.048576e6*1.);
}
else if (m == -18){
return (5.*std::sqrt(98790012321./Pi)*z/r*(-3. + 41.*std::pow(z/r,2.))*std::pow((x-I*y)/r,18.))/(524288.*1.);
}
else if (m == -17){
return (15.*std::sqrt(2533077239./Pi)*(1. - 78.*std::pow(z/r,2.) + 533.*std::pow(z/r,4.))*std::pow((x-I*y)/r,17.))/(1.048576e6*1.);
}
else if (m == -16){
return (3.*std::sqrt(240642337705./(2.*Pi))*z/r*(5. - 130.*std::pow(z/r,2.) + 533.*std::pow(z/r,4.))*std::pow((x-I*y)/r,16.))/(524288.*1.);
}
else if (m == -15){
return (std::sqrt(19511540895./Pi)*(-5. + 555.*std::pow(z/r,2.) - 7215.*std::pow(z/r,4.) + 19721.*std::pow(z/r,6.))*std::pow((x-I*y)/r,15.))/(1.048576e6*1.);
}
else if (m == -14){
return (3.*std::sqrt(2787362985./Pi)*z/r*(-35. + 1295.*std::pow(z/r,2.) - 10101.*std::pow(z/r,4.) + 19721.*std::pow(z/r,6.))*std::pow((x-I*y)/r,14.))/(524288.*1.);
}
else if (m == -13){
return (15.*std::sqrt(3902308179./(2.*Pi))*(1. - 140.*std::pow(z/r,2.) + 2590.*std::pow(z/r,4.) - 13468.*std::pow(z/r,6.) + 19721.*std::pow(z/r,8.))*std::pow((x-I*y)/r,13.))/(1.048576e6*1.);
}
else if (m == -12){
return (5.*std::sqrt(66339239043./Pi)*z/r*(9. - 420.*std::pow(z/r,2.) + 4662.*std::pow(z/r,4.) - 17316.*std::pow(z/r,6.) + 19721.*std::pow(z/r,8.))*std::pow((x-I*y)/r,12.))/(1.048576e6*1.);
}
else if (m == -11){
return (3.*std::sqrt(10051399855./(2.*Pi))*(-3. + 495.*std::pow(z/r,2.) - 11550.*std::pow(z/r,4.) + 85470.*std::pow(z/r,6.) - 238095.*std::pow(z/r,8.) + 216931.*std::pow(z/r,10.))*std::pow((x-I*y)/r,11.))/(1.048576e6*1.);
}
else if (m == -10){
return (3.*std::sqrt(110565398405./Pi)*z/r*(-3. + 165.*std::pow(z/r,2.) - 2310.*std::pow(z/r,4.) + 12210.*std::pow(z/r,6.) - 26455.*std::pow(z/r,8.) + 19721.*std::pow(z/r,10.))*std::pow((x-I*y)/r,10.))/(262144.*1.);
}
else if (m == -9){
return (std::sqrt(10699877265./Pi)*(3. - 558.*std::pow(z/r,2.) + 15345.*std::pow(z/r,4.) - 143220.*std::pow(z/r,6.) + 567765.*std::pow(z/r,8.) - 984126.*std::pow(z/r,10.) + 611351.*std::pow(z/r,12.))*std::pow((x-I*y)/r,9.))/(524288.*1.);
}
else if (m == -8){
return (15.*std::sqrt(9273226963./(2.*Pi))*z/r*(3. - 186.*std::pow(z/r,2.) + 3069.*std::pow(z/r,4.) - 20460.*std::pow(z/r,6.) + 63085.*std::pow(z/r,8.) - 89466.*std::pow(z/r,10.) + 47027.*std::pow(z/r,12.))*std::pow((x-I*y)/r,8.))/(262144.*1.);
}
else if (m == -7){
return (15.*std::sqrt(45680921./Pi)*(-3. + 609.*std::pow(z/r,2.) - 18879.*std::pow(z/r,4.) + 207669.*std::pow(z/r,6.) - 1038345.*std::pow(z/r,8.) + 2561251.*std::pow(z/r,10.) - 3026933.*std::pow(z/r,12.) + 1363783.*std::pow(z/r,14.))*std::pow((x-I*y)/r,7.))/(524288.*1.);
}
else if (m == -6){
return (std::sqrt(4796496705./Pi)*z/r*(-45. + 3045.*std::pow(z/r,2.) - 56637.*std::pow(z/r,4.) + 445005.*std::pow(z/r,6.) - 1730575.*std::pow(z/r,8.) + 3492615.*std::pow(z/r,10.) - 3492615.*std::pow(z/r,12.) + 1363783.*std::pow(z/r,14.))*std::pow((x-I*y)/r,6.))/(262144.*1.);
}
else if (m == -5){
return (3.*std::sqrt(1598832235./Pi)*(5. - 1080.*std::pow(z/r,2.) + 36540.*std::pow(z/r,4.) - 453096.*std::pow(z/r,6.) + 2670030.*std::pow(z/r,8.) - 8306760.*std::pow(z/r,10.) + 13970460.*std::pow(z/r,12.) - 11974680.*std::pow(z/r,14.) + 4091349.*std::pow(z/r,16.))*std::pow((x-I*y)/r,5.))/(1.048576e6*1.);
}
else if (m == -4){
return (3.*std::sqrt(7234535./(2.*Pi))*z/r*(1105. - 79560.*std::pow(z/r,2.) + 1615068.*std::pow(z/r,4.) - 14304888.*std::pow(z/r,6.) + 65564070.*std::pow(z/r,8.) - 166890360.*std::pow(z/r,10.) + 237497820.*std::pow(z/r,12.) - 176426952.*std::pow(z/r,14.) + 53187537.*std::pow(z/r,16.))*std::pow((x-I*y)/r,4.))/(524288.*1.);
}
else if (m == -3){
return (std::sqrt(7234535./Pi)*(-221. + 49725.*std::pow(z/r,2.) - 1790100.*std::pow(z/r,4.) + 24226020.*std::pow(z/r,6.) - 160929990.*std::pow(z/r,8.) + 590076630.*std::pow(z/r,10.) - 1251677700.*std::pow(z/r,12.) + 1526771700.*std::pow(z/r,14.) - 992401605.*std::pow(z/r,16.) + 265937685.*std::pow(z/r,18.))*std::pow((x-I*y)/r,3.))/(1.048576e6*1.);
}
else if (m == -2){
return (std::sqrt(1142295./(2.*Pi))*z/r*(-4199. + 314925.*std::pow(z/r,2.) - 6802380.*std::pow(z/r,4.) + 65756340.*std::pow(z/r,6.) - 339741090.*std::pow(z/r,8.) + 1019223270.*std::pow(z/r,10.) - 1829375100.*std::pow(z/r,12.) + 1933910820.*std::pow(z/r,14.) - 1109154735.*std::pow(z/r,16.) + 265937685.*std::pow(z/r,18.))*std::pow((x-I*y)/r,2.))/(262144.*1.);
}
else if (m == -1){
return (std::sqrt(9933./(2.*Pi))*(4199. - 965770.*std::pow(z/r,2.) + 36216375.*std::pow(z/r,4.) - 521515800.*std::pow(z/r,6.) + 3780989550.*std::pow(z/r,8.) - 15628090140.*std::pow(z/r,10.) + 39070225350.*std::pow(z/r,12.) - 60108039000.*std::pow(z/r,14.) + 55599936075.*std::pow(z/r,16.) - 28345065450.*std::pow(z/r,18.) + 6116566755.*std::pow(z/r,20.))*(x-I*y)/r)/(524288.*1.);
}
else if (m == 0){
return (std::sqrt(43./Pi)*(969969.*z/r - 74364290.*std::pow(z/r,3.) + 1673196525.*std::pow(z/r,5.) - 17210021400.*std::pow(z/r,7.) + 97045398450.*std::pow(z/r,9.) - 328189892940.*std::pow(z/r,11.) + 694247850450.*std::pow(z/r,13.) - 925663800600.*std::pow(z/r,15.) + 755505013725.*std::pow(z/r,17.) - 344616322050.*std::pow(z/r,19.) + 67282234305.*std::pow(z/r,21.)))/524288.;
}
else if (m == 1){
return -1.9073486328125e-6*(1.*std::sqrt(9933./(2.*Pi))*(4199. - 965770.*std::pow(z/r,2.) + 36216375.*std::pow(z/r,4.) - 521515800.*std::pow(z/r,6.) + 3780989550.*std::pow(z/r,8.) - 15628090140.*std::pow(z/r,10.) + 39070225350.*std::pow(z/r,12.) - 60108039000.*std::pow(z/r,14.) + 55599936075.*std::pow(z/r,16.) - 28345065450.*std::pow(z/r,18.) + 6116566755.*std::pow(z/r,20.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(1142295./(2.*Pi))*z/r*(-4199. + 314925.*std::pow(z/r,2.) - 6802380.*std::pow(z/r,4.) + 65756340.*std::pow(z/r,6.) - 339741090.*std::pow(z/r,8.) + 1019223270.*std::pow(z/r,10.) - 1829375100.*std::pow(z/r,12.) + 1933910820.*std::pow(z/r,14.) - 1109154735.*std::pow(z/r,16.) + 265937685.*std::pow(z/r,18.))*std::pow((x+I*y)/r,2.))/262144.;
}
else if (m == 3){
return -9.5367431640625e-7*(1.*std::sqrt(7234535./Pi)*(-221. + 49725.*std::pow(z/r,2.) - 1790100.*std::pow(z/r,4.) + 24226020.*std::pow(z/r,6.) - 160929990.*std::pow(z/r,8.) + 590076630.*std::pow(z/r,10.) - 1251677700.*std::pow(z/r,12.) + 1526771700.*std::pow(z/r,14.) - 992401605.*std::pow(z/r,16.) + 265937685.*std::pow(z/r,18.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(7234535./(2.*Pi))*z/r*(1105. - 79560.*std::pow(z/r,2.) + 1615068.*std::pow(z/r,4.) - 14304888.*std::pow(z/r,6.) + 65564070.*std::pow(z/r,8.) - 166890360.*std::pow(z/r,10.) + 237497820.*std::pow(z/r,12.) - 176426952.*std::pow(z/r,14.) + 53187537.*std::pow(z/r,16.))*std::pow((x+I*y)/r,4.))/524288.;
}
else if (m == 5){
return (-3.*1.*std::sqrt(1598832235./Pi)*(5. - 1080.*std::pow(z/r,2.) + 36540.*std::pow(z/r,4.) - 453096.*std::pow(z/r,6.) + 2670030.*std::pow(z/r,8.) - 8306760.*std::pow(z/r,10.) + 13970460.*std::pow(z/r,12.) - 11974680.*std::pow(z/r,14.) + 4091349.*std::pow(z/r,16.))*std::pow((x+I*y)/r,5.))/1.048576e6;
}
else if (m == 6){
return (1.*std::sqrt(4796496705./Pi)*z/r*(-45. + 3045.*std::pow(z/r,2.) - 56637.*std::pow(z/r,4.) + 445005.*std::pow(z/r,6.) - 1730575.*std::pow(z/r,8.) + 3492615.*std::pow(z/r,10.) - 3492615.*std::pow(z/r,12.) + 1363783.*std::pow(z/r,14.))*std::pow((x+I*y)/r,6.))/262144.;
}
else if (m == 7){
return (-15.*1.*std::sqrt(45680921./Pi)*(-3. + 609.*std::pow(z/r,2.) - 18879.*std::pow(z/r,4.) + 207669.*std::pow(z/r,6.) - 1038345.*std::pow(z/r,8.) + 2561251.*std::pow(z/r,10.) - 3026933.*std::pow(z/r,12.) + 1363783.*std::pow(z/r,14.))*std::pow((x+I*y)/r,7.))/524288.;
}
else if (m == 8){
return (15.*1.*std::sqrt(9273226963./(2.*Pi))*z/r*(3. - 186.*std::pow(z/r,2.) + 3069.*std::pow(z/r,4.) - 20460.*std::pow(z/r,6.) + 63085.*std::pow(z/r,8.) - 89466.*std::pow(z/r,10.) + 47027.*std::pow(z/r,12.))*std::pow((x+I*y)/r,8.))/262144.;
}
else if (m == 9){
return -1.9073486328125e-6*(1.*std::sqrt(10699877265./Pi)*(3. - 558.*std::pow(z/r,2.) + 15345.*std::pow(z/r,4.) - 143220.*std::pow(z/r,6.) + 567765.*std::pow(z/r,8.) - 984126.*std::pow(z/r,10.) + 611351.*std::pow(z/r,12.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (3.*1.*std::sqrt(110565398405./Pi)*z/r*(-3. + 165.*std::pow(z/r,2.) - 2310.*std::pow(z/r,4.) + 12210.*std::pow(z/r,6.) - 26455.*std::pow(z/r,8.) + 19721.*std::pow(z/r,10.))*std::pow((x+I*y)/r,10.))/262144.;
}
else if (m == 11){
return (-3.*1.*std::sqrt(10051399855./(2.*Pi))*(-3. + 495.*std::pow(z/r,2.) - 11550.*std::pow(z/r,4.) + 85470.*std::pow(z/r,6.) - 238095.*std::pow(z/r,8.) + 216931.*std::pow(z/r,10.))*std::pow((x+I*y)/r,11.))/1.048576e6;
}
else if (m == 12){
return (5.*1.*std::sqrt(66339239043./Pi)*z/r*(9. - 420.*std::pow(z/r,2.) + 4662.*std::pow(z/r,4.) - 17316.*std::pow(z/r,6.) + 19721.*std::pow(z/r,8.))*std::pow((x+I*y)/r,12.))/1.048576e6;
}
else if (m == 13){
return (-15.*1.*std::sqrt(3902308179./(2.*Pi))*(1. - 140.*std::pow(z/r,2.) + 2590.*std::pow(z/r,4.) - 13468.*std::pow(z/r,6.) + 19721.*std::pow(z/r,8.))*std::pow((x+I*y)/r,13.))/1.048576e6;
}
else if (m == 14){
return (3.*1.*std::sqrt(2787362985./Pi)*z/r*(-35. + 1295.*std::pow(z/r,2.) - 10101.*std::pow(z/r,4.) + 19721.*std::pow(z/r,6.))*std::pow((x+I*y)/r,14.))/524288.;
}
else if (m == 15){
return -9.5367431640625e-7*(1.*std::sqrt(19511540895./Pi)*(-5. + 555.*std::pow(z/r,2.) - 7215.*std::pow(z/r,4.) + 19721.*std::pow(z/r,6.))*std::pow((x+I*y)/r,15.));
}
else if (m == 16){
return (3.*1.*std::sqrt(240642337705./(2.*Pi))*z/r*(5. - 130.*std::pow(z/r,2.) + 533.*std::pow(z/r,4.))*std::pow((x+I*y)/r,16.))/524288.;
}
else if (m == 17){
return (-15.*1.*std::sqrt(2533077239./Pi)*(1. - 78.*std::pow(z/r,2.) + 533.*std::pow(z/r,4.))*std::pow((x+I*y)/r,17.))/1.048576e6;
}
else if (m == 18){
return (5.*1.*std::sqrt(98790012321./Pi)*z/r*(-3. + 41.*std::pow(z/r,2.))*std::pow((x+I*y)/r,18.))/524288.;
}
else if (m == 19){
return (-3.*1.*std::sqrt(164650020535./(2.*Pi))*(-1. + 41.*std::pow(z/r,2.))*std::pow((x+I*y)/r,19.))/1.048576e6;
}
else if (m == 20){
return (3.*1.*std::sqrt(6750650841935./Pi)*z/r*std::pow((x+I*y)/r,20.))/1.048576e6;
}
else if (m == 21){
return -9.5367431640625e-7*(1.*std::sqrt(2893136075115./(2.*Pi))*std::pow((x+I*y)/r,21.));
}
else{return 0.;}
}

else if (l == 22){
if(m == -22){
return (15.*std::sqrt(52602474093./(2.*Pi))*std::pow((x-I*y)/r,22.))/(2.097152e6*1.);
}
else if (m == -21){
return (15.*std::sqrt(578627215023./(2.*Pi))*z/r*std::pow((x-I*y)/r,21.))/(1.048576e6*1.);
}
else if (m == -20){
return (15.*std::sqrt(13456446861./Pi)*(-1. + 43.*std::pow(z/r,2.))*std::pow((x-I*y)/r,20.))/(2.097152e6*1.);
}
else if (m == -19){
return (15.*std::sqrt(94195128027./(2.*Pi))*z/r*(-3. + 43.*std::pow(z/r,2.))*std::pow((x-I*y)/r,19.))/(1.048576e6*1.);
}
else if (m == -18){
return (15.*std::sqrt(2297442147./(2.*Pi))*(3. - 246.*std::pow(z/r,2.) + 1763.*std::pow(z/r,4.))*std::pow((x-I*y)/r,18.))/(2.097152e6*1.);
}
else if (m == -17){
return (15.*std::sqrt(2297442147./Pi)*z/r*(15. - 410.*std::pow(z/r,2.) + 1763.*std::pow(z/r,4.))*std::pow((x-I*y)/r,17.))/(1.048576e6*1.);
}
else if (m == -16){
return (15.*std::sqrt(176726319./(2.*Pi))*(-5. + 585.*std::pow(z/r,2.) - 7995.*std::pow(z/r,4.) + 22919.*std::pow(z/r,6.))*std::pow((x-I*y)/r,16.))/(1.048576e6*1.);
}
else if (m == -15){
return (15.*std::sqrt(479685723./Pi)*z/r*(-35. + 1365.*std::pow(z/r,2.) - 11193.*std::pow(z/r,4.) + 22919.*std::pow(z/r,6.))*std::pow((x-I*y)/r,15.))/(1.048576e6*1.);
}
else if (m == -14){
return (15.*std::sqrt(12964479./(2.*Pi))*(35. - 5180.*std::pow(z/r,2.) + 101010.*std::pow(z/r,4.) - 552188.*std::pow(z/r,6.) + 848003.*std::pow(z/r,8.))*std::pow((x-I*y)/r,14.))/(2.097152e6*1.);
}
else if (m == -13){
return (15.*std::sqrt(12964479./(2.*Pi))*z/r*(315. - 15540.*std::pow(z/r,2.) + 181818.*std::pow(z/r,4.) - 709956.*std::pow(z/r,6.) + 848003.*std::pow(z/r,8.))*std::pow((x-I*y)/r,13.))/(1.048576e6*1.);
}
else if (m == -12){
return (15.*std::sqrt(90751353./Pi)*(-9. + 1575.*std::pow(z/r,2.) - 38850.*std::pow(z/r,4.) + 303030.*std::pow(z/r,6.) - 887445.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x-I*y)/r,12.))/(2.097152e6*1.);
}
else if (m == -11){
return (15.*std::sqrt(140252091./(2.*Pi))*z/r*(-99. + 5775.*std::pow(z/r,2.) - 85470.*std::pow(z/r,4.) + 476190.*std::pow(z/r,6.) - 1084655.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x-I*y)/r,11.))/(1.048576e6*1.);
}
else if (m == -10){
return (15.*std::sqrt(1542773001./(2.*Pi))*(3. - 594.*std::pow(z/r,2.) + 17325.*std::pow(z/r,4.) - 170940.*std::pow(z/r,6.) + 714285.*std::pow(z/r,8.) - 1301586.*std::pow(z/r,10.) + 848003.*std::pow(z/r,12.))*std::pow((x-I*y)/r,10.))/(2.097152e6*1.);
}
else if (m == -9){
return (15.*std::sqrt(20056049013./Pi)*z/r*(3. - 198.*std::pow(z/r,2.) + 3465.*std::pow(z/r,4.) - 24420.*std::pow(z/r,6.) + 79365.*std::pow(z/r,8.) - 118326.*std::pow(z/r,10.) + 65231.*std::pow(z/r,12.))*std::pow((x-I*y)/r,9.))/(524288.*1.);
}
else if (m == -8){
return (15.*std::sqrt(92424189./(2.*Pi))*(-3. + 651.*std::pow(z/r,2.) - 21483.*std::pow(z/r,4.) + 250635.*std::pow(z/r,6.) - 1324785.*std::pow(z/r,8.) + 3444441.*std::pow(z/r,10.) - 4279457.*std::pow(z/r,12.) + 2022161.*std::pow(z/r,14.))*std::pow((x-I*y)/r,8.))/(524288.*1.);
}
else if (m == -7){
return (15.*std::sqrt(92424189./Pi)*z/r*(-45. + 3255.*std::pow(z/r,2.) - 64449.*std::pow(z/r,4.) + 537075.*std::pow(z/r,6.) - 2207975.*std::pow(z/r,8.) + 4696965.*std::pow(z/r,10.) - 4937835.*std::pow(z/r,12.) + 2022161.*std::pow(z/r,14.))*std::pow((x-I*y)/r,7.))/(524288.*1.);
}
else if (m == -6){
return (15.*std::sqrt(3187041./Pi)*(45. - 10440.*std::pow(z/r,2.) + 377580.*std::pow(z/r,4.) - 4984056.*std::pow(z/r,6.) + 31150350.*std::pow(z/r,8.) - 102450040.*std::pow(z/r,10.) + 181615980.*std::pow(z/r,12.) - 163653960.*std::pow(z/r,14.) + 58642669.*std::pow(z/r,16.))*std::pow((x-I*y)/r,6.))/(2.097152e6*1.);
}
else if (m == -5){
return (15.*std::sqrt(1312311./Pi)*z/r*(765. - 59160.*std::pow(z/r,2.) + 1283772.*std::pow(z/r,4.) - 12104136.*std::pow(z/r,6.) + 58839550.*std::pow(z/r,8.) - 158331880.*std::pow(z/r,10.) + 237497820.*std::pow(z/r,12.) - 185474488.*std::pow(z/r,14.) + 58642669.*std::pow(z/r,16.))*std::pow((x-I*y)/r,5.))/(1.048576e6*1.);
}
else if (m == -4){
return (15.*std::sqrt(437437./(2.*Pi))*(-85. + 20655.*std::pow(z/r,2.) - 798660.*std::pow(z/r,4.) + 11553948.*std::pow(z/r,6.) - 81702918.*std::pow(z/r,8.) + 317733570.*std::pow(z/r,10.) - 712493460.*std::pow(z/r,12.) + 916063020.*std::pow(z/r,14.) - 625976397.*std::pow(z/r,16.) + 175928007.*std::pow(z/r,18.))*std::pow((x-I*y)/r,4.))/(1.048576e6*1.);
}
else if (m == -3){
return (15.*std::sqrt(1771./Pi)*z/r*(-20995. + 1700595.*std::pow(z/r,2.) - 39453804.*std::pow(z/r,4.) + 407689308.*std::pow(z/r,6.) - 2242291194.*std::pow(z/r,8.) + 7134562890.*std::pow(z/r,10.) - 13537375740.*std::pow(z/r,12.) + 15084504396.*std::pow(z/r,14.) - 9095068827.*std::pow(z/r,16.) + 2287064091.*std::pow(z/r,18.))*std::pow((x-I*y)/r,3.))/(1.048576e6*1.);
}
else if (m == -2){
return (3.*std::sqrt(8855./Pi)*(4199. - 1049750.*std::pow(z/r,2.) + 42514875.*std::pow(z/r,4.) - 657563400.*std::pow(z/r,6.) + 5096116350.*std::pow(z/r,8.) - 22422911940.*std::pow(z/r,10.) + 59454690750.*std::pow(z/r,12.) - 96695541000.*std::pow(z/r,14.) + 94278152475.*std::pow(z/r,16.) - 50528160150.*std::pow(z/r,18.) + 11435320455.*std::pow(z/r,20.))*std::pow((x-I*y)/r,2.))/(2.097152e6*1.);
}
else if (m == -1){
return (3.*std::sqrt(1265./(2.*Pi))*(88179.*z/r - 7348250.*std::pow(z/r,3.) + 178562475.*std::pow(z/r,5.) - 1972690200.*std::pow(z/r,7.) + 11890938150.*std::pow(z/r,9.) - 42807377340.*std::pow(z/r,11.) + 96042192750.*std::pow(z/r,13.) - 135373757400.*std::pow(z/r,15.) + 116461247175.*std::pow(z/r,17.) - 55846913850.*std::pow(z/r,19.) + 11435320455.*std::pow(z/r,21.))*(x-I*y)/r)/(524288.*1.);
}
else if (m == 0){
return (3.*std::sqrt(5./Pi)*(-88179. + 22309287.*std::pow(z/r,2.) - 929553625.*std::pow(z/r,4.) + 15058768725.*std::pow(z/r,6.) - 124772655150.*std::pow(z/r,8.) + 601681470390.*std::pow(z/r,10.) - 1805044411170.*std::pow(z/r,12.) + 3471239252250.*std::pow(z/r,14.) - 4281195077775.*std::pow(z/r,16.) + 3273855059475.*std::pow(z/r,18.) - 1412926920405.*std::pow(z/r,20.) + 263012370465.*std::pow(z/r,22.)))/1.048576e6;
}
else if (m == 1){
return (-3.*1.*std::sqrt(1265./(2.*Pi))*(88179.*z/r - 7348250.*std::pow(z/r,3.) + 178562475.*std::pow(z/r,5.) - 1972690200.*std::pow(z/r,7.) + 11890938150.*std::pow(z/r,9.) - 42807377340.*std::pow(z/r,11.) + 96042192750.*std::pow(z/r,13.) - 135373757400.*std::pow(z/r,15.) + 116461247175.*std::pow(z/r,17.) - 55846913850.*std::pow(z/r,19.) + 11435320455.*std::pow(z/r,21.))*(x+I*y)/r)/524288.;
}
else if (m == 2){
return (3.*1.*std::sqrt(8855./Pi)*(4199. - 1049750.*std::pow(z/r,2.) + 42514875.*std::pow(z/r,4.) - 657563400.*std::pow(z/r,6.) + 5096116350.*std::pow(z/r,8.) - 22422911940.*std::pow(z/r,10.) + 59454690750.*std::pow(z/r,12.) - 96695541000.*std::pow(z/r,14.) + 94278152475.*std::pow(z/r,16.) - 50528160150.*std::pow(z/r,18.) + 11435320455.*std::pow(z/r,20.))*std::pow((x+I*y)/r,2.))/2.097152e6;
}
else if (m == 3){
return (-15.*1.*std::sqrt(1771./Pi)*z/r*(-20995. + 1700595.*std::pow(z/r,2.) - 39453804.*std::pow(z/r,4.) + 407689308.*std::pow(z/r,6.) - 2242291194.*std::pow(z/r,8.) + 7134562890.*std::pow(z/r,10.) - 13537375740.*std::pow(z/r,12.) + 15084504396.*std::pow(z/r,14.) - 9095068827.*std::pow(z/r,16.) + 2287064091.*std::pow(z/r,18.))*std::pow((x+I*y)/r,3.))/1.048576e6;
}
else if (m == 4){
return (15.*1.*std::sqrt(437437./(2.*Pi))*(-85. + 20655.*std::pow(z/r,2.) - 798660.*std::pow(z/r,4.) + 11553948.*std::pow(z/r,6.) - 81702918.*std::pow(z/r,8.) + 317733570.*std::pow(z/r,10.) - 712493460.*std::pow(z/r,12.) + 916063020.*std::pow(z/r,14.) - 625976397.*std::pow(z/r,16.) + 175928007.*std::pow(z/r,18.))*std::pow((x+I*y)/r,4.))/1.048576e6;
}
else if (m == 5){
return (-15.*1.*std::sqrt(1312311./Pi)*z/r*(765. - 59160.*std::pow(z/r,2.) + 1283772.*std::pow(z/r,4.) - 12104136.*std::pow(z/r,6.) + 58839550.*std::pow(z/r,8.) - 158331880.*std::pow(z/r,10.) + 237497820.*std::pow(z/r,12.) - 185474488.*std::pow(z/r,14.) + 58642669.*std::pow(z/r,16.))*std::pow((x+I*y)/r,5.))/1.048576e6;
}
else if (m == 6){
return (15.*1.*std::sqrt(3187041./Pi)*(45. - 10440.*std::pow(z/r,2.) + 377580.*std::pow(z/r,4.) - 4984056.*std::pow(z/r,6.) + 31150350.*std::pow(z/r,8.) - 102450040.*std::pow(z/r,10.) + 181615980.*std::pow(z/r,12.) - 163653960.*std::pow(z/r,14.) + 58642669.*std::pow(z/r,16.))*std::pow((x+I*y)/r,6.))/2.097152e6;
}
else if (m == 7){
return (-15.*1.*std::sqrt(92424189./Pi)*z/r*(-45. + 3255.*std::pow(z/r,2.) - 64449.*std::pow(z/r,4.) + 537075.*std::pow(z/r,6.) - 2207975.*std::pow(z/r,8.) + 4696965.*std::pow(z/r,10.) - 4937835.*std::pow(z/r,12.) + 2022161.*std::pow(z/r,14.))*std::pow((x+I*y)/r,7.))/524288.;
}
else if (m == 8){
return (15.*1.*std::sqrt(92424189./(2.*Pi))*(-3. + 651.*std::pow(z/r,2.) - 21483.*std::pow(z/r,4.) + 250635.*std::pow(z/r,6.) - 1324785.*std::pow(z/r,8.) + 3444441.*std::pow(z/r,10.) - 4279457.*std::pow(z/r,12.) + 2022161.*std::pow(z/r,14.))*std::pow((x+I*y)/r,8.))/524288.;
}
else if (m == 9){
return (-15.*1.*std::sqrt(20056049013./Pi)*z/r*(3. - 198.*std::pow(z/r,2.) + 3465.*std::pow(z/r,4.) - 24420.*std::pow(z/r,6.) + 79365.*std::pow(z/r,8.) - 118326.*std::pow(z/r,10.) + 65231.*std::pow(z/r,12.))*std::pow((x+I*y)/r,9.))/524288.;
}
else if (m == 10){
return (15.*1.*std::sqrt(1542773001./(2.*Pi))*(3. - 594.*std::pow(z/r,2.) + 17325.*std::pow(z/r,4.) - 170940.*std::pow(z/r,6.) + 714285.*std::pow(z/r,8.) - 1301586.*std::pow(z/r,10.) + 848003.*std::pow(z/r,12.))*std::pow((x+I*y)/r,10.))/2.097152e6;
}
else if (m == 11){
return (-15.*1.*std::sqrt(140252091./(2.*Pi))*z/r*(-99. + 5775.*std::pow(z/r,2.) - 85470.*std::pow(z/r,4.) + 476190.*std::pow(z/r,6.) - 1084655.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x+I*y)/r,11.))/1.048576e6;
}
else if (m == 12){
return (15.*1.*std::sqrt(90751353./Pi)*(-9. + 1575.*std::pow(z/r,2.) - 38850.*std::pow(z/r,4.) + 303030.*std::pow(z/r,6.) - 887445.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x+I*y)/r,12.))/2.097152e6;
}
else if (m == 13){
return (-15.*1.*std::sqrt(12964479./(2.*Pi))*z/r*(315. - 15540.*std::pow(z/r,2.) + 181818.*std::pow(z/r,4.) - 709956.*std::pow(z/r,6.) + 848003.*std::pow(z/r,8.))*std::pow((x+I*y)/r,13.))/1.048576e6;
}
else if (m == 14){
return (15.*1.*std::sqrt(12964479./(2.*Pi))*(35. - 5180.*std::pow(z/r,2.) + 101010.*std::pow(z/r,4.) - 552188.*std::pow(z/r,6.) + 848003.*std::pow(z/r,8.))*std::pow((x+I*y)/r,14.))/2.097152e6;
}
else if (m == 15){
return (-15.*1.*std::sqrt(479685723./Pi)*z/r*(-35. + 1365.*std::pow(z/r,2.) - 11193.*std::pow(z/r,4.) + 22919.*std::pow(z/r,6.))*std::pow((x+I*y)/r,15.))/1.048576e6;
}
else if (m == 16){
return (15.*1.*std::sqrt(176726319./(2.*Pi))*(-5. + 585.*std::pow(z/r,2.) - 7995.*std::pow(z/r,4.) + 22919.*std::pow(z/r,6.))*std::pow((x+I*y)/r,16.))/1.048576e6;
}
else if (m == 17){
return (-15.*1.*std::sqrt(2297442147./Pi)*z/r*(15. - 410.*std::pow(z/r,2.) + 1763.*std::pow(z/r,4.))*std::pow((x+I*y)/r,17.))/1.048576e6;
}
else if (m == 18){
return (15.*1.*std::sqrt(2297442147./(2.*Pi))*(3. - 246.*std::pow(z/r,2.) + 1763.*std::pow(z/r,4.))*std::pow((x+I*y)/r,18.))/2.097152e6;
}
else if (m == 19){
return (-15.*1.*std::sqrt(94195128027./(2.*Pi))*z/r*(-3. + 43.*std::pow(z/r,2.))*std::pow((x+I*y)/r,19.))/1.048576e6;
}
else if (m == 20){
return (15.*1.*std::sqrt(13456446861./Pi)*(-1. + 43.*std::pow(z/r,2.))*std::pow((x+I*y)/r,20.))/2.097152e6;
}
else if (m == 21){
return (-15.*1.*std::sqrt(578627215023./(2.*Pi))*z/r*std::pow((x+I*y)/r,21.))/1.048576e6;
}
else if (m == 22){
return (15.*1.*std::sqrt(52602474093./(2.*Pi))*std::pow((x+I*y)/r,22.))/2.097152e6;
}
else{return 0.;}
}

else if (l == 23){
if(m == -23){
return (15.*std::sqrt(107492012277./Pi)*std::pow((x-I*y)/r,23.))/(4.194304e6*1.);
}
else if (m == -22){
return (15.*std::sqrt(2472316282371./(2.*Pi))*z/r*std::pow((x-I*y)/r,22.))/(2.097152e6*1.);
}
else if (m == -21){
return (std::sqrt(12361581411855./Pi)*(-1. + 45.*std::pow(z/r,2.))*std::pow((x-I*y)/r,21.))/(4.194304e6*1.);
}
else if (m == -20){
return (3.*std::sqrt(45325798510135./Pi)*z/r*(-1. + 15.*std::pow(z/r,2.))*std::pow((x-I*y)/r,20.))/(2.097152e6*1.);
}
else if (m == -19){
return (3.*std::sqrt(1054088337445./Pi)*(1. - 86.*std::pow(z/r,2.) + 645.*std::pow(z/r,4.))*std::pow((x-I*y)/r,19.))/(4.194304e6*1.);
}
else if (m == -18){
return (5.*std::sqrt(4427171017269./(2.*Pi))*z/r*(3. - 86.*std::pow(z/r,2.) + 387.*std::pow(z/r,4.))*std::pow((x-I*y)/r,18.))/(2.097152e6*1.);
}
else if (m == -17){
return (15.*std::sqrt(35993260303./Pi)*(-1. + 123.*std::pow(z/r,2.) - 1763.*std::pow(z/r,4.) + 5289.*std::pow(z/r,6.))*std::pow((x-I*y)/r,17.))/(4.194304e6*1.);
}
else if (m == -16){
return (3.*std::sqrt(25709471645./(2.*Pi))*z/r*(-35. + 1435.*std::pow(z/r,2.) - 12341.*std::pow(z/r,4.) + 26445.*std::pow(z/r,6.))*std::pow((x-I*y)/r,16.))/(1.048576e6*1.);
}
else if (m == -15){
return (std::sqrt(5932954995./Pi)*(35. - 5460.*std::pow(z/r,2.) + 111930.*std::pow(z/r,4.) - 641732.*std::pow(z/r,6.) + 1031355.*std::pow(z/r,8.))*std::pow((x-I*y)/r,15.))/(4.194304e6*1.);
}
else if (m == -14){
return (3.*std::sqrt(112726144905./(2.*Pi))*z/r*(35. - 1820.*std::pow(z/r,2.) + 22386.*std::pow(z/r,4.) - 91676.*std::pow(z/r,6.) + 114595.*std::pow(z/r,8.))*std::pow((x-I*y)/r,14.))/(2.097152e6*1.);
}
else if (m == -13){
return (15.*std::sqrt(609330513./Pi)*(-7. + 1295.*std::pow(z/r,2.) - 33670.*std::pow(z/r,4.) + 276094.*std::pow(z/r,6.) - 848003.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x-I*y)/r,13.))/(4.194304e6*1.);
}
else if (m == -12){
return (5.*std::sqrt(55393683./Pi)*z/r*(-693. + 42735.*std::pow(z/r,2.) - 666666.*std::pow(z/r,4.) + 3904758.*std::pow(z/r,6.) - 9328033.*std::pow(z/r,8.) + 7632027.*std::pow(z/r,10.))*std::pow((x-I*y)/r,12.))/(2.097152e6*1.);
}
else if (m == -11){
return (3.*std::sqrt(646259635./Pi)*(33. - 6930.*std::pow(z/r,2.) + 213675.*std::pow(z/r,4.) - 2222220.*std::pow(z/r,6.) + 9761895.*std::pow(z/r,8.) - 18656066.*std::pow(z/r,10.) + 12720045.*std::pow(z/r,12.))*std::pow((x-I*y)/r,11.))/(4.194304e6*1.);
}
else if (m == -10){
return (3.*std::sqrt(142823379335./(2.*Pi))*z/r*(33. - 2310.*std::pow(z/r,2.) + 42735.*std::pow(z/r,4.) - 317460.*std::pow(z/r,6.) + 1084655.*std::pow(z/r,8.) - 1696006.*std::pow(z/r,10.) + 978465.*std::pow(z/r,12.))*std::pow((x-I*y)/r,10.))/(2.097152e6*1.);
}
else if (m == -9){
return (std::sqrt(673310216865./Pi)*(-3. + 693.*std::pow(z/r,2.) - 24255.*std::pow(z/r,4.) + 299145.*std::pow(z/r,6.) - 1666665.*std::pow(z/r,8.) + 4555551.*std::pow(z/r,10.) - 5936021.*std::pow(z/r,12.) + 2935395.*std::pow(z/r,14.))*std::pow((x-I*y)/r,9.))/(4.194304e6*1.);
}
else if (m == -8){
return (15.*std::sqrt(44887347791./(2.*Pi))*z/r*(-3. + 231.*std::pow(z/r,2.) - 4851.*std::pow(z/r,4.) + 42735.*std::pow(z/r,6.) - 185185.*std::pow(z/r,8.) + 414141.*std::pow(z/r,10.) - 456617.*std::pow(z/r,12.) + 195693.*std::pow(z/r,14.))*std::pow((x-I*y)/r,8.))/(524288.*1.);
}
else if (m == -7){
return (15.*std::sqrt(1447978961./(2.*Pi))*(3. - 744.*std::pow(z/r,2.) + 28644.*std::pow(z/r,4.) - 401016.*std::pow(z/r,6.) + 2649570.*std::pow(z/r,8.) - 9185176.*std::pow(z/r,10.) + 17117828.*std::pow(z/r,12.) - 16177288.*std::pow(z/r,14.) + 6066483.*std::pow(z/r,16.))*std::pow((x-I*y)/r,7.))/(2.097152e6*1.);
}
else if (m == -6){
return (std::sqrt(1277628495./Pi)*z/r*(765. - 63240.*std::pow(z/r,2.) + 1460844.*std::pow(z/r,4.) - 14608440.*std::pow(z/r,6.) + 75071150.*std::pow(z/r,8.) - 212929080.*std::pow(z/r,10.) + 335772780.*std::pow(z/r,12.) - 275013896.*std::pow(z/r,14.) + 90997245.*std::pow(z/r,16.))*std::pow((x-I*y)/r,6.))/(2.097152e6*1.);
}
else if (m == -5){
return (3.*std::sqrt(44056155./(2.*Pi))*(-85. + 22185.*std::pow(z/r,2.) - 916980.*std::pow(z/r,4.) + 14121492.*std::pow(z/r,6.) - 105911190.*std::pow(z/r,8.) + 435412670.*std::pow(z/r,10.) - 1029157220.*std::pow(z/r,12.) + 1391058660.*std::pow(z/r,14.) - 996925373.*std::pow(z/r,16.) + 293213345.*std::pow(z/r,18.))*std::pow((x-I*y)/r,5.))/(2.097152e6*1.);
}
else if (m == -4){
return (3.*std::sqrt(16231215./(2.*Pi))*z/r*(-1615. + 140505.*std::pow(z/r,2.) - 3484524.*std::pow(z/r,4.) + 38329764.*std::pow(z/r,6.) - 223590290.*std::pow(z/r,8.) + 752076430.*std::pow(z/r,10.) - 1504152860.*std::pow(z/r,12.) + 1762007636.*std::pow(z/r,14.) - 1114210711.*std::pow(z/r,16.) + 293213345.*std::pow(z/r,18.))*std::pow((x-I*y)/r,4.))/(1.048576e6*1.);
}
else if (m == -3){
return (5.*std::sqrt(1082081./(2.*Pi))*(323. - 87210.*std::pow(z/r,2.) + 3793635.*std::pow(z/r,4.) - 62721432.*std::pow(z/r,6.) + 517451814.*std::pow(z/r,8.) - 2414775132.*std::pow(z/r,10.) + 6768687870.*std::pow(z/r,12.) - 11603464920.*std::pow(z/r,14.) + 11893551543.*std::pow(z/r,16.) - 6685264266.*std::pow(z/r,18.) + 1583352063.*std::pow(z/r,20.))*std::pow((x-I*y)/r,3.))/(2.097152e6*1.);
}
else if (m == -2){
return (5.*std::sqrt(35673./Pi)*(29393.*z/r - 2645370.*std::pow(z/r,3.) + 69044157.*std::pow(z/r,5.) - 815378616.*std::pow(z/r,7.) + 5232012786.*std::pow(z/r,9.) - 19976776092.*std::pow(z/r,11.) + 47380815090.*std::pow(z/r,13.) - 70394353848.*std::pow(z/r,15.) + 63665481789.*std::pow(z/r,17.) - 32018897274.*std::pow(z/r,19.) + 6861192273.*std::pow(z/r,21.))*std::pow((x-I*y)/r,2.))/(2.097152e6*1.);
}
else if (m == -1){
return (std::sqrt(3243./(2.*Pi))*(-29393. + 8083075.*std::pow(z/r,2.) - 363738375.*std::pow(z/r,4.) + 6329047725.*std::pow(z/r,6.) - 56057279850.*std::pow(z/r,8.) + 287760703230.*std::pow(z/r,10.) - 915602237550.*std::pow(z/r,12.) + 1861389164250.*std::pow(z/r,14.) - 2419805913525.*std::pow(z/r,16.) + 1945334165775.*std::pow(z/r,18.) - 880519675035.*std::pow(z/r,20.) + 171529806825.*std::pow(z/r,22.))*(x-I*y)/r)/(2.097152e6*1.);
}
else if (m == 0){
return (std::sqrt(47./Pi)*(-2028117.*z/r + 185910725.*std::pow(z/r,3.) - 5019589575.*std::pow(z/r,5.) + 62386327575.*std::pow(z/r,7.) - 429772478850.*std::pow(z/r,9.) + 1805044411170.*std::pow(z/r,11.) - 4859734953150.*std::pow(z/r,13.) + 8562390155550.*std::pow(z/r,15.) - 9821565178425.*std::pow(z/r,17.) + 7064634602025.*std::pow(z/r,19.) - 2893136075115.*std::pow(z/r,21.) + 514589420475.*std::pow(z/r,23.)))/1.048576e6;
}
else if (m == 1){
return -4.76837158203125e-7*(1.*std::sqrt(3243./(2.*Pi))*(-29393. + 8083075.*std::pow(z/r,2.) - 363738375.*std::pow(z/r,4.) + 6329047725.*std::pow(z/r,6.) - 56057279850.*std::pow(z/r,8.) + 287760703230.*std::pow(z/r,10.) - 915602237550.*std::pow(z/r,12.) + 1861389164250.*std::pow(z/r,14.) - 2419805913525.*std::pow(z/r,16.) + 1945334165775.*std::pow(z/r,18.) - 880519675035.*std::pow(z/r,20.) + 171529806825.*std::pow(z/r,22.))*(x+I*y)/r);
}
else if (m == 2){
return (5.*1.*std::sqrt(35673./Pi)*(29393.*z/r - 2645370.*std::pow(z/r,3.) + 69044157.*std::pow(z/r,5.) - 815378616.*std::pow(z/r,7.) + 5232012786.*std::pow(z/r,9.) - 19976776092.*std::pow(z/r,11.) + 47380815090.*std::pow(z/r,13.) - 70394353848.*std::pow(z/r,15.) + 63665481789.*std::pow(z/r,17.) - 32018897274.*std::pow(z/r,19.) + 6861192273.*std::pow(z/r,21.))*std::pow((x+I*y)/r,2.))/2.097152e6;
}
else if (m == 3){
return (-5.*1.*std::sqrt(1082081./(2.*Pi))*(323. - 87210.*std::pow(z/r,2.) + 3793635.*std::pow(z/r,4.) - 62721432.*std::pow(z/r,6.) + 517451814.*std::pow(z/r,8.) - 2414775132.*std::pow(z/r,10.) + 6768687870.*std::pow(z/r,12.) - 11603464920.*std::pow(z/r,14.) + 11893551543.*std::pow(z/r,16.) - 6685264266.*std::pow(z/r,18.) + 1583352063.*std::pow(z/r,20.))*std::pow((x+I*y)/r,3.))/2.097152e6;
}
else if (m == 4){
return (3.*1.*std::sqrt(16231215./(2.*Pi))*z/r*(-1615. + 140505.*std::pow(z/r,2.) - 3484524.*std::pow(z/r,4.) + 38329764.*std::pow(z/r,6.) - 223590290.*std::pow(z/r,8.) + 752076430.*std::pow(z/r,10.) - 1504152860.*std::pow(z/r,12.) + 1762007636.*std::pow(z/r,14.) - 1114210711.*std::pow(z/r,16.) + 293213345.*std::pow(z/r,18.))*std::pow((x+I*y)/r,4.))/1.048576e6;
}
else if (m == 5){
return (-3.*1.*std::sqrt(44056155./(2.*Pi))*(-85. + 22185.*std::pow(z/r,2.) - 916980.*std::pow(z/r,4.) + 14121492.*std::pow(z/r,6.) - 105911190.*std::pow(z/r,8.) + 435412670.*std::pow(z/r,10.) - 1029157220.*std::pow(z/r,12.) + 1391058660.*std::pow(z/r,14.) - 996925373.*std::pow(z/r,16.) + 293213345.*std::pow(z/r,18.))*std::pow((x+I*y)/r,5.))/2.097152e6;
}
else if (m == 6){
return (1.*std::sqrt(1277628495./Pi)*z/r*(765. - 63240.*std::pow(z/r,2.) + 1460844.*std::pow(z/r,4.) - 14608440.*std::pow(z/r,6.) + 75071150.*std::pow(z/r,8.) - 212929080.*std::pow(z/r,10.) + 335772780.*std::pow(z/r,12.) - 275013896.*std::pow(z/r,14.) + 90997245.*std::pow(z/r,16.))*std::pow((x+I*y)/r,6.))/2.097152e6;
}
else if (m == 7){
return (-15.*1.*std::sqrt(1447978961./(2.*Pi))*(3. - 744.*std::pow(z/r,2.) + 28644.*std::pow(z/r,4.) - 401016.*std::pow(z/r,6.) + 2649570.*std::pow(z/r,8.) - 9185176.*std::pow(z/r,10.) + 17117828.*std::pow(z/r,12.) - 16177288.*std::pow(z/r,14.) + 6066483.*std::pow(z/r,16.))*std::pow((x+I*y)/r,7.))/2.097152e6;
}
else if (m == 8){
return (15.*1.*std::sqrt(44887347791./(2.*Pi))*z/r*(-3. + 231.*std::pow(z/r,2.) - 4851.*std::pow(z/r,4.) + 42735.*std::pow(z/r,6.) - 185185.*std::pow(z/r,8.) + 414141.*std::pow(z/r,10.) - 456617.*std::pow(z/r,12.) + 195693.*std::pow(z/r,14.))*std::pow((x+I*y)/r,8.))/524288.;
}
else if (m == 9){
return -2.384185791015625e-7*(1.*std::sqrt(673310216865./Pi)*(-3. + 693.*std::pow(z/r,2.) - 24255.*std::pow(z/r,4.) + 299145.*std::pow(z/r,6.) - 1666665.*std::pow(z/r,8.) + 4555551.*std::pow(z/r,10.) - 5936021.*std::pow(z/r,12.) + 2935395.*std::pow(z/r,14.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (3.*1.*std::sqrt(142823379335./(2.*Pi))*z/r*(33. - 2310.*std::pow(z/r,2.) + 42735.*std::pow(z/r,4.) - 317460.*std::pow(z/r,6.) + 1084655.*std::pow(z/r,8.) - 1696006.*std::pow(z/r,10.) + 978465.*std::pow(z/r,12.))*std::pow((x+I*y)/r,10.))/2.097152e6;
}
else if (m == 11){
return (-3.*1.*std::sqrt(646259635./Pi)*(33. - 6930.*std::pow(z/r,2.) + 213675.*std::pow(z/r,4.) - 2222220.*std::pow(z/r,6.) + 9761895.*std::pow(z/r,8.) - 18656066.*std::pow(z/r,10.) + 12720045.*std::pow(z/r,12.))*std::pow((x+I*y)/r,11.))/4.194304e6;
}
else if (m == 12){
return (5.*1.*std::sqrt(55393683./Pi)*z/r*(-693. + 42735.*std::pow(z/r,2.) - 666666.*std::pow(z/r,4.) + 3904758.*std::pow(z/r,6.) - 9328033.*std::pow(z/r,8.) + 7632027.*std::pow(z/r,10.))*std::pow((x+I*y)/r,12.))/2.097152e6;
}
else if (m == 13){
return (-15.*1.*std::sqrt(609330513./Pi)*(-7. + 1295.*std::pow(z/r,2.) - 33670.*std::pow(z/r,4.) + 276094.*std::pow(z/r,6.) - 848003.*std::pow(z/r,8.) + 848003.*std::pow(z/r,10.))*std::pow((x+I*y)/r,13.))/4.194304e6;
}
else if (m == 14){
return (3.*1.*std::sqrt(112726144905./(2.*Pi))*z/r*(35. - 1820.*std::pow(z/r,2.) + 22386.*std::pow(z/r,4.) - 91676.*std::pow(z/r,6.) + 114595.*std::pow(z/r,8.))*std::pow((x+I*y)/r,14.))/2.097152e6;
}
else if (m == 15){
return -2.384185791015625e-7*(1.*std::sqrt(5932954995./Pi)*(35. - 5460.*std::pow(z/r,2.) + 111930.*std::pow(z/r,4.) - 641732.*std::pow(z/r,6.) + 1031355.*std::pow(z/r,8.))*std::pow((x+I*y)/r,15.));
}
else if (m == 16){
return (3.*1.*std::sqrt(25709471645./(2.*Pi))*z/r*(-35. + 1435.*std::pow(z/r,2.) - 12341.*std::pow(z/r,4.) + 26445.*std::pow(z/r,6.))*std::pow((x+I*y)/r,16.))/1.048576e6;
}
else if (m == 17){
return (-15.*1.*std::sqrt(35993260303./Pi)*(-1. + 123.*std::pow(z/r,2.) - 1763.*std::pow(z/r,4.) + 5289.*std::pow(z/r,6.))*std::pow((x+I*y)/r,17.))/4.194304e6;
}
else if (m == 18){
return (5.*1.*std::sqrt(4427171017269./(2.*Pi))*z/r*(3. - 86.*std::pow(z/r,2.) + 387.*std::pow(z/r,4.))*std::pow((x+I*y)/r,18.))/2.097152e6;
}
else if (m == 19){
return (-3.*1.*std::sqrt(1054088337445./Pi)*(1. - 86.*std::pow(z/r,2.) + 645.*std::pow(z/r,4.))*std::pow((x+I*y)/r,19.))/4.194304e6;
}
else if (m == 20){
return (3.*1.*std::sqrt(45325798510135./Pi)*z/r*(-1. + 15.*std::pow(z/r,2.))*std::pow((x+I*y)/r,20.))/2.097152e6;
}
else if (m == 21){
return -2.384185791015625e-7*(1.*std::sqrt(12361581411855./Pi)*(-1. + 45.*std::pow(z/r,2.))*std::pow((x+I*y)/r,21.));
}
else if (m == 22){
return (15.*1.*std::sqrt(2472316282371./(2.*Pi))*z/r*std::pow((x+I*y)/r,22.))/2.097152e6;
}
else if (m == 23){
return (-15.*1.*std::sqrt(107492012277./Pi)*std::pow((x+I*y)/r,23.))/4.194304e6;
}
else{return 0.;}
}

else if (l == 24){
if(m == -24){
return (105.*std::sqrt(35830670759./Pi)*std::pow((x-I*y)/r,24.))/(1.6777216e7*1.);
}
else if (m == -23){
return (105.*std::sqrt(107492012277./Pi)*z/r*std::pow((x-I*y)/r,23.))/(4.194304e6*1.);
}
else if (m == -22){
return (105.*std::sqrt(2287064091./(2.*Pi))*(-1. + 47.*std::pow(z/r,2.))*std::pow((x-I*y)/r,22.))/(4.194304e6*1.);
}
else if (m == -21){
return (105.*std::sqrt(17534158031./Pi)*z/r*(-3. + 47.*std::pow(z/r,2.))*std::pow((x-I*y)/r,21.))/(4.194304e6*1.);
}
else if (m == -20){
return (21.*std::sqrt(87670790155./Pi)*(1. - 90.*std::pow(z/r,2.) + 705.*std::pow(z/r,4.))*std::pow((x-I*y)/r,20.))/(8.388608e6*1.);
}
else if (m == -19){
return (105.*std::sqrt(192875738341./Pi)*z/r*(1. - 30.*std::pow(z/r,2.) + 141.*std::pow(z/r,4.))*std::pow((x-I*y)/r,19.))/(4.194304e6*1.);
}
else if (m == -18){
return (35.*std::sqrt(13456446861./(2.*Pi))*(-1. + 129.*std::pow(z/r,2.) - 1935.*std::pow(z/r,4.) + 6063.*std::pow(z/r,6.))*std::pow((x-I*y)/r,18.))/(4.194304e6*1.);
}
else if (m == -17){
return (105.*std::sqrt(4485482287./Pi)*z/r*(-7. + 301.*std::pow(z/r,2.) - 2709.*std::pow(z/r,4.) + 6063.*std::pow(z/r,6.))*std::pow((x-I*y)/r,17.))/(4.194304e6*1.);
}
else if (m == -16){
return (105.*std::sqrt(109402007./(2.*Pi))*(7. - 1148.*std::pow(z/r,2.) + 24682.*std::pow(z/r,4.) - 148092.*std::pow(z/r,6.) + 248583.*std::pow(z/r,8.))*std::pow((x-I*y)/r,16.))/(8.388608e6*1.);
}
else if (m == -15){
return (21.*std::sqrt(547010035./Pi)*z/r*(105. - 5740.*std::pow(z/r,2.) + 74046.*std::pow(z/r,4.) - 317340.*std::pow(z/r,6.) + 414305.*std::pow(z/r,8.))*std::pow((x-I*y)/r,15.))/(4.194304e6*1.);
}
else if (m == -14){
return (105.*std::sqrt(25246617./(2.*Pi))*(-7. + 1365.*std::pow(z/r,2.) - 37310.*std::pow(z/r,4.) + 320866.*std::pow(z/r,6.) - 1031355.*std::pow(z/r,8.) + 1077193.*std::pow(z/r,10.))*std::pow((x-I*y)/r,14.))/(4.194304e6*1.);
}
else if (m == -13){
return (105.*std::sqrt(43607793./Pi)*z/r*(-77. + 5005.*std::pow(z/r,2.) - 82082.*std::pow(z/r,4.) + 504218.*std::pow(z/r,6.) - 1260545.*std::pow(z/r,8.) + 1077193.*std::pow(z/r,10.))*std::pow((x-I*y)/r,13.))/(4.194304e6*1.);
}
else if (m == -12){
return (105.*std::sqrt(392863./Pi)*(77. - 17094.*std::pow(z/r,2.) + 555555.*std::pow(z/r,4.) - 6074068.*std::pow(z/r,6.) + 27984099.*std::pow(z/r,8.) - 55968198.*std::pow(z/r,10.) + 39856141.*std::pow(z/r,12.))*std::pow((x-I*y)/r,12.))/(8.388608e6*1.);
}
else if (m == -11){
return (105.*std::sqrt(5107219./Pi)*z/r*(231. - 17094.*std::pow(z/r,2.) + 333333.*std::pow(z/r,4.) - 2603172.*std::pow(z/r,6.) + 9328033.*std::pow(z/r,8.) - 15264054.*std::pow(z/r,10.) + 9197571.*std::pow(z/r,12.))*std::pow((x-I*y)/r,11.))/(4.194304e6*1.);
}
else if (m == -10){
return (21.*std::sqrt(25536095./(2.*Pi))*(-33. + 8085.*std::pow(z/r,2.) - 299145.*std::pow(z/r,4.) + 3888885.*std::pow(z/r,6.) - 22777755.*std::pow(z/r,8.) + 65296231.*std::pow(z/r,10.) - 89040315.*std::pow(z/r,12.) + 45987855.*std::pow(z/r,14.))*std::pow((x-I*y)/r,10.))/(4.194304e6*1.);
}
else if (m == -9){
return (35.*std::sqrt(260468169./Pi)*z/r*(-99. + 8085.*std::pow(z/r,2.) - 179487.*std::pow(z/r,4.) + 1666665.*std::pow(z/r,6.) - 7592585.*std::pow(z/r,8.) + 17808063.*std::pow(z/r,10.) - 20547765.*std::pow(z/r,12.) + 9197571.*std::pow(z/r,14.))*std::pow((x-I*y)/r,9.))/(4.194304e6*1.);
}
else if (m == -8){
return (105.*std::sqrt(955049953./Pi)*(3. - 792.*std::pow(z/r,2.) + 32340.*std::pow(z/r,4.) - 478632.*std::pow(z/r,6.) + 3333330.*std::pow(z/r,8.) - 12148136.*std::pow(z/r,10.) + 23744084.*std::pow(z/r,12.) - 23483160.*std::pow(z/r,14.) + 9197571.*std::pow(z/r,16.))*std::pow((x-I*y)/r,8.))/(1.6777216e7*1.);
}
else if (m == -7){
return (105.*std::sqrt(56179409./(2.*Pi))*z/r*(51. - 4488.*std::pow(z/r,2.) + 109956.*std::pow(z/r,4.) - 1162392.*std::pow(z/r,6.) + 6296290.*std::pow(z/r,8.) - 18774392.*std::pow(z/r,10.) + 31049956.*std::pow(z/r,12.) - 26614248.*std::pow(z/r,14.) + 9197571.*std::pow(z/r,16.))*std::pow((x-I*y)/r,7.))/(2.097152e6*1.);
}
else if (m == -6){
return (105.*std::sqrt(1812239./Pi)*(-17. + 4743.*std::pow(z/r,2.) - 208692.*std::pow(z/r,4.) + 3408636.*std::pow(z/r,6.) - 27025614.*std::pow(z/r,8.) + 117110994.*std::pow(z/r,10.) - 291003076.*std::pow(z/r,12.) + 412520844.*std::pow(z/r,14.) - 309390633.*std::pow(z/r,16.) + 95041567.*std::pow(z/r,18.))*std::pow((x-I*y)/r,6.))/(4.194304e6*1.);
}
else if (m == -5){
return (21.*std::sqrt(1430715./(2.*Pi))*z/r*(-1615. + 150195.*std::pow(z/r,2.) - 3965148.*std::pow(z/r,4.) + 46260060.*std::pow(z/r,6.) - 285270370.*std::pow(z/r,8.) + 1011413130.*std::pow(z/r,10.) - 2126560940.*std::pow(z/r,12.) + 2612632012.*std::pow(z/r,14.) - 1728947655.*std::pow(z/r,16.) + 475207835.*std::pow(z/r,18.))*std::pow((x-I*y)/r,5.))/(2.097152e6*1.);
}
else if (m == -4){
return (105.*std::sqrt(9867./(2.*Pi))*(323. - 93670.*std::pow(z/r,2.) + 4355655.*std::pow(z/r,4.) - 76659528.*std::pow(z/r,6.) + 670770870.*std::pow(z/r,8.) - 3309136292.*std::pow(z/r,10.) + 9776993590.*std::pow(z/r,12.) - 17620076360.*std::pow(z/r,14.) + 18941582087.*std::pow(z/r,16.) - 11142107110.*std::pow(z/r,18.) + 2756205443.*std::pow(z/r,20.))*std::pow((x-I*y)/r,4.))/(4.194304e6*1.);
}
else if (m == -3){
return (105.*std::sqrt(3289./(2.*Pi))*(6783.*z/r - 655690.*std::pow(z/r,3.) + 18293751.*std::pow(z/r,5.) - 229978584.*std::pow(z/r,7.) + 1565132030.*std::pow(z/r,9.) - 6317442012.*std::pow(z/r,11.) + 15793605030.*std::pow(z/r,13.) - 24668106904.*std::pow(z/r,15.) + 23398424931.*std::pow(z/r,17.) - 12314960490.*std::pow(z/r,19.) + 2756205443.*std::pow(z/r,21.))*std::pow((x-I*y)/r,3.))/(2.097152e6*1.);
}
else if (m == -2){
return (35.*std::sqrt(897./Pi)*(-2261. + 671517.*std::pow(z/r,2.) - 32456655.*std::pow(z/r,4.) + 603693783.*std::pow(z/r,6.) - 5691969954.*std::pow(z/r,8.) + 30989614194.*std::pow(z/r,10.) - 104237793198.*std::pow(z/r,12.) + 223366699710.*std::pow(z/r,14.) - 305267822937.*std::pow(z/r,16.) + 257382674241.*std::pow(z/r,18.) - 121918108851.*std::pow(z/r,20.) + 24805848987.*std::pow(z/r,22.))*std::pow((x-I*y)/r,2.))/(4.194304e6*1.);
}
else if (m == -1){
return (35.*std::sqrt(3./(2.*Pi))*(-676039.*z/r + 66927861.*std::pow(z/r,3.) - 1940907969.*std::pow(z/r,5.) + 25786348731.*std::pow(z/r,7.) - 189099890694.*std::pow(z/r,9.) + 842354058546.*std::pow(z/r,11.) - 2397469243554.*std::pow(z/r,13.) + 4452442880886.*std::pow(z/r,15.) - 5369122297539.*std::pow(z/r,17.) + 4050390505161.*std::pow(z/r,19.) - 1735881645069.*std::pow(z/r,21.) + 322476036831.*std::pow(z/r,23.))*(x-I*y)/r)/(2.097152e6*1.);
}
else if (m == 0){
return (7.*(676039. - 202811700.*std::pow(z/r,2.) + 10039179150.*std::pow(z/r,4.) - 194090796900.*std::pow(z/r,6.) + 1933976154825.*std::pow(z/r,8.) - 11345993441640.*std::pow(z/r,10.) + 42117702927300.*std::pow(z/r,12.) - 102748681866600.*std::pow(z/r,14.) + 166966608033225.*std::pow(z/r,16.) - 178970743251300.*std::pow(z/r,18.) + 121511715154830.*std::pow(z/r,20.) - 47342226683700.*std::pow(z/r,22.) + 8061900920775.*std::pow(z/r,24.)))/(8.388608e6*std::sqrt(Pi));
}
else if (m == 1){
return (-35.*1.*std::sqrt(3./(2.*Pi))*(-676039.*z/r + 66927861.*std::pow(z/r,3.) - 1940907969.*std::pow(z/r,5.) + 25786348731.*std::pow(z/r,7.) - 189099890694.*std::pow(z/r,9.) + 842354058546.*std::pow(z/r,11.) - 2397469243554.*std::pow(z/r,13.) + 4452442880886.*std::pow(z/r,15.) - 5369122297539.*std::pow(z/r,17.) + 4050390505161.*std::pow(z/r,19.) - 1735881645069.*std::pow(z/r,21.) + 322476036831.*std::pow(z/r,23.))*(x+I*y)/r)/2.097152e6;
}
else if (m == 2){
return (35.*1.*std::sqrt(897./Pi)*(-2261. + 671517.*std::pow(z/r,2.) - 32456655.*std::pow(z/r,4.) + 603693783.*std::pow(z/r,6.) - 5691969954.*std::pow(z/r,8.) + 30989614194.*std::pow(z/r,10.) - 104237793198.*std::pow(z/r,12.) + 223366699710.*std::pow(z/r,14.) - 305267822937.*std::pow(z/r,16.) + 257382674241.*std::pow(z/r,18.) - 121918108851.*std::pow(z/r,20.) + 24805848987.*std::pow(z/r,22.))*std::pow((x+I*y)/r,2.))/4.194304e6;
}
else if (m == 3){
return (-105.*1.*std::sqrt(3289./(2.*Pi))*(6783.*z/r - 655690.*std::pow(z/r,3.) + 18293751.*std::pow(z/r,5.) - 229978584.*std::pow(z/r,7.) + 1565132030.*std::pow(z/r,9.) - 6317442012.*std::pow(z/r,11.) + 15793605030.*std::pow(z/r,13.) - 24668106904.*std::pow(z/r,15.) + 23398424931.*std::pow(z/r,17.) - 12314960490.*std::pow(z/r,19.) + 2756205443.*std::pow(z/r,21.))*std::pow((x+I*y)/r,3.))/2.097152e6;
}
else if (m == 4){
return (105.*1.*std::sqrt(9867./(2.*Pi))*(323. - 93670.*std::pow(z/r,2.) + 4355655.*std::pow(z/r,4.) - 76659528.*std::pow(z/r,6.) + 670770870.*std::pow(z/r,8.) - 3309136292.*std::pow(z/r,10.) + 9776993590.*std::pow(z/r,12.) - 17620076360.*std::pow(z/r,14.) + 18941582087.*std::pow(z/r,16.) - 11142107110.*std::pow(z/r,18.) + 2756205443.*std::pow(z/r,20.))*std::pow((x+I*y)/r,4.))/4.194304e6;
}
else if (m == 5){
return (-21.*1.*std::sqrt(1430715./(2.*Pi))*z/r*(-1615. + 150195.*std::pow(z/r,2.) - 3965148.*std::pow(z/r,4.) + 46260060.*std::pow(z/r,6.) - 285270370.*std::pow(z/r,8.) + 1011413130.*std::pow(z/r,10.) - 2126560940.*std::pow(z/r,12.) + 2612632012.*std::pow(z/r,14.) - 1728947655.*std::pow(z/r,16.) + 475207835.*std::pow(z/r,18.))*std::pow((x+I*y)/r,5.))/2.097152e6;
}
else if (m == 6){
return (105.*1.*std::sqrt(1812239./Pi)*(-17. + 4743.*std::pow(z/r,2.) - 208692.*std::pow(z/r,4.) + 3408636.*std::pow(z/r,6.) - 27025614.*std::pow(z/r,8.) + 117110994.*std::pow(z/r,10.) - 291003076.*std::pow(z/r,12.) + 412520844.*std::pow(z/r,14.) - 309390633.*std::pow(z/r,16.) + 95041567.*std::pow(z/r,18.))*std::pow((x+I*y)/r,6.))/4.194304e6;
}
else if (m == 7){
return (-105.*1.*std::sqrt(56179409./(2.*Pi))*z/r*(51. - 4488.*std::pow(z/r,2.) + 109956.*std::pow(z/r,4.) - 1162392.*std::pow(z/r,6.) + 6296290.*std::pow(z/r,8.) - 18774392.*std::pow(z/r,10.) + 31049956.*std::pow(z/r,12.) - 26614248.*std::pow(z/r,14.) + 9197571.*std::pow(z/r,16.))*std::pow((x+I*y)/r,7.))/2.097152e6;
}
else if (m == 8){
return (105.*1.*std::sqrt(955049953./Pi)*(3. - 792.*std::pow(z/r,2.) + 32340.*std::pow(z/r,4.) - 478632.*std::pow(z/r,6.) + 3333330.*std::pow(z/r,8.) - 12148136.*std::pow(z/r,10.) + 23744084.*std::pow(z/r,12.) - 23483160.*std::pow(z/r,14.) + 9197571.*std::pow(z/r,16.))*std::pow((x+I*y)/r,8.))/1.6777216e7;
}
else if (m == 9){
return (-35.*1.*std::sqrt(260468169./Pi)*z/r*(-99. + 8085.*std::pow(z/r,2.) - 179487.*std::pow(z/r,4.) + 1666665.*std::pow(z/r,6.) - 7592585.*std::pow(z/r,8.) + 17808063.*std::pow(z/r,10.) - 20547765.*std::pow(z/r,12.) + 9197571.*std::pow(z/r,14.))*std::pow((x+I*y)/r,9.))/4.194304e6;
}
else if (m == 10){
return (21.*1.*std::sqrt(25536095./(2.*Pi))*(-33. + 8085.*std::pow(z/r,2.) - 299145.*std::pow(z/r,4.) + 3888885.*std::pow(z/r,6.) - 22777755.*std::pow(z/r,8.) + 65296231.*std::pow(z/r,10.) - 89040315.*std::pow(z/r,12.) + 45987855.*std::pow(z/r,14.))*std::pow((x+I*y)/r,10.))/4.194304e6;
}
else if (m == 11){
return (-105.*1.*std::sqrt(5107219./Pi)*z/r*(231. - 17094.*std::pow(z/r,2.) + 333333.*std::pow(z/r,4.) - 2603172.*std::pow(z/r,6.) + 9328033.*std::pow(z/r,8.) - 15264054.*std::pow(z/r,10.) + 9197571.*std::pow(z/r,12.))*std::pow((x+I*y)/r,11.))/4.194304e6;
}
else if (m == 12){
return (105.*1.*std::sqrt(392863./Pi)*(77. - 17094.*std::pow(z/r,2.) + 555555.*std::pow(z/r,4.) - 6074068.*std::pow(z/r,6.) + 27984099.*std::pow(z/r,8.) - 55968198.*std::pow(z/r,10.) + 39856141.*std::pow(z/r,12.))*std::pow((x+I*y)/r,12.))/8.388608e6;
}
else if (m == 13){
return (-105.*1.*std::sqrt(43607793./Pi)*z/r*(-77. + 5005.*std::pow(z/r,2.) - 82082.*std::pow(z/r,4.) + 504218.*std::pow(z/r,6.) - 1260545.*std::pow(z/r,8.) + 1077193.*std::pow(z/r,10.))*std::pow((x+I*y)/r,13.))/4.194304e6;
}
else if (m == 14){
return (105.*1.*std::sqrt(25246617./(2.*Pi))*(-7. + 1365.*std::pow(z/r,2.) - 37310.*std::pow(z/r,4.) + 320866.*std::pow(z/r,6.) - 1031355.*std::pow(z/r,8.) + 1077193.*std::pow(z/r,10.))*std::pow((x+I*y)/r,14.))/4.194304e6;
}
else if (m == 15){
return (-21.*1.*std::sqrt(547010035./Pi)*z/r*(105. - 5740.*std::pow(z/r,2.) + 74046.*std::pow(z/r,4.) - 317340.*std::pow(z/r,6.) + 414305.*std::pow(z/r,8.))*std::pow((x+I*y)/r,15.))/4.194304e6;
}
else if (m == 16){
return (105.*1.*std::sqrt(109402007./(2.*Pi))*(7. - 1148.*std::pow(z/r,2.) + 24682.*std::pow(z/r,4.) - 148092.*std::pow(z/r,6.) + 248583.*std::pow(z/r,8.))*std::pow((x+I*y)/r,16.))/8.388608e6;
}
else if (m == 17){
return (-105.*1.*std::sqrt(4485482287./Pi)*z/r*(-7. + 301.*std::pow(z/r,2.) - 2709.*std::pow(z/r,4.) + 6063.*std::pow(z/r,6.))*std::pow((x+I*y)/r,17.))/4.194304e6;
}
else if (m == 18){
return (35.*1.*std::sqrt(13456446861./(2.*Pi))*(-1. + 129.*std::pow(z/r,2.) - 1935.*std::pow(z/r,4.) + 6063.*std::pow(z/r,6.))*std::pow((x+I*y)/r,18.))/4.194304e6;
}
else if (m == 19){
return (-105.*1.*std::sqrt(192875738341./Pi)*z/r*(1. - 30.*std::pow(z/r,2.) + 141.*std::pow(z/r,4.))*std::pow((x+I*y)/r,19.))/4.194304e6;
}
else if (m == 20){
return (21.*1.*std::sqrt(87670790155./Pi)*(1. - 90.*std::pow(z/r,2.) + 705.*std::pow(z/r,4.))*std::pow((x+I*y)/r,20.))/8.388608e6;
}
else if (m == 21){
return (-105.*1.*std::sqrt(17534158031./Pi)*z/r*(-3. + 47.*std::pow(z/r,2.))*std::pow((x+I*y)/r,21.))/4.194304e6;
}
else if (m == 22){
return (105.*1.*std::sqrt(2287064091./(2.*Pi))*(-1. + 47.*std::pow(z/r,2.))*std::pow((x+I*y)/r,22.))/4.194304e6;
}
else if (m == 23){
return (-105.*1.*std::sqrt(107492012277./Pi)*z/r*std::pow((x+I*y)/r,23.))/4.194304e6;
}
else if (m == 24){
return (105.*1.*std::sqrt(35830670759./Pi)*std::pow((x+I*y)/r,24.))/1.6777216e7;
}
else{return 0.;}
}

else if (l == 25){
if(m == -25){
return (21.*std::sqrt(1827364208709./(2.*Pi))*std::pow((x-I*y)/r,25.))/(1.6777216e7*1.);
}
else if (m == -24){
return (105.*std::sqrt(1827364208709./Pi)*z/r*std::pow((x-I*y)/r,24.))/(1.6777216e7*1.);
}
else if (m == -23){
return (15.*std::sqrt(1827364208709./(2.*Pi))*(-1. + 49.*std::pow(z/r,2.))*std::pow((x-I*y)/r,23.))/(1.6777216e7*1.);
}
else if (m == -22){
return (15.*std::sqrt(1827364208709./(2.*Pi))*z/r*(-3. + 49.*std::pow(z/r,2.))*std::pow((x-I*y)/r,22.))/(4.194304e6*1.);
}
else if (m == -21){
return (15.*std::sqrt(38880089547./(2.*Pi))*(3. - 282.*std::pow(z/r,2.) + 2303.*std::pow(z/r,4.))*std::pow((x-I*y)/r,21.))/(8.388608e6*1.);
}
else if (m == -20){
return (3.*std::sqrt(4471210297905./Pi)*z/r*(15. - 470.*std::pow(z/r,2.) + 2303.*std::pow(z/r,4.))*std::pow((x-I*y)/r,20.))/(8.388608e6*1.);
}
else if (m == -19){
return (15.*std::sqrt(298080686527./(2.*Pi))*(-1. + 135.*std::pow(z/r,2.) - 2115.*std::pow(z/r,4.) + 6909.*std::pow(z/r,6.))*std::pow((x-I*y)/r,19.))/(8.388608e6*1.);
}
else if (m == -18){
return (15.*std::sqrt(22952212862579./(2.*Pi))*z/r*(-1. + 45.*std::pow(z/r,2.) - 423.*std::pow(z/r,4.) + 987.*std::pow(z/r,6.))*std::pow((x-I*y)/r,18.))/(4.194304e6*1.);
}
else if (m == -17){
return (15.*std::sqrt(533772392153./Pi)*(1. - 172.*std::pow(z/r,2.) + 3870.*std::pow(z/r,4.) - 24252.*std::pow(z/r,6.) + 42441.*std::pow(z/r,8.))*std::pow((x-I*y)/r,17.))/(1.6777216e7*1.);
}
else if (m == -16){
return (15.*std::sqrt(228759596637./(2.*Pi))*z/r*(21. - 1204.*std::pow(z/r,2.) + 16254.*std::pow(z/r,4.) - 72756.*std::pow(z/r,6.) + 99029.*std::pow(z/r,8.))*std::pow((x-I*y)/r,16.))/(8.388608e6*1.);
}
else if (m == -15){
return (3.*std::sqrt(27897511785./Pi)*(-21. + 4305.*std::pow(z/r,2.) - 123410.*std::pow(z/r,4.) + 1110690.*std::pow(z/r,6.) - 3728745.*std::pow(z/r,8.) + 4060189.*std::pow(z/r,10.))*std::pow((x-I*y)/r,15.))/(1.6777216e7*1.);
}
else if (m == -14){
return (15.*std::sqrt(507227487./(2.*Pi))*z/r*(-231. + 15785.*std::pow(z/r,2.) - 271502.*std::pow(z/r,4.) + 1745370.*std::pow(z/r,6.) - 4557355.*std::pow(z/r,8.) + 4060189.*std::pow(z/r,10.))*std::pow((x-I*y)/r,14.))/(4.194304e6*1.);
}
else if (m == -13){
return (15.*std::sqrt(39017499./(2.*Pi))*(77. - 18018.*std::pow(z/r,2.) + 615615.*std::pow(z/r,4.) - 7059052.*std::pow(z/r,6.) + 34034715.*std::pow(z/r,8.) - 71094738.*std::pow(z/r,10.) + 52782457.*std::pow(z/r,12.))*std::pow((x-I*y)/r,13.))/(8.388608e6*1.);
}
else if (m == -12){
return (15.*std::sqrt(9637322253./Pi)*z/r*(77. - 6006.*std::pow(z/r,2.) + 123123.*std::pow(z/r,4.) - 1008436.*std::pow(z/r,6.) + 3781635.*std::pow(z/r,8.) - 6463158.*std::pow(z/r,10.) + 4060189.*std::pow(z/r,12.))*std::pow((x-I*y)/r,12.))/(8.388608e6*1.);
}
else if (m == -11){
return (15.*std::sqrt(1823277183./(2.*Pi))*(-11. + 2849.*std::pow(z/r,2.) - 111111.*std::pow(z/r,4.) + 1518517.*std::pow(z/r,6.) - 9328033.*std::pow(z/r,8.) + 27984099.*std::pow(z/r,10.) - 39856141.*std::pow(z/r,12.) + 21460999.*std::pow(z/r,14.))*std::pow((x-I*y)/r,11.))/(8.388608e6*1.);
}
else if (m == -10){
return (3.*std::sqrt(3038795305./(2.*Pi))*z/r*(-495. + 42735.*std::pow(z/r,2.) - 999999.*std::pow(z/r,4.) + 9761895.*std::pow(z/r,6.) - 46640165.*std::pow(z/r,8.) + 114480405.*std::pow(z/r,10.) - 137963565.*std::pow(z/r,12.) + 64382997.*std::pow(z/r,14.))*std::pow((x-I*y)/r,10.))/(4.194304e6*1.);
}
else if (m == -9){
return (15.*std::sqrt(86822723./(2.*Pi))*(99. - 27720.*std::pow(z/r,2.) + 1196580.*std::pow(z/r,4.) - 18666648.*std::pow(z/r,6.) + 136666530.*std::pow(z/r,8.) - 522369848.*std::pow(z/r,10.) + 1068483780.*std::pow(z/r,12.) - 1103708520.*std::pow(z/r,14.) + 450680979.*std::pow(z/r,16.))*std::pow((x-I*y)/r,9.))/(1.6777216e7*1.);
}
else if (m == -8){
return (15.*std::sqrt(86822723./Pi)*z/r*(1683. - 157080.*std::pow(z/r,2.) + 4068372.*std::pow(z/r,4.) - 45333288.*std::pow(z/r,6.) + 258147890.*std::pow(z/r,8.) - 807298856.*std::pow(z/r,10.) + 1397248020.*std::pow(z/r,12.) - 1250869656.*std::pow(z/r,14.) + 450680979.*std::pow(z/r,16.))*std::pow((x-I*y)/r,8.))/(1.6777216e7*1.);
}
else if (m == -7){
return (15.*std::sqrt(2865149859./(2.*Pi))*(-17. + 5049.*std::pow(z/r,2.) - 235620.*std::pow(z/r,4.) + 4068372.*std::pow(z/r,6.) - 33999966.*std::pow(z/r,8.) + 154888734.*std::pow(z/r,10.) - 403649428.*std::pow(z/r,12.) + 598820580.*std::pow(z/r,14.) - 469076121.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x-I*y)/r,7.))/(1.6777216e7*1.);
}
else if (m == -6){
return (15.*std::sqrt(150797361./Pi)*z/r*(-323. + 31977.*std::pow(z/r,2.) - 895356.*std::pow(z/r,4.) + 11042724.*std::pow(z/r,6.) - 71777706.*std::pow(z/r,8.) + 267535086.*std::pow(z/r,10.) - 589949164.*std::pow(z/r,12.) + 758506068.*std::pow(z/r,14.) - 524261547.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x-I*y)/r,6.))/(4.194304e6*1.);
}
else if (m == -5){
return (3.*std::sqrt(24322155./Pi)*(323. - 100130.*std::pow(z/r,2.) + 4956435.*std::pow(z/r,4.) - 92520120.*std::pow(z/r,6.) + 855811110.*std::pow(z/r,8.) - 4450217772.*std::pow(z/r,10.) + 13822646110.*std::pow(z/r,12.) - 26126320120.*std::pow(z/r,14.) + 29392110135.*std::pow(z/r,16.) - 18057897730.*std::pow(z/r,18.) + 4657036783.*std::pow(z/r,20.))*std::pow((x-I*y)/r,5.))/(8.388608e6*1.);
}
else if (m == -4){
return (15.*std::sqrt(34051017./(2.*Pi))*(969.*z/r - 100130.*std::pow(z/r,3.) + 2973861.*std::pow(z/r,5.) - 39651480.*std::pow(z/r,7.) + 285270370.*std::pow(z/r,9.) - 1213695756.*std::pow(z/r,11.) + 3189841410.*std::pow(z/r,13.) - 5225264024.*std::pow(z/r,15.) + 5186842965.*std::pow(z/r,17.) - 2851247010.*std::pow(z/r,19.) + 665290969.*std::pow(z/r,21.))*std::pow((x-I*y)/r,4.))/(4.194304e6*1.);
}
else if (m == -3){
return (15.*std::sqrt(106743./Pi)*(-969. + 309111.*std::pow(z/r,2.) - 15970735.*std::pow(z/r,4.) + 316220553.*std::pow(z/r,6.) - 3162205530.*std::pow(z/r,8.) + 18200249606.*std::pow(z/r,10.) - 64528157694.*std::pow(z/r,12.) + 145365629970.*std::pow(z/r,14.) - 208357402957.*std::pow(z/r,16.) + 183844767315.*std::pow(z/r,18.) - 90954779619.*std::pow(z/r,20.) + 19293438101.*std::pow(z/r,22.))*std::pow((x-I*y)/r,3.))/(8.388608e6*1.);
}
else if (m == -2){
return (15.*std::sqrt(663./Pi)*(-156009.*z/r + 16588957.*std::pow(z/r,3.) - 514257667.*std::pow(z/r,5.) + 7273072719.*std::pow(z/r,7.) - 56568343370.*std::pow(z/r,9.) + 266385471506.*std::pow(z/r,11.) - 799156414518.*std::pow(z/r,13.) + 1560257761678.*std::pow(z/r,15.) - 1973267169181.*std::pow(z/r,17.) + 1557842501985.*std::pow(z/r,19.) - 697319977079.*std::pow(z/r,21.) + 135054066707.*std::pow(z/r,23.))*std::pow((x-I*y)/r,2.))/(4.194304e6*1.);
}
else if (m == -1){
return (5.*std::sqrt(663./(2.*Pi))*(52003. - 16848972.*std::pow(z/r,2.) + 895803678.*std::pow(z/r,4.) - 18513276012.*std::pow(z/r,6.) + 196372963413.*std::pow(z/r,8.) - 1221876216792.*std::pow(z/r,10.) + 4794938487108.*std::pow(z/r,12.) - 12329841823992.*std::pow(z/r,14.) + 21063479782653.*std::pow(z/r,16.) - 23679206030172.*std::pow(z/r,18.) + 16824699021438.*std::pow(z/r,20.) - 6846414320412.*std::pow(z/r,22.) + 1215486600363.*std::pow(z/r,24.))*(x-I*y)/r)/(8.388608e6*1.);
}
else if (m == 0){
return (std::sqrt(51./Pi)*(16900975.*z/r - 1825305300.*std::pow(z/r,3.) + 58227239070.*std::pow(z/r,5.) - 859544957700.*std::pow(z/r,7.) + 7091245901025.*std::pow(z/r,9.) - 36100888223400.*std::pow(z/r,11.) + 119873462177700.*std::pow(z/r,13.) - 267146572853160.*std::pow(z/r,15.) + 402684172315425.*std::pow(z/r,17.) - 405039050516100.*std::pow(z/r,19.) + 260382246760350.*std::pow(z/r,21.) - 96742811049300.*std::pow(z/r,23.) + 15801325804719.*std::pow(z/r,25.)))/8.388608e6;
}
else if (m == 1){
return (-5.*1.*std::sqrt(663./(2.*Pi))*(52003. - 16848972.*std::pow(z/r,2.) + 895803678.*std::pow(z/r,4.) - 18513276012.*std::pow(z/r,6.) + 196372963413.*std::pow(z/r,8.) - 1221876216792.*std::pow(z/r,10.) + 4794938487108.*std::pow(z/r,12.) - 12329841823992.*std::pow(z/r,14.) + 21063479782653.*std::pow(z/r,16.) - 23679206030172.*std::pow(z/r,18.) + 16824699021438.*std::pow(z/r,20.) - 6846414320412.*std::pow(z/r,22.) + 1215486600363.*std::pow(z/r,24.))*(x+I*y)/r)/8.388608e6;
}
else if (m == 2){
return (15.*1.*std::sqrt(663./Pi)*(-156009.*z/r + 16588957.*std::pow(z/r,3.) - 514257667.*std::pow(z/r,5.) + 7273072719.*std::pow(z/r,7.) - 56568343370.*std::pow(z/r,9.) + 266385471506.*std::pow(z/r,11.) - 799156414518.*std::pow(z/r,13.) + 1560257761678.*std::pow(z/r,15.) - 1973267169181.*std::pow(z/r,17.) + 1557842501985.*std::pow(z/r,19.) - 697319977079.*std::pow(z/r,21.) + 135054066707.*std::pow(z/r,23.))*std::pow((x+I*y)/r,2.))/4.194304e6;
}
else if (m == 3){
return (-15.*1.*std::sqrt(106743./Pi)*(-969. + 309111.*std::pow(z/r,2.) - 15970735.*std::pow(z/r,4.) + 316220553.*std::pow(z/r,6.) - 3162205530.*std::pow(z/r,8.) + 18200249606.*std::pow(z/r,10.) - 64528157694.*std::pow(z/r,12.) + 145365629970.*std::pow(z/r,14.) - 208357402957.*std::pow(z/r,16.) + 183844767315.*std::pow(z/r,18.) - 90954779619.*std::pow(z/r,20.) + 19293438101.*std::pow(z/r,22.))*std::pow((x+I*y)/r,3.))/8.388608e6;
}
else if (m == 4){
return (15.*1.*std::sqrt(34051017./(2.*Pi))*(969.*z/r - 100130.*std::pow(z/r,3.) + 2973861.*std::pow(z/r,5.) - 39651480.*std::pow(z/r,7.) + 285270370.*std::pow(z/r,9.) - 1213695756.*std::pow(z/r,11.) + 3189841410.*std::pow(z/r,13.) - 5225264024.*std::pow(z/r,15.) + 5186842965.*std::pow(z/r,17.) - 2851247010.*std::pow(z/r,19.) + 665290969.*std::pow(z/r,21.))*std::pow((x+I*y)/r,4.))/4.194304e6;
}
else if (m == 5){
return (-3.*1.*std::sqrt(24322155./Pi)*(323. - 100130.*std::pow(z/r,2.) + 4956435.*std::pow(z/r,4.) - 92520120.*std::pow(z/r,6.) + 855811110.*std::pow(z/r,8.) - 4450217772.*std::pow(z/r,10.) + 13822646110.*std::pow(z/r,12.) - 26126320120.*std::pow(z/r,14.) + 29392110135.*std::pow(z/r,16.) - 18057897730.*std::pow(z/r,18.) + 4657036783.*std::pow(z/r,20.))*std::pow((x+I*y)/r,5.))/8.388608e6;
}
else if (m == 6){
return (15.*1.*std::sqrt(150797361./Pi)*z/r*(-323. + 31977.*std::pow(z/r,2.) - 895356.*std::pow(z/r,4.) + 11042724.*std::pow(z/r,6.) - 71777706.*std::pow(z/r,8.) + 267535086.*std::pow(z/r,10.) - 589949164.*std::pow(z/r,12.) + 758506068.*std::pow(z/r,14.) - 524261547.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x+I*y)/r,6.))/4.194304e6;
}
else if (m == 7){
return (-15.*1.*std::sqrt(2865149859./(2.*Pi))*(-17. + 5049.*std::pow(z/r,2.) - 235620.*std::pow(z/r,4.) + 4068372.*std::pow(z/r,6.) - 33999966.*std::pow(z/r,8.) + 154888734.*std::pow(z/r,10.) - 403649428.*std::pow(z/r,12.) + 598820580.*std::pow(z/r,14.) - 469076121.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x+I*y)/r,7.))/1.6777216e7;
}
else if (m == 8){
return (15.*1.*std::sqrt(86822723./Pi)*z/r*(1683. - 157080.*std::pow(z/r,2.) + 4068372.*std::pow(z/r,4.) - 45333288.*std::pow(z/r,6.) + 258147890.*std::pow(z/r,8.) - 807298856.*std::pow(z/r,10.) + 1397248020.*std::pow(z/r,12.) - 1250869656.*std::pow(z/r,14.) + 450680979.*std::pow(z/r,16.))*std::pow((x+I*y)/r,8.))/1.6777216e7;
}
else if (m == 9){
return (-15.*1.*std::sqrt(86822723./(2.*Pi))*(99. - 27720.*std::pow(z/r,2.) + 1196580.*std::pow(z/r,4.) - 18666648.*std::pow(z/r,6.) + 136666530.*std::pow(z/r,8.) - 522369848.*std::pow(z/r,10.) + 1068483780.*std::pow(z/r,12.) - 1103708520.*std::pow(z/r,14.) + 450680979.*std::pow(z/r,16.))*std::pow((x+I*y)/r,9.))/1.6777216e7;
}
else if (m == 10){
return (3.*1.*std::sqrt(3038795305./(2.*Pi))*z/r*(-495. + 42735.*std::pow(z/r,2.) - 999999.*std::pow(z/r,4.) + 9761895.*std::pow(z/r,6.) - 46640165.*std::pow(z/r,8.) + 114480405.*std::pow(z/r,10.) - 137963565.*std::pow(z/r,12.) + 64382997.*std::pow(z/r,14.))*std::pow((x+I*y)/r,10.))/4.194304e6;
}
else if (m == 11){
return (-15.*1.*std::sqrt(1823277183./(2.*Pi))*(-11. + 2849.*std::pow(z/r,2.) - 111111.*std::pow(z/r,4.) + 1518517.*std::pow(z/r,6.) - 9328033.*std::pow(z/r,8.) + 27984099.*std::pow(z/r,10.) - 39856141.*std::pow(z/r,12.) + 21460999.*std::pow(z/r,14.))*std::pow((x+I*y)/r,11.))/8.388608e6;
}
else if (m == 12){
return (15.*1.*std::sqrt(9637322253./Pi)*z/r*(77. - 6006.*std::pow(z/r,2.) + 123123.*std::pow(z/r,4.) - 1008436.*std::pow(z/r,6.) + 3781635.*std::pow(z/r,8.) - 6463158.*std::pow(z/r,10.) + 4060189.*std::pow(z/r,12.))*std::pow((x+I*y)/r,12.))/8.388608e6;
}
else if (m == 13){
return (-15.*1.*std::sqrt(39017499./(2.*Pi))*(77. - 18018.*std::pow(z/r,2.) + 615615.*std::pow(z/r,4.) - 7059052.*std::pow(z/r,6.) + 34034715.*std::pow(z/r,8.) - 71094738.*std::pow(z/r,10.) + 52782457.*std::pow(z/r,12.))*std::pow((x+I*y)/r,13.))/8.388608e6;
}
else if (m == 14){
return (15.*1.*std::sqrt(507227487./(2.*Pi))*z/r*(-231. + 15785.*std::pow(z/r,2.) - 271502.*std::pow(z/r,4.) + 1745370.*std::pow(z/r,6.) - 4557355.*std::pow(z/r,8.) + 4060189.*std::pow(z/r,10.))*std::pow((x+I*y)/r,14.))/4.194304e6;
}
else if (m == 15){
return (-3.*1.*std::sqrt(27897511785./Pi)*(-21. + 4305.*std::pow(z/r,2.) - 123410.*std::pow(z/r,4.) + 1110690.*std::pow(z/r,6.) - 3728745.*std::pow(z/r,8.) + 4060189.*std::pow(z/r,10.))*std::pow((x+I*y)/r,15.))/1.6777216e7;
}
else if (m == 16){
return (15.*1.*std::sqrt(228759596637./(2.*Pi))*z/r*(21. - 1204.*std::pow(z/r,2.) + 16254.*std::pow(z/r,4.) - 72756.*std::pow(z/r,6.) + 99029.*std::pow(z/r,8.))*std::pow((x+I*y)/r,16.))/8.388608e6;
}
else if (m == 17){
return (-15.*1.*std::sqrt(533772392153./Pi)*(1. - 172.*std::pow(z/r,2.) + 3870.*std::pow(z/r,4.) - 24252.*std::pow(z/r,6.) + 42441.*std::pow(z/r,8.))*std::pow((x+I*y)/r,17.))/1.6777216e7;
}
else if (m == 18){
return (15.*1.*std::sqrt(22952212862579./(2.*Pi))*z/r*(-1. + 45.*std::pow(z/r,2.) - 423.*std::pow(z/r,4.) + 987.*std::pow(z/r,6.))*std::pow((x+I*y)/r,18.))/4.194304e6;
}
else if (m == 19){
return (-15.*1.*std::sqrt(298080686527./(2.*Pi))*(-1. + 135.*std::pow(z/r,2.) - 2115.*std::pow(z/r,4.) + 6909.*std::pow(z/r,6.))*std::pow((x+I*y)/r,19.))/8.388608e6;
}
else if (m == 20){
return (3.*1.*std::sqrt(4471210297905./Pi)*z/r*(15. - 470.*std::pow(z/r,2.) + 2303.*std::pow(z/r,4.))*std::pow((x+I*y)/r,20.))/8.388608e6;
}
else if (m == 21){
return (-15.*1.*std::sqrt(38880089547./(2.*Pi))*(3. - 282.*std::pow(z/r,2.) + 2303.*std::pow(z/r,4.))*std::pow((x+I*y)/r,21.))/8.388608e6;
}
else if (m == 22){
return (15.*1.*std::sqrt(1827364208709./(2.*Pi))*z/r*(-3. + 49.*std::pow(z/r,2.))*std::pow((x+I*y)/r,22.))/4.194304e6;
}
else if (m == 23){
return (-15.*1.*std::sqrt(1827364208709./(2.*Pi))*(-1. + 49.*std::pow(z/r,2.))*std::pow((x+I*y)/r,23.))/1.6777216e7;
}
else if (m == 24){
return (105.*1.*std::sqrt(1827364208709./Pi)*z/r*std::pow((x+I*y)/r,24.))/1.6777216e7;
}
else if (m == 25){
return (-21.*1.*std::sqrt(1827364208709./(2.*Pi))*std::pow((x+I*y)/r,25.))/1.6777216e7;
}
else{return 0.;}
}

else if (l == 26){
if(m == -26){
return (21.*std::sqrt(7450023312429./(2.*Pi))*std::pow((x-I*y)/r,26.))/(3.3554432e7*1.);
}
else if (m == -25){
return (21.*std::sqrt(96850303061577./(2.*Pi))*z/r*std::pow((x-I*y)/r,25.))/(1.6777216e7*1.);
}
else if (m == -24){
return (21.*std::sqrt(1899025550227./Pi)*(-1. + 51.*std::pow(z/r,2.))*std::pow((x-I*y)/r,24.))/(3.3554432e7*1.);
}
else if (m == -23){
return (105.*std::sqrt(5697076650681./(2.*Pi))*z/r*(-1. + 17.*std::pow(z/r,2.))*std::pow((x-I*y)/r,23.))/(1.6777216e7*1.);
}
else if (m == -22){
return (15.*std::sqrt(5697076650681./(2.*Pi))*(1. - 98.*std::pow(z/r,2.) + 833.*std::pow(z/r,4.))*std::pow((x-I*y)/r,22.))/(3.3554432e7*1.);
}
else if (m == -21){
return (3.*std::sqrt(9495127751135./(2.*Pi))*z/r*(15. - 490.*std::pow(z/r,2.) + 2499.*std::pow(z/r,4.))*std::pow((x-I*y)/r,21.))/(8.388608e6*1.);
}
else if (m == -20){
return (3.*std::sqrt(606071984115./Pi)*(-5. + 705.*std::pow(z/r,2.) - 11515.*std::pow(z/r,4.) + 39151.*std::pow(z/r,6.))*std::pow((x-I*y)/r,20.))/(1.6777216e7*1.);
}
else if (m == -19){
return (3.*std::sqrt(97577589442515./(2.*Pi))*z/r*(-5. + 235.*std::pow(z/r,2.) - 2303.*std::pow(z/r,4.) + 5593.*std::pow(z/r,6.))*std::pow((x-I*y)/r,19.))/(8.388608e6*1.);
}
else if (m == -18){
return (5.*std::sqrt(19515517888503./Pi)*(1. - 180.*std::pow(z/r,2.) + 4230.*std::pow(z/r,4.) - 27636.*std::pow(z/r,6.) + 50337.*std::pow(z/r,8.))*std::pow((x-I*y)/r,18.))/(3.3554432e7*1.);
}
else if (m == -17){
return (15.*std::sqrt(214670696773533./Pi)*z/r*(1. - 60.*std::pow(z/r,2.) + 846.*std::pow(z/r,4.) - 3948.*std::pow(z/r,6.) + 5593.*std::pow(z/r,8.))*std::pow((x-I*y)/r,17.))/(1.6777216e7*1.);
}
else if (m == -16){
return (3.*std::sqrt(24961708927155./(2.*Pi))*(-1. + 215.*std::pow(z/r,2.) - 6450.*std::pow(z/r,4.) + 60630.*std::pow(z/r,6.) - 212205.*std::pow(z/r,8.) + 240499.*std::pow(z/r,10.))*std::pow((x-I*y)/r,16.))/(1.6777216e7*1.);
}
else if (m == -15){
return (3.*std::sqrt(108059346005./Pi)*z/r*(-231. + 16555.*std::pow(z/r,2.) - 297990.*std::pow(z/r,4.) + 2000790.*std::pow(z/r,6.) - 5446595.*std::pow(z/r,8.) + 5050479.*std::pow(z/r,10.))*std::pow((x-I*y)/r,15.))/(1.6777216e7*1.);
}
else if (m == -14){
return (3.*std::sqrt(7906781415./Pi)*(77. - 18942.*std::pow(z/r,2.) + 678755.*std::pow(z/r,4.) - 8145060.*std::pow(z/r,6.) + 41016195.*std::pow(z/r,8.) - 89324158.*std::pow(z/r,10.) + 69023213.*std::pow(z/r,12.))*std::pow((x-I*y)/r,14.))/(3.3554432e7*1.);
}
else if (m == -13){
return (15.*std::sqrt(121642791./(2.*Pi))*z/r*(1001. - 82082.*std::pow(z/r,2.) + 1764763.*std::pow(z/r,4.) - 15126540.*std::pow(z/r,6.) + 59245615.*std::pow(z/r,8.) - 105564914.*std::pow(z/r,10.) + 69023213.*std::pow(z/r,12.))*std::pow((x-I*y)/r,13.))/(8.388608e6*1.);
}
else if (m == -12){
return (15.*std::sqrt(3689831327./Pi)*(-11. + 3003.*std::pow(z/r,2.) - 123123.*std::pow(z/r,4.) + 1764763.*std::pow(z/r,6.) - 11344905.*std::pow(z/r,8.) + 35547369.*std::pow(z/r,10.) - 52782457.*std::pow(z/r,12.) + 29581377.*std::pow(z/r,14.))*std::pow((x-I*y)/r,12.))/(1.6777216e7*1.);
}
else if (m == -11){
return (3.*std::sqrt(1051601928195./(2.*Pi))*z/r*(-55. + 5005.*std::pow(z/r,2.) - 123123.*std::pow(z/r,4.) + 1260545.*std::pow(z/r,6.) - 6302725.*std::pow(z/r,8.) + 16157895.*std::pow(z/r,10.) - 20300945.*std::pow(z/r,12.) + 9860459.*std::pow(z/r,14.))*std::pow((x-I*y)/r,11.))/(8.388608e6*1.);
}
else if (m == -10){
return (3.*std::sqrt(28421673735./(2.*Pi))*(55. - 16280.*std::pow(z/r,2.) + 740740.*std::pow(z/r,4.) - 12148136.*std::pow(z/r,6.) + 93280330.*std::pow(z/r,8.) - 373121320.*std::pow(z/r,10.) + 797122820.*std::pow(z/r,12.) - 858439960.*std::pow(z/r,14.) + 364836983.*std::pow(z/r,16.))*std::pow((x-I*y)/r,10.))/(3.3554432e7*1.);
}
else if (m == -9){
return (std::sqrt(483168453495./(2.*Pi))*z/r*(495. - 48840.*std::pow(z/r,2.) + 1333332.*std::pow(z/r,4.) - 15619032.*std::pow(z/r,6.) + 93280330.*std::pow(z/r,8.) - 305281080.*std::pow(z/r,10.) + 551854260.*std::pow(z/r,12.) - 515063976.*std::pow(z/r,14.) + 193148991.*std::pow(z/r,16.))*std::pow((x-I*y)/r,9.))/(1.6777216e7*1.);
}
else if (m == -8){
return (15.*std::sqrt(13804812957./Pi)*(-11. + 3465.*std::pow(z/r,2.) - 170940.*std::pow(z/r,4.) + 3111108.*std::pow(z/r,6.) - 27333306.*std::pow(z/r,8.) + 130592462.*std::pow(z/r,10.) - 356161260.*std::pow(z/r,12.) + 551854260.*std::pow(z/r,14.) - 450680979.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x-I*y)/r,8.))/(3.3554432e7*1.);
}
else if (m == -7){
return (15.*std::sqrt(42739359./(2.*Pi))*z/r*(-3553. + 373065.*std::pow(z/r,2.) - 11042724.*std::pow(z/r,4.) + 143555412.*std::pow(z/r,6.) - 980961982.*std::pow(z/r,8.) + 3834669566.*std::pow(z/r,10.) - 8849237460.*std::pow(z/r,12.) + 11883261732.*std::pow(z/r,14.) - 8562938601.*std::pow(z/r,16.) + 2553858881.*std::pow(z/r,18.))*std::pow((x-I*y)/r,7.))/(1.6777216e7*1.);
}
else if (m == -6){
return (3.*std::sqrt(783554915./(2.*Pi))*(323. - 106590.*std::pow(z/r,2.) + 5595975.*std::pow(z/r,4.) - 110427240.*std::pow(z/r,6.) + 1076665590.*std::pow(z/r,8.) - 5885771892.*std::pow(z/r,10.) + 19173347830.*std::pow(z/r,12.) - 37925303400.*std::pow(z/r,14.) + 44562231495.*std::pow(z/r,16.) - 28543128670.*std::pow(z/r,18.) + 7661576643.*std::pow(z/r,20.))*std::pow((x-I*y)/r,6.))/(3.3554432e7*1.);
}
else if (m == -5){
return (3.*std::sqrt(16454653215./Pi)*(323.*z/r - 35530.*std::pow(z/r,3.) + 1119195.*std::pow(z/r,5.) - 15775320.*std::pow(z/r,7.) + 119629510.*std::pow(z/r,9.) - 535070172.*std::pow(z/r,11.) + 1474872910.*std::pow(z/r,13.) - 2528353560.*std::pow(z/r,15.) + 2621307735.*std::pow(z/r,17.) - 1502269930.*std::pow(z/r,19.) + 364836983.*std::pow(z/r,21.))*std::pow((x-I*y)/r,5.))/(8.388608e6*1.);
}
else if (m == -4){
return (3.*std::sqrt(48254115./(2.*Pi))*(-323. + 110143.*std::pow(z/r,2.) - 6057865.*std::pow(z/r,4.) + 127215165.*std::pow(z/r,6.) - 1344846030.*std::pow(z/r,8.) + 8158732582.*std::pow(z/r,10.) - 30409821442.*std::pow(z/r,12.) + 71847380330.*std::pow(z/r,14.) - 107771070495.*std::pow(z/r,16.) + 99318437515.*std::pow(z/r,18.) - 51227404613.*std::pow(z/r,20.) + 11309946473.*std::pow(z/r,22.))*std::pow((x-I*y)/r,4.))/(8.388608e6*1.);
}
else if (m == -3){
return (15.*std::sqrt(139867./Pi)*(-22287.*z/r + 2533289.*std::pow(z/r,3.) - 83598537.*std::pow(z/r,5.) + 1253978055.*std::pow(z/r,7.) - 10310486230.*std::pow(z/r,9.) + 51177504378.*std::pow(z/r,11.) - 161405975346.*std::pow(z/r,13.) + 330497949518.*std::pow(z/r,15.) - 437423756715.*std::pow(z/r,17.) + 360682746765.*std::pow(z/r,19.) - 168318615157.*std::pow(z/r,21.) + 33929839419.*std::pow(z/r,23.))*std::pow((x-I*y)/r,3.))/(8.388608e6*1.);
}
else if (m == -2){
return (15.*std::sqrt(14469./(2.*Pi))*(7429. - 2585292.*std::pow(z/r,2.) + 146930762.*std::pow(z/r,4.) - 3232476764.*std::pow(z/r,6.) + 36365363595.*std::pow(z/r,8.) - 239203280536.*std::pow(z/r,10.) + 989431751308.*std::pow(z/r,12.) - 2674727591448.*std::pow(z/r,14.) + 4792220268011.*std::pow(z/r,16.) - 5637906197660.*std::pow(z/r,18.) + 4183919862474.*std::pow(z/r,20.) - 1774996305292.*std::pow(z/r,22.) + 327988447717.*std::pow(z/r,24.))*std::pow((x-I*y)/r,2.))/(1.6777216e7*1.);
}
else if (m == -1){
return (3.*std::sqrt(2067./(2.*Pi))*(1300075.*z/r - 150808700.*std::pow(z/r,3.) + 5142576670.*std::pow(z/r,5.) - 80811919100.*std::pow(z/r,7.) + 707104292125.*std::pow(z/r,9.) - 3805506735800.*std::pow(z/r,11.) + 13319273575300.*std::pow(z/r,13.) - 31205155233560.*std::pow(z/r,15.) + 49331679229525.*std::pow(z/r,17.) - 51928083399500.*std::pow(z/r,19.) + 34865998853950.*std::pow(z/r,21.) - 13505406670700.*std::pow(z/r,23.) + 2295919134019.*std::pow(z/r,25.))*(x-I*y)/r)/(8.388608e6*1.);
}
else if (m == 0){
return (std::sqrt(53./Pi)*(-1300075. + 456326325.*std::pow(z/r,2.) - 26466926850.*std::pow(z/r,4.) + 601681470390.*std::pow(z/r,6.) - 7091245901025.*std::pow(z/r,8.) + 49638721307175.*std::pow(z/r,10.) - 222622144044300.*std::pow(z/r,12.) + 667866432132900.*std::pow(z/r,14.) - 1369126185872445.*std::pow(z/r,16.) + 1923935489951475.*std::pow(z/r,18.) - 1822675727322450.*std::pow(z/r,20.) + 1112542327066950.*std::pow(z/r,22.) - 395033145117975.*std::pow(z/r,24.) + 61989816618513.*std::pow(z/r,26.)))/1.6777216e7;
}
else if (m == 1){
return (-3.*1.*std::sqrt(2067./(2.*Pi))*(1300075.*z/r - 150808700.*std::pow(z/r,3.) + 5142576670.*std::pow(z/r,5.) - 80811919100.*std::pow(z/r,7.) + 707104292125.*std::pow(z/r,9.) - 3805506735800.*std::pow(z/r,11.) + 13319273575300.*std::pow(z/r,13.) - 31205155233560.*std::pow(z/r,15.) + 49331679229525.*std::pow(z/r,17.) - 51928083399500.*std::pow(z/r,19.) + 34865998853950.*std::pow(z/r,21.) - 13505406670700.*std::pow(z/r,23.) + 2295919134019.*std::pow(z/r,25.))*(x+I*y)/r)/8.388608e6;
}
else if (m == 2){
return (15.*1.*std::sqrt(14469./(2.*Pi))*(7429. - 2585292.*std::pow(z/r,2.) + 146930762.*std::pow(z/r,4.) - 3232476764.*std::pow(z/r,6.) + 36365363595.*std::pow(z/r,8.) - 239203280536.*std::pow(z/r,10.) + 989431751308.*std::pow(z/r,12.) - 2674727591448.*std::pow(z/r,14.) + 4792220268011.*std::pow(z/r,16.) - 5637906197660.*std::pow(z/r,18.) + 4183919862474.*std::pow(z/r,20.) - 1774996305292.*std::pow(z/r,22.) + 327988447717.*std::pow(z/r,24.))*std::pow((x+I*y)/r,2.))/1.6777216e7;
}
else if (m == 3){
return (-15.*1.*std::sqrt(139867./Pi)*(-22287.*z/r + 2533289.*std::pow(z/r,3.) - 83598537.*std::pow(z/r,5.) + 1253978055.*std::pow(z/r,7.) - 10310486230.*std::pow(z/r,9.) + 51177504378.*std::pow(z/r,11.) - 161405975346.*std::pow(z/r,13.) + 330497949518.*std::pow(z/r,15.) - 437423756715.*std::pow(z/r,17.) + 360682746765.*std::pow(z/r,19.) - 168318615157.*std::pow(z/r,21.) + 33929839419.*std::pow(z/r,23.))*std::pow((x+I*y)/r,3.))/8.388608e6;
}
else if (m == 4){
return (3.*1.*std::sqrt(48254115./(2.*Pi))*(-323. + 110143.*std::pow(z/r,2.) - 6057865.*std::pow(z/r,4.) + 127215165.*std::pow(z/r,6.) - 1344846030.*std::pow(z/r,8.) + 8158732582.*std::pow(z/r,10.) - 30409821442.*std::pow(z/r,12.) + 71847380330.*std::pow(z/r,14.) - 107771070495.*std::pow(z/r,16.) + 99318437515.*std::pow(z/r,18.) - 51227404613.*std::pow(z/r,20.) + 11309946473.*std::pow(z/r,22.))*std::pow((x+I*y)/r,4.))/8.388608e6;
}
else if (m == 5){
return (-3.*1.*std::sqrt(16454653215./Pi)*(323.*z/r - 35530.*std::pow(z/r,3.) + 1119195.*std::pow(z/r,5.) - 15775320.*std::pow(z/r,7.) + 119629510.*std::pow(z/r,9.) - 535070172.*std::pow(z/r,11.) + 1474872910.*std::pow(z/r,13.) - 2528353560.*std::pow(z/r,15.) + 2621307735.*std::pow(z/r,17.) - 1502269930.*std::pow(z/r,19.) + 364836983.*std::pow(z/r,21.))*std::pow((x+I*y)/r,5.))/8.388608e6;
}
else if (m == 6){
return (3.*1.*std::sqrt(783554915./(2.*Pi))*(323. - 106590.*std::pow(z/r,2.) + 5595975.*std::pow(z/r,4.) - 110427240.*std::pow(z/r,6.) + 1076665590.*std::pow(z/r,8.) - 5885771892.*std::pow(z/r,10.) + 19173347830.*std::pow(z/r,12.) - 37925303400.*std::pow(z/r,14.) + 44562231495.*std::pow(z/r,16.) - 28543128670.*std::pow(z/r,18.) + 7661576643.*std::pow(z/r,20.))*std::pow((x+I*y)/r,6.))/3.3554432e7;
}
else if (m == 7){
return (-15.*1.*std::sqrt(42739359./(2.*Pi))*z/r*(-3553. + 373065.*std::pow(z/r,2.) - 11042724.*std::pow(z/r,4.) + 143555412.*std::pow(z/r,6.) - 980961982.*std::pow(z/r,8.) + 3834669566.*std::pow(z/r,10.) - 8849237460.*std::pow(z/r,12.) + 11883261732.*std::pow(z/r,14.) - 8562938601.*std::pow(z/r,16.) + 2553858881.*std::pow(z/r,18.))*std::pow((x+I*y)/r,7.))/1.6777216e7;
}
else if (m == 8){
return (15.*1.*std::sqrt(13804812957./Pi)*(-11. + 3465.*std::pow(z/r,2.) - 170940.*std::pow(z/r,4.) + 3111108.*std::pow(z/r,6.) - 27333306.*std::pow(z/r,8.) + 130592462.*std::pow(z/r,10.) - 356161260.*std::pow(z/r,12.) + 551854260.*std::pow(z/r,14.) - 450680979.*std::pow(z/r,16.) + 150226993.*std::pow(z/r,18.))*std::pow((x+I*y)/r,8.))/3.3554432e7;
}
else if (m == 9){
return -5.960464477539063e-8*(1.*std::sqrt(483168453495./(2.*Pi))*z/r*(495. - 48840.*std::pow(z/r,2.) + 1333332.*std::pow(z/r,4.) - 15619032.*std::pow(z/r,6.) + 93280330.*std::pow(z/r,8.) - 305281080.*std::pow(z/r,10.) + 551854260.*std::pow(z/r,12.) - 515063976.*std::pow(z/r,14.) + 193148991.*std::pow(z/r,16.))*std::pow((x+I*y)/r,9.));
}
else if (m == 10){
return (3.*1.*std::sqrt(28421673735./(2.*Pi))*(55. - 16280.*std::pow(z/r,2.) + 740740.*std::pow(z/r,4.) - 12148136.*std::pow(z/r,6.) + 93280330.*std::pow(z/r,8.) - 373121320.*std::pow(z/r,10.) + 797122820.*std::pow(z/r,12.) - 858439960.*std::pow(z/r,14.) + 364836983.*std::pow(z/r,16.))*std::pow((x+I*y)/r,10.))/3.3554432e7;
}
else if (m == 11){
return (-3.*1.*std::sqrt(1051601928195./(2.*Pi))*z/r*(-55. + 5005.*std::pow(z/r,2.) - 123123.*std::pow(z/r,4.) + 1260545.*std::pow(z/r,6.) - 6302725.*std::pow(z/r,8.) + 16157895.*std::pow(z/r,10.) - 20300945.*std::pow(z/r,12.) + 9860459.*std::pow(z/r,14.))*std::pow((x+I*y)/r,11.))/8.388608e6;
}
else if (m == 12){
return (15.*1.*std::sqrt(3689831327./Pi)*(-11. + 3003.*std::pow(z/r,2.) - 123123.*std::pow(z/r,4.) + 1764763.*std::pow(z/r,6.) - 11344905.*std::pow(z/r,8.) + 35547369.*std::pow(z/r,10.) - 52782457.*std::pow(z/r,12.) + 29581377.*std::pow(z/r,14.))*std::pow((x+I*y)/r,12.))/1.6777216e7;
}
else if (m == 13){
return (-15.*1.*std::sqrt(121642791./(2.*Pi))*z/r*(1001. - 82082.*std::pow(z/r,2.) + 1764763.*std::pow(z/r,4.) - 15126540.*std::pow(z/r,6.) + 59245615.*std::pow(z/r,8.) - 105564914.*std::pow(z/r,10.) + 69023213.*std::pow(z/r,12.))*std::pow((x+I*y)/r,13.))/8.388608e6;
}
else if (m == 14){
return (3.*1.*std::sqrt(7906781415./Pi)*(77. - 18942.*std::pow(z/r,2.) + 678755.*std::pow(z/r,4.) - 8145060.*std::pow(z/r,6.) + 41016195.*std::pow(z/r,8.) - 89324158.*std::pow(z/r,10.) + 69023213.*std::pow(z/r,12.))*std::pow((x+I*y)/r,14.))/3.3554432e7;
}
else if (m == 15){
return (-3.*1.*std::sqrt(108059346005./Pi)*z/r*(-231. + 16555.*std::pow(z/r,2.) - 297990.*std::pow(z/r,4.) + 2000790.*std::pow(z/r,6.) - 5446595.*std::pow(z/r,8.) + 5050479.*std::pow(z/r,10.))*std::pow((x+I*y)/r,15.))/1.6777216e7;
}
else if (m == 16){
return (3.*1.*std::sqrt(24961708927155./(2.*Pi))*(-1. + 215.*std::pow(z/r,2.) - 6450.*std::pow(z/r,4.) + 60630.*std::pow(z/r,6.) - 212205.*std::pow(z/r,8.) + 240499.*std::pow(z/r,10.))*std::pow((x+I*y)/r,16.))/1.6777216e7;
}
else if (m == 17){
return (-15.*1.*std::sqrt(214670696773533./Pi)*z/r*(1. - 60.*std::pow(z/r,2.) + 846.*std::pow(z/r,4.) - 3948.*std::pow(z/r,6.) + 5593.*std::pow(z/r,8.))*std::pow((x+I*y)/r,17.))/1.6777216e7;
}
else if (m == 18){
return (5.*1.*std::sqrt(19515517888503./Pi)*(1. - 180.*std::pow(z/r,2.) + 4230.*std::pow(z/r,4.) - 27636.*std::pow(z/r,6.) + 50337.*std::pow(z/r,8.))*std::pow((x+I*y)/r,18.))/3.3554432e7;
}
else if (m == 19){
return (-3.*1.*std::sqrt(97577589442515./(2.*Pi))*z/r*(-5. + 235.*std::pow(z/r,2.) - 2303.*std::pow(z/r,4.) + 5593.*std::pow(z/r,6.))*std::pow((x+I*y)/r,19.))/8.388608e6;
}
else if (m == 20){
return (3.*1.*std::sqrt(606071984115./Pi)*(-5. + 705.*std::pow(z/r,2.) - 11515.*std::pow(z/r,4.) + 39151.*std::pow(z/r,6.))*std::pow((x+I*y)/r,20.))/1.6777216e7;
}
else if (m == 21){
return (-3.*1.*std::sqrt(9495127751135./(2.*Pi))*z/r*(15. - 490.*std::pow(z/r,2.) + 2499.*std::pow(z/r,4.))*std::pow((x+I*y)/r,21.))/8.388608e6;
}
else if (m == 22){
return (15.*1.*std::sqrt(5697076650681./(2.*Pi))*(1. - 98.*std::pow(z/r,2.) + 833.*std::pow(z/r,4.))*std::pow((x+I*y)/r,22.))/3.3554432e7;
}
else if (m == 23){
return (-105.*1.*std::sqrt(5697076650681./(2.*Pi))*z/r*(-1. + 17.*std::pow(z/r,2.))*std::pow((x+I*y)/r,23.))/1.6777216e7;
}
else if (m == 24){
return (21.*1.*std::sqrt(1899025550227./Pi)*(-1. + 51.*std::pow(z/r,2.))*std::pow((x+I*y)/r,24.))/3.3554432e7;
}
else if (m == 25){
return (-21.*1.*std::sqrt(96850303061577./(2.*Pi))*z/r*std::pow((x+I*y)/r,25.))/1.6777216e7;
}
else if (m == 26){
return (21.*1.*std::sqrt(7450023312429./(2.*Pi))*std::pow((x+I*y)/r,26.))/3.3554432e7;
}
else{return 0.;}
}

else if (l == 27){
if(m == -27){
return (7.*std::sqrt(136583760727865./Pi)*std::pow((x-I*y)/r,27.))/(6.7108864e7*1.);
}
else if (m == -26){
return (21.*std::sqrt(409751282183595./(2.*Pi))*z/r*std::pow((x-I*y)/r,26.))/(3.3554432e7*1.);
}
else if (m == -25){
return (21.*std::sqrt(7731156267615./Pi)*(-1. + 53.*std::pow(z/r,2.))*std::pow((x-I*y)/r,25.))/(6.7108864e7*1.);
}
else if (m == -24){
return (21.*std::sqrt(33501677159665./Pi)*z/r*(-3. + 53.*std::pow(z/r,2.))*std::pow((x-I*y)/r,24.))/(3.3554432e7*1.);
}
else if (m == -23){
return (21.*std::sqrt(5912060675235./Pi)*(1. - 102.*std::pow(z/r,2.) + 901.*std::pow(z/r,4.))*std::pow((x-I*y)/r,23.))/(6.7108864e7*1.);
}
else if (m == -22){
return (105.*std::sqrt(1182412135047./(2.*Pi))*z/r*(5. - 170.*std::pow(z/r,2.) + 901.*std::pow(z/r,4.))*std::pow((x-I*y)/r,22.))/(3.3554432e7*1.);
}
else if (m == -21){
return (15.*std::sqrt(394137378349./Pi)*(-5. + 735.*std::pow(z/r,2.) - 12495.*std::pow(z/r,4.) + 44149.*std::pow(z/r,6.))*std::pow((x-I*y)/r,21.))/(6.7108864e7*1.);
}
else if (m == -20){
return (15.*std::sqrt(8276884945329./Pi)*z/r*(-5. + 245.*std::pow(z/r,2.) - 2499.*std::pow(z/r,4.) + 6307.*std::pow(z/r,6.))*std::pow((x-I*y)/r,20.))/(1.6777216e7*1.);
}
else if (m == -19){
return (15.*std::sqrt(176103935007./(2.*Pi))*(5. - 940.*std::pow(z/r,2.) + 23030.*std::pow(z/r,4.) - 156604.*std::pow(z/r,6.) + 296429.*std::pow(z/r,8.))*std::pow((x-I*y)/r,19.))/(3.3554432e7*1.);
}
else if (m == -18){
return (5.*std::sqrt(4050390505161./Pi)*z/r*(45. - 2820.*std::pow(z/r,2.) + 41454.*std::pow(z/r,4.) - 201348.*std::pow(z/r,6.) + 296429.*std::pow(z/r,8.))*std::pow((x-I*y)/r,18.))/(3.3554432e7*1.);
}
else if (m == -17){
return (15.*std::sqrt(4050390505161./(2.*Pi))*(-1. + 225.*std::pow(z/r,2.) - 7050.*std::pow(z/r,4.) + 69090.*std::pow(z/r,6.) - 251685.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x-I*y)/r,17.))/(3.3554432e7*1.);
}
else if (m == -16){
return (15.*std::sqrt(4050390505161./(2.*Pi))*z/r*(-11. + 825.*std::pow(z/r,2.) - 15510.*std::pow(z/r,4.) + 108570.*std::pow(z/r,6.) - 307615.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x-I*y)/r,16.))/(1.6777216e7*1.);
}
else if (m == -15){
return (15.*std::sqrt(31398376009./(2.*Pi))*(11. - 2838.*std::pow(z/r,2.) + 106425.*std::pow(z/r,4.) - 1333860.*std::pow(z/r,6.) + 7002765.*std::pow(z/r,8.) - 15872934.*std::pow(z/r,10.) + 12746447.*std::pow(z/r,12.))*std::pow((x-I*y)/r,15.))/(3.3554432e7*1.);
}
else if (m == -14){
return (15.*std::sqrt(1035111297./Pi)*z/r*(1001. - 86086.*std::pow(z/r,2.) + 1936935.*std::pow(z/r,4.) - 17340180.*std::pow(z/r,6.) + 70805735.*std::pow(z/r,8.) - 131312454.*std::pow(z/r,10.) + 89225129.*std::pow(z/r,12.))*std::pow((x-I*y)/r,14.))/(3.3554432e7*1.);
}
else if (m == -13){
return (15.*std::sqrt(176726319./(2.*Pi))*(-143. + 41041.*std::pow(z/r,2.) - 1764763.*std::pow(z/r,4.) + 26471445.*std::pow(z/r,6.) - 177736845.*std::pow(z/r,8.) + 580607027.*std::pow(z/r,10.) - 897301769.*std::pow(z/r,12.) + 522604327.*std::pow(z/r,14.))*std::pow((x-I*y)/r,13.))/(3.3554432e7*1.);
}
else if (m == -12){
return (15.*std::sqrt(58908773./Pi)*z/r*(-2145. + 205205.*std::pow(z/r,2.) - 5294289.*std::pow(z/r,4.) + 56724525.*std::pow(z/r,6.) - 296228075.*std::pow(z/r,8.) + 791736855.*std::pow(z/r,10.) - 1035348195.*std::pow(z/r,12.) + 522604327.*std::pow(z/r,14.))*std::pow((x-I*y)/r,12.))/(1.6777216e7*1.);
}
else if (m == -11){
return (15.*std::sqrt(2297442147./Pi)*(55. - 17160.*std::pow(z/r,2.) + 820820.*std::pow(z/r,4.) - 14118104.*std::pow(z/r,6.) + 113449050.*std::pow(z/r,8.) - 473964920.*std::pow(z/r,10.) + 1055649140.*std::pow(z/r,12.) - 1183255080.*std::pow(z/r,14.) + 522604327.*std::pow(z/r,16.))*std::pow((x-I*y)/r,11.))/(6.7108864e7*1.);
}
else if (m == -10){
return (15.*std::sqrt(742073813481./(2.*Pi))*z/r*(55. - 5720.*std::pow(z/r,2.) + 164164.*std::pow(z/r,4.) - 2016872.*std::pow(z/r,6.) + 12605450.*std::pow(z/r,8.) - 43087720.*std::pow(z/r,10.) + 81203780.*std::pow(z/r,12.) - 78883672.*std::pow(z/r,14.) + 30741431.*std::pow(z/r,16.))*std::pow((x-I*y)/r,10.))/(3.3554432e7*1.);
}
else if (m == -9){
return (5.*std::sqrt(20056049013./Pi)*(-55. + 18315.*std::pow(z/r,2.) - 952380.*std::pow(z/r,4.) + 18222204.*std::pow(z/r,6.) - 167904594.*std::pow(z/r,8.) + 839522970.*std::pow(z/r,10.) - 2391368460.*std::pow(z/r,12.) + 3862979820.*std::pow(z/r,14.) - 3283532847.*std::pow(z/r,16.) + 1137432947.*std::pow(z/r,18.))*std::pow((x-I*y)/r,9.))/(6.7108864e7*1.);
}
else if (m == -8){
return (15.*std::sqrt(1055581527./Pi)*z/r*(-1045. + 115995.*std::pow(z/r,2.) - 3619044.*std::pow(z/r,4.) + 49460268.*std::pow(z/r,6.) - 354465254.*std::pow(z/r,8.) + 1450085130.*std::pow(z/r,10.) - 3495076980.*std::pow(z/r,12.) + 4893107772.*std::pow(z/r,14.) - 3669830829.*std::pow(z/r,16.) + 1137432947.*std::pow(z/r,18.))*std::pow((x-I*y)/r,8.))/(3.3554432e7*1.);
}
else if (m == -7){
return (15.*std::sqrt(150797361./Pi)*(209. - 73150.*std::pow(z/r,2.) + 4059825.*std::pow(z/r,4.) - 84444360.*std::pow(z/r,6.) + 865554690.*std::pow(z/r,8.) - 4962513556.*std::pow(z/r,10.) + 16917659850.*std::pow(z/r,12.) - 34950769800.*std::pow(z/r,14.) + 42814693005.*std::pow(z/r,16.) - 28543128670.*std::pow(z/r,18.) + 7962030629.*std::pow(z/r,20.))*std::pow((x-I*y)/r,7.))/(6.7108864e7*1.);
}
else if (m == -6){
return (15.*std::sqrt(20697677./(2.*Pi))*(10659.*z/r - 1243550.*std::pow(z/r,3.) + 41410215.*std::pow(z/r,5.) - 615237480.*std::pow(z/r,7.) + 4904809910.*std::pow(z/r,9.) - 23008017396.*std::pow(z/r,11.) + 66369280950.*std::pow(z/r,13.) - 118832617320.*std::pow(z/r,15.) + 128444079015.*std::pow(z/r,17.) - 76615766430.*std::pow(z/r,19.) + 19336360099.*std::pow(z/r,21.))*std::pow((x-I*y)/r,6.))/(3.3554432e7*1.);
}
else if (m == -5){
return (15.*std::sqrt(62093031./Pi)*(-323. + 117249.*std::pow(z/r,2.) - 6839525.*std::pow(z/r,4.) + 151837455.*std::pow(z/r,6.) - 1691903070.*std::pow(z/r,8.) + 10790581802.*std::pow(z/r,10.) - 42181365226.*std::pow(z/r,12.) + 104294584350.*std::pow(z/r,14.) - 163394848815.*std::pow(z/r,16.) + 156987207685.*std::pow(z/r,18.) - 84277343073.*std::pow(z/r,20.) + 19336360099.*std::pow(z/r,22.))*std::pow((x-I*y)/r,5.))/(6.7108864e7*1.);
}
else if (m == -4){
return (15.*std::sqrt(2699697./(2.*Pi))*(-7429.*z/r + 898909.*std::pow(z/r,3.) - 31461815.*std::pow(z/r,5.) + 498894495.*std::pow(z/r,7.) - 4323752290.*std::pow(z/r,9.) + 22562125586.*std::pow(z/r,11.) - 74628569246.*std::pow(z/r,13.) + 159918362670.*std::pow(z/r,15.) - 221063618985.*std::pow(z/r,17.) + 190037146145.*std::pow(z/r,19.) - 92303756699.*std::pow(z/r,21.) + 19336360099.*std::pow(z/r,23.))*std::pow((x-I*y)/r,4.))/(8.388608e6*1.);
}
else if (m == -3){
return (15.*std::sqrt(29029./Pi)*(7429. - 2763588.*std::pow(z/r,2.) + 167197074.*std::pow(z/r,4.) - 3901265060.*std::pow(z/r,6.) + 46397188035.*std::pow(z/r,8.) - 321687170376.*std::pow(z/r,10.) + 1398851786332.*std::pow(z/r,12.) - 3965975394216.*std::pow(z/r,14.) + 7436203864155.*std::pow(z/r,16.) - 9137296251380.*std::pow(z/r,18.) + 7069381836594.*std::pow(z/r,20.) - 3121545226548.*std::pow(z/r,22.) + 599427163069.*std::pow(z/r,24.))*std::pow((x-I*y)/r,3.))/(3.3554432e7*1.);
}
else if (m == -2){
return (3.*std::sqrt(435435./(2.*Pi))*(185725.*z/r - 23029900.*std::pow(z/r,3.) + 835985370.*std::pow(z/r,5.) - 13933089500.*std::pow(z/r,7.) + 128881077875.*std::pow(z/r,9.) - 731107205400.*std::pow(z/r,11.) + 2690099589100.*std::pow(z/r,13.) - 6609958990360.*std::pow(z/r,15.) + 10935593917875.*std::pow(z/r,17.) - 12022758225500.*std::pow(z/r,19.) + 8415930757850.*std::pow(z/r,21.) - 3392983941900.*std::pow(z/r,23.) + 599427163069.*std::pow(z/r,25.))*std::pow((x-I*y)/r,2.))/(1.6777216e7*1.);
}
else if (m == -1){
return (3.*std::sqrt(1155./Pi)*(-185725. + 70018325.*std::pow(z/r,2.) - 4341136150.*std::pow(z/r,4.) + 105055494830.*std::pow(z/r,6.) - 1313193685375.*std::pow(z/r,8.) + 9717633271775.*std::pow(z/r,10.) - 45937902739300.*std::pow(z/r,12.) + 144881077870100.*std::pow(z/r,14.) - 311494317420715.*std::pow(z/r,16.) + 458079878559875.*std::pow(z/r,18.) - 453257985101350.*std::pow(z/r,20.) + 288436899609950.*std::pow(z/r,22.) - 106596245508025.*std::pow(z/r,24.) + 17383387729001.*std::pow(z/r,26.))*(x-I*y)/r)/(3.3554432e7*1.);
}
else if (m == 0){
return (std::sqrt(55./Pi)*(-35102025.*z/r + 4411154475.*std::pow(z/r,3.) - 164094946470.*std::pow(z/r,5.) + 2836498360410.*std::pow(z/r,7.) - 27577067392875.*std::pow(z/r,9.) + 166966608033225.*std::pow(z/r,11.) - 667866432132900.*std::pow(z/r,13.) + 1825501581163260.*std::pow(z/r,15.) - 3463083881912655.*std::pow(z/r,17.) + 4556689318306125.*std::pow(z/r,19.) - 4079321865912150.*std::pow(z/r,21.) + 2370198870707850.*std::pow(z/r,23.) - 805867616040669.*std::pow(z/r,25.) + 121683714103007.*std::pow(z/r,27.)))/1.6777216e7;
}
else if (m == 1){
return (-3.*1.*std::sqrt(1155./Pi)*(-185725. + 70018325.*std::pow(z/r,2.) - 4341136150.*std::pow(z/r,4.) + 105055494830.*std::pow(z/r,6.) - 1313193685375.*std::pow(z/r,8.) + 9717633271775.*std::pow(z/r,10.) - 45937902739300.*std::pow(z/r,12.) + 144881077870100.*std::pow(z/r,14.) - 311494317420715.*std::pow(z/r,16.) + 458079878559875.*std::pow(z/r,18.) - 453257985101350.*std::pow(z/r,20.) + 288436899609950.*std::pow(z/r,22.) - 106596245508025.*std::pow(z/r,24.) + 17383387729001.*std::pow(z/r,26.))*(x+I*y)/r)/3.3554432e7;
}
else if (m == 2){
return (3.*1.*std::sqrt(435435./(2.*Pi))*(185725.*z/r - 23029900.*std::pow(z/r,3.) + 835985370.*std::pow(z/r,5.) - 13933089500.*std::pow(z/r,7.) + 128881077875.*std::pow(z/r,9.) - 731107205400.*std::pow(z/r,11.) + 2690099589100.*std::pow(z/r,13.) - 6609958990360.*std::pow(z/r,15.) + 10935593917875.*std::pow(z/r,17.) - 12022758225500.*std::pow(z/r,19.) + 8415930757850.*std::pow(z/r,21.) - 3392983941900.*std::pow(z/r,23.) + 599427163069.*std::pow(z/r,25.))*std::pow((x+I*y)/r,2.))/1.6777216e7;
}
else if (m == 3){
return (-15.*1.*std::sqrt(29029./Pi)*(7429. - 2763588.*std::pow(z/r,2.) + 167197074.*std::pow(z/r,4.) - 3901265060.*std::pow(z/r,6.) + 46397188035.*std::pow(z/r,8.) - 321687170376.*std::pow(z/r,10.) + 1398851786332.*std::pow(z/r,12.) - 3965975394216.*std::pow(z/r,14.) + 7436203864155.*std::pow(z/r,16.) - 9137296251380.*std::pow(z/r,18.) + 7069381836594.*std::pow(z/r,20.) - 3121545226548.*std::pow(z/r,22.) + 599427163069.*std::pow(z/r,24.))*std::pow((x+I*y)/r,3.))/3.3554432e7;
}
else if (m == 4){
return (15.*1.*std::sqrt(2699697./(2.*Pi))*(-7429.*z/r + 898909.*std::pow(z/r,3.) - 31461815.*std::pow(z/r,5.) + 498894495.*std::pow(z/r,7.) - 4323752290.*std::pow(z/r,9.) + 22562125586.*std::pow(z/r,11.) - 74628569246.*std::pow(z/r,13.) + 159918362670.*std::pow(z/r,15.) - 221063618985.*std::pow(z/r,17.) + 190037146145.*std::pow(z/r,19.) - 92303756699.*std::pow(z/r,21.) + 19336360099.*std::pow(z/r,23.))*std::pow((x+I*y)/r,4.))/8.388608e6;
}
else if (m == 5){
return (-15.*1.*std::sqrt(62093031./Pi)*(-323. + 117249.*std::pow(z/r,2.) - 6839525.*std::pow(z/r,4.) + 151837455.*std::pow(z/r,6.) - 1691903070.*std::pow(z/r,8.) + 10790581802.*std::pow(z/r,10.) - 42181365226.*std::pow(z/r,12.) + 104294584350.*std::pow(z/r,14.) - 163394848815.*std::pow(z/r,16.) + 156987207685.*std::pow(z/r,18.) - 84277343073.*std::pow(z/r,20.) + 19336360099.*std::pow(z/r,22.))*std::pow((x+I*y)/r,5.))/6.7108864e7;
}
else if (m == 6){
return (15.*1.*std::sqrt(20697677./(2.*Pi))*(10659.*z/r - 1243550.*std::pow(z/r,3.) + 41410215.*std::pow(z/r,5.) - 615237480.*std::pow(z/r,7.) + 4904809910.*std::pow(z/r,9.) - 23008017396.*std::pow(z/r,11.) + 66369280950.*std::pow(z/r,13.) - 118832617320.*std::pow(z/r,15.) + 128444079015.*std::pow(z/r,17.) - 76615766430.*std::pow(z/r,19.) + 19336360099.*std::pow(z/r,21.))*std::pow((x+I*y)/r,6.))/3.3554432e7;
}
else if (m == 7){
return (-15.*1.*std::sqrt(150797361./Pi)*(209. - 73150.*std::pow(z/r,2.) + 4059825.*std::pow(z/r,4.) - 84444360.*std::pow(z/r,6.) + 865554690.*std::pow(z/r,8.) - 4962513556.*std::pow(z/r,10.) + 16917659850.*std::pow(z/r,12.) - 34950769800.*std::pow(z/r,14.) + 42814693005.*std::pow(z/r,16.) - 28543128670.*std::pow(z/r,18.) + 7962030629.*std::pow(z/r,20.))*std::pow((x+I*y)/r,7.))/6.7108864e7;
}
else if (m == 8){
return (15.*1.*std::sqrt(1055581527./Pi)*z/r*(-1045. + 115995.*std::pow(z/r,2.) - 3619044.*std::pow(z/r,4.) + 49460268.*std::pow(z/r,6.) - 354465254.*std::pow(z/r,8.) + 1450085130.*std::pow(z/r,10.) - 3495076980.*std::pow(z/r,12.) + 4893107772.*std::pow(z/r,14.) - 3669830829.*std::pow(z/r,16.) + 1137432947.*std::pow(z/r,18.))*std::pow((x+I*y)/r,8.))/3.3554432e7;
}
else if (m == 9){
return (-5.*1.*std::sqrt(20056049013./Pi)*(-55. + 18315.*std::pow(z/r,2.) - 952380.*std::pow(z/r,4.) + 18222204.*std::pow(z/r,6.) - 167904594.*std::pow(z/r,8.) + 839522970.*std::pow(z/r,10.) - 2391368460.*std::pow(z/r,12.) + 3862979820.*std::pow(z/r,14.) - 3283532847.*std::pow(z/r,16.) + 1137432947.*std::pow(z/r,18.))*std::pow((x+I*y)/r,9.))/6.7108864e7;
}
else if (m == 10){
return (15.*1.*std::sqrt(742073813481./(2.*Pi))*z/r*(55. - 5720.*std::pow(z/r,2.) + 164164.*std::pow(z/r,4.) - 2016872.*std::pow(z/r,6.) + 12605450.*std::pow(z/r,8.) - 43087720.*std::pow(z/r,10.) + 81203780.*std::pow(z/r,12.) - 78883672.*std::pow(z/r,14.) + 30741431.*std::pow(z/r,16.))*std::pow((x+I*y)/r,10.))/3.3554432e7;
}
else if (m == 11){
return (-15.*1.*std::sqrt(2297442147./Pi)*(55. - 17160.*std::pow(z/r,2.) + 820820.*std::pow(z/r,4.) - 14118104.*std::pow(z/r,6.) + 113449050.*std::pow(z/r,8.) - 473964920.*std::pow(z/r,10.) + 1055649140.*std::pow(z/r,12.) - 1183255080.*std::pow(z/r,14.) + 522604327.*std::pow(z/r,16.))*std::pow((x+I*y)/r,11.))/6.7108864e7;
}
else if (m == 12){
return (15.*1.*std::sqrt(58908773./Pi)*z/r*(-2145. + 205205.*std::pow(z/r,2.) - 5294289.*std::pow(z/r,4.) + 56724525.*std::pow(z/r,6.) - 296228075.*std::pow(z/r,8.) + 791736855.*std::pow(z/r,10.) - 1035348195.*std::pow(z/r,12.) + 522604327.*std::pow(z/r,14.))*std::pow((x+I*y)/r,12.))/1.6777216e7;
}
else if (m == 13){
return (-15.*1.*std::sqrt(176726319./(2.*Pi))*(-143. + 41041.*std::pow(z/r,2.) - 1764763.*std::pow(z/r,4.) + 26471445.*std::pow(z/r,6.) - 177736845.*std::pow(z/r,8.) + 580607027.*std::pow(z/r,10.) - 897301769.*std::pow(z/r,12.) + 522604327.*std::pow(z/r,14.))*std::pow((x+I*y)/r,13.))/3.3554432e7;
}
else if (m == 14){
return (15.*1.*std::sqrt(1035111297./Pi)*z/r*(1001. - 86086.*std::pow(z/r,2.) + 1936935.*std::pow(z/r,4.) - 17340180.*std::pow(z/r,6.) + 70805735.*std::pow(z/r,8.) - 131312454.*std::pow(z/r,10.) + 89225129.*std::pow(z/r,12.))*std::pow((x+I*y)/r,14.))/3.3554432e7;
}
else if (m == 15){
return (-15.*1.*std::sqrt(31398376009./(2.*Pi))*(11. - 2838.*std::pow(z/r,2.) + 106425.*std::pow(z/r,4.) - 1333860.*std::pow(z/r,6.) + 7002765.*std::pow(z/r,8.) - 15872934.*std::pow(z/r,10.) + 12746447.*std::pow(z/r,12.))*std::pow((x+I*y)/r,15.))/3.3554432e7;
}
else if (m == 16){
return (15.*1.*std::sqrt(4050390505161./(2.*Pi))*z/r*(-11. + 825.*std::pow(z/r,2.) - 15510.*std::pow(z/r,4.) + 108570.*std::pow(z/r,6.) - 307615.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x+I*y)/r,16.))/1.6777216e7;
}
else if (m == 17){
return (-15.*1.*std::sqrt(4050390505161./(2.*Pi))*(-1. + 225.*std::pow(z/r,2.) - 7050.*std::pow(z/r,4.) + 69090.*std::pow(z/r,6.) - 251685.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x+I*y)/r,17.))/3.3554432e7;
}
else if (m == 18){
return (5.*1.*std::sqrt(4050390505161./Pi)*z/r*(45. - 2820.*std::pow(z/r,2.) + 41454.*std::pow(z/r,4.) - 201348.*std::pow(z/r,6.) + 296429.*std::pow(z/r,8.))*std::pow((x+I*y)/r,18.))/3.3554432e7;
}
else if (m == 19){
return (-15.*1.*std::sqrt(176103935007./(2.*Pi))*(5. - 940.*std::pow(z/r,2.) + 23030.*std::pow(z/r,4.) - 156604.*std::pow(z/r,6.) + 296429.*std::pow(z/r,8.))*std::pow((x+I*y)/r,19.))/3.3554432e7;
}
else if (m == 20){
return (15.*1.*std::sqrt(8276884945329./Pi)*z/r*(-5. + 245.*std::pow(z/r,2.) - 2499.*std::pow(z/r,4.) + 6307.*std::pow(z/r,6.))*std::pow((x+I*y)/r,20.))/1.6777216e7;
}
else if (m == 21){
return (-15.*1.*std::sqrt(394137378349./Pi)*(-5. + 735.*std::pow(z/r,2.) - 12495.*std::pow(z/r,4.) + 44149.*std::pow(z/r,6.))*std::pow((x+I*y)/r,21.))/6.7108864e7;
}
else if (m == 22){
return (105.*1.*std::sqrt(1182412135047./(2.*Pi))*z/r*(5. - 170.*std::pow(z/r,2.) + 901.*std::pow(z/r,4.))*std::pow((x+I*y)/r,22.))/3.3554432e7;
}
else if (m == 23){
return (-21.*1.*std::sqrt(5912060675235./Pi)*(1. - 102.*std::pow(z/r,2.) + 901.*std::pow(z/r,4.))*std::pow((x+I*y)/r,23.))/6.7108864e7;
}
else if (m == 24){
return (21.*1.*std::sqrt(33501677159665./Pi)*z/r*(-3. + 53.*std::pow(z/r,2.))*std::pow((x+I*y)/r,24.))/3.3554432e7;
}
else if (m == 25){
return (-21.*1.*std::sqrt(7731156267615./Pi)*(-1. + 53.*std::pow(z/r,2.))*std::pow((x+I*y)/r,25.))/6.7108864e7;
}
else if (m == 26){
return (21.*1.*std::sqrt(409751282183595./(2.*Pi))*z/r*std::pow((x+I*y)/r,26.))/3.3554432e7;
}
else if (m == 27){
return (-7.*1.*std::sqrt(136583760727865./Pi)*std::pow((x+I*y)/r,27.))/6.7108864e7;
}
else{return 0.;}
}

else if (l == 28){
if(m == -28){
return (std::sqrt(54496920530418135./(2.*Pi))*std::pow((x-I*y)/r,28.))/(1.34217728e8*1.);
}
else if (m == -27){
return (7.*std::sqrt(7785274361488305./Pi)*z/r*std::pow((x-I*y)/r,27.))/(6.7108864e7*1.);
}
else if (m == -26){
return (7.*std::sqrt(141550442936151./(2.*Pi))*(-1. + 55.*std::pow(z/r,2.))*std::pow((x-I*y)/r,26.))/(6.7108864e7*1.);
}
else if (m == -25){
return (21.*std::sqrt(141550442936151./Pi)*z/r*(-3. + 55.*std::pow(z/r,2.))*std::pow((x-I*y)/r,25.))/(6.7108864e7*1.);
}
else if (m == -24){
return (21.*std::sqrt(2670763074267./Pi)*(3. - 318.*std::pow(z/r,2.) + 2915.*std::pow(z/r,4.))*std::pow((x-I*y)/r,24.))/(1.34217728e8*1.);
}
else if (m == -23){
return (21.*std::sqrt(173599599827355./Pi)*z/r*(3. - 106.*std::pow(z/r,2.) + 583.*std::pow(z/r,4.))*std::pow((x-I*y)/r,23.))/(6.7108864e7*1.);
}
else if (m == -22){
return (21.*std::sqrt(10211741166315./(2.*Pi))*(-1. + 153.*std::pow(z/r,2.) - 2703.*std::pow(z/r,4.) + 9911.*std::pow(z/r,6.))*std::pow((x-I*y)/r,22.))/(6.7108864e7*1.);
}
else if (m == -21){
return (3.*std::sqrt(71482188164205./Pi)*z/r*(-35. + 1785.*std::pow(z/r,2.) - 18921.*std::pow(z/r,4.) + 49555.*std::pow(z/r,6.))*std::pow((x-I*y)/r,21.))/(6.7108864e7*1.);
}
else if (m == -20){
return (3.*std::sqrt(71482188164205./(2.*Pi))*(5. - 980.*std::pow(z/r,2.) + 24990.*std::pow(z/r,4.) - 176596.*std::pow(z/r,6.) + 346885.*std::pow(z/r,8.))*std::pow((x-I*y)/r,20.))/(1.34217728e8*1.);
}
else if (m == -19){
return (3.*std::sqrt(23827396054735./(2.*Pi))*z/r*(45. - 2940.*std::pow(z/r,2.) + 44982.*std::pow(z/r,4.) - 227052.*std::pow(z/r,6.) + 346885.*std::pow(z/r,8.))*std::pow((x-I*y)/r,19.))/(3.3554432e7*1.);
}
else if (m == -18){
return (15.*std::sqrt(101393174701./Pi)*(-9. + 2115.*std::pow(z/r,2.) - 69090.*std::pow(z/r,4.) + 704718.*std::pow(z/r,6.) - 2667861.*std::pow(z/r,8.) + 3260719.*std::pow(z/r,10.))*std::pow((x-I*y)/r,18.))/(6.7108864e7*1.);
}
else if (m == -17){
return (15.*std::sqrt(25652473199353./(2.*Pi))*z/r*(-9. + 705.*std::pow(z/r,2.) - 13818.*std::pow(z/r,4.) + 100674.*std::pow(z/r,6.) - 296429.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x-I*y)/r,17.))/(3.3554432e7*1.);
}
else if (m == -16){
return (3.*std::sqrt(384787097990295./(2.*Pi))*(1. - 270.*std::pow(z/r,2.) + 10575.*std::pow(z/r,4.) - 138180.*std::pow(z/r,6.) + 755055.*std::pow(z/r,8.) - 1778574.*std::pow(z/r,10.) + 1482145.*std::pow(z/r,12.))*std::pow((x-I*y)/r,16.))/(6.7108864e7*1.);
}
else if (m == -15){
return (3.*std::sqrt(2690818867065./(2.*Pi))*z/r*(143. - 12870.*std::pow(z/r,2.) + 302445.*std::pow(z/r,4.) - 2822820.*std::pow(z/r,6.) + 11996985.*std::pow(z/r,8.) - 23121462.*std::pow(z/r,10.) + 16303595.*std::pow(z/r,12.))*std::pow((x-I*y)/r,15.))/(3.3554432e7*1.);
}
else if (m == -14){
return (3.*std::sqrt(8939597565./Pi)*(-143. + 43043.*std::pow(z/r,2.) - 1936935.*std::pow(z/r,4.) + 30345315.*std::pow(z/r,6.) - 212417205.*std::pow(z/r,8.) + 722218497.*std::pow(z/r,10.) - 1159926677.*std::pow(z/r,12.) + 701054585.*std::pow(z/r,14.))*std::pow((x-I*y)/r,14.))/(6.7108864e7*1.);
}
else if (m == -13){
return (15.*std::sqrt(12515436591./(2.*Pi))*z/r*(-429. + 43043.*std::pow(z/r,2.) - 1162161.*std::pow(z/r,4.) + 13005135.*std::pow(z/r,6.) - 70805735.*std::pow(z/r,8.) + 196968681.*std::pow(z/r,10.) - 267675387.*std::pow(z/r,12.) + 140210917.*std::pow(z/r,14.))*std::pow((x-I*y)/r,13.))/(3.3554432e7*1.);
}
else if (m == -12){
return (15.*std::sqrt(305254551./(2.*Pi))*(429. - 140712.*std::pow(z/r,2.) + 7059052.*std::pow(z/r,4.) - 127062936.*std::pow(z/r,6.) + 1066421070.*std::pow(z/r,8.) - 4644856216.*std::pow(z/r,10.) + 10767621228.*std::pow(z/r,12.) - 12542503848.*std::pow(z/r,14.) + 5748647597.*std::pow(z/r,16.))*std::pow((x-I*y)/r,12.))/(1.34217728e8*1.);
}
else if (m == -11){
return (3.*std::sqrt(25946636835./Pi)*z/r*(2145. - 234520.*std::pow(z/r,2.) + 7059052.*std::pow(z/r,4.) - 90759240.*std::pow(z/r,6.) + 592456150.*std::pow(z/r,8.) - 2111298280.*std::pow(z/r,10.) + 4141392780.*std::pow(z/r,12.) - 4180834616.*std::pow(z/r,14.) + 1690778705.*std::pow(z/r,16.))*std::pow((x-I*y)/r,11.))/(6.7108864e7*1.);
}
else if (m == -10){
return (3.*std::sqrt(112435426285./(2.*Pi))*(-55. + 19305.*std::pow(z/r,2.) - 1055340.*std::pow(z/r,4.) + 21177156.*std::pow(z/r,6.) - 204208290.*std::pow(z/r,8.) + 1066421070.*std::pow(z/r,10.) - 3166947420.*std::pow(z/r,12.) + 5324647860.*std::pow(z/r,14.) - 4703438943.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x-I*y)/r,10.))/(6.7108864e7*1.);
}
else if (m == -9){
return (3.*std::sqrt(112435426285./Pi)*z/r*(-1045. + 122265.*std::pow(z/r,2.) - 4010292.*std::pow(z/r,4.) + 57480852.*std::pow(z/r,6.) - 431106390.*std::pow(z/r,8.) + 1842000030.*std::pow(z/r,10.) - 4628615460.*std::pow(z/r,12.) + 6744553956.*std::pow(z/r,14.) - 5256784701.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x-I*y)/r,9.))/(6.7108864e7*1.);
}
else if (m == -8){
return (15.*std::sqrt(607759061./Pi)*(209. - 77330.*std::pow(z/r,2.) + 4523805.*std::pow(z/r,4.) - 98920536.*std::pow(z/r,6.) + 1063395762.*std::pow(z/r,8.) - 6380374572.*std::pow(z/r,10.) + 22718000370.*std::pow(z/r,12.) - 48931077720.*std::pow(z/r,14.) + 62387124093.*std::pow(z/r,16.) - 43222451986.*std::pow(z/r,18.) + 12511762417.*std::pow(z/r,20.))*std::pow((x-I*y)/r,8.))/(1.34217728e8*1.);
}
else if (m == -7){
return (15.*std::sqrt(260468169./Pi)*(4389.*z/r - 541310.*std::pow(z/r,3.) + 18999981.*std::pow(z/r,5.) - 296761608.*std::pow(z/r,7.) + 2481256778.*std::pow(z/r,9.) - 12180715092.*std::pow(z/r,11.) + 36698308290.*std::pow(z/r,13.) - 68503508808.*std::pow(z/r,15.) + 77066447409.*std::pow(z/r,17.) - 47772183774.*std::pow(z/r,19.) + 12511762417.*std::pow(z/r,21.))*std::pow((x-I*y)/r,7.))/(6.7108864e7*1.);
}
else if (m == -6){
return (3.*std::sqrt(100280245065./(2.*Pi))*(-57. + 21945.*std::pow(z/r,2.) - 1353275.*std::pow(z/r,4.) + 31666635.*std::pow(z/r,6.) - 370952010.*std::pow(z/r,8.) + 2481256778.*std::pow(z/r,10.) - 10150595910.*std::pow(z/r,12.) + 26213077350.*std::pow(z/r,14.) - 42814693005.*std::pow(z/r,16.) + 42814693005.*std::pow(z/r,18.) - 23886091887.*std::pow(z/r,20.) + 5687164735.*std::pow(z/r,22.))*std::pow((x-I*y)/r,6.))/(6.7108864e7*1.);
}
else if (m == -5){
return (3.*std::sqrt(256471215./Pi)*(-22287.*z/r + 2860165.*std::pow(z/r,3.) - 105826105.*std::pow(z/r,5.) + 1768807755.*std::pow(z/r,7.) - 16115803990.*std::pow(z/r,9.) + 88197400018.*std::pow(z/r,11.) - 305298692370.*std::pow(z/r,13.) + 683287549590.*std::pow(z/r,15.) - 984737939115.*std::pow(z/r,17.) + 881081313945.*std::pow(z/r,19.) - 444736282277.*std::pow(z/r,21.) + 96681800495.*std::pow(z/r,23.))*std::pow((x-I*y)/r,5.))/(6.7108864e7*1.);
}
else if (m == -4){
return (3.*std::sqrt(23315565./(2.*Pi))*(7429. - 2941884.*std::pow(z/r,2.) + 188770890.*std::pow(z/r,4.) - 4656348620.*std::pow(z/r,6.) + 58370655915.*std::pow(z/r,8.) - 425457225336.*std::pow(z/r,10.) + 1940342800396.*std::pow(z/r,12.) - 5757061056120.*std::pow(z/r,14.) + 11274244568235.*std::pow(z/r,16.) - 14442823107020.*std::pow(z/r,18.) + 11630273344074.*std::pow(z/r,20.) - 5336835387324.*std::pow(z/r,22.) + 1063499805445.*std::pow(z/r,24.))*std::pow((x-I*y)/r,4.))/(1.34217728e8*1.);
}
else if (m == -3){
return (3.*std::sqrt(23315565./Pi)*(37145.*z/r - 4903140.*std::pow(z/r,3.) + 188770890.*std::pow(z/r,5.) - 3325963300.*std::pow(z/r,7.) + 32428142175.*std::pow(z/r,9.) - 193389647880.*std::pow(z/r,11.) + 746285692460.*std::pow(z/r,13.) - 1919020352040.*std::pow(z/r,15.) + 3315954284775.*std::pow(z/r,17.) - 3800742922900.*std::pow(z/r,19.) + 2769112700970.*std::pow(z/r,21.) - 1160181605940.*std::pow(z/r,23.) + 212699961089.*std::pow(z/r,25.))*std::pow((x-I*y)/r,3.))/(3.3554432e7*1.);
}
else if (m == -2){
return (3.*std::sqrt(57855./(2.*Pi))*(-37145. + 14969435.*std::pow(z/r,2.) - 987982710.*std::pow(z/r,4.) + 25358222890.*std::pow(z/r,6.) - 335090802475.*std::pow(z/r,8.) + 2613708259305.*std::pow(z/r,10.) - 12989338015940.*std::pow(z/r,12.) + 42964733437340.*std::pow(z/r,14.) - 96670650234015.*std::pow(z/r,16.) + 148481064084925.*std::pow(z/r,18.) - 153169939792870.*std::pow(z/r,20.) + 101450219862810.*std::pow(z/r,22.) - 38962765599485.*std::pow(z/r,24.) + 6593698793759.*std::pow(z/r,26.))*std::pow((x-I*y)/r,2.))/(3.3554432e7*1.);
}
else if (m == -1){
return (std::sqrt(11571./Pi)*(-5014575.*z/r + 673624575.*std::pow(z/r,3.) - 26675533170.*std::pow(z/r,5.) + 489051441450.*std::pow(z/r,7.) - 5026362037125.*std::pow(z/r,9.) + 32077328636925.*std::pow(z/r,11.) - 134889279396300.*std::pow(z/r,13.) + 386682600936060.*std::pow(z/r,15.) - 767678693034825.*std::pow(z/r,17.) + 1054997034287625.*std::pow(z/r,19.) - 984663898668450.*std::pow(z/r,21.) + 595468681803450.*std::pow(z/r,23.) - 210398934237219.*std::pow(z/r,25.) + 32968493968795.*std::pow(z/r,27.))*(x-I*y)/r)/(3.3554432e7*1.);
}
else if (m == 0){
return (std::sqrt(57./Pi)*(5014575. - 2035917450.*std::pow(z/r,2.) + 136745788725.*std::pow(z/r,4.) - 3610088822340.*std::pow(z/r,6.) + 49638721307175.*std::pow(z/r,8.) - 408140597414550.*std::pow(z/r,10.) + 2170565904431925.*std::pow(z/r,12.) - 7823578204985400.*std::pow(z/r,14.) + 19624141997505045.*std::pow(z/r,16.) - 34630838819126550.*std::pow(z/r,18.) + 42832879592077575.*std::pow(z/r,20.) - 36343049350853700.*std::pow(z/r,22.) + 20146690401016725.*std::pow(z/r,24.) - 6570920561562378.*std::pow(z/r,26.) + 956086325095055.*std::pow(z/r,28.)))/6.7108864e7;
}
else if (m == 1){
return -2.9802322387695312e-8*(1.*std::sqrt(11571./Pi)*(-5014575.*z/r + 673624575.*std::pow(z/r,3.) - 26675533170.*std::pow(z/r,5.) + 489051441450.*std::pow(z/r,7.) - 5026362037125.*std::pow(z/r,9.) + 32077328636925.*std::pow(z/r,11.) - 134889279396300.*std::pow(z/r,13.) + 386682600936060.*std::pow(z/r,15.) - 767678693034825.*std::pow(z/r,17.) + 1054997034287625.*std::pow(z/r,19.) - 984663898668450.*std::pow(z/r,21.) + 595468681803450.*std::pow(z/r,23.) - 210398934237219.*std::pow(z/r,25.) + 32968493968795.*std::pow(z/r,27.))*(x+I*y)/r);
}
else if (m == 2){
return (3.*1.*std::sqrt(57855./(2.*Pi))*(-37145. + 14969435.*std::pow(z/r,2.) - 987982710.*std::pow(z/r,4.) + 25358222890.*std::pow(z/r,6.) - 335090802475.*std::pow(z/r,8.) + 2613708259305.*std::pow(z/r,10.) - 12989338015940.*std::pow(z/r,12.) + 42964733437340.*std::pow(z/r,14.) - 96670650234015.*std::pow(z/r,16.) + 148481064084925.*std::pow(z/r,18.) - 153169939792870.*std::pow(z/r,20.) + 101450219862810.*std::pow(z/r,22.) - 38962765599485.*std::pow(z/r,24.) + 6593698793759.*std::pow(z/r,26.))*std::pow((x+I*y)/r,2.))/3.3554432e7;
}
else if (m == 3){
return (-3.*1.*std::sqrt(23315565./Pi)*(37145.*z/r - 4903140.*std::pow(z/r,3.) + 188770890.*std::pow(z/r,5.) - 3325963300.*std::pow(z/r,7.) + 32428142175.*std::pow(z/r,9.) - 193389647880.*std::pow(z/r,11.) + 746285692460.*std::pow(z/r,13.) - 1919020352040.*std::pow(z/r,15.) + 3315954284775.*std::pow(z/r,17.) - 3800742922900.*std::pow(z/r,19.) + 2769112700970.*std::pow(z/r,21.) - 1160181605940.*std::pow(z/r,23.) + 212699961089.*std::pow(z/r,25.))*std::pow((x+I*y)/r,3.))/3.3554432e7;
}
else if (m == 4){
return (3.*1.*std::sqrt(23315565./(2.*Pi))*(7429. - 2941884.*std::pow(z/r,2.) + 188770890.*std::pow(z/r,4.) - 4656348620.*std::pow(z/r,6.) + 58370655915.*std::pow(z/r,8.) - 425457225336.*std::pow(z/r,10.) + 1940342800396.*std::pow(z/r,12.) - 5757061056120.*std::pow(z/r,14.) + 11274244568235.*std::pow(z/r,16.) - 14442823107020.*std::pow(z/r,18.) + 11630273344074.*std::pow(z/r,20.) - 5336835387324.*std::pow(z/r,22.) + 1063499805445.*std::pow(z/r,24.))*std::pow((x+I*y)/r,4.))/1.34217728e8;
}
else if (m == 5){
return (-3.*1.*std::sqrt(256471215./Pi)*(-22287.*z/r + 2860165.*std::pow(z/r,3.) - 105826105.*std::pow(z/r,5.) + 1768807755.*std::pow(z/r,7.) - 16115803990.*std::pow(z/r,9.) + 88197400018.*std::pow(z/r,11.) - 305298692370.*std::pow(z/r,13.) + 683287549590.*std::pow(z/r,15.) - 984737939115.*std::pow(z/r,17.) + 881081313945.*std::pow(z/r,19.) - 444736282277.*std::pow(z/r,21.) + 96681800495.*std::pow(z/r,23.))*std::pow((x+I*y)/r,5.))/6.7108864e7;
}
else if (m == 6){
return (3.*1.*std::sqrt(100280245065./(2.*Pi))*(-57. + 21945.*std::pow(z/r,2.) - 1353275.*std::pow(z/r,4.) + 31666635.*std::pow(z/r,6.) - 370952010.*std::pow(z/r,8.) + 2481256778.*std::pow(z/r,10.) - 10150595910.*std::pow(z/r,12.) + 26213077350.*std::pow(z/r,14.) - 42814693005.*std::pow(z/r,16.) + 42814693005.*std::pow(z/r,18.) - 23886091887.*std::pow(z/r,20.) + 5687164735.*std::pow(z/r,22.))*std::pow((x+I*y)/r,6.))/6.7108864e7;
}
else if (m == 7){
return (-15.*1.*std::sqrt(260468169./Pi)*(4389.*z/r - 541310.*std::pow(z/r,3.) + 18999981.*std::pow(z/r,5.) - 296761608.*std::pow(z/r,7.) + 2481256778.*std::pow(z/r,9.) - 12180715092.*std::pow(z/r,11.) + 36698308290.*std::pow(z/r,13.) - 68503508808.*std::pow(z/r,15.) + 77066447409.*std::pow(z/r,17.) - 47772183774.*std::pow(z/r,19.) + 12511762417.*std::pow(z/r,21.))*std::pow((x+I*y)/r,7.))/6.7108864e7;
}
else if (m == 8){
return (15.*1.*std::sqrt(607759061./Pi)*(209. - 77330.*std::pow(z/r,2.) + 4523805.*std::pow(z/r,4.) - 98920536.*std::pow(z/r,6.) + 1063395762.*std::pow(z/r,8.) - 6380374572.*std::pow(z/r,10.) + 22718000370.*std::pow(z/r,12.) - 48931077720.*std::pow(z/r,14.) + 62387124093.*std::pow(z/r,16.) - 43222451986.*std::pow(z/r,18.) + 12511762417.*std::pow(z/r,20.))*std::pow((x+I*y)/r,8.))/1.34217728e8;
}
else if (m == 9){
return (-3.*1.*std::sqrt(112435426285./Pi)*z/r*(-1045. + 122265.*std::pow(z/r,2.) - 4010292.*std::pow(z/r,4.) + 57480852.*std::pow(z/r,6.) - 431106390.*std::pow(z/r,8.) + 1842000030.*std::pow(z/r,10.) - 4628615460.*std::pow(z/r,12.) + 6744553956.*std::pow(z/r,14.) - 5256784701.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x+I*y)/r,9.))/6.7108864e7;
}
else if (m == 10){
return (3.*1.*std::sqrt(112435426285./(2.*Pi))*(-55. + 19305.*std::pow(z/r,2.) - 1055340.*std::pow(z/r,4.) + 21177156.*std::pow(z/r,6.) - 204208290.*std::pow(z/r,8.) + 1066421070.*std::pow(z/r,10.) - 3166947420.*std::pow(z/r,12.) + 5324647860.*std::pow(z/r,14.) - 4703438943.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x+I*y)/r,10.))/6.7108864e7;
}
else if (m == 11){
return (-3.*1.*std::sqrt(25946636835./Pi)*z/r*(2145. - 234520.*std::pow(z/r,2.) + 7059052.*std::pow(z/r,4.) - 90759240.*std::pow(z/r,6.) + 592456150.*std::pow(z/r,8.) - 2111298280.*std::pow(z/r,10.) + 4141392780.*std::pow(z/r,12.) - 4180834616.*std::pow(z/r,14.) + 1690778705.*std::pow(z/r,16.))*std::pow((x+I*y)/r,11.))/6.7108864e7;
}
else if (m == 12){
return (15.*1.*std::sqrt(305254551./(2.*Pi))*(429. - 140712.*std::pow(z/r,2.) + 7059052.*std::pow(z/r,4.) - 127062936.*std::pow(z/r,6.) + 1066421070.*std::pow(z/r,8.) - 4644856216.*std::pow(z/r,10.) + 10767621228.*std::pow(z/r,12.) - 12542503848.*std::pow(z/r,14.) + 5748647597.*std::pow(z/r,16.))*std::pow((x+I*y)/r,12.))/1.34217728e8;
}
else if (m == 13){
return (-15.*1.*std::sqrt(12515436591./(2.*Pi))*z/r*(-429. + 43043.*std::pow(z/r,2.) - 1162161.*std::pow(z/r,4.) + 13005135.*std::pow(z/r,6.) - 70805735.*std::pow(z/r,8.) + 196968681.*std::pow(z/r,10.) - 267675387.*std::pow(z/r,12.) + 140210917.*std::pow(z/r,14.))*std::pow((x+I*y)/r,13.))/3.3554432e7;
}
else if (m == 14){
return (3.*1.*std::sqrt(8939597565./Pi)*(-143. + 43043.*std::pow(z/r,2.) - 1936935.*std::pow(z/r,4.) + 30345315.*std::pow(z/r,6.) - 212417205.*std::pow(z/r,8.) + 722218497.*std::pow(z/r,10.) - 1159926677.*std::pow(z/r,12.) + 701054585.*std::pow(z/r,14.))*std::pow((x+I*y)/r,14.))/6.7108864e7;
}
else if (m == 15){
return (-3.*1.*std::sqrt(2690818867065./(2.*Pi))*z/r*(143. - 12870.*std::pow(z/r,2.) + 302445.*std::pow(z/r,4.) - 2822820.*std::pow(z/r,6.) + 11996985.*std::pow(z/r,8.) - 23121462.*std::pow(z/r,10.) + 16303595.*std::pow(z/r,12.))*std::pow((x+I*y)/r,15.))/3.3554432e7;
}
else if (m == 16){
return (3.*1.*std::sqrt(384787097990295./(2.*Pi))*(1. - 270.*std::pow(z/r,2.) + 10575.*std::pow(z/r,4.) - 138180.*std::pow(z/r,6.) + 755055.*std::pow(z/r,8.) - 1778574.*std::pow(z/r,10.) + 1482145.*std::pow(z/r,12.))*std::pow((x+I*y)/r,16.))/6.7108864e7;
}
else if (m == 17){
return (-15.*1.*std::sqrt(25652473199353./(2.*Pi))*z/r*(-9. + 705.*std::pow(z/r,2.) - 13818.*std::pow(z/r,4.) + 100674.*std::pow(z/r,6.) - 296429.*std::pow(z/r,8.) + 296429.*std::pow(z/r,10.))*std::pow((x+I*y)/r,17.))/3.3554432e7;
}
else if (m == 18){
return (15.*1.*std::sqrt(101393174701./Pi)*(-9. + 2115.*std::pow(z/r,2.) - 69090.*std::pow(z/r,4.) + 704718.*std::pow(z/r,6.) - 2667861.*std::pow(z/r,8.) + 3260719.*std::pow(z/r,10.))*std::pow((x+I*y)/r,18.))/6.7108864e7;
}
else if (m == 19){
return (-3.*1.*std::sqrt(23827396054735./(2.*Pi))*z/r*(45. - 2940.*std::pow(z/r,2.) + 44982.*std::pow(z/r,4.) - 227052.*std::pow(z/r,6.) + 346885.*std::pow(z/r,8.))*std::pow((x+I*y)/r,19.))/3.3554432e7;
}
else if (m == 20){
return (3.*1.*std::sqrt(71482188164205./(2.*Pi))*(5. - 980.*std::pow(z/r,2.) + 24990.*std::pow(z/r,4.) - 176596.*std::pow(z/r,6.) + 346885.*std::pow(z/r,8.))*std::pow((x+I*y)/r,20.))/1.34217728e8;
}
else if (m == 21){
return (-3.*1.*std::sqrt(71482188164205./Pi)*z/r*(-35. + 1785.*std::pow(z/r,2.) - 18921.*std::pow(z/r,4.) + 49555.*std::pow(z/r,6.))*std::pow((x+I*y)/r,21.))/6.7108864e7;
}
else if (m == 22){
return (21.*1.*std::sqrt(10211741166315./(2.*Pi))*(-1. + 153.*std::pow(z/r,2.) - 2703.*std::pow(z/r,4.) + 9911.*std::pow(z/r,6.))*std::pow((x+I*y)/r,22.))/6.7108864e7;
}
else if (m == 23){
return (-21.*1.*std::sqrt(173599599827355./Pi)*z/r*(3. - 106.*std::pow(z/r,2.) + 583.*std::pow(z/r,4.))*std::pow((x+I*y)/r,23.))/6.7108864e7;
}
else if (m == 24){
return (21.*1.*std::sqrt(2670763074267./Pi)*(3. - 318.*std::pow(z/r,2.) + 2915.*std::pow(z/r,4.))*std::pow((x+I*y)/r,24.))/1.34217728e8;
}
else if (m == 25){
return (-21.*1.*std::sqrt(141550442936151./Pi)*z/r*(-3. + 55.*std::pow(z/r,2.))*std::pow((x+I*y)/r,25.))/6.7108864e7;
}
else if (m == 26){
return (7.*1.*std::sqrt(141550442936151./(2.*Pi))*(-1. + 55.*std::pow(z/r,2.))*std::pow((x+I*y)/r,26.))/6.7108864e7;
}
else if (m == 27){
return (-7.*1.*std::sqrt(7785274361488305./Pi)*z/r*std::pow((x+I*y)/r,27.))/6.7108864e7;
}
else if (m == 28){
return (1.*std::sqrt(54496920530418135./(2.*Pi))*std::pow((x+I*y)/r,28.))/1.34217728e8;
}
else{return 0.;}
}

else if (l == 29){
if(m == -29){
return (std::sqrt(110873045217057585./Pi)*std::pow((x-I*y)/r,29.))/(2.68435456e8*1.);
}
else if (m == -28){
return (std::sqrt(3215318311294669965./(2.*Pi))*z/r*std::pow((x-I*y)/r,28.))/(1.34217728e8*1.);
}
else if (m == -27){
return (std::sqrt(56409093180608245./Pi)*(-1. + 57.*std::pow(z/r,2.))*std::pow((x-I*y)/r,27.))/(2.68435456e8*1.);
}
else if (m == -26){
return (7.*std::sqrt(24175325648832105./(2.*Pi))*z/r*(-1. + 19.*std::pow(z/r,2.))*std::pow((x-I*y)/r,26.))/(6.7108864e7*1.);
}
else if (m == -25){
return (7.*std::sqrt(439551375433311./(2.*Pi))*(1. - 110.*std::pow(z/r,2.) + 1045.*std::pow(z/r,4.))*std::pow((x-I*y)/r,25.))/(1.34217728e8*1.);
}
else if (m == -24){
return (21.*std::sqrt(732585625722185./Pi)*z/r*(3. - 110.*std::pow(z/r,2.) + 627.*std::pow(z/r,4.))*std::pow((x-I*y)/r,24.))/(1.34217728e8*1.);
}
else if (m == -23){
return (21.*std::sqrt(41467110889935./(2.*Pi))*(-1. + 159.*std::pow(z/r,2.) - 2915.*std::pow(z/r,4.) + 11077.*std::pow(z/r,6.))*std::pow((x-I*y)/r,23.))/(1.34217728e8*1.);
}
else if (m == -22){
return (3.*std::sqrt(3773507090984085./(2.*Pi))*z/r*(-7. + 371.*std::pow(z/r,2.) - 4081.*std::pow(z/r,4.) + 11077.*std::pow(z/r,6.))*std::pow((x-I*y)/r,22.))/(6.7108864e7*1.);
}
else if (m == -21){
return (3.*std::sqrt(73990335117335./Pi)*(7. - 1428.*std::pow(z/r,2.) + 37842.*std::pow(z/r,4.) - 277508.*std::pow(z/r,6.) + 564927.*std::pow(z/r,8.))*std::pow((x-I*y)/r,21.))/(2.68435456e8*1.);
}
else if (m == -20){
return (3.*std::sqrt(73990335117335./(2.*Pi))*z/r*(105. - 7140.*std::pow(z/r,2.) + 113526.*std::pow(z/r,4.) - 594660.*std::pow(z/r,6.) + 941545.*std::pow(z/r,8.))*std::pow((x-I*y)/r,20.))/(1.34217728e8*1.);
}
else if (m == -19){
return (15.*std::sqrt(14798067023467./Pi)*(-3. + 735.*std::pow(z/r,2.) - 24990.*std::pow(z/r,4.) + 264894.*std::pow(z/r,6.) - 1040655.*std::pow(z/r,8.) + 1318163.*std::pow(z/r,10.))*std::pow((x-I*y)/r,19.))/(2.68435456e8*1.);
}
else if (m == -18){
return (5.*std::sqrt(488336211774411./Pi)*z/r*(-9. + 735.*std::pow(z/r,2.) - 14994.*std::pow(z/r,4.) + 113526.*std::pow(z/r,6.) - 346885.*std::pow(z/r,8.) + 359499.*std::pow(z/r,10.))*std::pow((x-I*y)/r,18.))/(6.7108864e7*1.);
}
else if (m == -17){
return (15.*std::sqrt(3463377388471./Pi)*(3. - 846.*std::pow(z/r,2.) + 34545.*std::pow(z/r,4.) - 469812.*std::pow(z/r,6.) + 2667861.*std::pow(z/r,8.) - 6521438.*std::pow(z/r,10.) + 5632151.*std::pow(z/r,12.))*std::pow((x-I*y)/r,17.))/(1.34217728e8*1.);
}
else if (m == -16){
return (15.*std::sqrt(6127513841141./(2.*Pi))*z/r*(39. - 3666.*std::pow(z/r,2.) + 89817.*std::pow(z/r,4.) - 872508.*std::pow(z/r,6.) + 3853577.*std::pow(z/r,8.) - 7707154.*std::pow(z/r,10.) + 5632151.*std::pow(z/r,12.))*std::pow((x-I*y)/r,16.))/(6.7108864e7*1.);
}
else if (m == -15){
return (3.*std::sqrt(4376795600815./Pi)*(-13. + 4095.*std::pow(z/r,2.) - 192465.*std::pow(z/r,4.) + 3143595.*std::pow(z/r,6.) - 22903335.*std::pow(z/r,8.) + 80925117.*std::pow(z/r,10.) - 134875195.*std::pow(z/r,12.) + 84482265.*std::pow(z/r,14.))*std::pow((x-I*y)/r,15.))/(1.34217728e8*1.);
}
else if (m == -14){
return (15.*std::sqrt(238734305499./Pi)*z/r*(-143. + 15015.*std::pow(z/r,2.) - 423423.*std::pow(z/r,4.) + 4939935.*std::pow(z/r,6.) - 27992965.*std::pow(z/r,8.) + 80925117.*std::pow(z/r,10.) - 114125165.*std::pow(z/r,12.) + 61953661.*std::pow(z/r,14.))*std::pow((x-I*y)/r,14.))/(6.7108864e7*1.);
}
else if (m == -13){
return (15.*std::sqrt(5551960593./Pi)*(143. - 49192.*std::pow(z/r,2.) + 2582580.*std::pow(z/r,4.) - 48552504.*std::pow(z/r,6.) + 424834410.*std::pow(z/r,8.) - 1925915992.*std::pow(z/r,10.) + 4639706708.*std::pow(z/r,12.) - 5608436680.*std::pow(z/r,14.) + 2664007423.*std::pow(z/r,16.))*std::pow((x-I*y)/r,13.))/(2.68435456e8*1.);
}
else if (m == -12){
return (15.*std::sqrt(220227770189./(2.*Pi))*z/r*(429. - 49192.*std::pow(z/r,2.) + 1549548.*std::pow(z/r,4.) - 20808216.*std::pow(z/r,6.) + 141611470.*std::pow(z/r,8.) - 525249816.*std::pow(z/r,10.) + 1070701548.*std::pow(z/r,12.) - 1121687336.*std::pow(z/r,14.) + 470118957.*std::pow(z/r,16.))*std::pow((x-I*y)/r,12.))/(1.34217728e8*1.);
}
else if (m == -11){
return (15.*std::sqrt(5371409029./Pi)*(-143. + 52767.*std::pow(z/r,2.) - 3025308.*std::pow(z/r,4.) + 63531468.*std::pow(z/r,6.) - 639852642.*std::pow(z/r,8.) + 3483642162.*std::pow(z/r,10.) - 10767621228.*std::pow(z/r,12.) + 18813755772.*std::pow(z/r,14.) - 17245942791.*std::pow(z/r,16.) + 6424959079.*std::pow(z/r,18.))*std::pow((x-I*y)/r,11.))/(2.68435456e8*1.);
}
else if (m == -10){
return (3.*std::sqrt(510283857755./(2.*Pi))*z/r*(-715. + 87945.*std::pow(z/r,2.) - 3025308.*std::pow(z/r,4.) + 45379620.*std::pow(z/r,6.) - 355473690.*std::pow(z/r,8.) + 1583473710.*std::pow(z/r,10.) - 4141392780.*std::pow(z/r,12.) + 6271251924.*std::pow(z/r,14.) - 5072336115.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x-I*y)/r,10.))/(6.7108864e7*1.);
}
else if (m == -9){
return (5.*std::sqrt(3980214090489./(2.*Pi))*(11. - 4290.*std::pow(z/r,2.) + 263835.*std::pow(z/r,4.) - 6050616.*std::pow(z/r,6.) + 68069430.*std::pow(z/r,8.) - 426568428.*std::pow(z/r,10.) + 1583473710.*std::pow(z/r,12.) - 3549765240.*std::pow(z/r,14.) + 4703438943.*std::pow(z/r,16.) - 3381557410.*std::pow(z/r,18.) + 1014467223.*std::pow(z/r,20.))*std::pow((x-I*y)/r,9.))/(1.34217728e8*1.);
}
else if (m == -8){
return (15.*std::sqrt(9975473911./Pi)*(1463.*z/r - 190190.*std::pow(z/r,3.) + 7018011.*std::pow(z/r,5.) - 114961704.*std::pow(z/r,7.) + 1005914910.*std::pow(z/r,9.) - 5157600084.*std::pow(z/r,11.) + 16200154110.*std::pow(z/r,13.) - 31474585128.*std::pow(z/r,15.) + 36797492907.*std::pow(z/r,17.) - 23670901870.*std::pow(z/r,19.) + 6424959079.*std::pow(z/r,21.))*std::pow((x-I*y)/r,8.))/(1.34217728e8*1.);
}
else if (m == -7){
return (15.*std::sqrt(2965681433./(2.*Pi))*(-133. + 54131.*std::pow(z/r,2.) - 3518515.*std::pow(z/r,4.) + 86555469.*std::pow(z/r,6.) - 1063395762.*std::pow(z/r,8.) + 7443770334.*std::pow(z/r,10.) - 31805200518.*std::pow(z/r,12.) + 85629386010.*std::pow(z/r,14.) - 145569956217.*std::pow(z/r,16.) + 151278581951.*std::pow(z/r,18.) - 87582336919.*std::pow(z/r,20.) + 21611225993.*std::pow(z/r,22.))*std::pow((x-I*y)/r,7.))/(1.34217728e8*1.);
}
else if (m == -6){
return (15.*std::sqrt(128942671./(2.*Pi))*(-9177.*z/r + 1245013.*std::pow(z/r,3.) - 48555507.*std::pow(z/r,5.) + 853189623.*std::pow(z/r,7.) - 8152700842.*std::pow(z/r,9.) + 46692741186.*std::pow(z/r,11.) - 168812218134.*std::pow(z/r,13.) + 393895175646.*std::pow(z/r,15.) - 590842763469.*std::pow(z/r,17.) + 549380113401.*std::pow(z/r,19.) - 287770535591.*std::pow(z/r,21.) + 64833677979.*std::pow(z/r,23.))*std::pow((x-I*y)/r,6.))/(6.7108864e7*1.);
}
else if (m == -5){
return (3.*std::sqrt(13538980455./Pi)*(437. - 183540.*std::pow(z/r,2.) + 12450130.*std::pow(z/r,4.) - 323703380.*std::pow(z/r,6.) + 4265948115.*std::pow(z/r,8.) - 32610803368.*std::pow(z/r,10.) + 155642470620.*std::pow(z/r,12.) - 482320623240.*std::pow(z/r,14.) + 984737939115.*std::pow(z/r,16.) - 1312983918820.*std::pow(z/r,18.) + 1098760226802.*std::pow(z/r,20.) - 523219155620.*std::pow(z/r,22.) + 108056129965.*std::pow(z/r,24.))*std::pow((x-I*y)/r,5.))/(2.68435456e8*1.);
}
else if (m == -4){
return (3.*std::sqrt(796410615./(2.*Pi))*(37145.*z/r - 5200300.*std::pow(z/r,3.) + 211652210.*std::pow(z/r,5.) - 3930683900.*std::pow(z/r,7.) + 40289509975.*std::pow(z/r,9.) - 251992571480.*std::pow(z/r,11.) + 1017662307900.*std::pow(z/r,13.) - 2733150198360.*std::pow(z/r,15.) + 4923689695575.*std::pow(z/r,17.) - 5873875426300.*std::pow(z/r,19.) + 4447362822770.*std::pow(z/r,21.) - 1933636009900.*std::pow(z/r,23.) + 367390841881.*std::pow(z/r,25.))*std::pow((x-I*y)/r,4.))/(1.34217728e8*1.);
}
else if (m == -3){
return (3.*std::sqrt(1856435./Pi)*(-37145. + 15935205.*std::pow(z/r,2.) - 1115464350.*std::pow(z/r,4.) + 30266266030.*std::pow(z/r,6.) - 421565848275.*std::pow(z/r,8.) + 3456839955855.*std::pow(z/r,10.) - 18017468860820.*std::pow(z/r,12.) + 62368161441300.*std::pow(z/r,14.) - 146565179387055.*std::pow(z/r,16.) + 234695875489075.*std::pow(z/r,18.) - 251989255788270.*std::pow(z/r,20.) + 173447150088030.*std::pow(z/r,22.) - 69127487353925.*std::pow(z/r,24.) + 12123897782073.*std::pow(z/r,26.))*std::pow((x-I*y)/r,3.))/(2.68435456e8*1.);
}
else if (m == -2){
return (std::sqrt(5569305./(2.*Pi))*(-334305.*z/r + 47805615.*std::pow(z/r,3.) - 2007835830.*std::pow(z/r,5.) + 38913770610.*std::pow(z/r,7.) - 421565848275.*std::pow(z/r,9.) + 2828323600245.*std::pow(z/r,11.) - 12473632288260.*std::pow(z/r,13.) + 37420896864780.*std::pow(z/r,15.) - 77593330263735.*std::pow(z/r,17.) + 111171730494825.*std::pow(z/r,19.) - 107995395337830.*std::pow(z/r,21.) + 67870623947490.*std::pow(z/r,23.) - 24885895447413.*std::pow(z/r,25.) + 4041299260691.*std::pow(z/r,27.))*std::pow((x-I*y)/r,2.))/(3.3554432e7*1.);
}
else if (m == -1){
return (std::sqrt(25665./(2.*Pi))*(334305. - 145088370.*std::pow(z/r,2.) + 10373818455.*std::pow(z/r,4.) - 290466916740.*std::pow(z/r,6.) + 4222144111185.*std::pow(z/r,8.) - 36591915630270.*std::pow(z/r,10.) + 204582073751055.*std::pow(z/r,12.) - 773365201872120.*std::pow(z/r,14.) + 2030083654914315.*std::pow(z/r,16.) - 3741722814940110.*std::pow(z/r,18.) + 4824853103475405.*std::pow(z/r,20.) - 4260909234238020.*std::pow(z/r,22.) + 2454654232767555.*std::pow(z/r,24.) - 830806048013634.*std::pow(z/r,26.) + 125280277081421.*std::pow(z/r,28.))*(x-I*y)/r)/(6.7108864e7*1.);
}
else if (m == 0){
return (std::sqrt(59./Pi)*(145422675.*z/r - 21037813650.*std::pow(z/r,3.) + 902522205585.*std::pow(z/r,5.) - 18050444111700.*std::pow(z/r,7.) + 204070298707275.*std::pow(z/r,9.) - 1447043936287950.*std::pow(z/r,11.) + 6845630929362225.*std::pow(z/r,13.) - 22427590854291480.*std::pow(z/r,15.) + 51946258228689825.*std::pow(z/r,17.) - 85665759184155150.*std::pow(z/r,19.) + 99943385714847675.*std::pow(z/r,21.) - 80586761604066900.*std::pow(z/r,23.) + 42710983650155457.*std::pow(z/r,25.) - 13385208551330770.*std::pow(z/r,27.) + 1879204156221315.*std::pow(z/r,29.)))/6.7108864e7;
}
else if (m == 1){
return -1.4901161193847656e-8*(1.*std::sqrt(25665./(2.*Pi))*(334305. - 145088370.*std::pow(z/r,2.) + 10373818455.*std::pow(z/r,4.) - 290466916740.*std::pow(z/r,6.) + 4222144111185.*std::pow(z/r,8.) - 36591915630270.*std::pow(z/r,10.) + 204582073751055.*std::pow(z/r,12.) - 773365201872120.*std::pow(z/r,14.) + 2030083654914315.*std::pow(z/r,16.) - 3741722814940110.*std::pow(z/r,18.) + 4824853103475405.*std::pow(z/r,20.) - 4260909234238020.*std::pow(z/r,22.) + 2454654232767555.*std::pow(z/r,24.) - 830806048013634.*std::pow(z/r,26.) + 125280277081421.*std::pow(z/r,28.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(5569305./(2.*Pi))*(-334305.*z/r + 47805615.*std::pow(z/r,3.) - 2007835830.*std::pow(z/r,5.) + 38913770610.*std::pow(z/r,7.) - 421565848275.*std::pow(z/r,9.) + 2828323600245.*std::pow(z/r,11.) - 12473632288260.*std::pow(z/r,13.) + 37420896864780.*std::pow(z/r,15.) - 77593330263735.*std::pow(z/r,17.) + 111171730494825.*std::pow(z/r,19.) - 107995395337830.*std::pow(z/r,21.) + 67870623947490.*std::pow(z/r,23.) - 24885895447413.*std::pow(z/r,25.) + 4041299260691.*std::pow(z/r,27.))*std::pow((x+I*y)/r,2.))/3.3554432e7;
}
else if (m == 3){
return (-3.*1.*std::sqrt(1856435./Pi)*(-37145. + 15935205.*std::pow(z/r,2.) - 1115464350.*std::pow(z/r,4.) + 30266266030.*std::pow(z/r,6.) - 421565848275.*std::pow(z/r,8.) + 3456839955855.*std::pow(z/r,10.) - 18017468860820.*std::pow(z/r,12.) + 62368161441300.*std::pow(z/r,14.) - 146565179387055.*std::pow(z/r,16.) + 234695875489075.*std::pow(z/r,18.) - 251989255788270.*std::pow(z/r,20.) + 173447150088030.*std::pow(z/r,22.) - 69127487353925.*std::pow(z/r,24.) + 12123897782073.*std::pow(z/r,26.))*std::pow((x+I*y)/r,3.))/2.68435456e8;
}
else if (m == 4){
return (3.*1.*std::sqrt(796410615./(2.*Pi))*(37145.*z/r - 5200300.*std::pow(z/r,3.) + 211652210.*std::pow(z/r,5.) - 3930683900.*std::pow(z/r,7.) + 40289509975.*std::pow(z/r,9.) - 251992571480.*std::pow(z/r,11.) + 1017662307900.*std::pow(z/r,13.) - 2733150198360.*std::pow(z/r,15.) + 4923689695575.*std::pow(z/r,17.) - 5873875426300.*std::pow(z/r,19.) + 4447362822770.*std::pow(z/r,21.) - 1933636009900.*std::pow(z/r,23.) + 367390841881.*std::pow(z/r,25.))*std::pow((x+I*y)/r,4.))/1.34217728e8;
}
else if (m == 5){
return (-3.*1.*std::sqrt(13538980455./Pi)*(437. - 183540.*std::pow(z/r,2.) + 12450130.*std::pow(z/r,4.) - 323703380.*std::pow(z/r,6.) + 4265948115.*std::pow(z/r,8.) - 32610803368.*std::pow(z/r,10.) + 155642470620.*std::pow(z/r,12.) - 482320623240.*std::pow(z/r,14.) + 984737939115.*std::pow(z/r,16.) - 1312983918820.*std::pow(z/r,18.) + 1098760226802.*std::pow(z/r,20.) - 523219155620.*std::pow(z/r,22.) + 108056129965.*std::pow(z/r,24.))*std::pow((x+I*y)/r,5.))/2.68435456e8;
}
else if (m == 6){
return (15.*1.*std::sqrt(128942671./(2.*Pi))*(-9177.*z/r + 1245013.*std::pow(z/r,3.) - 48555507.*std::pow(z/r,5.) + 853189623.*std::pow(z/r,7.) - 8152700842.*std::pow(z/r,9.) + 46692741186.*std::pow(z/r,11.) - 168812218134.*std::pow(z/r,13.) + 393895175646.*std::pow(z/r,15.) - 590842763469.*std::pow(z/r,17.) + 549380113401.*std::pow(z/r,19.) - 287770535591.*std::pow(z/r,21.) + 64833677979.*std::pow(z/r,23.))*std::pow((x+I*y)/r,6.))/6.7108864e7;
}
else if (m == 7){
return (-15.*1.*std::sqrt(2965681433./(2.*Pi))*(-133. + 54131.*std::pow(z/r,2.) - 3518515.*std::pow(z/r,4.) + 86555469.*std::pow(z/r,6.) - 1063395762.*std::pow(z/r,8.) + 7443770334.*std::pow(z/r,10.) - 31805200518.*std::pow(z/r,12.) + 85629386010.*std::pow(z/r,14.) - 145569956217.*std::pow(z/r,16.) + 151278581951.*std::pow(z/r,18.) - 87582336919.*std::pow(z/r,20.) + 21611225993.*std::pow(z/r,22.))*std::pow((x+I*y)/r,7.))/1.34217728e8;
}
else if (m == 8){
return (15.*1.*std::sqrt(9975473911./Pi)*(1463.*z/r - 190190.*std::pow(z/r,3.) + 7018011.*std::pow(z/r,5.) - 114961704.*std::pow(z/r,7.) + 1005914910.*std::pow(z/r,9.) - 5157600084.*std::pow(z/r,11.) + 16200154110.*std::pow(z/r,13.) - 31474585128.*std::pow(z/r,15.) + 36797492907.*std::pow(z/r,17.) - 23670901870.*std::pow(z/r,19.) + 6424959079.*std::pow(z/r,21.))*std::pow((x+I*y)/r,8.))/1.34217728e8;
}
else if (m == 9){
return (-5.*1.*std::sqrt(3980214090489./(2.*Pi))*(11. - 4290.*std::pow(z/r,2.) + 263835.*std::pow(z/r,4.) - 6050616.*std::pow(z/r,6.) + 68069430.*std::pow(z/r,8.) - 426568428.*std::pow(z/r,10.) + 1583473710.*std::pow(z/r,12.) - 3549765240.*std::pow(z/r,14.) + 4703438943.*std::pow(z/r,16.) - 3381557410.*std::pow(z/r,18.) + 1014467223.*std::pow(z/r,20.))*std::pow((x+I*y)/r,9.))/1.34217728e8;
}
else if (m == 10){
return (3.*1.*std::sqrt(510283857755./(2.*Pi))*z/r*(-715. + 87945.*std::pow(z/r,2.) - 3025308.*std::pow(z/r,4.) + 45379620.*std::pow(z/r,6.) - 355473690.*std::pow(z/r,8.) + 1583473710.*std::pow(z/r,10.) - 4141392780.*std::pow(z/r,12.) + 6271251924.*std::pow(z/r,14.) - 5072336115.*std::pow(z/r,16.) + 1690778705.*std::pow(z/r,18.))*std::pow((x+I*y)/r,10.))/6.7108864e7;
}
else if (m == 11){
return (-15.*1.*std::sqrt(5371409029./Pi)*(-143. + 52767.*std::pow(z/r,2.) - 3025308.*std::pow(z/r,4.) + 63531468.*std::pow(z/r,6.) - 639852642.*std::pow(z/r,8.) + 3483642162.*std::pow(z/r,10.) - 10767621228.*std::pow(z/r,12.) + 18813755772.*std::pow(z/r,14.) - 17245942791.*std::pow(z/r,16.) + 6424959079.*std::pow(z/r,18.))*std::pow((x+I*y)/r,11.))/2.68435456e8;
}
else if (m == 12){
return (15.*1.*std::sqrt(220227770189./(2.*Pi))*z/r*(429. - 49192.*std::pow(z/r,2.) + 1549548.*std::pow(z/r,4.) - 20808216.*std::pow(z/r,6.) + 141611470.*std::pow(z/r,8.) - 525249816.*std::pow(z/r,10.) + 1070701548.*std::pow(z/r,12.) - 1121687336.*std::pow(z/r,14.) + 470118957.*std::pow(z/r,16.))*std::pow((x+I*y)/r,12.))/1.34217728e8;
}
else if (m == 13){
return (-15.*1.*std::sqrt(5551960593./Pi)*(143. - 49192.*std::pow(z/r,2.) + 2582580.*std::pow(z/r,4.) - 48552504.*std::pow(z/r,6.) + 424834410.*std::pow(z/r,8.) - 1925915992.*std::pow(z/r,10.) + 4639706708.*std::pow(z/r,12.) - 5608436680.*std::pow(z/r,14.) + 2664007423.*std::pow(z/r,16.))*std::pow((x+I*y)/r,13.))/2.68435456e8;
}
else if (m == 14){
return (15.*1.*std::sqrt(238734305499./Pi)*z/r*(-143. + 15015.*std::pow(z/r,2.) - 423423.*std::pow(z/r,4.) + 4939935.*std::pow(z/r,6.) - 27992965.*std::pow(z/r,8.) + 80925117.*std::pow(z/r,10.) - 114125165.*std::pow(z/r,12.) + 61953661.*std::pow(z/r,14.))*std::pow((x+I*y)/r,14.))/6.7108864e7;
}
else if (m == 15){
return (-3.*1.*std::sqrt(4376795600815./Pi)*(-13. + 4095.*std::pow(z/r,2.) - 192465.*std::pow(z/r,4.) + 3143595.*std::pow(z/r,6.) - 22903335.*std::pow(z/r,8.) + 80925117.*std::pow(z/r,10.) - 134875195.*std::pow(z/r,12.) + 84482265.*std::pow(z/r,14.))*std::pow((x+I*y)/r,15.))/1.34217728e8;
}
else if (m == 16){
return (15.*1.*std::sqrt(6127513841141./(2.*Pi))*z/r*(39. - 3666.*std::pow(z/r,2.) + 89817.*std::pow(z/r,4.) - 872508.*std::pow(z/r,6.) + 3853577.*std::pow(z/r,8.) - 7707154.*std::pow(z/r,10.) + 5632151.*std::pow(z/r,12.))*std::pow((x+I*y)/r,16.))/6.7108864e7;
}
else if (m == 17){
return (-15.*1.*std::sqrt(3463377388471./Pi)*(3. - 846.*std::pow(z/r,2.) + 34545.*std::pow(z/r,4.) - 469812.*std::pow(z/r,6.) + 2667861.*std::pow(z/r,8.) - 6521438.*std::pow(z/r,10.) + 5632151.*std::pow(z/r,12.))*std::pow((x+I*y)/r,17.))/1.34217728e8;
}
else if (m == 18){
return (5.*1.*std::sqrt(488336211774411./Pi)*z/r*(-9. + 735.*std::pow(z/r,2.) - 14994.*std::pow(z/r,4.) + 113526.*std::pow(z/r,6.) - 346885.*std::pow(z/r,8.) + 359499.*std::pow(z/r,10.))*std::pow((x+I*y)/r,18.))/6.7108864e7;
}
else if (m == 19){
return (-15.*1.*std::sqrt(14798067023467./Pi)*(-3. + 735.*std::pow(z/r,2.) - 24990.*std::pow(z/r,4.) + 264894.*std::pow(z/r,6.) - 1040655.*std::pow(z/r,8.) + 1318163.*std::pow(z/r,10.))*std::pow((x+I*y)/r,19.))/2.68435456e8;
}
else if (m == 20){
return (3.*1.*std::sqrt(73990335117335./(2.*Pi))*z/r*(105. - 7140.*std::pow(z/r,2.) + 113526.*std::pow(z/r,4.) - 594660.*std::pow(z/r,6.) + 941545.*std::pow(z/r,8.))*std::pow((x+I*y)/r,20.))/1.34217728e8;
}
else if (m == 21){
return (-3.*1.*std::sqrt(73990335117335./Pi)*(7. - 1428.*std::pow(z/r,2.) + 37842.*std::pow(z/r,4.) - 277508.*std::pow(z/r,6.) + 564927.*std::pow(z/r,8.))*std::pow((x+I*y)/r,21.))/2.68435456e8;
}
else if (m == 22){
return (3.*1.*std::sqrt(3773507090984085./(2.*Pi))*z/r*(-7. + 371.*std::pow(z/r,2.) - 4081.*std::pow(z/r,4.) + 11077.*std::pow(z/r,6.))*std::pow((x+I*y)/r,22.))/6.7108864e7;
}
else if (m == 23){
return (-21.*1.*std::sqrt(41467110889935./(2.*Pi))*(-1. + 159.*std::pow(z/r,2.) - 2915.*std::pow(z/r,4.) + 11077.*std::pow(z/r,6.))*std::pow((x+I*y)/r,23.))/1.34217728e8;
}
else if (m == 24){
return (21.*1.*std::sqrt(732585625722185./Pi)*z/r*(3. - 110.*std::pow(z/r,2.) + 627.*std::pow(z/r,4.))*std::pow((x+I*y)/r,24.))/1.34217728e8;
}
else if (m == 25){
return (-7.*1.*std::sqrt(439551375433311./(2.*Pi))*(1. - 110.*std::pow(z/r,2.) + 1045.*std::pow(z/r,4.))*std::pow((x+I*y)/r,25.))/1.34217728e8;
}
else if (m == 26){
return (7.*1.*std::sqrt(24175325648832105./(2.*Pi))*z/r*(-1. + 19.*std::pow(z/r,2.))*std::pow((x+I*y)/r,26.))/6.7108864e7;
}
else if (m == 27){
return -3.725290298461914e-9*(1.*std::sqrt(56409093180608245./Pi)*(-1. + 57.*std::pow(z/r,2.))*std::pow((x+I*y)/r,27.));
}
else if (m == 28){
return (1.*std::sqrt(3215318311294669965./(2.*Pi))*z/r*std::pow((x+I*y)/r,28.))/1.34217728e8;
}
else if (m == 29){
return -3.725290298461914e-9*(1.*std::sqrt(110873045217057585./Pi)*std::pow((x+I*y)/r,29.));
}
else{return 0.;}
}

else if (l == 30){
if(m == -30){
return (std::sqrt(450883717216034179./Pi)*std::pow((x-I*y)/r,30.))/(5.36870912e8*1.);
}
else if (m == -29){
return (std::sqrt(6763255758240512685./Pi)*z/r*std::pow((x-I*y)/r,29.))/(2.68435456e8*1.);
}
else if (m == -28){
return (std::sqrt(114631453529500215./(2.*Pi))*(-1. + 59.*std::pow(z/r,2.))*std::pow((x-I*y)/r,28.))/(2.68435456e8*1.);
}
else if (m == -27){
return (std::sqrt(1108104050785168745./Pi)*z/r*(-3. + 59.*std::pow(z/r,2.))*std::pow((x-I*y)/r,27.))/(2.68435456e8*1.);
}
else if (m == -26){
return (std::sqrt(174963797492395065./Pi)*(1. - 114.*std::pow(z/r,2.) + 1121.*std::pow(z/r,4.))*std::pow((x-I*y)/r,26.))/(5.36870912e8*1.);
}
else if (m == -25){
return (7.*std::sqrt(4998965642639859./(2.*Pi))*z/r*(5. - 190.*std::pow(z/r,2.) + 1121.*std::pow(z/r,4.))*std::pow((x-I*y)/r,25.))/(1.34217728e8*1.);
}
else if (m == -24){
return (7.*std::sqrt(757419036763615./Pi)*(-1. + 165.*std::pow(z/r,2.) - 3135.*std::pow(z/r,4.) + 12331.*std::pow(z/r,6.))*std::pow((x-I*y)/r,24.))/(2.68435456e8*1.);
}
else if (m == -23){
return (3.*std::sqrt(15905799772035915./(2.*Pi))*z/r*(-7. + 385.*std::pow(z/r,2.) - 4389.*std::pow(z/r,4.) + 12331.*std::pow(z/r,6.))*std::pow((x-I*y)/r,23.))/(1.34217728e8*1.);
}
else if (m == -22){
return (3.*std::sqrt(300109429661055./Pi)*(7. - 1484.*std::pow(z/r,2.) + 40810.*std::pow(z/r,4.) - 310156.*std::pow(z/r,6.) + 653543.*std::pow(z/r,8.))*std::pow((x-I*y)/r,22.))/(5.36870912e8*1.);
}
else if (m == -21){
return (std::sqrt(3901422585593715./Pi)*z/r*(63. - 4452.*std::pow(z/r,2.) + 73458.*std::pow(z/r,4.) - 398772.*std::pow(z/r,6.) + 653543.*std::pow(z/r,8.))*std::pow((x-I*y)/r,21.))/(2.68435456e8*1.);
}
else if (m == -20){
return (3.*std::sqrt(15299696414093./(2.*Pi))*(-21. + 5355.*std::pow(z/r,2.) - 189210.*std::pow(z/r,4.) + 2081310.*std::pow(z/r,6.) - 8473905.*std::pow(z/r,8.) + 11110231.*std::pow(z/r,10.))*std::pow((x-I*y)/r,20.))/(2.68435456e8*1.);
}
else if (m == -19){
return (15.*std::sqrt(168296660555023./Pi)*z/r*(-21. + 1785.*std::pow(z/r,2.) - 37842.*std::pow(z/r,4.) + 297330.*std::pow(z/r,6.) - 941545.*std::pow(z/r,8.) + 1010021.*std::pow(z/r,10.))*std::pow((x-I*y)/r,19.))/(2.68435456e8*1.);
}
else if (m == -18){
return (5.*std::sqrt(504889981665069./Pi)*(3. - 882.*std::pow(z/r,2.) + 37485.*std::pow(z/r,4.) - 529788.*std::pow(z/r,6.) + 3121965.*std::pow(z/r,8.) - 7908978.*std::pow(z/r,10.) + 7070147.*std::pow(z/r,12.))*std::pow((x-I*y)/r,18.))/(5.36870912e8*1.);
}
else if (m == -17){
return (15.*std::sqrt(12945896965771./Pi)*z/r*(39. - 3822.*std::pow(z/r,2.) + 97461.*std::pow(z/r,4.) - 983892.*std::pow(z/r,6.) + 4509505.*std::pow(z/r,8.) - 9346974.*std::pow(z/r,10.) + 7070147.*std::pow(z/r,12.))*std::pow((x-I*y)/r,17.))/(1.34217728e8*1.);
}
else if (m == -16){
return (15.*std::sqrt(39349230899./(2.*Pi))*(-39. + 12831.*std::pow(z/r,2.) - 628719.*std::pow(z/r,4.) + 10688223.*std::pow(z/r,6.) - 80925117.*std::pow(z/r,8.) + 296725429.*std::pow(z/r,10.) - 512525741.*std::pow(z/r,12.) + 332296909.*std::pow(z/r,14.))*std::pow((x-I*y)/r,16.))/(1.34217728e8*1.);
}
else if (m == -15){
return (std::sqrt(13575484660155./Pi)*z/r*(-585. + 64155.*std::pow(z/r,2.) - 1886157.*std::pow(z/r,4.) + 22903335.*std::pow(z/r,6.) - 134875195.*std::pow(z/r,8.) + 404625585.*std::pow(z/r,10.) - 591375855.*std::pow(z/r,12.) + 332296909.*std::pow(z/r,14.))*std::pow((x-I*y)/r,15.))/(1.34217728e8*1.);
}
else if (m == -14){
return (15.*std::sqrt(2715096932031./Pi)*(13. - 4680.*std::pow(z/r,2.) + 256620.*std::pow(z/r,4.) - 5029752.*std::pow(z/r,6.) + 45806670.*std::pow(z/r,8.) - 215800312.*std::pow(z/r,10.) + 539500780.*std::pow(z/r,12.) - 675858120.*std::pow(z/r,14.) + 332296909.*std::pow(z/r,16.))*std::pow((x-I*y)/r,14.))/(5.36870912e8*1.);
}
else if (m == -13){
return (15.*std::sqrt(4196058894957./Pi)*z/r*(143. - 17160.*std::pow(z/r,2.) + 564564.*std::pow(z/r,4.) - 7903896.*std::pow(z/r,6.) + 55985930.*std::pow(z/r,8.) - 215800312.*std::pow(z/r,10.) + 456500660.*std::pow(z/r,12.) - 495629288.*std::pow(z/r,14.) + 215015647.*std::pow(z/r,16.))*std::pow((x-I*y)/r,13.))/(2.68435456e8*1.);
}
else if (m == -12){
return (5.*std::sqrt(97582764999./(2.*Pi))*(-143. + 55341.*std::pow(z/r,2.) - 3320460.*std::pow(z/r,4.) + 72828756.*std::pow(z/r,6.) - 764701938.*std::pow(z/r,8.) + 4333310982.*std::pow(z/r,10.) - 13919120124.*std::pow(z/r,12.) + 25237965060.*std::pow(z/r,14.) - 23976066807.*std::pow(z/r,16.) + 9245672821.*std::pow(z/r,18.))*std::pow((x-I*y)/r,12.))/(2.68435456e8*1.);
}
else if (m == -11){
return (15.*std::sqrt(4326169248289./Pi)*z/r*(-143. + 18447.*std::pow(z/r,2.) - 664092.*std::pow(z/r,4.) + 10404108.*std::pow(z/r,6.) - 84966882.*std::pow(z/r,8.) + 393937362.*std::pow(z/r,10.) - 1070701548.*std::pow(z/r,12.) + 1682531004.*std::pow(z/r,14.) - 1410356871.*std::pow(z/r,16.) + 486614359.*std::pow(z/r,18.))*std::pow((x-I*y)/r,11.))/(2.68435456e8*1.);
}
else if (m == -10){
return (3.*std::sqrt(527581615645./Pi)*(143. - 58630.*std::pow(z/r,2.) + 3781635.*std::pow(z/r,4.) - 90759240.*std::pow(z/r,6.) + 1066421070.*std::pow(z/r,8.) - 6967284324.*std::pow(z/r,10.) + 26919053070.*std::pow(z/r,12.) - 62712519240.*std::pow(z/r,14.) + 86229713955.*std::pow(z/r,16.) - 64249590790.*std::pow(z/r,18.) + 19951188719.*std::pow(z/r,20.))*std::pow((x-I*y)/r,10.))/(5.36870912e8*1.);
}
else if (m == -9){
return (5.*std::sqrt(45221281341./(2.*Pi))*(3003.*z/r - 410410.*std::pow(z/r,3.) + 15882867.*std::pow(z/r,5.) - 272277720.*std::pow(z/r,7.) + 2488315830.*std::pow(z/r,9.) - 13301179164.*std::pow(z/r,11.) + 43484624190.*std::pow(z/r,13.) - 87797526936.*std::pow(z/r,15.) + 106519058415.*std::pow(z/r,17.) - 71012705610.*std::pow(z/r,19.) + 19951188719.*std::pow(z/r,21.))*std::pow((x-I*y)/r,9.))/(1.34217728e8*1.);
}
else if (m == -8){
return (15.*std::sqrt(2155547743921./Pi)*(-7. + 3003.*std::pow(z/r,2.) - 205205.*std::pow(z/r,4.) + 5294289.*std::pow(z/r,6.) - 68069430.*std::pow(z/r,8.) + 497663166.*std::pow(z/r,10.) - 2216863194.*std::pow(z/r,12.) + 6212089170.*std::pow(z/r,14.) - 10974690867.*std::pow(z/r,16.) + 11835450935.*std::pow(z/r,18.) - 7101270561.*std::pow(z/r,20.) + 1813744429.*std::pow(z/r,22.))*std::pow((x-I*y)/r,8.))/(2.68435456e8*1.);
}
else if (m == -7){
return (15.*std::sqrt(4932603533./(2.*Pi))*(-3059.*z/r + 437437.*std::pow(z/r,3.) - 17934917.*std::pow(z/r,5.) + 330514899.*std::pow(z/r,7.) - 3305148990.*std::pow(z/r,9.) + 19770800322.*std::pow(z/r,11.) - 74520708906.*std::pow(z/r,13.) + 180978864486.*std::pow(z/r,15.) - 282114112287.*std::pow(z/r,17.) + 272215371505.*std::pow(z/r,19.) - 147774058817.*std::pow(z/r,21.) + 34461144151.*std::pow(z/r,23.))*std::pow((x-I*y)/r,7.))/(1.34217728e8*1.);
}
else if (m == -6){
return (5.*std::sqrt(399940827./Pi)*(3059. - 1358196.*std::pow(z/r,2.) + 97111014.*std::pow(z/r,4.) - 2654367716.*std::pow(z/r,6.) + 36687153789.*std::pow(z/r,8.) - 293497230312.*std::pow(z/r,10.) + 1463039223828.*std::pow(z/r,12.) - 4726742107752.*std::pow(z/r,14.) + 10044326978973.*std::pow(z/r,16.) - 13917629539492.*std::pow(z/r,18.) + 12086362494822.*std::pow(z/r,20.) - 5964698374068.*std::pow(z/r,22.) + 1275062333587.*std::pow(z/r,24.))*std::pow((x-I*y)/r,6.))/(5.36870912e8*1.);
}
else if (m == -5){
return (3.*std::sqrt(399940827./Pi)*(76475.*z/r - 11318300.*std::pow(z/r,3.) + 485555070.*std::pow(z/r,5.) - 9479884700.*std::pow(z/r,7.) + 101908760525.*std::pow(z/r,9.) - 667039159800.*std::pow(z/r,11.) + 2813536968900.*std::pow(z/r,13.) - 7877903512920.*std::pow(z/r,15.) + 14771069086725.*std::pow(z/r,17.) - 18312670446700.*std::pow(z/r,19.) + 14388526779550.*std::pow(z/r,21.) - 6483367797900.*std::pow(z/r,23.) + 1275062333587.*std::pow(z/r,25.))*std::pow((x-I*y)/r,5.))/(2.68435456e8*1.);
}
else if (m == -4){
return (3.*std::sqrt(1076763765./(2.*Pi))*(-2185. + 994175.*std::pow(z/r,2.) - 73568950.*std::pow(z/r,4.) + 2104071970.*std::pow(z/r,6.) - 30809625275.*std::pow(z/r,8.) + 264962777365.*std::pow(z/r,10.) - 1445251512900.*std::pow(z/r,12.) + 5225140085100.*std::pow(z/r,14.) - 12801593208495.*std::pow(z/r,16.) + 21335988680825.*std::pow(z/r,18.) - 23806471580710.*std::pow(z/r,20.) + 17004622557650.*std::pow(z/r,22.) - 7023648447725.*std::pow(z/r,24.) + 1275062333587.*std::pow(z/r,26.))*std::pow((x-I*y)/r,4.))/(2.68435456e8*1.);
}
else if (m == -3){
return (std::sqrt(21113015./Pi)*(-1002915.*z/r + 152108775.*std::pow(z/r,3.) - 6753629610.*std::pow(z/r,5.) + 137967004890.*std::pow(z/r,7.) - 1571290889025.*std::pow(z/r,9.) + 11056174073685.*std::pow(z/r,11.) - 51028495724700.*std::pow(z/r,13.) + 159889286604060.*std::pow(z/r,15.) - 345643016629365.*std::pow(z/r,17.) + 515432568657825.*std::pow(z/r,19.) - 520341450264090.*std::pow(z/r,21.) + 339353119737450.*std::pow(z/r,23.) - 128954185500231.*std::pow(z/r,25.) + 21676059670979.*std::pow(z/r,27.))*std::pow((x-I*y)/r,3.))/(2.68435456e8*1.);
}
else if (m == -2){
return (std::sqrt(822585./Pi)*(334305. - 154448910.*std::pow(z/r,2.) + 11712375675.*std::pow(z/r,4.) - 346686319980.*std::pow(z/r,6.) + 5311729688265.*std::pow(z/r,8.) - 48395759381970.*std::pow(z/r,10.) + 283775134557915.*std::pow(z/r,12.) - 1122626905943400.*std::pow(z/r,14.) + 3077868767128155.*std::pow(z/r,16.) - 5914336062324690.*std::pow(z/r,18.) + 7937661557330505.*std::pow(z/r,20.) - 7284780303697260.*std::pow(z/r,22.) + 4355031703297275.*std::pow(z/r,24.) - 1527611120541198.*std::pow(z/r,26.) + 238436656380769.*std::pow(z/r,28.))*std::pow((x-I*y)/r,2.))/(5.36870912e8*1.);
}
else if (m == -1){
return (std::sqrt(28365./(2.*Pi))*(9694845.*z/r - 1493006130.*std::pow(z/r,3.) + 67931778915.*std::pow(z/r,5.) - 1436271897060.*std::pow(z/r,7.) + 17115573439965.*std::pow(z/r,9.) - 127588820188830.*std::pow(z/r,11.) + 633036838629195.*std::pow(z/r,13.) - 2170412018157240.*std::pow(z/r,15.) + 5250482014512735.*std::pow(z/r,17.) - 9027144516179790.*std::pow(z/r,19.) + 10961532626789745.*std::pow(z/r,21.) - 9185157774226980.*std::pow(z/r,23.) + 5051836775824839.*std::pow(z/r,25.) - 1640767499840546.*std::pow(z/r,27.) + 238436656380769.*std::pow(z/r,29.))*(x-I*y)/r)/(6.7108864e7*1.);
}
else if (m == 0){
return (std::sqrt(61./Pi)*(-9694845. + 4508102925.*std::pow(z/r,2.) - 347123925225.*std::pow(z/r,4.) + 10529425731825.*std::pow(z/r,6.) - 166966608033225.*std::pow(z/r,8.) + 1591748329916745.*std::pow(z/r,10.) - 9888133564634325.*std::pow(z/r,12.) + 42051732851796525.*std::pow(z/r,14.) - 126155198555389575.*std::pow(z/r,16.) + 271274904083157975.*std::pow(z/r,18.) - 419762220002360235.*std::pow(z/r,20.) + 463373879223384675.*std::pow(z/r,22.) - 355924863751295475.*std::pow(z/r,24.) + 180700315442965395.*std::pow(z/r,26.) - 54496920530418135.*std::pow(z/r,28.) + 7391536347803839.*std::pow(z/r,30.)))/1.34217728e8;
}
else if (m == 1){
return -1.4901161193847656e-8*(1.*std::sqrt(28365./(2.*Pi))*(9694845.*z/r - 1493006130.*std::pow(z/r,3.) + 67931778915.*std::pow(z/r,5.) - 1436271897060.*std::pow(z/r,7.) + 17115573439965.*std::pow(z/r,9.) - 127588820188830.*std::pow(z/r,11.) + 633036838629195.*std::pow(z/r,13.) - 2170412018157240.*std::pow(z/r,15.) + 5250482014512735.*std::pow(z/r,17.) - 9027144516179790.*std::pow(z/r,19.) + 10961532626789745.*std::pow(z/r,21.) - 9185157774226980.*std::pow(z/r,23.) + 5051836775824839.*std::pow(z/r,25.) - 1640767499840546.*std::pow(z/r,27.) + 238436656380769.*std::pow(z/r,29.))*(x+I*y)/r);
}
else if (m == 2){
return (1.*std::sqrt(822585./Pi)*(334305. - 154448910.*std::pow(z/r,2.) + 11712375675.*std::pow(z/r,4.) - 346686319980.*std::pow(z/r,6.) + 5311729688265.*std::pow(z/r,8.) - 48395759381970.*std::pow(z/r,10.) + 283775134557915.*std::pow(z/r,12.) - 1122626905943400.*std::pow(z/r,14.) + 3077868767128155.*std::pow(z/r,16.) - 5914336062324690.*std::pow(z/r,18.) + 7937661557330505.*std::pow(z/r,20.) - 7284780303697260.*std::pow(z/r,22.) + 4355031703297275.*std::pow(z/r,24.) - 1527611120541198.*std::pow(z/r,26.) + 238436656380769.*std::pow(z/r,28.))*std::pow((x+I*y)/r,2.))/5.36870912e8;
}
else if (m == 3){
return -3.725290298461914e-9*(1.*std::sqrt(21113015./Pi)*(-1002915.*z/r + 152108775.*std::pow(z/r,3.) - 6753629610.*std::pow(z/r,5.) + 137967004890.*std::pow(z/r,7.) - 1571290889025.*std::pow(z/r,9.) + 11056174073685.*std::pow(z/r,11.) - 51028495724700.*std::pow(z/r,13.) + 159889286604060.*std::pow(z/r,15.) - 345643016629365.*std::pow(z/r,17.) + 515432568657825.*std::pow(z/r,19.) - 520341450264090.*std::pow(z/r,21.) + 339353119737450.*std::pow(z/r,23.) - 128954185500231.*std::pow(z/r,25.) + 21676059670979.*std::pow(z/r,27.))*std::pow((x+I*y)/r,3.));
}
else if (m == 4){
return (3.*1.*std::sqrt(1076763765./(2.*Pi))*(-2185. + 994175.*std::pow(z/r,2.) - 73568950.*std::pow(z/r,4.) + 2104071970.*std::pow(z/r,6.) - 30809625275.*std::pow(z/r,8.) + 264962777365.*std::pow(z/r,10.) - 1445251512900.*std::pow(z/r,12.) + 5225140085100.*std::pow(z/r,14.) - 12801593208495.*std::pow(z/r,16.) + 21335988680825.*std::pow(z/r,18.) - 23806471580710.*std::pow(z/r,20.) + 17004622557650.*std::pow(z/r,22.) - 7023648447725.*std::pow(z/r,24.) + 1275062333587.*std::pow(z/r,26.))*std::pow((x+I*y)/r,4.))/2.68435456e8;
}
else if (m == 5){
return (-3.*1.*std::sqrt(399940827./Pi)*(76475.*z/r - 11318300.*std::pow(z/r,3.) + 485555070.*std::pow(z/r,5.) - 9479884700.*std::pow(z/r,7.) + 101908760525.*std::pow(z/r,9.) - 667039159800.*std::pow(z/r,11.) + 2813536968900.*std::pow(z/r,13.) - 7877903512920.*std::pow(z/r,15.) + 14771069086725.*std::pow(z/r,17.) - 18312670446700.*std::pow(z/r,19.) + 14388526779550.*std::pow(z/r,21.) - 6483367797900.*std::pow(z/r,23.) + 1275062333587.*std::pow(z/r,25.))*std::pow((x+I*y)/r,5.))/2.68435456e8;
}
else if (m == 6){
return (5.*1.*std::sqrt(399940827./Pi)*(3059. - 1358196.*std::pow(z/r,2.) + 97111014.*std::pow(z/r,4.) - 2654367716.*std::pow(z/r,6.) + 36687153789.*std::pow(z/r,8.) - 293497230312.*std::pow(z/r,10.) + 1463039223828.*std::pow(z/r,12.) - 4726742107752.*std::pow(z/r,14.) + 10044326978973.*std::pow(z/r,16.) - 13917629539492.*std::pow(z/r,18.) + 12086362494822.*std::pow(z/r,20.) - 5964698374068.*std::pow(z/r,22.) + 1275062333587.*std::pow(z/r,24.))*std::pow((x+I*y)/r,6.))/5.36870912e8;
}
else if (m == 7){
return (-15.*1.*std::sqrt(4932603533./(2.*Pi))*(-3059.*z/r + 437437.*std::pow(z/r,3.) - 17934917.*std::pow(z/r,5.) + 330514899.*std::pow(z/r,7.) - 3305148990.*std::pow(z/r,9.) + 19770800322.*std::pow(z/r,11.) - 74520708906.*std::pow(z/r,13.) + 180978864486.*std::pow(z/r,15.) - 282114112287.*std::pow(z/r,17.) + 272215371505.*std::pow(z/r,19.) - 147774058817.*std::pow(z/r,21.) + 34461144151.*std::pow(z/r,23.))*std::pow((x+I*y)/r,7.))/1.34217728e8;
}
else if (m == 8){
return (15.*1.*std::sqrt(2155547743921./Pi)*(-7. + 3003.*std::pow(z/r,2.) - 205205.*std::pow(z/r,4.) + 5294289.*std::pow(z/r,6.) - 68069430.*std::pow(z/r,8.) + 497663166.*std::pow(z/r,10.) - 2216863194.*std::pow(z/r,12.) + 6212089170.*std::pow(z/r,14.) - 10974690867.*std::pow(z/r,16.) + 11835450935.*std::pow(z/r,18.) - 7101270561.*std::pow(z/r,20.) + 1813744429.*std::pow(z/r,22.))*std::pow((x+I*y)/r,8.))/2.68435456e8;
}
else if (m == 9){
return (-5.*1.*std::sqrt(45221281341./(2.*Pi))*(3003.*z/r - 410410.*std::pow(z/r,3.) + 15882867.*std::pow(z/r,5.) - 272277720.*std::pow(z/r,7.) + 2488315830.*std::pow(z/r,9.) - 13301179164.*std::pow(z/r,11.) + 43484624190.*std::pow(z/r,13.) - 87797526936.*std::pow(z/r,15.) + 106519058415.*std::pow(z/r,17.) - 71012705610.*std::pow(z/r,19.) + 19951188719.*std::pow(z/r,21.))*std::pow((x+I*y)/r,9.))/1.34217728e8;
}
else if (m == 10){
return (3.*1.*std::sqrt(527581615645./Pi)*(143. - 58630.*std::pow(z/r,2.) + 3781635.*std::pow(z/r,4.) - 90759240.*std::pow(z/r,6.) + 1066421070.*std::pow(z/r,8.) - 6967284324.*std::pow(z/r,10.) + 26919053070.*std::pow(z/r,12.) - 62712519240.*std::pow(z/r,14.) + 86229713955.*std::pow(z/r,16.) - 64249590790.*std::pow(z/r,18.) + 19951188719.*std::pow(z/r,20.))*std::pow((x+I*y)/r,10.))/5.36870912e8;
}
else if (m == 11){
return (-15.*1.*std::sqrt(4326169248289./Pi)*z/r*(-143. + 18447.*std::pow(z/r,2.) - 664092.*std::pow(z/r,4.) + 10404108.*std::pow(z/r,6.) - 84966882.*std::pow(z/r,8.) + 393937362.*std::pow(z/r,10.) - 1070701548.*std::pow(z/r,12.) + 1682531004.*std::pow(z/r,14.) - 1410356871.*std::pow(z/r,16.) + 486614359.*std::pow(z/r,18.))*std::pow((x+I*y)/r,11.))/2.68435456e8;
}
else if (m == 12){
return (5.*1.*std::sqrt(97582764999./(2.*Pi))*(-143. + 55341.*std::pow(z/r,2.) - 3320460.*std::pow(z/r,4.) + 72828756.*std::pow(z/r,6.) - 764701938.*std::pow(z/r,8.) + 4333310982.*std::pow(z/r,10.) - 13919120124.*std::pow(z/r,12.) + 25237965060.*std::pow(z/r,14.) - 23976066807.*std::pow(z/r,16.) + 9245672821.*std::pow(z/r,18.))*std::pow((x+I*y)/r,12.))/2.68435456e8;
}
else if (m == 13){
return (-15.*1.*std::sqrt(4196058894957./Pi)*z/r*(143. - 17160.*std::pow(z/r,2.) + 564564.*std::pow(z/r,4.) - 7903896.*std::pow(z/r,6.) + 55985930.*std::pow(z/r,8.) - 215800312.*std::pow(z/r,10.) + 456500660.*std::pow(z/r,12.) - 495629288.*std::pow(z/r,14.) + 215015647.*std::pow(z/r,16.))*std::pow((x+I*y)/r,13.))/2.68435456e8;
}
else if (m == 14){
return (15.*1.*std::sqrt(2715096932031./Pi)*(13. - 4680.*std::pow(z/r,2.) + 256620.*std::pow(z/r,4.) - 5029752.*std::pow(z/r,6.) + 45806670.*std::pow(z/r,8.) - 215800312.*std::pow(z/r,10.) + 539500780.*std::pow(z/r,12.) - 675858120.*std::pow(z/r,14.) + 332296909.*std::pow(z/r,16.))*std::pow((x+I*y)/r,14.))/5.36870912e8;
}
else if (m == 15){
return -7.450580596923828e-9*(1.*std::sqrt(13575484660155./Pi)*z/r*(-585. + 64155.*std::pow(z/r,2.) - 1886157.*std::pow(z/r,4.) + 22903335.*std::pow(z/r,6.) - 134875195.*std::pow(z/r,8.) + 404625585.*std::pow(z/r,10.) - 591375855.*std::pow(z/r,12.) + 332296909.*std::pow(z/r,14.))*std::pow((x+I*y)/r,15.));
}
else if (m == 16){
return (15.*1.*std::sqrt(39349230899./(2.*Pi))*(-39. + 12831.*std::pow(z/r,2.) - 628719.*std::pow(z/r,4.) + 10688223.*std::pow(z/r,6.) - 80925117.*std::pow(z/r,8.) + 296725429.*std::pow(z/r,10.) - 512525741.*std::pow(z/r,12.) + 332296909.*std::pow(z/r,14.))*std::pow((x+I*y)/r,16.))/1.34217728e8;
}
else if (m == 17){
return (-15.*1.*std::sqrt(12945896965771./Pi)*z/r*(39. - 3822.*std::pow(z/r,2.) + 97461.*std::pow(z/r,4.) - 983892.*std::pow(z/r,6.) + 4509505.*std::pow(z/r,8.) - 9346974.*std::pow(z/r,10.) + 7070147.*std::pow(z/r,12.))*std::pow((x+I*y)/r,17.))/1.34217728e8;
}
else if (m == 18){
return (5.*1.*std::sqrt(504889981665069./Pi)*(3. - 882.*std::pow(z/r,2.) + 37485.*std::pow(z/r,4.) - 529788.*std::pow(z/r,6.) + 3121965.*std::pow(z/r,8.) - 7908978.*std::pow(z/r,10.) + 7070147.*std::pow(z/r,12.))*std::pow((x+I*y)/r,18.))/5.36870912e8;
}
else if (m == 19){
return (-15.*1.*std::sqrt(168296660555023./Pi)*z/r*(-21. + 1785.*std::pow(z/r,2.) - 37842.*std::pow(z/r,4.) + 297330.*std::pow(z/r,6.) - 941545.*std::pow(z/r,8.) + 1010021.*std::pow(z/r,10.))*std::pow((x+I*y)/r,19.))/2.68435456e8;
}
else if (m == 20){
return (3.*1.*std::sqrt(15299696414093./(2.*Pi))*(-21. + 5355.*std::pow(z/r,2.) - 189210.*std::pow(z/r,4.) + 2081310.*std::pow(z/r,6.) - 8473905.*std::pow(z/r,8.) + 11110231.*std::pow(z/r,10.))*std::pow((x+I*y)/r,20.))/2.68435456e8;
}
else if (m == 21){
return -3.725290298461914e-9*(1.*std::sqrt(3901422585593715./Pi)*z/r*(63. - 4452.*std::pow(z/r,2.) + 73458.*std::pow(z/r,4.) - 398772.*std::pow(z/r,6.) + 653543.*std::pow(z/r,8.))*std::pow((x+I*y)/r,21.));
}
else if (m == 22){
return (3.*1.*std::sqrt(300109429661055./Pi)*(7. - 1484.*std::pow(z/r,2.) + 40810.*std::pow(z/r,4.) - 310156.*std::pow(z/r,6.) + 653543.*std::pow(z/r,8.))*std::pow((x+I*y)/r,22.))/5.36870912e8;
}
else if (m == 23){
return (-3.*1.*std::sqrt(15905799772035915./(2.*Pi))*z/r*(-7. + 385.*std::pow(z/r,2.) - 4389.*std::pow(z/r,4.) + 12331.*std::pow(z/r,6.))*std::pow((x+I*y)/r,23.))/1.34217728e8;
}
else if (m == 24){
return (7.*1.*std::sqrt(757419036763615./Pi)*(-1. + 165.*std::pow(z/r,2.) - 3135.*std::pow(z/r,4.) + 12331.*std::pow(z/r,6.))*std::pow((x+I*y)/r,24.))/2.68435456e8;
}
else if (m == 25){
return (-7.*1.*std::sqrt(4998965642639859./(2.*Pi))*z/r*(5. - 190.*std::pow(z/r,2.) + 1121.*std::pow(z/r,4.))*std::pow((x+I*y)/r,25.))/1.34217728e8;
}
else if (m == 26){
return (1.*std::sqrt(174963797492395065./Pi)*(1. - 114.*std::pow(z/r,2.) + 1121.*std::pow(z/r,4.))*std::pow((x+I*y)/r,26.))/5.36870912e8;
}
else if (m == 27){
return -3.725290298461914e-9*(1.*std::sqrt(1108104050785168745./Pi)*z/r*(-3. + 59.*std::pow(z/r,2.))*std::pow((x+I*y)/r,27.));
}
else if (m == 28){
return (1.*std::sqrt(114631453529500215./(2.*Pi))*(-1. + 59.*std::pow(z/r,2.))*std::pow((x+I*y)/r,28.))/2.68435456e8;
}
else if (m == 29){
return -3.725290298461914e-9*(1.*std::sqrt(6763255758240512685./Pi)*z/r*std::pow((x+I*y)/r,29.));
}
else if (m == 30){
return (1.*std::sqrt(450883717216034179./Pi)*std::pow((x+I*y)/r,30.))/5.36870912e8;
}
else{return 0.;}
}

else{return 0.;}
}