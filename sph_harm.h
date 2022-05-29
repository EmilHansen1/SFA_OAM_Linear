#include <cmath>
#include <complex>
//#include <iostream>

using dcmplx = std::complex<double>;
constexpr dcmplx I(0., 1.);
constexpr double Pi = 3.14159265358979323846;

dcmplx sph_harm(dcmplx x, dcmplx y, dcmplx z, dcmplx r, int l, int m){
if(std::abs(m) > l){ return 0.; }
if(l < 0){ return 0.0; }


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

else{return 0.;}
}


/*
int main(){

double r; 
double theta; 
double phi; 
double x, y, z;
int l, m;

r = 1.1;  
theta = 0.8; 
phi = 0.69; 

l = 10; 
m = -9; 

x = r * std::sin(theta) * std::cos(phi);
y = r * std::sin(theta) * std::sin(phi);
z = r * std::cos(theta);
 

std::cout << sph_harm(x, y, z, r, l, m) << "\n"; 
return 0; 
}
*/

