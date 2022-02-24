#ifndef XDIFFEQ_HPP
#define XDIFFEQ_HPP

#include "XNumeric.hpp"
#include <utility>

typedef XVecD XForceFunc(XNum, const XVecD&);

extern XVecD XForceCache;

XVecD XRungeKuttaStep(XNum t, XNum h, const XVecD &x, XForceFunc *func);
XNum XRungeKuttaInterval(XNum t, XNum h, XNum delta, const XVecD &x, XForceFunc *func);
XVecD XNumerovStep(int n, XNum t, XNum h, XNum x0, XNum x1, XOneArguFunc *func);
XNum XNumerovInit(XNum t, XNum h, XNum x0, XNum v0, XOneArguFunc *func);
XVecD XVerletStep(int n, XNum t, XNum h, XNum x0, XNum x1, XOneArguFunc *func);
std::pair<XVecD, XVecD> XVerletMVStep(XNum t, XNum h, const XVecD &x, const XVecD &v, XForceFunc *func);
void XVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, XVecD &theta, XVecD &vtheta, XNum beta, XNum m, XVecD &Q, XForceFunc *func);
void XVerletNHChain3(XNum t, XNum h, XVecD &x, XVecD &v, XVecD &theta, XVecD &vtheta, XNum beta, XNum m, XVecD &Q, XForceFunc *func);
XVecD NHForce(XNum beta, XNum m, int f, XVecD &Q, XVecD &v, XVecD &vtheta);

#endif