#include "XDiffEq.hpp"
#include <cmath>

XVecD XRungeKuttaStep(XNum t, XNum h, const XVecD &x, XForceFunc *func) {
  XVecD res = x;
  XVecD k1 = func(t, x);
  XVecD::size_type i, n = x.size();
  for (i = 0; i < n; ++i) {
    k1[i] *= h;
    res[i] += k1[i] / 6;
    k1[i] = x[i] + k1[i] / 2;
  }
  XVecD k2 = func(t + h / 2, k1);
  for (i = 0; i < n; ++i) {
    k2[i] *= h;
    res[i] += k2[i] / 3;
    k2[i] = x[i] + k2[i] / 2;
  }
  XVecD k3 = func(t + h / 2, k2);
  for (i = 0; i < n; ++i) {
    k3[i] *= h;
    res[i] += k3[i] / 3;
    k3[i] = x[i] + k3[i];
  }
  XVecD k4 = func(t + h, k3);
  for (i = 0; i < n; ++i)
    res[i] += k4[i] * h / 6;
  return res;
}

XNum XRungeKuttaInterval(XNum t, XNum h, XNum delta, const XVecD &x, XForceFunc *func) {
  XVecD x1 = XRungeKuttaStep(t, h, x, func);
  XVecD x2 = XRungeKuttaStep(t, h / 2, x, func);
  x2 = XRungeKuttaStep(t + h / 2, h / 2, x2, func);
  XNum d = 0;
  XVecD::size_type i, n = x.size();
  for (i = 0; i < n; ++i)
    d += (x1[i] - x2[i]) * (x1[i] - x2[i]);
  d = std::sqrt(d);
  return (15.0 / 16.0) * h * std::pow(delta / d, 0.2);
}

XVecD XNumerovStep(int n, XNum t, XNum h, XNum x0, XNum x1, XOneArguFunc *func) {
  XNum w0 = (1 - h * h * func(t) / 12) * x0;
  XNum w1 = (1 - h * h * func(t + h) / 12) * x1;
  int i;
  XVecD res;
  res.push_back(x0);
  res.push_back(x1);
  for (i = 0; i < n - 1; ++i) {
    XNum w2 = h * h * func(t + h) * x1 + 2 * w1 - w0;
    XNum x2 = w2 / (1 - h * h * func(t + 2 * h) / 12);
    res.push_back(x2);
    w0 = w1;
    w1 = w2;
    x0 = x1;
    x1 = x2;
    t += h;
  }
  return res;
}

XNum XNumerovInit(XNum t, XNum h, XNum x0, XNum v0, XOneArguFunc *func) {
  XNum top = (2 + 5 * h * h * func(t) / 6) * (1 - h * h * func(t - h) / 12) * x0;
  top += 2 * h * v0 * (1 - h * h * func(t - h) / 6);
  XNum bottom = (1 - h * h * func(t + h) / 12) * (1 - h * h * func(t - h) / 6);
  bottom += (1 - h * h * func(t - h) / 12) * (1 - h * h * func(t + h) / 6);
  return top / bottom;
}

XVecD XVerletStep(int n, XNum t, XNum h, XNum x0, XNum x1, XOneArguFunc *func) {
  int i;
  XVecD res;
  res.push_back(x0);
  res.push_back(x1);
  for (i = 0; i < n - 1; ++i) {
    XNum x2 = 2 * x1 - x0 + h * h * func(t + h);
    res.push_back(x2);
    x0 = x1;
    x1 = x2;
    t += h;
  }
  return res;
}

std::pair<XVecD, XVecD> XVerletMVStep(XNum t, XNum h, const XVecD &x, const XVecD &v, XForceFunc *func) {
  XVecD f = func(t, x);
  XVecD resx;
  int i;
  for (i = 0; i < x.size(); ++i)
    resx.push_back(x[i] + h * v[i] + h * h * f[i] / 2);
  XVecD f2 = func(t + h, resx);
  XVecD resv;
  for (i = 0; i < v.size(); ++i)
    resv.push_back(v[i] + h * (f2[i] + f[i]) / 2);
  return std::pair<XVecD, XVecD>(resx, resv);
}

XVecD NHForce(XNum beta, XNum m, int f, XVecD &Q, XVecD &v, XVecD &vtheta) {
  XVecD res;
  XNum sum = 0;
  int i;
  for (i = 0; i < f; ++i)
    sum += m*v[i]*v[i];
  res.push_back((sum-f*(1/beta))/Q[0]);
  int M = vtheta.size();
  for (i = 1; i < M; ++i)
    res.push_back((Q[i-1]*vtheta[i-1]*vtheta[i-1]-(1/beta))/Q[i]);
  return res;
}

static XNum GjPlus(XNum h, int j, XVecD &vtheta) {
  if (j < 0 || j >= vtheta.size()) return 1;
  return 1 + h*vtheta[j]/2;
}

static XNum GjMinus(XNum h, int j, XVecD &vtheta) {
  if (j < 0 || j >= vtheta.size()) return 1;
  return 1 - h*vtheta[j]/2;
}

static XNum gjPlus(XNum h, int j, XVecD &vtheta) {
  if (j < 0 || j >= vtheta.size()) return 1;
  return 1 + h*vtheta[j]/4;
}

static XNum gjMinus(XNum h, int j, XVecD &vtheta) {
  if (j < 0 || j >= vtheta.size()) return 1;
  return 1 - h*vtheta[j]/4;
}

XVecD XForceCache;

void XVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, XVecD &theta, XVecD &vtheta, XNum beta, XNum m, XVecD &Q, XForceFunc *func) {
  int i;
  int N = v.size();
  int M = vtheta.size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  XVecD f2 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 0; i < N; ++i)
    v[i] = (v[i]*gjMinus(h,0,vtheta)+h*XForceCache[i]/2)/gjPlus(h,0,vtheta);
  int M2 = M / 2;
  for (i = 1; i <= M2; ++i)
    theta[2*i-2] = theta[2*i-2]+h*vtheta[2*i-2]/2;
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-1] = (vtheta[2*i-1]*gjMinus(h,2*i,vtheta)+h*f2[2*i-1]/2)/gjPlus(h,2*i,vtheta);
  for (i = 0; i < N; ++i)
    x[i] = x[i] + h*v[i];
  for (i = 1; i <= M2; ++i)
    theta[2*i-1] = theta[2*i-1]+h*vtheta[2*i-1];
  XVecD f3 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-2] = (vtheta[2*i-2]*GjMinus(h,2*i-1,vtheta)+h*f3[2*i-2])/GjPlus(h,2*i-1,vtheta);
  XForceCache = func(t+h, x);
  for (i = 0; i < N; ++i)
    v[i] = (v[i]*gjMinus(h,0,vtheta)+h*XForceCache[i]/2)/gjPlus(h,0,vtheta);
  for (i = 1; i <= M2; ++i)
    theta[2*i-2] = theta[2*i-2]+h*vtheta[2*i-2]/2;
  f3 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-1] = (vtheta[2*i-1]*gjMinus(h,2*i,vtheta)+h*f3[2*i-1]/2)/gjPlus(h,2*i,vtheta);
}

void XVerletNHChain3(XNum t, XNum h, XVecD &x, XVecD &v, XVecD &theta, XVecD &vtheta, XNum beta, XNum m, XVecD &Q, XForceFunc *func) {
  int i;
  int N = v.size();
  int M = vtheta.size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  XVecD f2 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 0; i < N; ++i)
    v[i] = v[i]*std::exp(-0.5*h*vtheta[0])+0.5*h*XForceCache[i]*std::exp(-0.25*h*vtheta[0]);
  int M2 = M / 2;
  for (i = 1; i <= M2; ++i)
    theta[2*i-2] = theta[2*i-2]+h*vtheta[2*i-2]/2;
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-1] = vtheta[2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[2*i]));
  for (i = 0; i < N; ++i)
    x[i] = x[i] + h*v[i];
  for (i = 1; i <= M2; ++i)
    theta[2*i-1] = theta[2*i-1]+h*vtheta[2*i-1];
  XVecD f3 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-2] = vtheta[2*i-2]*std::exp(-h*vtheta[2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[2*i-1]);
  XForceCache = func(t+h, x);
  for (i = 0; i < N; ++i)
    v[i] = v[i]*std::exp(-0.5*h*vtheta[0])+0.5*h*XForceCache[i]*std::exp(-0.25*h*vtheta[0]);
  for (i = 1; i <= M2; ++i)
    theta[2*i-2] = theta[2*i-2]+h*vtheta[2*i-2]/2;
  f3 = NHForce(beta, m, N, Q, v, vtheta);
  for (i = 1; i <= M2; ++i)
    vtheta[2*i-1] = vtheta[2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[2*i]));
}