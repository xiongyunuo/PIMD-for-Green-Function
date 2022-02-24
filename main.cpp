#include "XNumeric.hpp"
#include "XRandom.hpp"
#include "XDiffEq.hpp"
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

int N = 8;
int P = 64;
int d = 2;
XNum m = 1;
XNum kB = 1;
XNum T;
XNum beta;
XNum hBar = 1.0;
XNum omgP;
XNum L = 3.0;
XNum omg0 = 0.0 / hBar;
XNum s = 0.5;
XNum g = 3.0;
int cor_num = 60;
XNum incre = 0.5*std::sqrt(2)*L / cor_num;
int J0 = 2;
int J = 20;

XNum minimum_image(XNum d) {
  if (std::abs(d) > L / 2)
    return L - std::abs(d);
  return d;
}

XNum minimum_image2(XNum d) {
  if (std::abs(d) > L / 2) {
    if (d < 0)
      return L - std::abs(d);
    else
      return -(L - std::abs(d));
  }
  return d;
}

XNum distance(int i, int j, const std::vector<XNum> &Rs) {
  XNum a = minimum_image(Rs[d*i]-Rs[d*j]);
  XNum b = minimum_image(Rs[d*i+1]-Rs[d*j+1]);
  return a*a+b*b;
}

XNum ENk(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l) {
    if (l == N) {
      for (j = 1; j <= P-J; ++j) {
        int index = (l - 1) * P + j - 1;
        int index2 = (l - 1) * P + j;
        if (j == P-J) {
          if (l == N2)
            index2 = (N2 - k) * P;
          else
            index2 = l * P;
        }
        if (j == J0)
          index2 = -1;
        if (index2 != -1)
          res += distance(index2, index, Rs);
      }
      continue;
    }
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      res += distance(index2, index, Rs);
    }
  }
  return 0.5*m*omgP*omgP*res;
}

std::vector<XVecD> ENkCache;
std::vector<XNum> VBCache;

std::vector<XNum> expVB(int N, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum sum = 0;
    XNum tmp = 0.5*(ENkCache[N2-1][N2-1]+res[N2-1]);
    for (k = 1; k <= N2; ++k)
      sum += std::exp(-beta*(ENkCache[N2-1][k-1]+res[N2-k]-tmp));
    res.push_back(tmp-(1/beta)*std::log(sum / N2));
  }
  return res;
}

std::vector<XNum> gradient(int N2, int k, int l, int j, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  res.push_back(0);
  if (l >= N2 - k + 1 && l <= N2) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      if (l == N && j == J0)
        index2 = -1;
      if (l == N && j == P-J)
        index2 = (N2 - k) * P;
      int index3 = (l - 1) * P + j - 2;
      if (j == 1) {
        if (l == N2 - k + 1) {
          if (N2 == N)
            index3 = (N2 - 1) * P + P-J - 1;
          else
            index3 = (N2 - 1) * P + P - 1;
        }
        else
          index3 = (l - 2) * P + P - 1;
      }
      if (l == N && j == J0+1)
        index3 = -1;
      XNum a, b;
      if (index2 == -1)
        a = 0;
      else
        a = Rs[d*index] - Rs[d*index2];
      if (index3 == -1)
        b = 0;
      else
        b = Rs[d*index] - Rs[d*index3];
      res[0] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
      if (index2 == -1)
        a = 0;
      else
        a = Rs[d*index+1] - Rs[d*index2+1];
      if (index3 == -1)
        b = 0;
      else
        b = Rs[d*index+1] - Rs[d*index3+1];
      res[1] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
  }
  return res;
}

std::vector<XNum> XForceVBCache;

void forceVB(int N, const std::vector<XNum> &Rs) {
  int N2, k, l, j;
  XForceVBCache.clear();
  std::vector<XNum> tmp;
  for (l = 1; l <= N; ++l) {
    int num = P;
    if (l == N)
      num = P-J;
    for (j = 1; j <= num; ++j) {
      XVecD res;
      res.push_back(0);
      res.push_back(0);
      for (N2 = 1; N2 <= N; ++N2) {
        int i;
        XNum sum2 = 0;
        XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
        for (k = 1; k <= N2; ++k)
          sum2 += exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
        for (i = 0; i < d; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            tmp = gradient(N2, k, l, j, Rs);
            sum += (tmp[i] + res[d*(N2-k)+i])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
          }
          res.push_back(sum / sum2);
        }
      }
      XForceVBCache.push_back(res[d*N]);
      XForceVBCache.push_back(res[d*N+1]);
    }
  }
}

XVecD force(XNum t, const XVecD &Rs) {
  XVecD res;
  std::vector<XNum> tmp;
  int l, j;
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rs);
  VBCache = expVB(N, Rs);
  forceVB(N, Rs);
  for (l = 1; l <= N; ++l) {
    int num = P;
    if (l == N)
      num = P-J;
    for (j = 1; j <= num; ++j) {
      int index = (l - 1) * P + j - 1;
      XNum a = Rs[d*index] - L/2;
      res.push_back(-XForceVBCache[d*index]/m-minimum_image2(a)*omg0*omg0/P);
      a = Rs[d*index+1] - L/2;
      res.push_back(-XForceVBCache[d*index+1]/m-minimum_image2(a)*omg0*omg0/P);
    }
  }
  int k;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      if (l == N && j > P-J) continue;
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        if (k == N && j > P-J) continue;
        int index2 = (k - 1) * P + j - 1;
        if (k == N && j > J0)
          index = (l - 1) * P + j - 1 + J;
        if (l == N && j > J0)
          index2 = (k - 1) * P + j - 1 + J;
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += (2*minimum_image2(a)/(s*s))*inter/P;
      }
    }
  return res;
}

XNum GauEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int j, k, l;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      if (l == N && j > P-J) continue;
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        if (k == N && j > P-J) continue;
        int index2 = (k - 1) * P + j - 1;
        if (k == N && j > J0)
          index = (l - 1) * P + j - 1 + J;
        if (l == N && j > J0)
          index2 = (k - 1) * P + j - 1 + J;
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        res += 0.5*inter;
      }
    }
  return res/P;
}

XNum energy(const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
    XNum sum2 = 0;
    for (k = 1; k <= N2; ++k)
      sum2 += exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
    XNum sum = 0;
    for (k = 1; k <= N2; ++k)
      sum += (res[N2-k]-ENkCache[N2-1][k-1])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
    res.push_back(sum / sum2);
  }
  return res[N];
}

void period_boundary(std::vector<XNum> &Rs) {
  int i, j, k;
  for (i = 0; i < N; ++i) {
    int num = P;
    if (i == N-1)
      num = P-J;
    for (j = 0; j < num; ++j) {
      int index = i * P + j;
      for (k = 0; k < d; ++k) {
        if (Rs[d*index+k] < 0) {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] += (n + 1) * L;
        }
        else {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] -= n * L;
        }
      }
    }
  }
}

void initial(std::vector<XNum> &Rs, std::vector<XNum> &Vs) {
  int i, j;
  XNum v1 = 0, v2 = 0;
  for (i = 0; i < N; ++i) {
    int num = P;
    if (i == N-1)
      num = P-J;
    for (j = 0; j < num; ++j) {
      XNum tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v1 += tmp;
      Vs.push_back(tmp);
      tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v2 += tmp;
      Vs.push_back(tmp);
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
    }
  }
  v1 /= N * P - J;
  v2 /= N * P - J;
  for (i = 0; i < N; ++i) {
    int num = P;
    if (i == N-1)
      num = P-J;
    for (j = 0; j < num; ++j) {
      int index = i * P + j;
      Vs[d*index] -= v1;
      Vs[d*index+1] -= v2;
    }
  }
}

void logR(const std::vector<XNum> &Rs, std::ostream &out) {
  out << "{";
  for (int i = 0; i < N * P - J; ++i) {
    out << "{" << Rs[d*i] << "," << Rs[d*i+1] << "}";
    if (i != N * P - J - 1)
      out << ",";
  }
  out << "}";
}

void velocity_rescale(std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P - J; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1];
  XNum lam = std::sqrt((N*P-J)*d*(1/beta)/(m*sum));
  for (i = 0; i < N * P - J; ++i) {
    Vs[d*i] *= lam;
    Vs[d*i+1] *= lam;
  }
}

XNum temperature(const std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P - J; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1];
  return m*sum/(d*(N*P-J));
}

XNum trapEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int i;
  for (i = 0; i < N * P - J; ++i) {
    XNum a = Rs[d*i] - L/2;
    res += 0.5 * m * omg0 * omg0 * minimum_image(a) * minimum_image(a);
    a = Rs[d*i+1] - L/2;
    res += 0.5 * m * omg0 * omg0 * minimum_image(a) * minimum_image(a);
  }
  return res/P;
}

void XMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < N2*P-J; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < N2*P-J; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void XMMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int d2 = d;
  int d = 1;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < (N2*P-J)*d2; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < (N2*P-J)*d2; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void pair_correlation(const std::vector<XNum> &Rs, std::vector<XNum> &corr) {
  int n = corr.size();
  int in = (N-1)*P+J0-1;
  int in2 = (N-1)*P+J0;
  XNum a = Rs[d*in]-Rs[d*in2];
  XNum b = Rs[d*in+1]-Rs[d*in2+1];
  XNum dis = std::sqrt(minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b));
  int index = int(dis / incre);
  if (index < 0)
    index = 0;
  if (index >= n)
    index = n - 1;
  corr[index] += 1;
}

void pair_correlation2(const std::vector<XNum> &Rs, std::vector<XNum> &corr) {
  int n = corr.size();
  int in = (N-1)*P+J0;
  XNum a = Rs[d*in]-L/2;
  XNum b = Rs[d*in+1]-L/2;
  XNum dis = std::sqrt(minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b));
  int index = int(dis / incre);
  if (index < 0)
    index = 0;
  if (index >= n)
    index = n - 1;
  corr[index] += 1;
}

std::ofstream out("data.txt");

XNum RealT;
bool printInfo = true;

XNum simulation() {
  XForceCache.clear();
  ENkCache.clear();
  VBCache.clear();
  std::vector<XNum> Rs, Vs;
  std::vector<XNum> corr;
  int i, j;
  initial(Rs, Vs);
  velocity_rescale(Vs);
  period_boundary(Rs);
  XNum t = 0, h = 0.005;
  int skip = 50000;
  for (i = 0; i < cor_num; ++i)
    corr.push_back(0);
  std::vector<XNum> corr2 = corr;
  std::vector<XVecD> theta, vtheta, Q;
  for (i = 0; i < d*(N*P-J); ++i) {
    XVecD tmp;
    for (j = 0; j < 4; ++j) {
      tmp.push_back(0);
      tmp.push_back(0);
    }
    theta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1);
      tmp.push_back(1);
    }
    vtheta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1.0*(1+0.5*j));
      tmp.push_back(1.0*(1+0.5*j));
    }
    Q.push_back(tmp);
  }
  //std::cout << "{";
  for (i = 0; i < skip; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    t += h;
    /*if (i % 50 == 0) {
      logR(Rs);
      std::cout << ",";
      std::cout << "{" << i << "," << temperature(Vs) << "},";
    }*/
  }
  //std::cout << "}" << std::endl;
  if (printInfo) {
    logR(Rs, std::cout);
    logR(Rs, out);
    std::cout << std::endl;
    out << std::endl;
  }
  int count = 0;
  int steps = 10000000;
  XNum e = 0;
  XNum trape = 0;
  XNum ge = 0;
  XNum temp = 0;
  for (i = 0; i < steps; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    if (i % 1 == 0) {
      e += energy(Rs);
      temp += temperature(Vs);
      trape += trapEnergy(Rs);
      ge += GauEnergy(Rs);
      pair_correlation(Rs, corr);
      pair_correlation2(Rs, corr2);
      ++count;
    }
    if (i % 10000 == 0 && printInfo) {
      std::cout << "i=" << i << std::endl;
      out << "i=" << i << std::endl;
    }
    t += h;
  }
  e /= count;
  e += P * d * N / (2 * beta);
  //std::cout << "Ek=" << e << std::endl;
  //out << "Ek=" << e << std::endl;
  trape /= count;
  e += trape;
  ge /= count;
  //std::cout << "Eg=" << ge << std::endl;
  //out << "Eg=" << ge << std::endl;
  e += ge;
  temp /= count;
  //std::cout << temp << std::endl;
  //out << temp << std::endl;
  RealT = temp;
  XNum norm = 0;
  for (i = 0; i < cor_num; ++i)
    norm += corr[i]*incre;
  for (i = 0; i < cor_num; ++i)
    corr[i] /= norm;
  for (i = 0; i < cor_num; ++i)
    corr[i] /= 2*M_PI*(i+0.5)*incre;
  XNum den = N/(L*L);
  XNum mult = den/corr[0];
  for (i = 0; i < cor_num; ++i)
    corr[i] *= mult;
  norm = 0;
  for (i = 0; i < cor_num; ++i)
    norm += corr2[i]*incre;
  for (i = 0; i < cor_num; ++i)
    corr2[i] /= norm;
  for (i = 0; i < cor_num; ++i)
    corr2[i] /= 2*M_PI*(i+0.5)*incre;
  mult = den/corr2[0];
  for (i = 0; i < cor_num; ++i)
    corr2[i] *= mult;
  for (i = 0; i < cor_num; ++i)
    std::cout << "{" << i*incre << "," << corr[i] << "},";
  std::cout << std::endl;
  for (i = 0; i < cor_num; ++i)
    out << "{" << i*incre << "," << corr[i] << "},";
  out << std::endl;
  /*for (i = 0; i < cor_num; ++i)
    std::cout << "{" << i*incre << "," << corr2[i] << "},";
  std::cout << std::endl;
  for (i = 0; i < cor_num; ++i)
    out << "{" << i*incre << "," << corr2[i] << "},";
  out << std::endl;*/
  norm = 0;
  for (i = 0; i < cor_num; ++i)
    norm += corr[i]*incre*(i+0.5)*incre;
  //std::cout << norm << std::endl;
  return e;
}

int main() {
  clock_t t;
  t = clock();
  XSetRandSeed(5);
  for (T = 1.625/1.0; T <= 1.625/1.0; T += 1) {
    P = (int)(72/(T*6.0));
    J = (int)(P/3.0);
    beta = 1 / (kB * T);
    omgP = std::sqrt(P) / (beta * hBar);
    XNum e = simulation();
    std::cout << "{" << T << ", " << RealT << "}, ";
    out << "{" << T << ", " << RealT << "}, ";
    std::cout << std::endl;
    out << std::endl;
  }
  std::cout << std::endl;
  out << std::endl;
  t = clock() - t;
  std::cout << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  out << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  return 0;
}
