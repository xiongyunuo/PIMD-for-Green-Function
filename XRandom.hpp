#ifndef XRANDOM_HPP
#define XRANDOM_HPP

typedef unsigned int XRandUInt;
typedef float XRandF;

XRandUInt XRandInteger();
void XShiftRegisterInit();
void XShiftRegisterDestroy();
XRandUInt XShiftRegisterRand();
void XSetRandSeed(XRandUInt seed);
XRandUInt XGetRandSeed();
XRandF XRandFloat();
XRandF XRandGauss();

#endif