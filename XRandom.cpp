#include "XRandom.hpp"
#include <utility>
#include <cmath>

static XRandUInt XSeed = 1;
static XRandUInt XMod = 2147483647;
static XRandUInt XA = 16807;
static XRandUInt XQ = 127766;
static XRandUInt XR = 120485;
static int XNP = 521;
static int XNQ = 32;

namespace {
  struct XRandStackItem {
    XRandUInt num;
    XRandStackItem *next;
  };
}

static XRandStackItem *XFrontItem = nullptr;
static XRandStackItem *XTrailItem = nullptr;
static XRandStackItem *XPItem = nullptr;
static XRandStackItem *XQItem = nullptr;

void XShiftRegisterInit() {
  int i;
  XFrontItem = new XRandStackItem();
  XFrontItem->num = XRandInteger();
  XFrontItem->next = nullptr;
  XTrailItem = XFrontItem;
  XQItem = XFrontItem;
  for (i = 1; i <= XNP - XNQ; ++i) {
    XRandStackItem *temp = new XRandStackItem();
    temp->num = XRandInteger();
    temp->next = nullptr;
    XTrailItem->next = temp;
    XTrailItem = temp;
  }
  XPItem = XTrailItem;
  for (i = XNP - XNQ + 1; i < XNP; ++i) {
    XRandStackItem *temp = new XRandStackItem();
    temp->num = XRandInteger();
    temp->next = nullptr;
    XTrailItem->next = temp;
    XTrailItem = temp;
  }
}

XRandUInt XShiftRegisterRand() {
  XRandUInt res = XQItem->num ^ XPItem->num;
  XRandStackItem *temp = new XRandStackItem();
  temp->num = res;
  temp->next = nullptr;
  XTrailItem->next = temp;
  XTrailItem = temp;
  XRandStackItem *old = XFrontItem;
  XFrontItem = XFrontItem->next;
  XQItem = XQItem->next;
  XPItem = XPItem->next;
  delete old;
  return res;
}

void XShiftRegisterDestroy() {
  XRandStackItem *cur = XFrontItem, *old;
  while (cur->next) {
    old = cur;
    cur = cur->next;
    delete old;
  }
  delete cur;
  XFrontItem = nullptr;
  XTrailItem = nullptr;
  XQItem = nullptr;
  XPItem = nullptr;
}

XRandUInt XRandInteger() {
  int res = (int)(XA * (XSeed % XQ)) - (int)(XR * (XSeed / XQ));
  if (res < 0) res += XMod;
  XSeed = res;
  return XSeed;
}

void XSetRandSeed(XRandUInt seed) {
  XSeed = seed;
}

XRandUInt XGetRandSeed() {
  return XSeed;
}

XRandF XRandFloat() {
  return (XRandF)(XRandInteger()) / XMod;
}

XRandF XRandGauss() {
  static XRandF first, second;
  static bool has = false;
  if (has) {
    has = false;
    return second;
  }
  has = true;
  XRandF r = std::sqrt(-2 * std::log(XRandFloat()));
  XRandF angle = 2 * M_PI * XRandFloat();
  first = r * cos(angle);
  second = r * sin(angle);
  return first;
}