# Elliptical ElGamal encryption
# This is a program for analyzing the resistance of elliptical ElGamal encryption to timing attacks.
# Author: HIRATA Tomonori, Graduate School of Informatics Nagoya University, 2020

import math
import sys
import random
import sympy
import inspect
import numpy as np
import matplotlib.pyplot as plt

# =============================
# global constants
# =============================
INF = -1
MODULO_P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f # 256bit
A = 0x0000000000000000000000000000000000000000000000000000000000000000
B = 0x0000000000000000000000000000000000000000000000000000000000000007
N = 100
LOOP = 1000

# =============================
# variables for analysis
# =============================
all_mult = 0
all_reduction = 0
mult_list = []
reduction_list = []

# =============================
# classes
# =============================

# Point is the class that represents coordinates
class Point:
  def __init__(self, x: int, y: int):
    self.x = x
    self.y = y

  # check if P is a infinity point
  def is_inf(self):
    return (self.x == self.y == INF)
  
  # inverse P
  def inv(self):
    self.y = (-self.y) + MODULO_P

# Operation is the class that performs multiplication, modulo, and division. 
class Operation:
  def __init__(self):
    self.mult = 0
    self.reduction = 0
  
  def multiply(self, x: int, y: int):
    self.mult += 1
    return x*y
  
  def modulo(self, x: int, y: int):
    if x >= y:
      self.reduction += 1
    return x % y
  
  def div(self, x: int, y: int):
    self.reduction += 1
    return x // y

# opr is the global variable. 
opr = Operation()

# =============================
# sub function
# =============================

# Extended Euclidean algorithm
### Measured value: 大体453.2496669540005
def egcd(a: int, b: int):
    global opr
    (x, lastx) = (0, 1)
    (y, lasty) = (1, 0)
    while b != 0:
        q = opr.div(a, b)
        (a, b) = (b, a - opr.multiply(q, b))
        (x, lastx) = (lastx - opr.multiply(q, x), x)
        (y, lasty) = (lasty - opr.multiply(q, y), y)
    return (lastx, lasty, a)

# ax ≡ 1 (mod m)
def modinv(a: int, m: int):
    global opr
    (inv, _, _) = egcd(a, m)
    return opr.modulo(inv, m) # not perform reductions.

# return true if called from ElGamalDec
def calledFromElGamalDec():
  stack = inspect.stack()
  for i in range(len(stack)-1, -1, -1):
    if stack[i][0].f_code.co_name == "ElGamalDec":
      return True
  return False

# =============================
# main function
# =============================

# Point addition
# double
### measured value: 459.129860112433
# add
### measured value: 456.066388308977
def PointAdd(P: Point, Q: Point):
  global opr
  if P.is_inf():
    return Q
  elif Q.is_inf():
    return P
  elif P.x == Q.x:
    if ((P.y + Q.y) == MODULO_P):
      return Point(INF, INF)
    else:
      tmp1 = opr.multiply(opr.multiply(3, P.x), P.x) + A
      tmp2 = opr.multiply(2, P.y)
      inv_tmp2 = modinv(tmp2, MODULO_P)
      lmd = opr.modulo(opr.multiply(tmp1, inv_tmp2), MODULO_P) # perform reductions almost all. 
      x_3 = opr.modulo(opr.multiply(lmd, lmd) - P.x - Q.x, MODULO_P)  # perform reductions almost all. 
      y_3 = opr.modulo(opr.multiply(lmd, (P.x - x_3)) - P.y, MODULO_P)  # perform reductions with probability 1/2. 
  else:
    tmp1 = (Q.y - P.y)
    tmp2 = (Q.x - P.x)
    inv_tmp2 = modinv(tmp2, MODULO_P)
    lmd = opr.modulo(opr.multiply(tmp1, inv_tmp2), MODULO_P) # perform reductions with probability 1/2. 
    x_3 = opr.modulo(opr.multiply(lmd, lmd) - P.x - Q.x, MODULO_P) # perform reductions almost all. 
    y_3 = opr.modulo(opr.multiply(lmd, (P.x - x_3)) - P.y, MODULO_P) # perform reductions with probability 1/2. 
  return Point(x_3, y_3)

# Binary method
# Q = dP
def Binary(d: int, P: Point):
  global opr
  Q = Point(INF, INF)
  d_bit_seq = bin(d)[2:]
  for i in d_bit_seq:
    Q = PointAdd(Q, Q)  # measured value (multiplication): 457.4243529411765
    if int(i) == 1:
      Q = PointAdd(P, Q)  # measured value (multiplication): 452.6232283464567
  return Q

# convert message to point
def MsgToPoint(m: int):
  for i in range(0, N):
    x = 100*m+i
    z = x**3
    z %= MODULO_P
    z += A*x + B
    if sympy.legendre_symbol(z, MODULO_P) == 1:
      y = sympy.sqrt_mod(z, MODULO_P)
      return Point(x, y)

# convert point to message
def PointToMsg(P: Point):
  return P.x // 100

# Encryption
def ElGamalEnc(M: Point, GP: Point, PUBLIC_KEY: Point):
  r = random.randint(1, MODULO_P)
  return (Binary(r, GP), PointAdd(M, Binary(r, PUBLIC_KEY)))

# Decryption
def ElGamalDec(C1: Point, C2: Point, d: int):
  C1.inv()
  B = Binary(d, C1)
  ans = PointAdd(C2, B)
  return ans


# =============================
# main function
# =============================
if __name__ == '__main__':
  GP = MsgToPoint(0xffffffff222ffff9dcbbac55eeefffffffff)

  # generate secret key and public key
  d = 0x49be667ef9dcbbac55b06295ce870b0702900cdb2dce28d959f2815b16f81798
  d_bit_seq = bin(d)[2:]
  PUBLIC_KEY = Binary(d, GP)

  correct = 0
  for i in range(LOOP):
    print("i = {0}".format(i))
    plaintext = random.randint(1, MODULO_P//100)
    M = MsgToPoint(plaintext)

    # core logic
    C1, C2 = ElGamalEnc(M, GP, PUBLIC_KEY)
    opr.mult = 0
    opr.reduction = 0
    DecM = ElGamalDec(C1, C2, d)
    all_mult += opr.mult
    all_reduction += opr.reduction
    mult_list.append(opr.mult)
    reduction_list.append(opr.reduction)

    decPlaintext = PointToMsg(DecM)
    if plaintext == decPlaintext:
      correct += 1

  print("========== Result ==========")
  print("correctness: {0}/{1}".format(correct, LOOP))
  print("GP = {0}".format(vars(GP)))
  print("PUBLIC_KEY = {0}".format(vars(PUBLIC_KEY)))
  print('d = {0} ({1}bits)'.format(d, len(d_bit_seq)))
  print("Hamming weight of d: {0}".format(d_bit_seq.count("1")))
  print("Average multiplication = {0}".format(all_mult/LOOP))
  print("Average reduction = {0}".format(all_reduction/LOOP))

  # display histgram
  plt.hist(mult_list, bins=LOOP)
  plt.show()
  plt.hist(reduction_list, bins=LOOP)
  plt.show()


