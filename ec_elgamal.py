# Elliptical ElGamal encryption
# This is a program for analyzing the resistance of elliptical ElGamal encryption to timing attacks.
# Author: HIRATA Tomonori, Graduate School of Informatics Nagoya University, 2020

import math
import sys
import random
import sympy
import inspect
from tqdm import tqdm
import numpy as np
from scipy.stats import norm
import statistics as stat
import scipy.stats as stats
import matplotlib.pyplot as plt
import yaml
from argparse import ArgumentParser

# =============================
# option analysis
# =============================
def get_option(ecparams, loop, tlength, weight):
  argparser = ArgumentParser()
  argparser.add_argument('-ecp', '--ecparams', type=str,
                          default=ecparams,
                          help='Specify type of ec-parameter')
  argparser.add_argument('-l', '--loop', type=int,
                          default=loop,
                          help='Specify number of loops')
  argparser.add_argument('-t', '--tlength', type=int,
                          default=tlength,
                          help='Specify number of t-length')
  argparser.add_argument('-w', '--weight', type=int,
                          default=weight,
                          help='Specify number of hamming weight')
  return argparser.parse_args()


# =============================
# yaml import
# =============================
args = get_option('config1', 1000, 256, 128)
with open('./parameters/'+args.ecparams+'.yml', 'r') as yml:
  config = yaml.load(yml, Loader=yaml.FullLoader)

# =============================
# global constants
# =============================
INF = 999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
MODULO_P = int(config['ec_params']['modulo_p'])
A = int(config['ec_params']['a'])
B = int(config['ec_params']['b'])
BIT = len(bin(MODULO_P)[2:])
N = 100 # N is used when convert Message to Point, and Point to Message
P_DOUBLE = 5/6
P_ADD = 2/3
LOOP = args.loop
ALPHA = (12*math.log(2)*math.log(2**BIT))/(math.pi**2) + 1.467
SIGMA_BASE = config['ec_params']['sigma_base'] # depend on MODULO_P

# =============================
# variables for analysis
# =============================
egcd_count_list = []
mult_list = []
reduction_list = []
total_list = []

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
  global opr, egcd_count_list
  (x, lastx) = (0, 1)
  (y, lasty) = (1, 0)

  egcd_count = 0
  if calledFromElGamalDec():
    while b != 0:
      egcd_count += 1
      q = opr.div(a, b)
      (a, b) = (b, a - opr.multiply(q, b))
      (x, lastx) = (lastx - opr.multiply(q, x), x)
      (y, lasty) = (lasty - opr.multiply(q, y), y)
    egcd_count_list.append(egcd_count)
  else:
    while b != 0:
      q = opr.div(a, b)
      (a, b) = (b, a - opr.multiply(q, b))
      (x, lastx) = (lastx - opr.multiply(q, x), x)
      (y, lasty) = (lasty - opr.multiply(q, y), y)

  # while b != 0:
  #   q = opr.div(a, b)
  #   (a, b) = (b, a - opr.multiply(q, b))
  #   (x, lastx) = (lastx - opr.multiply(q, x), x)
  #   (y, lasty) = (lasty - opr.multiply(q, y), y)
  return (lastx, lasty, a)

# ax ≡ 1 (mod m)
def modinv(a: int, m: int):
  global opr
  (inv, _, _) = egcd(a, m)
  return opr.modulo(inv, m) # not perform reductions.

def generate_d(t: int, w: int):
  # generate a list with the correct number of 1's
  x = [1]*(w-1)+[0]*(t-w)
  random.shuffle(x)
  x = [1]+x
  # convert back to a number
  return int(''.join(map(str, x)), 2)

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
  GP = MsgToPoint(generate_d(args.tlength-1, args.weight))

  # generate secret key and public key
  # d = 0x49be667ef9dcbbac55b06295ce870b0702900cdb2dce28d959f9425b16f81798
  d = generate_d(args.tlength, args.weight)
  d_bit_seq = bin(d)[2:]
  t = len(d_bit_seq)
  w = d_bit_seq.count("1")

  if BIT < args.tlength:
    print('tlength must be least than BIT(length of MODULO_P)')
    sys.exit()

  # analysis of distributions
  iteration = (t+w-1)
  mu_mult = (t-1)*(3*ALPHA+6)+w*(3*ALPHA+3)
  mu_reduction = ALPHA*(t+w-1)+3*(t-1)*P_DOUBLE+3*w*P_ADD
  sigma_mult = 3*SIGMA_BASE*math.sqrt(iteration)
  sigma_reduction = SIGMA_BASE*math.sqrt(iteration)
  mu = mu_mult + mu_reduction
  sigma = sigma_mult + sigma_reduction

  PUBLIC_KEY = Binary(d, GP)

  correct = 0
  for i in tqdm(range(LOOP)):
    plaintext = random.randint(1, MODULO_P//100)
    M = MsgToPoint(plaintext)

    # core logic
    C1, C2 = ElGamalEnc(M, GP, PUBLIC_KEY)
    opr.mult = 0
    opr.reduction = 0
    DecM = ElGamalDec(C1, C2, d)
    mult_list.append(opr.mult)
    reduction_list.append(opr.reduction)
    total_list.append(opr.mult+opr.reduction)

    decPlaintext = PointToMsg(DecM)
    if plaintext == decPlaintext:
      correct += 1

  print("========== Result ==========")
  print("correctness: {0}/{1}".format(correct, LOOP))
  print("=================== global constants ==================")
  print("GP = {0}".format(vars(GP)))
  print("PUBLIC_KEY = {0}".format(vars(PUBLIC_KEY)))
  print('d = {0} ({1}bits)'.format(d, t))
  print("Hamming weight of d: {0}".format(w))
  print("MODULO_P: {0} ({1} bits)".format(MODULO_P, BIT))
  print("A: {0}".format(A))
  print("B: {0}".format(B))
  print("=================== egcd ==================")
  print("Average of egcd count (theoretical, namely ALPHA): {0}".format(ALPHA))
  print("Average of egcd count (experimental) = {0}".format(stat.mean(egcd_count_list)))
  print("Standard deviation of egcd count (theoretical, namely SIGMA_BASE): {0}".format(SIGMA_BASE))
  print("Standard deviation of egcd count (experimental) = {0}".format(stat.pstdev(egcd_count_list)))
  print("Variance value of egcd count (theoretical): {0}".format(SIGMA_BASE**2))
  print("Variance value of egcd count (experimental) = {0}".format(stat.pvariance(egcd_count_list)))
  print("================== multiplication ===================")
  print("Average multiplication (theoretical) = {0}".format(mu_mult))
  print("Average multiplication (experimental) = {0}".format(stat.mean(mult_list)))
  print("Standard deviation of multiplication (theoretical) = {0}".format(sigma_mult))
  print("Standard deviation of multiplication (experimental) = {0}".format(stat.pstdev(mult_list)))
  print("Variance value of multiplication (theoretical) = {0}".format(sigma_mult**2))
  print("Variance value of multiplication (experimental) = {0}".format(stat.pvariance(mult_list)))
  print("================== reduction ===================")
  print("Average reduction (theoretical) = {0}".format(mu_reduction))
  print("Average reduction (experimental) = {0}".format(stat.mean(reduction_list)))
  print("Standard deviation of reduction (theoretical) = {0}".format(sigma_reduction))
  print("Standard deviation of reduction (experimental) = {0}".format(stat.pstdev(reduction_list)))
  print("Variance value of reduction (theoretical) = {0}".format(sigma_reduction**2))
  print("Variance value of reduction (experimental) = {0}".format(stat.pvariance(reduction_list)))
  print("================= total count ====================")
  print("Average computation (theoretical) = {0}".format(mu))
  print("Average computation (experimental) = {0}".format(stat.mean(total_list)))
  print("Standard deviation of computation (theoretical) = {0}".format(sigma))
  print("Standard deviation of computation (experimental) = {0}".format(stat.pstdev(total_list)))
  print("Variance value of computation (theoretical) = {0}".format(sigma**2))
  print("Variance value of computation (experimental) = {0}".format(stat.pvariance(total_list)))
  print("================== test whether normal distribution (p > 0.05) ===================")
  if LOOP <= 5000:
    print("Shapiro-Wilk test")
    print("multiply: {0}".format(stats.shapiro(mult_list)))
    print("reduction: {0}".format(stats.shapiro(reduction_list)))
    print("total count: {0}".format(stats.shapiro(total_list)))
  else:
    print("Kolmogorov–Smirnov test")
    print("multiply: {0}".format(stats.kstest(mult_list, "norm")))
    print("reduction: {0}".format(stats.kstest(reduction_list, "norm")))
    print("total count: {0}".format(stats.kstest(total_list, "norm")))


  # ====================
  # display histgram
  # ====================

  # multiplication (172000 ~ 176000)
  plt.hist(mult_list, bins=len(set(mult_list)), density=True)
  X = np.arange(mu_mult-sigma_mult*5, mu_mult+sigma_mult*5, 1)
  Y = norm.pdf(X, loc=mu_mult, scale=sigma_mult)
  plt.plot(X, Y, 'r-')
  plt.title("Multiplication ($\mu = {0} \ \ \sigma = {1}$)".format(round(mu_mult, 2), round(sigma_mult, 2)))
  plt.savefig("picture/ec_elgamal_mult.png")
  plt.show()

  # reduction (57600 ~ 59200)
  plt.hist(reduction_list, bins=len(set(reduction_list)), density=True)
  X = np.arange(mu_reduction-sigma_reduction*5, mu_reduction+sigma_reduction*5, 1)
  Y = norm.pdf(X, loc=mu_reduction, scale=sigma_reduction)
  plt.plot(X, Y, 'b-')
  plt.title("Reduction ($\mu = {0} \ \ \sigma = {1}$)".format(
      round(mu_reduction, 2), round(sigma_reduction, 2)))
  plt.savefig("picture/ec_elgamal_reduction.png")
  plt.show()

  # total computation
  # mu = stat.mean(total_list)
  # sigma = stat.pstdev(total_list)
  plt.hist(total_list, bins=len(set(total_list)), density=True)
  X = np.arange(mu-sigma*5, mu+sigma*5, 1)
  Y = norm.pdf(X, loc=mu, scale=sigma)
  plt.plot(X, Y, 'g-')
  plt.title("Total ($\mu = {0} \ \  \sigma = {1}$)".format(round(mu, 2), round(sigma, 2)))
  plt.savefig("picture/ec_elgamal_total.png")
  plt.show()

  # total computation
  mu = stat.mean(total_list)
  sigma = stat.pstdev(total_list)
  plt.hist(total_list, bins=len(set(total_list)), density=True)
  X = np.arange(mu-sigma*5, mu+sigma*5, 1)
  Y = norm.pdf(X, loc=mu, scale=sigma)
  plt.plot(X, Y, 'g-')
  plt.savefig("picture/ec_elgamal_total_sum.png")
  plt.show()

  # egcd_count
  plt.hist(egcd_count_list, bins=len(set(egcd_count_list)), density=True)
  X = np.arange(60, 200, 0.1)
  Y = norm.pdf(X, loc=ALPHA, scale=SIGMA_BASE)
  plt.plot(X, Y, 'r-')
  plt.title("Egcd Iteration ($\mu = {0} \ \ \sigma = {1}$)".format(round(ALPHA, 2), round(SIGMA_BASE, 2)))
  plt.savefig("picture/ec_elgamal_egcd_count.png")
  plt.show()
