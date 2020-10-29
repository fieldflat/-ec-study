# 楕円ElGamal暗号
# alpha = (12*math.log(2)*math.log(2**256))/(math.pi**2)+1.467

import math
import sys
import random
import sympy
import inspect
import numpy as np
import matplotlib.pyplot as plt

# =============================
# グローバル定数
# =============================
INF = -1
MODULO_P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f # 256bit
A = 0x0000000000000000000000000000000000000000000000000000000000000000
B = 0x0000000000000000000000000000000000000000000000000000000000000007
N = 100
LOOP = 100
NR1 = 0
NR1_SUM = 0
NR2 = 0
NR2_SUM = 0
tmp_mult_sum_list = []

# =============================
# class
# =============================

# 座標を表すclass
class Point:
  def __init__(self, x: int, y: int):
    self.x = x
    self.y = y

  # 点Pが無限遠点であるかチェックする
  def is_inf(self):
    return (self.x == self.y == INF)
  
  def inv(self):
    self.y = (-self.y) + MODULO_P

# 分析クラス
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
    return x // y

opr = Operation()

# =============================
# サブ関数 (主要関数内で使用されるサブ関数)
# =============================

# 拡張ユークリッド互除法
# multiply: 3*alpha  (alpha ... 平均操作回数)
### 実測値: 大体453.2496669540005
# reduction: alpha
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
# multiply: 3*alpha  (alpha ... 平均操作回数)
# reduction: alpha
def modinv(a: int, m: int):
    global opr
    (inv, _, _) = egcd(a, m)
    return opr.modulo(inv, m) # ほぼmodulo発生しない！

# ElGamalDecからの関数呼び出しであればTrueを返す．
def calledFromElGamalDec():
  stack = inspect.stack()
  for i in range(len(stack)-1, -1, -1):
    if stack[i][0].f_code.co_name == "ElGamalDec":
      return True
  return False

# =============================
# 主要関数
# =============================

# 楕円曲線上の点演算
# 2倍算
### 実測値: 459.129860112433
# 加算
### 実測値: 456.066388308977
def PointAdd(P: Point, Q: Point):
  global opr
  if P.is_inf():
    # 一回の復号あたり1回実行される
    return Q
  elif Q.is_inf():
    # 一回の復号あたり1回実行される
    return P
  elif P.x == Q.x:
    if ((P.y + Q.y) == MODULO_P):
      return Point(INF, INF)
    else:
      tmp1 = opr.multiply(opr.multiply(3, P.x), P.x) + A
      tmp2 = opr.multiply(2, P.y)
      inv_tmp2 = modinv(tmp2, MODULO_P)
      lmd = opr.modulo(opr.multiply(tmp1, inv_tmp2), MODULO_P) # reductionはほぼ発生する．
      x_3 = opr.modulo(opr.multiply(lmd, lmd) - P.x - Q.x, MODULO_P)  # reductionはほぼ発生する。
      y_3 = opr.modulo(opr.multiply(lmd, (P.x - x_3)) - P.y, MODULO_P)  # reduction発生しない時もある。
  else:
    tmp1 = (Q.y - P.y)
    tmp2 = (Q.x - P.x)
    inv_tmp2 = modinv(tmp2, MODULO_P)
    lmd = opr.modulo(opr.multiply(tmp1, inv_tmp2), MODULO_P) # reduction発生しない時もある。
    x_3 = opr.modulo(opr.multiply(lmd, lmd) - P.x - Q.x, MODULO_P) # reductionはほぼ発生する。
    y_3 = opr.modulo(opr.multiply(lmd, (P.x - x_3)) - P.y, MODULO_P) # reduction発生しない時もある。
  return Point(x_3, y_3)

# 楕円曲線上のバイナリ法
# Q = dPを求める
def Binary(d: int, P: Point):
  global opr
  Q = Point(INF, INF)
  d_bit_seq = bin(d)[2:]
  for i in d_bit_seq:
    Q = PointAdd(Q, Q)  # 実測値(multiplication): 457.4243529411765
    if int(i) == 1:
      Q = PointAdd(P, Q)  # 実測値(multiplication): 452.6232283464567
  return Q

# MsgToPoint
def MsgToPoint(m: int):
  for i in range(0, N):
    x = 100*m+i
    # z = pow(x, 3, MODULO_P)
    z = x**3
    z %= MODULO_P
    z += A*x + B
    if sympy.legendre_symbol(z, MODULO_P) == 1:
      y = sympy.sqrt_mod(z, MODULO_P)
      return Point(x, y)

# PointToMsg
def PointToMsg(P: Point):
  return P.x // 100

# 楕円ElGamal暗号のEncryption
def ElGamalEnc(M: Point, GP: Point, PUBLIC_KEY: Point):
  r = random.randint(1, MODULO_P)
  return (Binary(r, GP), PointAdd(M, Binary(r, PUBLIC_KEY)))

# 楕円ElGamal暗号のDecryption
def ElGamalDec(C1: Point, C2: Point, d: int):
  C1.inv()
  B = Binary(d, C1)
  ans = PointAdd(C2, B)
  return ans


# =============================
# main function
# =============================
if __name__ == '__main__':
  # 生成元
  x = 0xffffffff222ffff9dcbbac55eeefffffffff
  GP = MsgToPoint(x)

  # 秘密鍵と公開鍵の生成
  d = 0x49be667ef9dcbbac55b06295ce870b0702900cdb2dce28d959f2815b16f81798
  d_bit_seq = bin(d)[2:]
  print("ハミング重み: {0}".format(d_bit_seq.count("1")))
  PUBLIC_KEY = Binary(d, GP)

  # 分析用
  all_mult = 0
  all_reduction = 0
  mult_list = []
  reduction_list = []

  # =========================

  correct = 0
  for i in range(LOOP):
    print("i = {0}".format(i))
    # 平文の生成
    plaintext = random.randint(1, MODULO_P//100)

    # 平文 → Point
    M = MsgToPoint(plaintext)

    # 暗号化・復号
    # コアロジック
    C1, C2 = ElGamalEnc(M, GP, PUBLIC_KEY) # 攻撃者も暗号文を自由に作成することができる．

    opr.mult = 0
    opr.reduction = 0
    DecM = ElGamalDec(C1, C2, d)
    print(vars(opr))
    all_mult += opr.mult
    all_reduction += opr.reduction
    print("平均 mult = {0}".format(all_mult/(i+1)))
    print("平均 reduction = {0}".format(all_reduction/(i+1)))
    mult_list.append(opr.mult)
    reduction_list.append(opr.reduction)

    # Point → 平文
    decPlaintext = PointToMsg(DecM)
    if plaintext == decPlaintext:
      correct += 1

  print("========== 結果 ==========")
  print(correct)
  print("ハミング重み: {0}".format(d_bit_seq.count("1")))
  print("平均 mult = {0}".format(all_mult/LOOP))
  print("平均 reduction = {0}".format(all_reduction/LOOP))
  print("平均 NR1 = {0}".format(NR1_SUM/LOOP))
  print("平均 NR2 = {0}".format(NR2_SUM/LOOP))
  print("tmp_mult_sum / LOOP = {0}".format(sum(tmp_mult_sum_list)/len(tmp_mult_sum_list)))

  plt.hist(mult_list, bins=LOOP)
  plt.show()
  plt.hist(reduction_list, bins=LOOP)
  plt.show()


