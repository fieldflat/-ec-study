# 楕円ElGamal暗号

import math
import sys
import random
import sympy

# =============================
# グローバル定数
# =============================
INF = -1
MODULO_P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f # 256bit
A = 0x0000000000000000000000000000000000000000000000000000000000000000
B = 0x0000000000000000000000000000000000000000000000000000000000000007
N = 100
LOOP = 100


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

# =============================
# サブ関数 (主要関数内で使用されるサブ関数)
# =============================

# ユークリッド互除法
def egcd(a: int, b: int):
    (x, lastx) = (0, 1)
    (y, lasty) = (1, 0)
    while b != 0:
        q = a // b
        (a, b) = (b, a % b)
        (x, lastx) = (lastx - q * x, x)
        (y, lasty) = (lasty - q * y, y)
    return (lastx, lasty, a)

# ax ≡ 1 (mod m)
def modinv(a: int, m: int):
    (inv, _, _) = egcd(a, m)
    return inv % m

# =============================
# 主要関数
# =============================

# 楕円曲線上の点演算
def PointAdd(P: Point, Q: Point):
  if P.is_inf():
    return Q
  elif Q.is_inf():
    return P
  elif P.x == Q.x:
    if ((P.y + Q.y) % MODULO_P == 0):
      return Point(INF, INF)
    else:
      tmp1 = (3*(P.x)**2 + A)
      tmp2 = 2*P.y
      inv_tmp2 = modinv(tmp2, MODULO_P)
      lmd = (tmp1 * inv_tmp2) % MODULO_P
  else:
    tmp1 = (Q.y - P.y)
    tmp2 = (Q.x - P.x)
    inv_tmp2 = modinv(tmp2, MODULO_P)
    lmd = (tmp1 * inv_tmp2) % MODULO_P
  x_3 = (lmd**2 - P.x - Q.x) % MODULO_P
  y_3 = (lmd*(P.x - x_3) - P.y) % MODULO_P
  return Point(x_3, y_3)

# 楕円曲線上のバイナリ法
# Q = dPを求める
def Binary(d: int, P: Point):
  Q = Point(INF, INF)
  d_bit_seq = bin(d)[2:]
  for i in d_bit_seq:
    Q = PointAdd(Q, Q)
    if int(i) == 1:
      Q = PointAdd(Q, P)
  return Q


# MsgToPoint
def MsgToPoint(m: int):
  for i in range(0, N):
    x = 100*m+i
    z = pow(x, 3, MODULO_P)
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
def ElGamalDec(C1, C2, d):
  C1.inv()
  return PointAdd(C2, Binary(d, C1))


# =============================
# main function
# =============================
if __name__ == '__main__':
  # 生成元
  x = 0xffffffffffffffffffffffffeeefffffffff
  GP = MsgToPoint(x)

  # 秘密鍵と公開鍵の生成
  d = 0x49be667ef9dcbbac55b06295ce870b07029b3cdb2dce28d959f2815b16f81798
  PUBLIC_KEY = Binary(d, GP)

  # =========================

  correct = 0
  for i in range(LOOP):
    # 平文の生成
    plaintext = random.randint(1, MODULO_P//100)

    # 平文 → Point
    M = MsgToPoint(plaintext)

    # 暗号化・復号
    C1, C2 = ElGamalEnc(M, GP, PUBLIC_KEY)
    DecM = ElGamalDec(C1, C2, d)

    # Point → 平文
    decPlaintext = PointToMsg(DecM)

    if plaintext == decPlaintext:
      correct += 1
  
  print(correct)


