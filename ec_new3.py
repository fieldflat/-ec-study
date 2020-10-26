import math
import sys

# =============================
# グローバル定数
# =============================
INF = -1
MODULO_P = 17
A = -3
B = 3


# =============================
# サブ関数 (主要関数内で使用されるサブ関数)
# =============================

# 座標を表すclass
class Point:
  def __init__(self, x: int, y: int):
    self.x = x
    self.y = y

  # 点Pが無限遠点であるかチェックする
  def is_inf(self):
    return (self.x == self.y == INF)

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

if __name__ == '__main__':
  P = Point(7, 6)
  Q = Binary(22, P)
