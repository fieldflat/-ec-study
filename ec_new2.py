import math

def egcd(a, b):
    (x, lastx) = (0, 1)
    (y, lasty) = (1, 0)
    while b != 0:
        q = a // b
        (a, b) = (b, a % b)
        (x, lastx) = (lastx - q * x, x)
        (y, lasty) = (lasty - q * y, y)
    return (lastx, lasty, a)

# ax ≡ 1 (mod m)
def modinv(a, m):
    (inv, q, gcd_val) = egcd(a, m)
    return inv % m

# 点加算 (2点が異なる場合の足し算)
def ec1_plus(point_p, point_q, p):
  x_1 = point_p[0]
  y_1 = point_p[1]
  x_2 = point_q[0]
  y_2 = point_q[1]
  margin_y = y_2 - y_1
  margin_x = x_2 - x_1
  inv_margin_x = modinv(margin_x, p)
  s = (margin_y * inv_margin_x) % p
  x_3 = ((s ** 2) - x_1 - x_2) % p
  y_3 = (-(s * (x_3 - x_1) + y_1)) % p
  print('点加算 ({0}, {1}) + ({2}, {3}) = ({4}, {5})'.format(x_1,
                                                          y_1, x_2, y_2, x_3, y_3))
  return (x_3, y_3)

# 点2倍加算 (2点が同じ場合の足し算)
def ec2_plus(point_p, point_q, p, a):
  x_1 = point_p[0]
  y_1 = point_p[1]
  x_2 = point_q[0]
  y_2 = point_q[1]
  s_1 = 3*(x_1**2) + a
  s_2 = 2 * y_1
  inv_s_2 = modinv(s_2, p)
  s = (s_1 * inv_s_2) % p
  x_3 = ((s ** 2) - x_1 - x_2) % p
  y_3 = (-(s * (x_3 - x_1) + y_1)) % p
  print('点2倍加算 ({0}, {1}) + ({2}, {3}) = ({4}, {5})'.format(x_1,
                                                            y_1, x_2, y_2, x_3, y_3))
  return (x_3, y_3)

# 加算
def ec_plus(point_p, point_q, p, a):
  x_1 = point_p[0]
  y_1 = point_p[1]
  x_2 = point_q[0]
  y_2 = point_q[1]
  if point_p == point_q:
    s_1 = 3*(x_1**2) + a
    s_2 = 2 * y_1
    inv_s_2 = modinv(s_2, p)
    s = (s_1 * inv_s_2) % p
  else:
    margin_y = y_2 - y_1
    margin_x = x_2 - x_1
    inv_margin_x = modinv(margin_x, p)
    s = (margin_y * inv_margin_x) % p
  x_3 = ((s ** 2) - x_1 - x_2) % p
  y_3 = (-(s * (x_3 - x_1) + y_1)) % p
  print('({0}, {1}) + ({2}, {3}) = ({4}, {5})'.format(x_1, y_1, x_2, y_2, x_3, y_3))
  return (x_3, y_3)


if __name__ == '__main__':
  point_p = (7, 6)
  point_q = (10, 2)
  p = 17
  point_r = ec1_plus(point_p, point_q, p)
  # print(point_r)
  a = -3
  point_r = ec2_plus(point_p, point_p, p, a)

  tmp_point_p = point_p
  for i in range(2, 23 + 1):
    print("{0}P = ".format(i), end="")
    tmp_point_p = ec_plus(point_p, tmp_point_p, p, a)
