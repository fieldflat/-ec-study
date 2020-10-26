from hashlib import sha256

class ECDSA:
    def __init__(self):
        #各種定数
        self.Pcur = 2**256 - 2**32 - 2**9 - 2**8 - 2**7 - 2**6 - 2**4 - 1
        self.A = 0
        self.B = 7
        G_x = 55066263022277343669578718895168534326250603453777594175500187360389116729240
        G_y = 32670510020758816978083085130507043184471273380659243275938904335757337482424
        self.G = (G_x, G_y)
        self.secret_key = None
        self.public_key = None

    def key_pair(self, secret):  # 文字列をバイナリにしてec_multを呼びdouble and add使った秘密鍵回の加算を行いで秘密鍵と公開鍵のペアを得る
        secret_hash = sha256(secret.encode('utf-8')).hexdigest()
        self.secret_key = int(secret_hash, 16)
        public_point = self.ec_multi(self.secret_key)
        self.public_key = int("04" + "%064x" %
                              public_point[0] + "%064x" % public_point[1], 16)
        return (self.public_key, self.secret_key)

    def ec_add(self, p, q):  # 楕円曲線上での加算
        tmp = ((q[1]-p[1]) * self.inv_mod(q[0]-p[0])) % self.Pcur
        x = (tmp ** 2 - p[0] - q[0]) % self.Pcur
        y = (tmp * (p[0] - x) - p[1]) % self.Pcur
        return (x, y)

    def ec_double(self, P):  # 楕円曲線上での倍算
        tmp = ((3 * (P[0] ** 2) + self.A) * self.inv_mod(2 * P[1])) % self.Pcur
        x = (tmp ** 2 - 2 * P[0]) % self.Pcur
        y = (tmp * (P[0] - x) - P[1]) % self.Pcur
        return (x, y)

    def ec_multi(self, scalar):  # double and addで加算と倍算を繰り返し秘密鍵から公開鍵を作る
        if scalar == 0:
            pass  # エラー処理しないとマズいが適当にpass
        scalar_bin = str(bin(scalar))[2:]
        G_p = self.G
        for i in range(1, len(scalar_bin)):
            G_p = self.ec_double(G_p)
            if scalar_bin[i] == "1":
                G_p = self.ec_add(G_p, self.G)
        return (G_p)

    def inv_mod(self, a):  # 拡張ユークリッド互除法
        low, high = a % self.Pcur, self.Pcur
        c0, c1 = 1, 0
        while low > 1:
            r = high // low
            c2 = c1 - c0 * r
            new = high - low * r
            high = low
            low = new
            c1 = c0
            c0 = c2
        return c0 % self.Pcur


if __name__ == "__main__":
    ecdsa = ECDSA()
    # Satoshi Nakamotoをサンプルとして秘密鍵としている例が世の中には多い(当然これで作られた鍵ペアは使っちゃダメですよ
    ecdsa.key_pair("Satoshi Nakamoto")
    print(ecdsa.secret_key, ecdsa.public_key)
