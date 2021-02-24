import math

INF_POINT = None


class EllipticCurve:
	def __init__(self, a, b, p):
		self.a = a
		self.b = b
		self.p = p
		self.points = []
		self.define_points()


	def define_points(self):
		self.points.append(INF_POINT)
		for x in range(self.p):
			for y in range(self.p):
				if self.equal_modp(y * y, x * x * x + self.a * x + self.b):
					self.points.append((x,y))	


	def addition(self, P1, P2):
		if P1 == INF_POINT:
			return P2
		if P2 == INF_POINT:
			return P1

		x1 = P1[0]
		y1 = P1[1]
		x2 = P2[0]
		y2 = P2[1]

		if self.equal_modp(x1, x2) and self.equal_modp(y1, -y2):
			return INF_POINT

		if self.equal_modp(x1, x2) and self.equal_modp(y1, y2):
			u = self.reduce_modp((3 * x1 * x1 + self.a) * self.inverse_modp(2 * y1))
		else:
			u = self.reduce_modp((y1 - y2) * self.inverse_modp(x1 - x2))

		v = self.reduce_modp(y1 - u * x1)
		x3 = self.reduce_modp(u * u - x1 - x2)
		y3 = self.reduce_modp(-u * x3 - v)
		return (x3, y3)


	def test_associativity(self):
		n = len(self.points)
		for i in range (n):
			for j in range(n):
				for k in range(n):
					P = self.addition(self.points[i], self.addition(self.points[j], self.points[k]))
					Q = self.addition(self.addition(self.points[i], self.points[j]), self.points[k])
					if P != Q:
						return False
		return True


	def number_points(self):
		return len(self.points)


	def discriminant(self):
		D = -16 *(4 * self.a * self.a * self.a + 27 * self.b * self.b)
		return self.reduce_modp(D)


	def print_points(self):
		print(self.points)


	def verify_hasse(self):
		return abs(self.number_points() - (self.p + 1)) <= 2 * math.sqrt(p)


	def lower_hasse_bound(self):
		return math.ceil(self.p + 1 - 2 * math.sqrt(self.p))


	def upper_hasse_bound(self):
		return math.floor(self.p + 1 + 2 * math.sqrt(self.p))

	# helper functions

	def reduce_modp(self, x):
		return x % self.p


	def equal_modp(self, x, y):
		return self.reduce_modp(x - y) == 0


	def inverse_modp(self, x):
		for y in range(self.p):
			if self.equal_modp(x * y, 1):
				return y
		return None



p = 101

epsilon = 2 * math.sqrt(p)
lower_hasse_bound = math.ceil((p + 1) - epsilon)
upper_hasse_bound = math.floor((p + 1) + epsilon)

group_orders = []

for a in range(0, p):
	for b in range(0, p):
		ec = EllipticCurve(a, b, p)
		if ec.discriminant() == 0:
			continue
		group_orders.append(ec.number_points())
		if ec.verify_hasse() == False:
			# this case should not happen
			print("THE HASSE THEOREM FAILED???")
			break


print("Hasse lower bound =", lower_hasse_bound)
print("Minimum group order =", min(group_orders))

print("Hasse upper bound =", upper_hasse_bound)
print("Maximum group order =", max(group_orders))


# test to see that all numbers in the "Hasse interval" are
# equal to the order of some elliptic curve

not_group_orders = []
for i in range(lower_hasse_bound, upper_hasse_bound + 1):
	if i not in group_orders:
		not_group_orders.append(i)

if len(not_group_orders) == 0:
	print("Every integer in the possibe range is a group order.")
else:
	# this case should not happen
	print("The following are not the orders of any elliptic curve.")
	print(not_group_orders)

# average group order
print("Average group order =", sum(group_orders) / len(group_orders))
print("p + 1 =", p + 1)

