import numpy as np
import matplotlib.pyplot as plt


def f(x):
	return 1 / (1 + x ** 2)


def linear_spline(nodes, values, x):
	""" 线性样条插值 """
	n = len(nodes) - 1
	for i in range(n):
		if nodes[i] <= x <= nodes[i + 1]:
			return values[i] + (values[i + 1] - values[i]) / (nodes[i + 1] - nodes[i]) * (x - nodes[i])
	return None


def natural_cubic_spline(nodes, values, x):
	""" 自然三次样条插值 """
	n = len(nodes) - 1
	h = np.diff(nodes)
	alpha = np.zeros(n)
	for i in range(1, n):
		alpha[i] = (3 / h[i]) * (values[i + 1] - values[i]) - (3 / h[i - 1]) * (values[i] - values[i - 1])

	l = np.ones(n + 1)
	mu = np.zeros(n)
	z = np.zeros(n + 1)

	for i in range(1, n):
		l[i] = 2 * (nodes[i + 1] - nodes[i - 1]) - h[i - 1] * mu[i - 1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

	b = np.zeros(n)
	c = np.zeros(n + 1)
	d = np.zeros(n)

	for j in range(n - 1, -1, -1):
		c[j] = z[j] - mu[j] * c[j + 1]
		b[j] = (values[j + 1] - values[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
		d[j] = (c[j + 1] - c[j]) / (3 * h[j])

	for i in range(n):
		if nodes[i] <= x <= nodes[i + 1]:
			dx = x - nodes[i]
			return values[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
	return None


def hermite_spline(nodes, values, derivatives, x):
	""" 三次 Hermite 样条插值 """
	n = len(nodes) - 1
	for i in range(n):
		if nodes[i] <= x <= nodes[i + 1]:
			h = nodes[i + 1] - nodes[i]
			t = (x - nodes[i]) / h
			h00 = (1 + 2 * t) * (1 - t) ** 2
			h10 = t * (1 - t) ** 2
			h01 = t ** 2 * (3 - 2 * t)
			h11 = t ** 2 * (t - 1)
			return h00 * values[i] + h10 * h * derivatives[i] + h01 * values[i + 1] + h11 * h * derivatives[i + 1]
	return None


def compute_error(nodes, method, *args):
	""" 计算误差 """
	midpoints = (nodes[:-1] + nodes[1:]) / 2
	errors = [abs(f(x) - method(nodes, *args, x)) for x in midpoints]
	return max(errors)


# _nodes = np.array([0, 5 / 3, 10 / 3, 5])
node_num = np.array([4, 11, 21, 41, 1601, 3201])
lin_err = np.zeros(node_num.size)
cub_err = np.zeros(node_num.size)
her_err = np.zeros(node_num.size)

for i in range(node_num.size):
	_nodes = np.linspace(0, 5, node_num[i])
	_values = f(_nodes)
	_derivatives = -2 * _nodes / (1 + _nodes ** 2) ** 2  # 计算导数

	grid = np.linspace(0, 5, 100)
	lin_err[i] = compute_error(_nodes, linear_spline, _values)
	cub_err[i] = compute_error(_nodes, natural_cubic_spline, _values)
	her_err[i] = compute_error(_nodes, hermite_spline, _values, _derivatives)

lin_ord = np.log(lin_err[0:-2] / lin_err[1:-1]) / np.log(node_num[1:-1] / node_num[0:-2])
cub_ord = np.log(cub_err[0:-2] / cub_err[1:-1]) / np.log(node_num[1:-1] / node_num[0:-2])
her_ord = np.log(her_err[0:-2] / her_err[1:-1]) / np.log(node_num[1:-1] / node_num[0:-2])

print(lin_ord)
print(cub_ord)
print(her_ord)
