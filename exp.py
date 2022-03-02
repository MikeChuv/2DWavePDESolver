import numpy as np
from opensimplex import OpenSimplex

def testU(x, y, t, a, b):
	n = 1 # количество периодов по x - должно быть 1 (и больше) или сдвигать фазу(брать косинус) т.к. у синуса 0 в 0
	m = 1 # количество периодов по y - может быть меньше 1
	nu = 2 * np.pi * n / a
	mu = 2 * np.pi * m / b
	# lambda = nu + mu
	k = np.sqrt(nu ** 2 + mu ** 2)
	return np.sin(nu * x) * np.cos(mu * y) * np.sin(k * t) # тестовый пример 1
	# return np.sin(nu * x) * np.cos(mu * y) * np.sin(k * t) + 0.5 * np.sin(2 * nu * x) * np.cos(2 * mu * y) * np.sin(2 * k * t) # тестовый пример 2



x = np.linspace(-1, 1, 10)
y = np.linspace(-0.5, 0.5, 5)


x1, y1 = np.meshgrid(x, y)


u = testU(x1, y1, 0.1, 2, 1)

# print(u.shape)

verts = np.zeros((x1.shape[0]*x1.shape[1],)+(3,))

# print(x1.shape)
print(verts.shape)

for i in range(x1.shape[0]):
	for j in range(x1.shape[1]):
		verts[j + i * x1.shape[1], 0] = x1[i, j]
		verts[j + i * x1.shape[1], 1] = y1[i, j]
		verts[j + i * x1.shape[1], 2] = 0

# print(verts[0:10], )

ypoints = range(-20, 22, 1)
xpoints = range(-20, 22, 1)

faces = []
colors = []
for m in range(x1.shape[0] - 1):
	yoff = m * x1.shape[1]
	for n in range(x1.shape[1] - 1):
		faces.append([yoff + n, yoff + n + x1.shape[1], yoff + n + x1.shape[1] + 1])
		faces.append([yoff + n, yoff + n + 1, yoff + n + x1.shape[1] + 1])
		colors.append([0, 0, 0, 0])
		colors.append([0, 0, 0, 0])


print(verts.shape)

print(max(faces))