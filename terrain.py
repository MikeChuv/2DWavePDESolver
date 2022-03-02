import numpy as np
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import sys
from opensimplex import OpenSimplex

# первое начальное условие из примера задачи u(x, y, 0): 
def phi(x, y, a, b):
	return np.arctan(np.cos(np.pi * x / a))

# второе начальное условие из примера задачи dudt(x, y, 0): 
def psi(x, y, a, b):
	return np.sin(2 * np.pi * x / a) * np.sin(np.pi * y / b)

# тестовая функция (первое начальное условие - u(x, y, 0))
def testphi(x, y, a, b):
	return np.zeros((x.shape[0], y.shape[1]))

# тестовая функция (второе начальное условие - dudt(x, y, 0))
def testpsi(x, y, a, b):
	n = m = 1
	nu = 2 * np.pi * n / a
	mu = 2 * np.pi * m / b
	# lambda = nu + mu
	k = np.sqrt(nu ** 2 + mu ** 2)
	# return np.sin(nu * x) * np.cos(mu * y) # тестовый пример 1
	return np.sin(nu * x) * np.cos(mu * y) + np.sin(2 * nu * x) * np.cos(2 * mu * y)
	
# тестовый пример - функция u(x, y, t)
def testU(x, y, t, a, b):
	n = m = 1
	nu = 2 * np.pi * n / a
	mu = 2 * np.pi * m / b
	# lambda = nu + mu
	k = np.sqrt(nu ** 2 + mu ** 2)
	# return (1 / k) * np.sin(nu * x) * np.cos(mu * y) * np.sin(k * t) # тестовый пример 1
	return (1 / k) * np.sin(nu * x) * np.cos(mu * y) * np.sin(k * t) + (1 / (2*k)) * np.sin(2 * nu * x) * np.cos(2 * mu * y) * np.sin(2 * k * t) # тестовый пример 2


class SolverVec:
	u = None
	def __init__(self, a, b, step, ut0, dudt0):
		self.a = float(a)
		self.b = float(b)
		self.dx = step
		self.dy = step
		self.nodes_x = int(self.a / self.dx)
		self.nodes_y = int(self.b / self.dy)
		self.x = np.linspace(-self.a/2, self.a/2, self.nodes_x)
		self.y = np.linspace(-self.b/2, self.b/2, self.nodes_y)
		self.u0 = np.zeros((self.nodes_y, self.nodes_x))
		self.u1 = np.zeros((self.nodes_y, self.nodes_x))
		self.I = ut0
		self.V = dudt0
		self.invdxsq = 1 / self.dx ** 2
		self.invdysq = 1 / self.dy ** 2
	
	def uxx_matr(self, first_step):
		u = self.u0 if first_step else self.u1
		uxx = np.zeros((u.shape[0], u.shape[1]))
		# uxx[:, 0] *= 0
		# uxx[:, -1] *= 0
		# uxx = (u[1:-1, :-2] - 2 * u[1:-1, 1:-1] + u[1:-1, 2:]) * self.invdxsq
		uxx[:, 1:-1] = (u[:, :-2] - 2 * u[:, 1:-1] + u[:, 2:]) * self.invdxsq
		#print(uxx.shape)
		return uxx

	def uyy_matr(self, first_step):
		u = self.u0 if first_step else self.u1
		uyy = np.zeros((u.shape[0], u.shape[1]))
		uyy[0, :] = (2 * u[1, :] - 2 * u[0, :]) * self.invdysq
		uyy[-1, :] = (2 * u[-2, :] - 2 * u[-1, :]) * self.invdysq
		uyy[1:-1, :] = (u[:-2, :] - 2 * u[1:-1, :] + u[2:, :]) * self.invdysq
		# (u[:-2, 1:-1] - 2 * u[1:-1, 1:-1] + u[2:, 1:-1]) * self.invdysq
		#print(uyy.shape)
		return uyy

	def first_step(self, dt):
		for i in range(self.nodes_y):
			for j in range(self.nodes_x):
				self.u0[i, j] = self.I(self.x[j], self.y[i], self.a, self.b)
		uxxf = self.uxx_matr(True)
		uyyf = self.uyy_matr(True)
		for i in range(self.nodes_y):
			for j in range(self.nodes_x):
				self.u1[i, j] = self.u0[i, j] + dt * self.V(self.x[j], self.y[i], self.a, self.b) + dt**2 * (uxxf[i, j] + uyyf[i, j]) / 2

	def advance(self, dt):
		if self.u is not None:
			self.u0 = self.u1
			self.u1 = self.u	
		self.u = 2*self.u1.copy() - self.u0.copy()
		self.u[:, 1:-1] += dt**2 * (self.uxx_matr(False)[:, 1:-1] + self.uyy_matr(False)[:, 1:-1] )
		self.u[:,0] *= 0
		self.u[:, -1] *= 0
		# self.u[0,:] *= 0
		# self.u[-1, :] *= 0
		return self.u





class Terrain(object):
	def __init__(self, *args):
		"""
		Initialize the graphics window and mesh
		"""

		# setup the view window
		self.app = QtGui.QApplication(sys.argv)
		self.w = gl.GLViewWidget()
		self.w.setGeometry(0, 110, 1280, 720)
		self.w.show()
		self.w.setWindowTitle('Terrain')
		self.w.setCameraPosition(distance=30, elevation=8)

		if len(args) == 5 : # решаем диффур
			self.mode = 'task'
			self.a, self.b, self.h, self.phi, self.psi = args
			self.dt = round(self.h*self.h/np.sqrt(self.h**2 + self.h**2), 4) - (3 * 1e-4)
			self.solverVec = SolverVec(self.a, self.b, self.h, self.phi, self.psi)
			self.solverVec.first_step(self.dt)
			self.meshSize = (int(self.b / self.h), int(self.a / self.h))
			self.x, self.y = np.meshgrid(self.solverVec.x, self.solverVec.y)
		else: # тестовый пример
			self.mode = 'test'
			self.a, self.b, self.h, self.U = args
			self.meshSize = (int(self.b / self.h), int(self.a / self.h))
			x = np.linspace(-self.a/2, self.a/2, self.meshSize[1])
			y = np.linspace(-self.b/2, self.b/2, self.meshSize[0])
			self.x, self.y = np.meshgrid(x, y)
			self.dt = round(self.h*self.h/np.sqrt(self.h**2 + self.h**2), 4) - (3 * 1e-4)
		#if self.mode == 'task': self.t = 2 * self.dt

		# create the veritices array
		verts = np.zeros((self.x.shape[0]*self.x.shape[1],)+(3,))

		for i in range(self.x.shape[0]):
			for j in range(self.x.shape[1]):
				verts[j + i * self.x.shape[1], 0] = self.x[i, j]
				verts[j + i * self.x.shape[1], 1] = self.y[i, j]
				verts[j + i * self.x.shape[1], 2] = 0

		# create the faces and colors arrays
		faces = []
		colors = []
		for m in range(self.meshSize[0] - 1):
			yoff = m * self.meshSize[1]
			for n in range(self.meshSize[1] - 1):
				faces.append([n + yoff, yoff + n + self.meshSize[1], yoff + n + self.meshSize[1] + 1])
				faces.append([n + yoff, yoff + n + 1, yoff + n + self.meshSize[1] + 1])
				colors.append([0, 0, 0, 0])
				colors.append([0, 0, 0, 0])

		faces = np.array(faces)
		colors = np.array(colors)

		# create the mesh item
		self.m1 = gl.GLMeshItem(
			vertexes=verts,
			faces=faces, faceColors=colors,
			smooth=False, drawEdges=False,
		)
		self.m1.setGLOptions('opaque') # additive opaque translucent
		self.w.addItem(self.m1)
		self.t = 0

	def update(self):
		"""
		update the mesh and shift the noise each time
		"""
		if self.mode == 'test':
			u = self.U(self.x, self.y, self.t, self.a, self.b)
			self.t += self.dt * 0.1
		else:
			u = self.solverVec.advance(self.dt)

		verts = np.zeros((self.x.shape[0]*self.x.shape[1],)+(3,))

		for i in range(self.x.shape[0]):
			for j in range(self.x.shape[1]):
				verts[j + i * self.x.shape[1], 0] = self.x[i, j]
				verts[j + i * self.x.shape[1], 1] = self.y[i, j]
				verts[j + i * self.x.shape[1], 2] = u[i, j]

		faces = []
		colors = []
		for m in range(self.meshSize[0] - 1):
			yoff = m * self.meshSize[1]
			for n in range(self.meshSize[1] - 1):
				faces.append([n + yoff, yoff + n + self.meshSize[1], yoff + n + self.meshSize[1] + 1])
				faces.append([n + yoff, yoff + n + 1, yoff + n + self.meshSize[1] + 1])
				colors.append([n / self.meshSize[1], 1 - n / self.meshSize[1], m / self.meshSize[0], 1])
				colors.append([n / self.meshSize[1], 1 - n / self.meshSize[1], m / self.meshSize[0], 1])

		faces = np.array(faces, dtype=np.uint32)
		colors = np.array(colors, dtype=np.float32)

		self.m1.setMeshData(
			vertexes=verts, faces=faces, faceColors=colors
		)
		

	def start(self):
		"""
		get the graphics window open and setup
		"""
		if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
			QtGui.QApplication.instance().exec_()

	def animation(self):
		"""
		calls the update method to run in a loop
		"""
		timer = QtCore.QTimer()
		timer.timeout.connect(self.update)
		timer.start(10)
		self.start()
		self.update()


if __name__ == '__main__':
	# (Task) args: a, b, h, phi, psi
	# (Test) args: a, b, h, U
	# t = Terrain(2, 1, 0.05, testphi, testpsi)
	t = Terrain(2, 1, 0.05, testU)
	t.animation()