from scipy.signal import convolve2d
import pywt
import numpy as np


def W(x, level=4):
    coef = pywt.wavedec2(x, 'db2', 'periodization', level=level)
    array, slices = pywt.coeffs_to_array(coef)
    return array, slices


def WH(array, slices):
    coef = pywt.array_to_coeffs(array, slices, output_format="wavedec2")
    return pywt.waverec2(coef, 'db2', 'periodization')


def DivergenceIm(p1, p2):
    z = p2[:, 1:-1] - p2[:, :-2]
    v = np.column_stack((p2[:, 0], z, -p2[:, -1]))

    z = p1[1:-1, :] - p1[:-2, :]
    u = np.row_stack((p1[0, :], z, -p1[-1, :]))

    divp = v + u
    return divp


def GradientIm(u):
    z = u[1:, :] - u[:-1, :]
    dux = np.row_stack((z, np.zeros(z.shape[1])))

    z = u[:, 1:] - u[:, :-1]
    duy = np.column_stack((z, np.zeros(z.shape[0])))

    return dux, duy


def chambolle_prox_TV_stop(g, alpha, MaxIter=8, tol=1e-4, px=np.array([0]), py=np.array([0])):
    if (px.shape != g.shape) or (py.shape != g.shape):
        px = np.zeros(g.shape)
        py = np.zeros(g.shape)
    tau = 0.249
    cont = True
    k = 0

    while cont:
        k = k + 1
        divp = DivergenceIm(px, py)
        u = divp - g/alpha
        upx, upy = GradientIm(u)
        tmp = np.sqrt(upx**2 + upy**2)
        err = np.sqrt(((tmp*px-upx)**2 + (tmp*py-upy)**2).sum())
        px = (px + tau * upx)/(1 + tau * tmp)
        py = (py + tau * upy)/(1 + tau * tmp)
        cont = (k < MaxIter) and (err > tol)

    f = g - alpha*DivergenceIm(px, py)
    return f, px, py


def denoiseTV(g, alpha, MaxIter=10, tol=1e-2, px=np.array([0]), py=np.array([0])):
    f, _, _ = chambolle_prox_TV_stop(g, alpha, MaxIter, tol, px, py)
    return f


def TVnorm(x):
    dh, dv = GradientIm(x)
    J = (np.sqrt(dh**2 + dv**2)).sum()
    return J


def TVnorm_ani(x):
    dh, dv = GradientIm(x)
    J = (np.sqrt(dh**2 + dv**2)).sum()
    return J


def MSE(m, m_):
    e = m-m_
    return np.sum(np.abs(e)**2)/np.shape(m)[0]/np.shape(m)[1]
