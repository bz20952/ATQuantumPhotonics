import numpy as np
from scipy import linalg
from vals import b, density, E, zeta
from matplotlib import pyplot as plt
import smoa


plt.rcParams.update({'font.size': 18})


def mass_mat(L):

    return np.array([
        [ 156, 22*L, 54, -13*L],
        [ 22*L, 4*L**2, 13*L, -3*L**2],
        [ 54, 13*L, 156, -22*L],
        [-13*L, -3*L**2, -22*L, 4*L**2]
    ])


def stiffness_mat(L):

    return np.array([
        [12, 6*L, -12 , 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]
    ])


def fe(beam: str):
     
    # Element 1
    h_1 = 0.18E-6
    L_1 = 7.5E-6
    I_1, A_1 = smoa.rect_beam(b, h_1)
    M_1 = ((density * A_1 * L_1)/420) * mass_mat(L_1)
    K_1 = ((E*I_1)/(L_1**3)) * stiffness_mat(L_1)

    # Element 2
    h_2 = 0.34E-6
    L_2 = 25E-6/2
    if beam == 'i':
        I_2, A_2 = smoa.i_beam(0.24E-6, 0.05E-6, 0.05E-6, 0.16E-6)
    elif beam == 't':
        I_2, A_2 = smoa.t_beam(0.29E-6, 0.05E-6, 0.05E-6, 0.16E-6)
    else:
        I_2, A_2 = smoa.rect_beam(b, h_2)
    I_0, _ = smoa.rect_beam(b, h_2)
    m_0 = density * b * h_2 * L_2 * 2
    print(I_2/I_0)
    print(((density * A_2 * L_2)*2)/m_0)
    M_2 = ((density * A_2 * L_2)/420) * mass_mat(L_2)
    K_2 = ((E*I_2)/(L_2**3)) * stiffness_mat(L_2)

    # Element 3
    A_3 = A_1
    L_3 = L_1
    I_3 = I_1
    M_3 = ((density * A_3 * L_3)/420) * mass_mat(L_3)
    K_3 = ((E*I_3)/(L_3**3)) * stiffness_mat(L_3)

    # Assemble matrices
    M = np.zeros((10,10))
    M[:4,:4] += M_1
    M[2:6,2:6] += M_2
    M[4:8,4:8] += M_2
    M[6:,6:] += M_3

    K = np.zeros((10,10))
    K[:4,:4] += K_1
    K[2:6,2:6] += K_2
    K[4:8,4:8] += K_2
    K[6:,6:] += K_3

    # print("Mass matrix:\n", M)
    # print("Stiffness matrix:\n", K)

    # Apply simply-supported BCs
    M = np.delete(M, 0, axis=0)
    M = np.delete(M, 0, axis=1)
    M = np.delete(M, -2, axis=0)
    M = np.delete(M, -2, axis=1)

    K = np.delete(K, 0, axis=0)
    K = np.delete(K, 0, axis=1)
    K = np.delete(K, -2, axis=0)
    K = np.delete(K, -2, axis=1)

    # # Damping (assumes Rayleigh damping)
    # a1 = np.array([[w[0]**2, w[0]],
    #                [w[1]**2, w[1]]])
    # a2 = np.array([1E-10, 1E-10]).reshape(2,1)
    # alpha_beta = 2*(np.linalg.inv(a1) @ a2)
    # C = alpha_beta[0,0]*M + alpha_beta[1,0]*K

    # Solve for eigenvalues
    eigenvalues, _ = linalg.eig(K, M)
    index = eigenvalues.argsort()
    w_squared = eigenvalues[index]

    # Sqrt and convert to MHz
    w = np.sqrt(w_squared)
    f = (np.sqrt(w_squared) / (2 * np.pi)) / 1E6
    # print(f)
    f_damped = f[0] * (1 - (zeta**2))
    print(f_damped)

    return M, K, w


def get_frequency_response(M, K, w, frequencies):

    """Analyse response amplitude over a range of frequencies using the FRF formula."""

    gains = []
    for w in frequencies:
        s = complex(0, w)

        # Gives response to impulse (delta function)
        frf = np.linalg.inv(K - (w**2)*M)

        tf_gain = abs(frf)
        gains.append(tf_gain)

    return gains


def plot_frf():

    """Generic FRF plotting."""

    fig, ax = plt.subplots(3, 1)

    frequencies = np.linspace(0, 1.2, 800)*1.5E7 

    for j, beam in enumerate(['i', 't', 'rect'][::-1]):
        M, K, w = fe(beam)
        gains = get_frequency_response(M, K, w, frequencies)

        max_gain = 0
        for i in [1, 3]:
            max_gain = max([gain[i, i] for gain in gains])

        for i in [1, 3]:
            ax[j].plot(frequencies / (2*np.pi) / 1E6, [gain[i, i] for gain in gains]/max_gain, label=f'N{(i+1)/2:.0f}')

        ax[j].set_yscale('log')
        ax[j].grid(True)

    # ax[-1].set_xticks(fontsize=18)
    # ax[-1].set_xlabel(r'$f$ [MHz]', fontsize=18)
    # plt.ylabel('Amplitude [m/N]', fontsize=18)
    # plt.yscale('log')
    ax[0].legend(fontsize=18)

    plt.show()


if __name__ == '__main__':
    plot_frf()