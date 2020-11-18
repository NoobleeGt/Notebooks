import matplotlib.pyplot as plt
plt.style.use('../my_params.mplstyle')


def solve_euler(system, t, x, u, h):
    dx_dt = system.x_prime(t, x, u)

    new_state = [state + h * dstate for state, dstate in zip(x, dx_dt)]

    return new_state


def run(system, t, x, u):
    h = t[1] - t[0]
    for i in range(len(t) - 1):
        x0 = [state[-1] for state in x]

        new_state = solve_euler(system, t[i], x0, u[i], h)
        
        for i, state in enumerate(new_state):
            x[i].append(state)


class Example1:
    def __init__(self):
        pass

    def x_prime(self, t, x, u):
        dx_dt = u - abs(x[0]) * x[0]

        return [dx_dt]


def example1(u):
    system = Example1()

    t_vec = [i/500 for i in range(10001)]

    x = [[0]]
    u_vec = [u if i <= 5 else 0 for i in t_vec]

    run(system, t_vec, x, u_vec)

    fig, ax = plt.subplots()
    ax.plot(t_vec, u_vec, '-', label='u')
    ax.plot(t_vec, x[0], '--', label='v')
    ax.set_xlim(0, 20)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('thrust u, velocity v')
    ax.legend()



