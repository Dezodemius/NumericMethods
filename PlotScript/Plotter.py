import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


data_path = "data/"


def main():
    def animate(i, *data):
        part1 = [float(x) for x in data[0][i]]
        part2 = [float(x) for x in data[1][i]]

        x1 = np.linspace(t0, t1, len(part1))
        x2 = np.linspace(t0, t1, len(part2))

        plt.clf()
        plt.plot(x1, part1, color='black')
        plt.plot(x2, part2, color='red')

        if (i % 10) == 0:
            print('{0}%'.format(i * 100 / len(exact_sol_data_parts)))
        if i == (len(exact_sol_data_parts) // 2):
            plt.savefig("pics/{}.png".format(sys.argv[5]))
        return line,

    def init():
        line.set_data([], [])
        return line,

    def save_pics(data_parts_, n):
        i = 0
        for part in data_parts_:
            if i % n == 0:
                part = [float(x) for x in part]
                plt.plot(np.linspace(t0, t1, len(part)), part)
                plt.savefig("pics/exact_sol_data{}.png".format(i))
                plt.clf()
            i += 1

    t0 = int(sys.argv[1])
    t1 = int(sys.argv[2])

    exact_sol_name = sys.argv[3]
    f = open(data_path + exact_sol_name, "r")
    exact_sol_data = f.read()
    f.close()

    numeric_sol_name = sys.argv[4]
    f = open(data_path + numeric_sol_name, "r")
    numeric_sol_data = f.read()
    f.close()

    exact_sol_data_parts = [x.split() for x in exact_sol_data.replace(",", ".").split('\n\n')]
    numeric_sol_data_parts = [x.split() for x in numeric_sol_data.replace(",", ".").split('\n\n')]

    fig = plt.figure()
    ax = plt.axes(xlim=(t0, t1), ylim=(0.0, 1.2))
    line, = ax.plot([], [], lw=3)

    anim = FuncAnimation(fig, animate, fargs=[exact_sol_data_parts, numeric_sol_data_parts], init_func=init,
                         frames=len(numeric_sol_data_parts), interval=1, blit=True, cache_frame_data=True)
    anim.save('gifs/{0}.gif'.format(sys.argv[5]), writer='pillow', fps=60)


main()
