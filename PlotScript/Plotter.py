import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Путь к папке с численными расчётами.
data_path = "data/"


def main():
    """
    Точка входа в программу.
    :return:
    """
    def animate(i, *data):
        """
        Функция анимации.
        :param i: Номер кадра.
        :param data: Наборы данных.
        :return: Данные для построения анимации.
        """
        part1 = [float(x) for x in data[0][i]]
        part2 = [float(x) for x in data[1][i]]

        x1 = np.linspace(t0, t1, len(part1))
        x2 = np.linspace(t0, t1, len(part2))

        plt.clf()
        plt.plot(x1, part1, color='black', label='точ')
        plt.plot(x2, part2, color='red', label='числ')
        plt.grid(True)

        plt.xlim(t0, t1)
        plt.legend(prop={'size': 20})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        if (i % 10) == 0:
            print('{0}%'.format(i * 100 / len(exact_sol_data_parts)))
        if i == 0 or i == 400 or i == 990:
            plt.savefig("pics/{1}/05_{0}_{1}.png".format(i, sys.argv[5]))
        return line,

    def init():
        """
        Функция инициализации данных для анимации.
        :return: Начальные данные для анимирования.
        """
        line.set_data([], [])
        return line,

    # Настройка папки для сохранения снимков графиков.
    if not os.path.exists("pics\\" + sys.argv[5]):
        os.makedirs("pics\\" + sys.argv[5])

    # Гранницы графика по оси абсцисс.
    t0 = int(sys.argv[1])
    t1 = int(sys.argv[2])

    # Чтение из файла и обработка данных аналитического решения.
    exact_sol_name = sys.argv[3]
    f = open(data_path + exact_sol_name, "r")
    exact_sol_data = f.read()
    f.close()
    exact_sol_data_parts = [x.split() for x in exact_sol_data.replace(",", ".").split('\n\n')]
    exact_sol_data_parts.remove([])

    # Чтение из файла и обработка данных численного решения.
    numeric_sol_name = sys.argv[4]
    f = open(data_path + numeric_sol_name, "r")
    numeric_sol_data = f.read()
    f.close()
    numeric_sol_data_parts = [x.split() for x in numeric_sol_data.replace(",", ".").split('\n\n')]

    # Настройка графического полотна и старт построения анимации.
    fig = plt.figure()
    ax = plt.axes(xlim=(t0, t1), ylim=(0.0, 1.2))
    line, = ax.plot([], [], lw=3)

    anim = FuncAnimation(fig, animate, fargs=[exact_sol_data_parts, numeric_sol_data_parts], init_func=init,
                         frames=len(numeric_sol_data_parts), interval=1, blit=True, cache_frame_data=True)
    anim.save('gifs/{0}.gif'.format(sys.argv[5]), writer='pillow', fps=60)


main()
