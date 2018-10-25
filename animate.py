from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np


def animate_imstack(imstack, figsize=(10, 10)):
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(xlim=(0, imstack.shape[2]), ylim=(imstack.shape[1], 0))
    image = ax.imshow(imstack[0])

    def init():
        image.set_array(imstack[0])
        return image,

    def animate(i):
        image.set_array(imstack[i])
        return image,

    return animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=imstack.shape[0], interval=30)
