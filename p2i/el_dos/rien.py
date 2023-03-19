import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()
rect = plt.Rectangle((0, 0), 0.5, 0.5, color='r')
ax.add_patch(rect)
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])

def update(frame):
    x, y = rect.get_xy()
    rect.set_xy([x+0.1, y+0.1])
    return rect,

ani = animation.FuncAnimation(fig, update, frames=100, interval=50, blit=True)
plt.show()
