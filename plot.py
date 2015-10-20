import matplotlib.pyplot as plt

def plot_signature(array, ylim=0.2, ax=None):
    color = ((0.196,0.714,0.863),)*16 + ((0.102,0.098,0.098),)*16 + ((0.816,0.180,0.192),)*16 + ((0.777,0.773,0.757),)*16 + ((0.604,0.777,0.408),)*16 + ((0.902,0.765,0.737),)*16
    color = list(color)

    width = max(array.shape)
    x = np.arange(width)
    if ax == None:
        f,ax = plt.subplots(1)
    bars = ax.bar(x, array)
    for h in range(len(x)):
        bars[h].set_color(color[h])
    plt.ylim(0, ylim)
    plt.xlim(0, width)
