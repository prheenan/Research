# see:
# stackoverflow.com/questions/41469207/how-to-draw-rectangle-outside-of-the-plot-frame-in-matplotlib/46656255#46656255

def run():
    import random
    import matplotlib.pyplot as pyplot

    x = random.sample(range(50), 50)
    y= random.sample(range(50), 50)
    
    fig = pyplot.figure()
    ax = pyplot.subplot(111)
    ax.scatter(x,y,label='a')
    ax.set_aspect('equal')
    ax.set_xlim(0,60)
    ax.set_ylim(0,60)
    ax.plot([0,60], [0, 60], color='k', linestyle='-', linewidth=1.25)
    
    ax.add_patch(pyplot.Rectangle((0,60),60, 10,facecolor='silver',
                                  clip_on=False,linewidth = 0))
    TITLE = ax.text(26,61, r'$\mathregular{Title}$',fontsize = 14,zorder = 5,
                    color = 'k')
    pyplot.show()

if __name__ == "__main__":
    run()
