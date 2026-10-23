"""A compact, fixed realization of the tree/occupancy coupling."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT=Path(__file__).resolve().parents[1]

def plot_growth_example():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
                         'pdf.fonttype':42,'savefig.bbox':'tight'})
    fig,(tree,draws)=plt.subplots(1,2,figsize=(3.45,1.23),
                                 gridspec_kw={'width_ratios':[1,1.9]})
    ink,blue,orange='#30383d','#176b89','#b75c23'
    tree.set(xlim=(-.2,1.8),ylim=(-.1,2.1),aspect='equal');tree.axis('off')
    pos={1:(.8,1.8),2:(.8,1.05),3:(.2,.3),4:(1.4,.3)}
    for a,b in ((1,2),(2,3),(2,4)):
        tree.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color=ink,lw=1)
    for i,(x,y) in pos.items():
        tree.add_patch(Circle((x,y),.2,fc='#f7dfc9' if i>2 else 'white',
                             ec=orange if i>2 else ink,lw=1,zorder=2))
        tree.text(x,y,str(i),ha='center',va='center',zorder=3)
    draws.set(xlim=(0,4.2),ylim=(0,2.1));draws.axis('off')
    draws.text(0,1.85,r'Arrival $t$')
    draws.text(0,1.12,r'Draw $X_t$')
    draws.text(0,.32,r'Parent $P_t$')
    for x,t,v,p in ((2.15,2,1,1),(2.83,3,2,2),(3.51,4,6,2)):
        draws.text(x,1.85,str(t),ha='center')
        draws.text(x,1.12,str(v),ha='center',color=orange if t==4 else ink)
        draws.text(x,.32,str(p),ha='center',color=blue if t==4 else ink)
    draws.annotate('',xy=(3.51,.55),xytext=(3.51,1.05),
                   arrowprops={'arrowstyle':'->','color':blue,'lw':1})
    fig.subplots_adjust(left=.01,right=.99,bottom=.04,top=.99,wspace=.12)
    (ROOT/'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures/growth_example.pdf')
    fig.savefig(ROOT/'figures/growth_example.png',dpi=220)
    plt.close(fig)

if __name__=='__main__':
    plot_growth_example()
