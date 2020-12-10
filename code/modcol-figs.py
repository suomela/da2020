#!/usr/bin/env python3

import itertools
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot

black = "#000000"
light = "#f0f0f0"
dark = "#808080"
blue = "#0088cc"
orange = "#f26924"

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = "Bernino Sans"
matplotlib.rcParams['font.weight'] = "light"

def get_fig():
    fig = matplotlib.pyplot.figure(figsize=(960/72, 540/72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim((0, 16))
    ax.set_ylim((0, 9))
    return fig, ax

def save_fig(fig, name, pdf=False):
    print(name)
    fig.savefig(f"figs/{name}.png", dpi=4*72)
    if pdf:
        fig.savefig(f"figs/{name}.pdf")
    matplotlib.pyplot.close(fig)

def pos(q, b):
    return 1j ** (1 - 4*b/q)

def angle(q, b):
    return 90 - 360*b/q

def draw_clock(ax, q, x, y, a, b, f, scale=1, highlight=False):
    r = 0.43 * scale
    c = x + y*1j
    ax.add_patch(matplotlib.patches.Circle(
        xy=(c.real, c.imag), radius=r,
        facecolor=light,
        edgecolor=blue if highlight else dark,
        linewidth=5 if highlight else 1,
    ))
    for i in range(q):
        p = pos(q, i)
        p0 = 1.0 * p * r + c
        p1 = 0.9 * p * r + c
        ax.plot(
            (p0.real, p1.real), (p0.imag, p1.imag),
            color=dark,
            linewidth=1,
        )

    ax.add_patch(matplotlib.patches.Arc(
        xy=(c.real, c.imag), width=0.6*r, height=0.6*r,
        theta1=angle(q, b+a), theta2=angle(q, b),
        color=black,
        linewidth=3*scale,
    ))

    p = pos(q, b+f*a)
    p0 = c
    p1 = 0.75 * p * r + c
    ax.plot(
        (p0.real, p1.real), (p0.imag, p1.imag),
        color=orange,
        linewidth=4*scale,
        solid_capstyle="round",
    )

def draw_grid(output):
    q = 7
    ff = 8
    for t in range(q):
        for f in range(ff+1):
            fig, ax = get_fig()
            for b in range(q):
                for a in range(q):
                    draw_clock(ax, q, 3+a, 0.5+q-b, a, (b + t*a) % q, f/ff)
                if output:
                    for a in range(1):
                        draw_clock(ax, q, 13+a, 0.5+q-b, a, (b + t*a) % q, f/ff)
            name = "clocks" if output else "input"
            save_fig(fig, f"{name}-{t}-{f}")

def draw_pairs():
    q = 7
    ff = 8
    for clocks in [
        [ (2,3), (2,5) ],
        [ (1,2), (3,4) ],
        [ (5,3), (4,1) ],
        [ (2,3), (0,2) ],
    ]:
        name = "-".join(f"{a}{b}" for a,b in clocks)
        for t in range(q):
            for f in range(ff+1):
                fig, ax = get_fig()
                for i,clock in enumerate(clocks):
                    a,b = clock
                    draw_clock(ax, q, 5.5+5*i, 4.5, a, (b + t*a) % q, f/ff, scale=3)
                if f == 0:
                    values = [ (b + t*a) % q for a,b in clocks ]
                elif f == ff:
                    values = [ (b + (t+1)*a) % q for a,b in clocks ]
                else:
                    values = None
                if values is not None and values[0] == values[1]:
                    ax.plot(
                        [7.0, 9.0], [4.5, 4.5],
                        color=blue,
                        linewidth=12,
                        solid_capstyle="round",
                    )
                save_fig(fig, f"pairs-{name}-{t}-{f}")

def draw_sim():
    q = 7
    ff = 8
    for init in [
        [(2, 3), (3, 2), (3, 6), (4, 6)],
        [(1, 1), (1, 4), (2, 4), (4, 1)],
        [(1, 2), (3, 0), (3, 2), (4, 0)],
    ]:
        print(init)
        name = "-".join(f"{a}{b}" for a,b in init)
        clocks = [ (a,b) for a,b in init ]
        t = 1
        while True:
            presteps = [ "a", "b", "c" ]
            for step in presteps + list(range(ff+1)):
                fig, ax = get_fig()
                f = 0 if step in presteps else step
                stepname = step if step in presteps else f"x{step}"

                conflicts = [ False for c in clocks ]
                highlights = [ False for c in clocks ]

                for i,j in itertools.combinations(range(len(clocks)), r=2):
                    conflict = False
                    if clocks[i][1] == clocks[j][1]:
                        conflict = True
                        conflicts[i] = True
                        conflicts[j] = True

                    x1,y1 = i % 2, i // 2
                    x2,y2 = j % 2, j // 2
                    def mix(s,t):
                        r = 0.3
                        return (1-r)*s + r*t, (1-r)*t + r*s
                    x1,x2 = mix(x1,x2)
                    y1,y2 = mix(y1,y2)
                    ax.plot(
                        [5.5 + 5*x1, 5.5 + 5*x2],
                        [2.0 + 5*y1, 2.0 + 5*y2],
                        color=blue if conflict and step == "b" else light,
                        zorder=1 if conflict and step == "b" else 0,
                        linewidth=12,
                        solid_capstyle="round",
                    )

                if step == "c":
                    for i in range(len(clocks)):
                        if not conflicts[i] and clocks[i][0] != 0:
                            clocks[i] = (0, clocks[i][1])
                            highlights[i] = True

                for i,clock in enumerate(clocks):
                    x,y = i % 2, i // 2
                    a,b = clock
                    draw_clock(ax, q, 5.5+5*x, 2.0+5*y, a, b, f/ff, scale=3, highlight=highlights[i])

                save_fig(fig, f"sim-{name}-{t}-{stepname}")

            t += 1
            clocks = [ (a, (a+b) % q) for a,b in clocks ]
            print(clocks)
            if all(a == 0 for a,b in clocks):
                break
        print()

def draw_one():
    q = 7

    fig, ax = get_fig()
    a, b = 2, 3
    draw_clock(ax, q, 13.5, 6.5, a, b, 0, scale=3)
    b = (b + a) % q
    draw_clock(ax, q, 13.5, 2.0, a, b, 0, scale=3)
    save_fig(fig, f"one-1", pdf=True)

    fig, ax = get_fig()
    a, b = 2, 3
    draw_clock(ax, q, 13.5, 6.5, a, b, 0, scale=3)
    a = 0
    draw_clock(ax, q, 13.5, 2.0, a, b, 0, scale=3)
    save_fig(fig, f"one-2", pdf=True)


def main():
    os.makedirs("figs", exist_ok=True)
    draw_grid(False)
    # draw_grid(True)
    # draw_pairs()
    # draw_sim()
    # draw_one()

main()
