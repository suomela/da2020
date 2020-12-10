#!/usr/bin/env python3

import itertools
import os
import random
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot

black = "#000000"
white = "#ffffff"
light = "#f0f0f0"
medium = "#b0b0b0"
dark = "#808080"
blue = "#0088cc"
orange = "#f26924"
blues = [ "#EAF3F9", "#D8E7F4", "#C5DBEE", "#B0CFE9", "#9CC2E3" ]

tokenradius = 0.11
noderadius = 0.08
arrowoffset = 1.5 * noderadius
edge_density = 0.5
scale = 1.2
randomness = 0.1

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
    # fig.savefig(f"figs/{name}.pdf")
    matplotlib.pyplot.close(fig)

def interpol(a, b, t0, t1, t):
    td = t1 - t0
    if td == 0:
        f = 0.5
    else:
        f = (t - t0) / td
    return a + f * (b - a)

class Node:
    def __init__(self, pos, i, j):
        self.pos = pos
        self.i = i
        self.j = j
        self.sneigh = []
        self.rneigh = []
        self.l_root = self
        self.l_parent = None
        self.l_new = None
        self.l_msgs = []


class Graph:
    def __init__(self):
        self.rng_pos = random.Random(1)
        self.rng_graph = random.Random(1)
        self.snodes = []
        self.snodemap = {}
        self.rnodes = []
        self.sedges = []
        self.sedgeset = set()
        self.redges = []
        self.redgeset = set()
        self.gen_nodes()
        self.gen_ids()
        self.gen_sedges()
        self.gen_redges()
        self.gen_ports()
        self.token = None
        self.leader_node = None

    def gen_ids(self):
        ids = list(range(1, len(self.rnodes) + 1))
        self.rng_graph.shuffle(ids)
        self.idmap = {}
        for i,v in enumerate(self.rnodes):
            v.id = ids[i]
            self.idmap[v.id] = v
        self.idorder = sorted(self.rnodes, key=lambda v: v.id)

    def gen_nodes(self):
        w = 16
        h = 9
        smargin = -1.5
        rmargin = 0.5
        radius = 3 * max(w, h)
        c1 = scale * 1j ** 0.17
        c2 = c1 * scale * 1j ** (4/6)
        for i in range(-radius, radius):
            for j in range(-radius, radius):
                offset = self.rng_pos.uniform(0, randomness) * 1j ** self.rng_pos.uniform(0, 4)
                p = i * c1 + j * c2 + offset
                v = Node(p, i, j)
                v.supp = smargin <= p.real <= w - smargin and smargin <= p.imag <= h - smargin
                v.real = rmargin <= p.real <= w - rmargin and rmargin <= p.imag <= h - rmargin
                if not v.supp:
                    continue
                v.sindex = len(self.snodes)
                self.snodes.append(v)
                self.snodemap[(i,j)] = v
                if v.real:
                    v.rindex = len(self.rnodes)
                    v.comp = v.rindex
                    self.rnodes.append(v)

    def gen_sedges(self):
        for i,j in self.snodemap.keys():
            for i2, j2 in [ (i+1,j), (i,j+1), (i+1,j-1) ]:
                if (i2,j2) in self.snodemap:
                    self.add_sedge(self.snodemap[(i,j)], self.snodemap[(i2,j2)])

    def gen_redges(self):
        n_comp = len(self.rnodes)
        while n_comp > 1:
            u,v = self.rng_graph.choice([ (u,v) for u,v in self.sedges if u.real and v.real and u.comp != v.comp ])
            self.add_redge(u, v)
            n_comp -= 1
        rest = [ (u,v) for u,v in self.sedges if u.real and v.real and (u.rindex, v.rindex) not in self.redgeset ]
        for u,v in self.rng_graph.sample(rest, round(len(rest) * edge_density)):
            self.add_redge(u, v)

    def gen_ports(self):
        for v in self.rnodes:
            self.rng_graph.shuffle(v.rneigh)

    def add_sedge(self, u, v):
        if u.sindex > v.sindex:
            u, v = v, u
        assert u.sindex < v.sindex
        assert (u.sindex, v.sindex) not in self.sedgeset
        self.sedgeset.add((u.sindex, v.sindex))
        self.sedges.append((u,v))
        u.sneigh.append(v)
        v.sneigh.append(u)

    def add_redge(self, u, v):
        assert u.real and v.real
        assert u.rindex < v.rindex
        assert (u.sindex, v.sindex) in self.sedgeset
        self.redgeset.add((u.rindex, v.rindex))
        self.redges.append((u,v))
        u.rneigh.append(v)
        v.rneigh.append(u)
        oldcomp = max(u.comp, v.comp)
        newcomp = min(u.comp, v.comp)
        for v in self.rnodes:
            if v.comp == oldcomp:
                v.comp = newcomp

    def draw_sedge(self, ax, u, v, state):
        ax.plot(
            (u.pos.real, v.pos.real), (u.pos.imag, v.pos.imag),
            color=light,
            linewidth=1,
            zorder=1,
        )

    def draw_redge(self, ax, u, v, state):
        color = black if state == "graph" else dark
        ax.plot(
            (u.pos.real, v.pos.real), (u.pos.imag, v.pos.imag),
            color=color,
            linewidth=1,
            zorder=2,
        )

    def draw_parent(self, ax, u, v, state, time):
        d = v.pos - u.pos
        d -= arrowoffset * d / abs(d)
        if state == "leader":
            if u.l_root == v.l_root:
                color = black
            else:
                color = medium
        elif state == "apsp" or state == "token":
            if u.visittime < time:
                color = medium
            else:
                color = black
        ax.arrow(
            u.pos.real, u.pos.imag, d.real, d.imag,
            color=color,
            linewidth=3,
            width=0.02,
            head_length=0.04,
            length_includes_head=True,
            zorder=3,
        )

    def draw_msg(self, ax, u, v, state, time):
        td = time - self.leader_time + 1
        assert 0 <= td <= 1
        if td == 0:
            return
        x = u.pos
        y = v.pos
        d = y - x
        start = x * td + y * (1 - td)
        end = y - arrowoffset * d / abs(d)
        delta = end - start
        color = blue
        ax.arrow(
            start.real, start.imag, delta.real, delta.imag,
            color=color,
            linewidth=4,
            width=0.02,
            head_length=0.04,
            length_includes_head=True,
            zorder=5,
        )

    def draw_node(self, ax, v, state, time):
        facecolor = light
        edgecolor = black
        if (state == "apsp" or state == "token") and v.visittime < time:
            facecolor = medium
            edgecolor = dark
        if (state == "leader" or state == "apsp" or state == "token") and v.l_root == v:
            facecolor = blue
        if (state == "wave") and self.leader_node == v:
            facecolor = blue
        ax.add_patch(matplotlib.patches.Circle(
            xy=(v.pos.real, v.pos.imag), radius=noderadius,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1,
            zorder=10,
        ))

    def draw_token(self, ax, state, time):
        smaller = [ (v,t) for v,t in self.tokentrail if t <= time ]
        larger = [ (v,t) for v,t in self.tokentrail if t > time ]
        v1, t1 = smaller[-1]
        if len(larger) == 0:
            v2, t2 = v1, time
        else:
            v2, t2 = larger[0]
        p = interpol(v1.pos, v2.pos, t1, t2, time)
        ax.add_patch(matplotlib.patches.Circle(
            xy=(p.real, p.imag), radius=tokenradius,
            facecolor=orange,
            edgecolor=orange,
            linewidth=1,
            zorder=15,
        ))

    def draw_wave(self, ax, v, state, time):
        deltas = [
            ( 0, +1),
            (+1,  0),
            (+1, -1),
            ( 0, -1),
            (-1,  0),
            (-1, +1),
        ]
        deltapairs = [ (deltas[i], deltas[(i+1) % len(deltas)]) for i in range(len(deltas)) ]
        for a in self.snodes:
            at = a.wavetimes[v.rindex]
            if at is None:
                continue
            for dd in deltapairs:
                i1, j1 = a.i+dd[0][0], a.j+dd[0][1]
                i2, j2 = a.i+dd[1][0], a.j+dd[1][1]
                b = self.snodemap.get((i1,j1), None)
                c = self.snodemap.get((i2,j2), None)
                if b is None or c is None:
                    continue
                bt = b.wavetimes[v.rindex]
                ct = c.wavetimes[v.rindex]
                if bt is None:
                    continue
                if ct is None:
                    continue
                if at <= time and (bt <= time or ct <= time):
                    continue
                if at > time and (bt > time or ct > time):
                    continue
                p1 = interpol(a.pos, b.pos, at, bt, time)
                p2 = interpol(a.pos, c.pos, at, ct, time)
                for i in range(5):
                    ax.plot(
                        (p1.real, p2.real), (p1.imag, p2.imag),
                        color=blues[i],
                        linewidth=20-4*i,
                        solid_capstyle='round',
                        zorder=-100+i,
                    )

    def draw(self, name, state, time=None):
        fig, ax = get_fig()
        # for u, v in self.sedges:
        #     self.draw_sedge(ax, u, v, state)
        for u, v in self.redges:
            self.draw_redge(ax, u, v, state)
        for v in self.rnodes:
            if v.l_parent is not None:
                self.draw_parent(ax, v, v.l_parent, state, time)
            if state == "leader":
                for u in v.l_msgs:
                    self.draw_msg(ax, v, u, state, time)
        for v in self.rnodes:
            self.draw_node(ax, v, state, time)
        if time is not None:
            if state == "apsp" or state == "token":
                self.draw_token(ax, state, time)
            if state != "token":
                for v in self.rnodes:
                    if v.wavetimes[v.rindex] is not None and v.wavetimes[v.rindex] <= time:
                        self.draw_wave(ax, v, state, time)
        save_fig(fig, name)

    def leader_init(self):
        self.leader_time = 0

    def run_leader(self):
        self.leader_init()
        while self.leader_step():
            pass

    def leader_step(self):
        self.leader_step_1()
        return self.leader_step_2()

    def leader_step_1(self):
        self.leader_time += 1
        for v in self.rnodes:
            v.l_new = None
            v.l_msgs = []
            for u in v.rneigh:
                if u.l_root != v.l_root:
                    self.leader_msg(u, v)

    def leader_step_2(self):
        progress = False
        self.leaders = []
        for v in self.rnodes:
            if v.l_new is not None:
                v.l_parent = v.l_new[0]
                v.l_root = v.l_new[1]
                progress = True
                v.l_new = None
            v.l_msgs = []
            if v.l_root == v:
                self.leaders.append(v)
        if len(self.leaders) == 1:
            self.leader_node = self.leaders[0]
        return progress

    def leader_msg(self, u, v):
        if u.l_root.id >= v.l_root.id:
            return
        assert u.l_root.id < v.l_root.id
        v.l_msgs.append(u)
        if v.l_new is None or v.l_new[1].id > u.l_root.id:
            v.l_new = (u, u.l_root)

    def init_waves(self):
        self.apsp_time = 0
        for v in self.snodes:
            v.wavetimes = [ None for v in self.rnodes ]
        for v in self.rnodes:
            v.to_visit = [ u for u in v.rneigh if u.l_parent == v ]
            v.to_visit.reverse()
            v.visited = False

    def run_apsp(self):
        self.init_waves()
        self.tokentrail = []
        v = self.leader_node
        while v is not None:
            self.tokentrail.append((v, self.apsp_time))
            if not v.visited:
                self.apsp_wave_start(v)
            self.apsp_time += 1
            self.tokentrail.append((v, self.apsp_time))
            self.apsp_time += 1
            if len(v.to_visit):
                v = v.to_visit.pop()
            else:
                v = v.l_parent
        return max(max(u.wavetimes) for u in self.snodes)

    def apsp_wave_start(self, v):
        assert not v.visited
        v.visited = True
        v.visittime = self.apsp_time
        d = self.apsp_time
        queue = [ v ]
        while len(queue) > 0:
            new = []
            for u in queue:
                if u.wavetimes[v.rindex] is None:
                    u.wavetimes[v.rindex] = d
                    new += u.rneigh + [ x for x in u.sneigh if not x.real ]
            queue = new
            d += 1

    def run_wave(self, v):
        self.leader_node = v
        self.init_waves()
        self.apsp_wave_start(v)
        return max(u.wavetimes[v.rindex] for u in self.snodes)

    def run_all(self):
        self.init_waves()
        for v in self.rnodes:
            self.apsp_wave_start(v)
        return max(max(u.wavetimes) for u in self.snodes)

    def run_leader_wave(self):
        tt = self.run_wave(self.idmap[1])
        return tt


def do_graph():
    graph = Graph()
    graph.draw("graph", state="graph")

def do_wave(part, total):
    graph = Graph()
    for node in [1, 2, 3]:
        tt = graph.run_wave(graph.idmap[node])
        print(tt)
        for t in range(tt+1):
            ff = 8
            for f in range(ff):
                if f % total == part:
                    graph.draw(f"wave-{node}-{t:03d}-{f}", state="wave", time=t + f/ff)

def do_all(part, total):
    graph = Graph()
    tt = graph.run_all()
    print(tt)
    for t in range(4):
        ff = 64
        for f in range(ff):
            if f % total == part:
                graph.draw(f"all-{t:03d}-{f:02d}", state="wave", time=t + f/ff)

def do_leader(part, total):
    graph = Graph()
    graph.run_leader_wave()
    graph.leader_init()
    while True:
        t = graph.leader_time
        ff = 32
        for f in range(ff):
            if f % total == part:
                graph.draw(f"leader-{t:03d}-{f:02d}", state="leader", time=t + f/ff)
            if f == 0:
                graph.leader_step_1()
        if not graph.leader_step_2():
            return

def do_apsp(part, total):
    graph = Graph()
    graph.run_leader()
    tt = graph.run_apsp()
    print(tt)
    for t in range(tt+1):
        ff = 8
        for f in range(ff):
            if f % total == part:
                graph.draw(f"apsp-{t:03d}-{f}", state="apsp", time=t + f/ff)

def do_slow(part, total):
    graph = Graph()
    graph.run_leader()
    tt = graph.run_apsp()
    print(tt)
    for t in range(tt+1):
        ff = 16
        for f in range(ff):
            if f % total == part:
                graph.draw(f"slow-{t:03d}-{f:02d}", state="apsp", time=t + f/ff)

def do_token(part, total):
    graph = Graph()
    graph.run_leader()
    tt = graph.run_apsp()
    print(tt)
    for t in range(tt+1):
        ff = 8
        for f in range(ff):
            if f % total == part:
                graph.draw(f"token-{t:03d}-{f}", state="token", time=t + f/ff)

def main():
    os.makedirs("figs", exist_ok=True)
    args = sys.argv[1:]
    what = args[0]
    if len(args) > 1:
        part = int(args[1])
        total = int(args[2])
    else:
        part = 0
        total = 1
    if what == "graph":
        do_graph()
    elif what == "wave":
        do_wave(part, total)
    elif what == "all":
        do_all(part, total)
    elif what == "leader":
        do_leader(part, total)
    elif what == "apsp":
        do_apsp(part, total)
    elif what == "slow":
        do_slow(part, total)
    elif what == "token":
        do_token(part, total)
    else:
        assert False, what

main()
