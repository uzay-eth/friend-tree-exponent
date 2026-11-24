#!/usr/bin/env python3
"""Physical-arrow controls for the sequential fixed-shape review."""
from fractions import Fraction as Q


class Tree:
    def __init__(self):
        # Focal o=0 starts with a live exterior core p=1 and real leaf x=2.
        self.adj = [{1,2}, {0,3}, {0}, {1}]
        self.up = [1,None,0,1]

    def send(self, sender, receiver):
        assert receiver in self.adj[sender]
        rate = Q(1,len(self.adj[sender]))
        x = len(self.adj)
        self.adj[receiver].add(x)
        self.adj.append({receiver})
        self.up.append(receiver)
        return x,rate

    def stock(self, leaf, core, n):
        assert self.adj[leaf]=={core}
        for _ in range(n):
            _, rate = self.send(leaf,core)
            assert rate==1

    def component(self, x, N):
        if len(self.adj[x])>N:
            return set()
        out, todo = {x},[x]
        while todo:
            i=todo.pop()
            for j in self.adj[i]-out:
                if len(self.adj[j])<=N:
                    todo.append(j)
                    out.add(j)
        return out


def nested_cuts():
    # Core path o--x--z--w--q. Grow hubs o,w, leaving x,z bare.
    # If w is cut, q is irrelevant even when q is itself a second hub.
    for N in (2,8,31,128):
        T=Tree()
        z,_=T.send(0,2)
        w,_=T.send(2,z)
        q,_=T.send(z,w)
        qleaf,_=T.send(w,q)
        wleaf,_=T.send(q,w)
        # qleaf remains a leaf attached to q; wleaf remains at w.
        T.stock(3,1,N+3)
        oleaf,_=T.send(1,0)
        T.stock(oleaf,0,N+3)
        T.stock(wleaf,w,N+3)
        T.stock(qleaf,q,N+3)
        assert T.component(2,N)=={2,z}
        assert (len(T.adj[2]), len(T.adj[z]))==(2,2)
        # A delayed send makes at least one endpoint gain. Both bounded
        # degrees cannot be frozen forever: this edge has rate 1/2 each way.
        leaf,rate=T.send(z,2)
        assert rate==Q(1,2) and len(T.adj[2])==3
        assert T.component(2,N)==(set() if N==2 else {2,z,leaf})


def delayed_and_multiple():
    for N in (8,31,128):
        T=Tree()
        T.stock(2,0,N+4)
        # x=2 and x'=4 each produce a core child. Promotions can occur
        # before the first sends from their own eventual hub neighbours.
        tops=(2,4)
        for x in tops:
            w,_=T.send(0,x)
            tip,_=T.send(x,w)
            T.stock(tip,w,N+4)
            assert T.component(x,N)=={x}
        assert all(T.component(x,N)=={x} for x in tops)
        leaf, rate=T.send(0,2)
        assert rate>0 and T.component(2,N)=={2,leaf}
        # A late send can be repaired: promote the newly added leaf and
        # grow that neighbour too. The final degree is 3, not 2; the minimum
        # must include this extra neighbour even though it was absent earlier.
        tip,_=T.send(2,leaf)
        assert T.component(2,N)=={2,leaf,tip}
        T.stock(tip,leaf,N+4)
        assert T.component(2,N)=={2} and len(T.adj[2])==3


def exact_controls():
    # Two positive masses in distinct local clocks. From parent-clock origin
    # the late neighbour's scale is e^(-tau) times its own-clock mass.
    # Use a delay log(4), parent mass 2, and late neighbour own mass 1.
    m=min(Q(2),Q(1,4))
    assert m==Q(1,4) and Q(2)/m==8
    # At mu=1/2 the integrated simultaneous gate equals sqrt(m)=1/2;
    # using only the parent gate or the wrong origin gives another answer.
    assert Q(1,2)**2==m and m!=min(Q(2),Q(1))
    barrier=lambda j:min(Q(1),Q(2,j+1))
    for j in range(1,65):
        for extra in (0,1,9):
            d=j+1+extra
            drift=j*(barrier(j+1)-barrier(j))+Q(j,d)*(barrier(j-1)-barrier(j))+Q(1,d)*(1-barrier(j))
            assert drift<=0
    # Uniform positivity after the setup: two first-send failure bounds,
    # and two zero-mass bounds 4/j+2/j via Gamma(j) and mean phantom loss 1.
    j=32
    assert 2*barrier(j)+2*Q(6,j)<1
    assert 1/(Q(2)-Q(1,2))==Q(2,3)


if __name__=='__main__':
    nested_cuts()
    delayed_and_multiple()
    exact_controls()
    print('PASS: 4 nested-cut old pairs; 3 multiple/delayed/repaired singleton histories')
    print('PASS: common-clock mass minimum, 192 barrier drifts, positive union-bound budget, virgin holding mean')
