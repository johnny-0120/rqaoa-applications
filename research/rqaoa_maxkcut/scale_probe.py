import sys, time
import numpy as np
from rqk import maxkcut_problem, rqaoa, make_transfer
from baselines import wifi_instance, bravyi_graph, cocolour, greedy, tabu, descent
k=3
for kind, n in [("wifi",30),("wifi",96),("bravyi",60)]:
    rng=np.random.default_rng(1)
    if kind=="wifi": E,W,_=wifi_instance(n,rng)
    else: E,W=bravyi_graph(n,4,rng)
    pr=maxkcut_problem(n,E,W,k)
    t=time.perf_counter(); rule,info=make_transfer(pr); x,c,tr=rqaoa(pr,nc=6,transfer=rule); tt=time.perf_counter()-t
    xd=descent(n,E,W,k,x,rng)
    tb=[tabu(n,E,W,k,time_limit=tt,rng=np.random.default_rng(s)) for s in range(3)]
    tb01=[tabu(n,E,W,k,time_limit=0.01*tt,rng=np.random.default_rng(s)) for s in range(3)]
    print(kind,n,len(E),"RQtransfer",round(c,3),"+descent",round(cocolour(xd,E,W),3),"time",round(tt,1),
          "greedy",round(cocolour(greedy(n,E,W,k),E,W),3),
          "tabu=T",[round(b[1],3) for b in tb],"tabu 1%",[round(b[1],3) for b in tb01], flush=True)
