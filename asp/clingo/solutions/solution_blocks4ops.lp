pred(1).
p_arity(1,0). % empty
pred(2).
p_arity(2,1). % holding(?x)
pred(3).
p_arity(3,1). % clear(?x)
pred(4).
p_arity(4,2). % ontable(?x,?x) / on(?x,?y)

% PUTDOWN(?x): 1p + 3e
action(1).
a_arity(1,1).
prec(1,(1,(0,0)),0) :- opt_req_prec > 0. % -empty
prec(1,(2,(1,1)),1). % holding(?x)
prec(1,(4,(1,1)),0) :- opt_req_prec > 0. % -ontable(?x,?x)
eff(1,(1,(0,0)),1).  % empty
eff(1,(2,(1,1)),0).  % -holding(?x)
eff(1,(4,(1,1)),1).  % ontable(?x,?x)

% PICK(?x): 3p + 3e
action(2).
a_arity(2,1).
prec(2,(1,(0,0)),1). % empty
prec(2,(2,(1,1)),0) :- opt_req_prec > 0. % -holding(?x)
prec(2,(3,(1,1)),1). % clear(?x)
prec(2,(4,(1,1)),1). % ontable(?x,?x)
eff(2,(1,(0,0)),0).  % -empty
eff(2,(2,(1,1)),1).  % holding(?x)
eff(2,(4,(1,1)),0).  % -ontable(?x,?x)

% UNSTACK(?x,?y): 3p + 4e
action(3).
a_arity(3,2).
prec(3,(1,(0,0)),1). % empty
prec(3,(2,(1,1)),0) :- opt_req_prec > 0. % -holding(?x)
prec(3,(3,(1,1)),1). % clear(?x)
prec(3,(4,(1,2)),1). % on(?x,?y)
eff(3,(1,(0,0)),0).  % -empty
eff(3,(2,(1,1)),1).  % holding(?x)
eff(3,(3,(2,2)),1).  % clear(?y)
eff(3,(4,(1,2)),0).  % -on(?x,?y)

% STACK(?x,?y): 2p + 4e
action(4).
a_arity(4,2).
prec(4,(1,(0,0)),0) :- opt_req_prec > 0. % -empty
prec(4,(2,(1,1)),1). % holding(?x)
prec(4,(3,(2,2)),1). % clear(?y)
prec(4,(4,(1,2)),0) :- opt_req_prec > 0. % -on(?x,?y)
eff(4,(1,(0,0)),1).  % empty
eff(4,(2,(1,1)),0).  % -holding(?x)
eff(4,(3,(2,2)),0).  % -clear(?y)
eff(4,(4,(1,2)),1).  % on(?x,?y)

