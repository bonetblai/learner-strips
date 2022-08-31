% Predicates
pred(1).      % atX(?p)
p_arity(1,1).
pred(2).      % atY(?p)
p_arity(2,1).
pred(3).      % connX(?x1,?x2)
p_arity(3,2).
p_static(3).
pred(4).      % connY(?y1,?y2)
p_arity(4,2).
p_static(4).

% LEFT(?p,?q): 2p + 2e
action(1).
a_arity(1,2).
prec(1,(1,(1,1)),1). % atX(?p)
prec(1,(1,(2,2)),0) :- opt_req_prec > 0. % -atX(?q)
prec(1,(3,(1,2)),1). % connX(?p,?q)
eff(1,(1,(1,1)),0).  % -atX(?p)
eff(1,(1,(2,2)),1).  % atX(?q)
unequal(1,(1,2)).

% RIGHT(?p,?q): 2p + 2e
action(4).
a_arity(4,2).
prec(4,(1,(1,1)),1). % atX(?p)
prec(4,(1,(2,2)),0) :- opt_req_prec > 0. % -atX(?q)
prec(4,(3,(2,1)),1). % connX(?q,?p)
eff(4,(1,(1,1)),0).  % -atX(?p)
eff(4,(1,(2,2)),1).  % atX(?q)
unequal(2,(1,2)).

% DOWN(?p,?q): 2p + 2e
action(2).
a_arity(2,2).
prec(2,(2,(1,1)),1). % atY(?p)
prec(2,(2,(2,2)),0) :- opt_req_prec > 0. % -atY(?q)
prec(2,(4,(1,2)),1). % connY(?p,?q)
eff(2,(2,(1,1)),0).  % -atY(?p)
eff(2,(2,(2,2)),1).  % atY(?q)
unequal(3,(1,2)).

% UP(?p,?q): 2p + 2e
action(3).
a_arity(3,2).
prec(3,(2,(1,1)),1). % atY(?p)
prec(3,(2,(2,2)),0) :- opt_req_prec > 0. % -atY(?q)
prec(3,(4,(2,1)),1). % connY(?q,?p)
eff(3,(2,(1,1)),0).  % -atY(?p)
eff(3,(2,(2,2)),1).  % atY(?q)
unequal(4,(1,2)).

