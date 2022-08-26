% Predicates
pred(1).
p_arity(1,1). % at(?x)
pred(2).
p_arity(2,1). % empty(?g)
pred(3).
p_arity(3,2). % ball(?b,?r)
pred(4).
p_arity(4,2). % holding(?g,?b)
pred(5).
p_arity(5,2). % connected(?x,?y)
p_static(5).

% MOVE(?x,?y): 2p + 2e
action(1).
a_arity(1,2).
prec(1,(1,(1,1)),1). % at(?y)
%%prec(1,(1,(2,2)),0). % -at(?x)
prec(1,(5,(1,2)),1). % connected(?y,?x)
eff(1,(1,(2,2)),1).  % at(?x)
eff(1,(1,(1,1)),0).  % -at(?y)
unequal(1,(1,2)).

% DROP(?g,?r,?b): 2p + 3e
action(2).
a_arity(2,3).
prec(2,(1,(2,2)),1). % at(?r)
%%prec(2,(2,(1,1)),0). % -empty(?g)
%%prec(2,(3,(3,2)),0). % -ball(?b,?r)
prec(2,(4,(1,3)),1). % holding(?g,?b)
eff(2,(2,(1,1)),1).  % empty(?g)
eff(2,(3,(3,2)),1).  % ball(?b,?r)
eff(2,(4,(1,3)),0).  % -holding(?g,?b)

% PICK(?r,?g,?b): 3p + 3e
action(3).
a_arity(3,3).
prec(3,(1,(2,2)),1). % at(?r)
prec(3,(2,(1,1)),1). % empty(?g)
prec(3,(3,(3,2)),1). % ball(?b,?r)
%%prec(3,(4,(1,3)),0). % -holding(?g,?b)
eff(3,(2,(1,1)),0).  % -empty(?g)
eff(3,(3,(3,2)),0).  % -ball(?b,?r)
eff(3,(4,(1,3)),1).  % holding(?g,?b)

