% (Break symmetries) Order predicates by arity
% CHECK: see constraints_javier.lp with psym=1 (redundant)
:- p_arity(P,X), p_arity(P+1,Y), X > Y, P=1..num_predicates-num_static-1.
:- not p_used(P), p_used(P+1), P=1..num_predicates-num_static-1.
:- p_arity(P,X), p_arity(P+1,Y), X > Y, P=num_predicates-num_static+1..num_predicates-1.
:- not p_used(P), p_used(P+1), P=num_predicates-num_static+1..num_predicates-1.

%*
% Simple small problems (e.g. packages without a package)
% may have fluents with equal denotation across all states

% Non-static predicates must differ at some state (if non-partial problem)
:- not partial(I,File) : filename(File); pred(P), not p_static(P), p_arity(P,0), root(I,R), val(I,(P,( 0, 0)),R,X), val(I,(P,(0,0)),S,X) : node(I,S), S != R.
:- not partial(I,File) : filename(File); pred(P), not p_static(P), p_arity(P,1), root(I,R), val(I,(P,(O1,O1)),S,X) : node(I,S), object(I,O1),               val(I,(P,(O1,O1)),R,X), S != R.
:- not partial(I,File) : filename(File); pred(P), not p_static(P), p_arity(P,2), root(I,R), val(I,(P,(O1,O2)),S,X) : node(I,S), object(I,O1), object(I,O2), val(I,(P,(O1,O2)),R,X), S != R.
*%

