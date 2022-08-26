%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

% Whether effect must modify precs (opt_req_prec=1) or whether
% this a preferred model (opt_req_prec=2); opt_req_prec=2 requires
% --heuristic=Domain
%#const opt_req_prec=1.

% Order predicates by arity and use (opt_psym>=1), and enforce
% the use of predicates in order (opt_psym=2).
%#const opt_psym=2.

% Symmetry breaking on objects controlled with opt_osym={0..3}.
% opt_osym=0 : no symmetry breaking on objects
% opt_osym=1 : complex formulation (see comments below)
% opt_osym=2 : simpler formulation (see comments below)
% opt_osym=3 : separating the levels (see comments below)
%#const opt_osym=2.

% Symmetry breaking on action variables controlled with opt_asym={0..3}.
% opt_asym=0 : no symmetry breaking on action variables
% opt_asym=1 : complex formulation (see comments below)
% opt_asym=2 : simpler formulation (see comments below)
% opt_asym=3 : separating the levels (see comments below)
%#const opt_asym=2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Every effect changes some precondition

:- eff(A,M,V), not prec(A,M,1-V),           opt_req_prec = 1, opt_synthesis = 1.

% For opt_req_prec=2, some heuristics may apply


%*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on predicates

p_next(P,P+1) :- pred(P), pred(P+1), not p_static(P), not p_static(P+1), opt_psym >= 1, opt_synthesis = 1.
p_next(P,P+1) :- pred(P), pred(P+1),     p_static(P),     p_static(P+1), opt_psym >= 1, opt_synthesis = 1.

% Order predicates by arity and by use
:- p_next(P,PP), p_arity(P,N), p_arity(PP,M), N > M,                     opt_psym >= 1, opt_synthesis = 1.
:- p_next(P,PP), not p_used(P), p_used(PP),                              opt_psym >= 1, opt_synthesis = 1.

a_pred(A,P) :-  eff(A,(P,T),V),                                          opt_psym >= 1, opt_synthesis = 1.
a_pred(A,P) :- prec(A,(P,T),V),                                          opt_psym >= 1, opt_synthesis = 1.

% Action A cannot use predicate PP if previous actions and
% A itself haven't used previous predicate P of same arity
:- p_next(P,PP), p_arity(P,N), p_arity(PP,N), a_pred(A,PP),
   #false : a_pred(AA,P), AA <= A;
   opt_synthesis = 1,
   opt_psym = 2.
*%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on objects

% Option opt_osym=1

% Consider consecutive objects O and O+1. A ground atom of
% the form P(O,OO) true at root node R yields a value of
% B**(2*(M-P) + 1) while an atom of the form P(OO,O) yields
% a value of B**(2*(M-P)), where B=1+num_objects. Similarly,
% for object O+1.
% Then, it cannot hold that the aggregated (sum) value yielded
% by object O is less than the aggregated value yielded by O+1.

:- object(I,O), object(I,O+1), M = num_predicates, root(I,R),
   1 #sum{ -((1+num_objects)**((M-P)*2 + 1)),P,OO : pred(P), val(I,(P,  (O,OO)),R,1);
           -((1+num_objects)**((M-P)*2 + 1)),P,OO : pred(P), val(I,(P,  (O,OO)),  1);
           -((1+num_objects)**((M-P)*2    )),P,OO : pred(P), val(I,(P,  (OO,O)),R,1);
           -((1+num_objects)**((M-P)*2    )),P,OO : pred(P), val(I,(P,  (OO,O)),  1);
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : pred(P), val(I,(P,(O+1,OO)),R,1);
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : pred(P), val(I,(P,(O+1,OO)),  1);
             (1+num_objects)**((M-P)*2     ),P,OO : pred(P), val(I,(P,(OO,O+1)),R,1);
             (1+num_objects)**((M-P)*2     ),P,OO : pred(P), val(I,(P,(OO,O+1)),  1) },
   opt_synthesis = 1,
   opt_osym = 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simpler symmetry breaking on objects

% Option opt_osym=2

% For two consecutive objects O and O+1, number of true atoms
% at root node R that involve O should be bigger than number
% of true atoms at R that involve O+1.

:- object(I,O), object(I,O+1), root(I,R),
   1 #sum{ -1,P,OO,1 : val(I,(P,  (O,OO)),R,1);
           -1,P,OO,2 : val(I,(P,  (OO,O)),R,1);
           -1,P,OO,1 : val(I,(P,  (O,OO)),  1);
           -1,P,OO,2 : val(I,(P,  (OO,O)),  1);
            1,P,OO,1 : val(I,(P,(O+1,OO)),R,1);
            1,P,OO,2 : val(I,(P,(OO,O+1)),R,1);
            1,P,OO,1 : val(I,(P,(O+1,OO)),  1);
            1,P,OO,2 : val(I,(P,(OO,O+1)),  1) },
   opt_synthesis = 1,
   opt_osym = 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on objects separating the levels

% Option opt_osym=3

% O occurs less or equal times than O+1 as first argument of P
lesseq(O,O+1,P,1) :- object(I,O), object(I,O+1), pred(P), root(I,R),
   0 #sum{ -1,OO : val(I,(P,  (O,OO)),  1);
           -1,OO : val(I,(P,  (O,OO)),R,1);
            1,OO : val(I,(P,(O+1,OO)),  1);
            1,OO : val(I,(P,(O+1,OO)),R,1) },
   opt_synthesis = 1,
   opt_osym = 3.

% O occurs less or equal times than O+1 as second argument of P
lesseq(O,O+1,P,2) :- object(I,O), object(I,O+1), pred(P), root(I,R),
   0 #sum{ -1,OO : val(I,(P,  (OO,O)),  1);
           -1,OO : val(I,(P,  (OO,O)),R,1);
            1,OO : val(I,(P,(OO,O+1)),  1);
            1,OO : val(I,(P,(OO,O+1)),R,1) },
   opt_synthesis = 1,
   opt_osym = 3.

% O occurs less times than O+1 as first argument of P
less(O,O+1,P,1) :- object(I,O), object(I,O+1), pred(P), root(I,R),
   1 #sum{ -1,OO : val(I,(P,  (O,OO)),  1);
           -1,OO : val(I,(P,  (O,OO)),R,1);
            1,OO : val(I,(P,(O+1,OO)),  1);
            1,OO : val(I,(P,(O+1,OO)),R,1) },
   opt_synthesis = 1,
   opt_osym = 3.

% O occurs less times than O+1 as second argument of P
less(O,O+1,P,2) :- object(I,O), object(I,O+1), pred(P), root(I,R),
   1 #sum{ -1,OO : val(I,(P,  (OO,O)),  1);
           -1,OO : val(I,(P,  (OO,O)),R,1);
            1,OO : val(I,(P,(OO,O+1)),  1);
            1,OO : val(I,(P,(OO,O+1)),R,1) },
   opt_synthesis = 1,
   opt_osym = 3.

% It cannot be that object O occurs less than O+1 as argument I of P,
% and it occurs less than or equal on the previous pairs (PP,II)
:- object(I,O), object(I,O+1),
   less(O,O+1,P,I),
   lesseq(O,O+1,PP,II) : pred(PP), II = 1..2, (PP,II) < (P,I);
   opt_synthesis = 1,
   opt_osym = 3.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on action variables

% Directive to avoid warnings when sreversible/1 isn't defined
#defined sreversible/1.

% Option opt_asym=1

:- V = 1..2, Vs = 3, action(A), M = num_predicates,
   V = 1 : sreversible(A);
   1 #sum{ -((1+Vs)**((M-P)*2 + 1)),I,P,VV : a_atom(I,A,(P,( V,  VV)));
           -((1+Vs)**((M-P)*2    )),I,P,VV : a_atom(I,A,(P,( VV,  V)));
             (1+Vs)**((M-P)*2 + 1 ),I,P,VV : a_atom(I,A,(P,(V+1, VV)));
             (1+Vs)**((M-P)*2     ),I,P,VV : a_atom(I,A,(P,( VV,V+1))) },
   opt_synthesis = 1,
   opt_asym = 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simpler symmetry breaking on action variables

% Option opt_asym=2

:- V = 1..2, action(A), not act_inverse(A),
   V = 1 : sreversible(A);
   1 #sum{ -1,I,P,VV,1 : a_atom(I,A,(P,(  V, VV)));
           -1,I,P,VV,2 : a_atom(I,A,(P,( VV,  V)));
            1,I,P,VV,1 : a_atom(I,A,(P,(V+1, VV)));
            1,I,P,VV,2 : a_atom(I,A,(P,( VV,V+1))) },
   opt_synthesis = 1,
   opt_asym = 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on action variables separating the levels

% Option opt_asym=3

% V occurs less or equal times than V+1 on the first argument of P
% (in a condition of type I of action A)
lesseq(V,V+1,A,I,P,1) :- V = 1..2, action(A), I = 1..2, pred(P),
   V = 1 : sreversible(A);
   0 #sum{ -1,VV : a_atom(I,A,(P,(  V,VV)));
            1,VV : a_atom(I,A,(P,(V+1,VV))) },
   opt_synthesis = 1,
   opt_asym = 3.

% V occurs less or equal times than V+1 on the second argument of P
lesseq(V,V+1,A,I,P,2) :- V = 1..2, action(A), I = 1..2, pred(P),
   V = 1 : sreversible(A);
   0 #sum{ -1,VV : a_atom(I,A,(P,(VV,  V)));
            1,VV : a_atom(I,A,(P,(VV,V+1)))},
   opt_synthesis = 1,
   opt_asym = 3.

% V occurs less times than V+1 on the first argument of P
less(V,V+1,A,I,P,1) :- V = 1..2, action(A), I = 1..2, pred(P),
   V = 1 : sreversible(A);
   1 #sum{ -1,VV : a_atom(I,A,(P,(  V,VV)));
            1,VV : a_atom(I,A,(P,(V+1,VV)))},
   opt_synthesis = 1,
   opt_asym = 3.

% V occurs less times than V+1 on the second argument of P
less(V,V+1,A,I,P,2) :- V = 1..2, action(A), I = 1..2, pred(P),
   V = 1 : sreversible(A);
   1 #sum{ -1,VV : a_atom(I,A,(P,(VV,  V)));
            1,VV : a_atom(I,A,(P,(VV,V+1)))},
   opt_synthesis = 1,
   opt_asym = 3.

% It cannot be that object O occurs less than O+1 as argument I of P,
% and it occurs less than or equal on the previous pairs (PP,II)
:- V = 1..2, action(A),
   V = 1 : sreversible(A);
   less(V,V+1,A,I,P,J),
   lesseq(V,V+1,A,II,PP,JJ) : II = 1..2, pred(PP), JJ = 1..2, (II,PP,JJ) < (I,P,J);
   opt_synthesis = 1,
   opt_asym = 3.


