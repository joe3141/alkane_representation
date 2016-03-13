add_branch_to_carbon(carb(c,h,h,c),N,Ans):-
                                          branch_name(N,B),
                                          Ans = carb(c,B,h,c).

add_branch_to_carbon(carb(c,H,h,c),N,Ans):-
                                          \+H = h,
                                          branch_name(N,B),
                                          Ans = carb(c,H,B,c).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branch_name(S,N):-

                   HS is S*2+1,
                   atomic_list_concat([c,S,h,HS],N).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
straight_chain_alkane(1,A):- A =[carb(h,h,h,h)].

straight_chain_alkane(N,A):-
                            N > 0,
                            N1 is N -1,
                            straight_chain_helper(N1,T),
                            A = [carb(h,h,h,c)|T].

straight_chain_helper(1,A):- A = [carb(c,h,h,h)].

straight_chain_helper(N,A):-
                            N > 0,
                            N1 is N -1,
                            straight_chain_helper(N1,T),
                            A = [carb(c,h,h,c)|T].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq(Low,High,Low):-
                       Low=<High.
seq(Low,High,Val) :-
                   Low=<High,
                   Now is Low + 1,
                   seq(Now,High,Val).
composition(0,[]).
composition(N,[H|T]) :-
                      seq(1,N,H),
                      M is N - H,
                      composition(M,T).
break_down(N,L2):-
                 findall(R,composition(N,R),R1),
                 sortlisted(R1,K),
                 sort(K,L),
                 display1(L,L2).
display1([],_):-false.
display1([_H|T],R):-
                   display1(T,R).
display1([H|_T],R):-
                   R = H.
sortlisted([],[]).
sortlisted([H|T],R):-
                     msort(H,X),sortlisted(T,R1),append([X],R1,R).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branched_alkane(N,BA):-
                       N1 is N -1,
                     %  N2 is (N+1)//2,
                      N2 is 3,
                        seq(N2,N1,R),
                        B is N - R,
                        break_down(B,Br),
                        straight_chain_alkane(R,SC),
                        add_branches(SC,Br,BA).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

add_branches(X,[],X).
add_branches(SC,[H|T],Sol):-
                            length(SC,R),
                            R1 is R -1,
                            seq(2,R1,Index),
                            nth1(Index,SC,Carb,Rest),
                            O is R - Index +1,
                            H < O,
                            %Maha
                            H <Index,
                             R2 is ( R + 1 )//2,
                            H < R2,
                            %EndMaha
                            add_branch_to_carbon(Carb,H,Carb1),
                            nth1(Index,Sol1,Carb1,Rest),
                            add_branches(Sol1,T,Sol).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reverse_alkane([H|T],Rev):-
                           append(Mid,[Last],T),
                           reverse(Mid,MidRev),
                           append(MidRev,[Last],Rev1),
                           append([H],Rev1,Rev).
isomers_helper([],[]).
isomers_helper([H|T],BA):-
                            reverse_alkane(H,RevH),
                            member(RevH,T),
                            isomers_helper(T,BA).
isomers_helper([H|T],[H|BA]):-
                            reverse_alkane(H,RevH),
                            \+ member(RevH,T),
                            isomers_helper(T,BA).
isomers(N,A):-
               straight_chain_alkane(N,SC),
               setof(R,branched_alkane(N,R),F),
               reverse(F,L),
               isomers_helper(L,BA),
               append([SC],BA,A).